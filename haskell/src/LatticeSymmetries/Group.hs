{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE OverloadedRecordDot #-}

module LatticeSymmetries.Group
  ( RepElement (..)
  , Representation (..)
  , AbelianPermutationGroup (..)
  , abelianSubgroup
  , abelianRepresentations
  , groupRepresentations

    -- * Helpers
  , Symmetry
  , mkSymmetry
  , fromGenerators
  , isRepElementReal
  , isRepresentationReal

    -- * Hm...

  -- , Symmetries (..)
  -- , areSymmetriesReal
  -- , emptySymmetries
  , emptyRepresentation
  , nullRepresentation
  -- , compileGroupRepresentation

    -- ** Low-level interface for FFI

  -- , newCpermutation_group
  -- , destroyCpermutation_group

    -- * Automorphisms

  -- , mkMultiplicationTable
  -- , abelianRepresentations
  -- , GroupElement (..)
  -- , shrinkMultiplicationTable
  -- , getGroupElements
  -- , selectPrimeElements
  -- , groupGenerators
  -- , groupElementCycle
  -- , groupElementOrder
  -- , Coset (..)
  -- , abelianization
  -- , commutatorSubgroup

    -- * Reference implementation
  , isRepresentativeSlow
  , stateInfoSlow
  )
where

import Control.Arrow (left)
import Control.Monad.ST.Strict (runST)
import Data.Aeson (FromJSON (..), ToJSON (..), object, withObject, (.:), (.=))
import Data.Complex
import Data.IntSet qualified as IntSet
import Data.List qualified
import Data.List.NonEmpty qualified
import Data.Map.Strict qualified as Map
import Data.Ratio
import Data.Validity
import Data.Vector qualified as B
import Data.Vector.Generic ((!))
import Data.Vector.Generic qualified as G
import Data.Vector.Unboxed qualified as U
import Data.Vector.Unboxed.Mutable qualified as UM
import GHC.Records (HasField (..))
import LatticeSymmetries.Automorphisms
import LatticeSymmetries.Dense
import LatticeSymmetries.Generator
import LatticeSymmetries.Permutation
import LatticeSymmetries.Utils (eitherToParser, sortVectorBy)
import Prelude hiding (group, identity, permutations, second, toList)

-- | Representation of a permutation
-- data Symmetry = Symmetry
--   { permutation :: !Permutation
--   , phase :: !(Ratio Int)
--   }
--   deriving stock (Show, Eq, Ord, Generic)
--   deriving anyclass (NFData)
type Symmetry = RepElement Permutation

instance HasPeriodicity a => HasField "sector" (RepElement a) Int where
  getField s
    | denominator sector == 1 = numerator sector
    | otherwise = error "this should not have happened: invalid RepElement"
    where
      sector = s.phase * fromIntegral (getPeriodicity s.element)

instance HasField "length" a Int => HasField "size" (RepElement a) Int where
  getField = (.length) . (.element)

phaseToCharacter :: Ratio Int -> Complex Double
phaseToCharacter φ
  | φ == 0 = 1 :+ 0
  | φ == 1 % 4 = 0 :+ (-1)
  | φ == 1 % 2 = (-1) :+ 0
  | φ == 3 % 4 = 0 :+ 1
  | 0 <= φ && φ < 1 % 4 = cos (2 * pi * realToFrac φ) :+ (-sin (2 * pi * realToFrac φ))
  | 1 % 4 <= φ && φ < 1 % 2 = (-cos (2 * pi * realToFrac (1 % 2 - φ))) :+ (-sin (2 * pi * realToFrac (1 % 2 - φ)))
  | 1 % 2 <= φ && φ < 3 % 4 = (-cos (2 * pi * realToFrac (φ - 1 % 2))) :+ sin (2 * pi * realToFrac (φ - 1 % 2))
  | 3 % 4 <= φ && φ < 1 = cos (2 * pi * realToFrac (1 - φ)) :+ sin (2 * pi * realToFrac (1 - φ))
  | otherwise = error "should never happen"

instance HasField "character" (RepElement a) (Complex Double) where
  getField ((.phase) -> φ) = phaseToCharacter φ

instance FromJSON Symmetry where
  parseJSON = withObject "Symmetry" $ \v ->
    eitherToParser =<< (mkSymmetry <$> v .: "permutation" <*> v .: "sector")

instance ToJSON Symmetry where
  toJSON s = object ["__type__" .= ("Symmetry" :: Text), "permutation" .= s.element, "sector" .= s.sector]

-- | Create a new 'Symmetry'.
mkSymmetry :: Permutation -> Int -> Either Text (RepElement Permutation)
mkSymmetry p k = Control.Arrow.left fromString . prettyValidate $ RepElement p (modOne (k % getPeriodicity p))

instance Validity (RepElement Permutation) where
  validate x =
    check (x.phase >= 0) "phase φ is non-negative"
      <> check (x.phase < 1) "phase φ is less than one"
      <> check (denominator (x.phase * fromIntegral (getPeriodicity x.element)) == 1) "phase is consistent with periodicity"

modOne :: Integral a => Ratio a -> Ratio a
modOne x
  | x >= 1 = x - fromIntegral (numerator x `div` denominator x)
  | x < 0 = modOne (x + fromIntegral (denominator x))
  | otherwise = x

-- instance Semigroup Symmetry where
--   (<>) (Symmetry pa λa) (Symmetry pb λb) = undefined -- uncurry Symmetry (combine (pa, λa) (pb, λb))

data RepElement a = RepElement {element :: !a, phase :: !(Ratio Int)}
  deriving stock (Show, Eq, Ord, Generic)
  deriving anyclass (NFData)

newtype ComparingElement a = ComparingElement (RepElement a)

instance Eq a => Eq (ComparingElement a) where
  (ComparingElement a) == (ComparingElement b) = a.element == b.element

instance Ord a => Ord (ComparingElement a) where
  compare (ComparingElement a) (ComparingElement b) = compare a.element b.element

instance Semigroup a => Semigroup (RepElement a) where
  (<>) (RepElement pa λa) (RepElement pb λb) = RepElement (pa <> pb) (modOne (λa + λb))

newtype Representation a = Representation {unRepresentation :: B.Vector (RepElement a)}
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

instance HasField "length" a Int => HasField "numberBits" (Representation a) (Maybe Int) where
  getField (unRepresentation -> r)
    | G.null r = Nothing
    | otherwise = Just (G.head r).size

unionWithConflictsM :: (Ord k, Eq a) => Map.Map k a -> Map.Map k a -> Either (k, a, a) (Map.Map k a)
unionWithConflictsM a0 b0 = go a0 (Map.toList b0)
  where
    go !acc ((k, e) : rest)
      | Just e' <- Map.lookup k acc = if e' == e then go acc rest else Left (k, e', e)
      | otherwise = go (Map.insert k e acc) rest
    go !acc [] = Right acc

differenceWithConflictsM :: (Ord k, Eq a) => Map.Map k a -> Map.Map k a -> Either (k, a, a) (Map.Map k a)
differenceWithConflictsM a0 b0 = go a0 (Map.toList b0)
  where
    go !acc ((k, e) : rest)
      | Just e' <- Map.lookup k acc = if e' == e then go (Map.delete k acc) rest else Left (k, e', e)
      | otherwise = go acc rest
    go !acc [] = Right acc

fromListWithConflictsM :: (Ord k, Eq a) => [(k, a)] -> Either (k, a, a) (Map.Map k a)
fromListWithConflictsM = go Map.empty
  where
    go !acc ((k, e) : rest)
      | Just e' <- Map.lookup k acc = if e' == e then go acc rest else Left (k, e', e)
      | otherwise = go (Map.insert k e acc) rest
    go !acc [] = Right acc

fromGenerators :: forall a. (Semigroup a, Ord a, Show a) => [RepElement a] -> Either Text (Representation a)
fromGenerators [] = Right . Representation $ G.empty
fromGenerators gs = Control.Arrow.left conflictToText $ go Map.empty =<< fromListWithConflictsM (to <$> gs)
  where
    to x = (x.element, x.phase)
    from (p, φ) = RepElement p φ
    conflictToText (p, φ1, φ2) = "conflicting phases for group element " <> show p <> ": " <> show φ1 <> " != " <> show φ2
    go !interior !boundary
      | Map.null boundary = Right . Representation . G.fromList . fmap from . Map.toList $ interior
      | otherwise = do
          interior' <- unionWithConflictsM interior boundary
          everything <- fromListWithConflictsM [to (from h <> g) | h <- Map.toList boundary, g <- gs]
          boundary' <- differenceWithConflictsM everything interior'
          go interior' boundary'

isRepElementReal :: RepElement a -> Bool
isRepElementReal x = x.phase == 0 || x.phase == 1 % 2

isRepresentationReal :: Representation a -> Bool
isRepresentationReal = G.all isRepElementReal . unRepresentation

-- compileGroupRepresentation :: HasCallStack => B.Vector Symmetry -> Symmetries
-- compileGroupRepresentation symmetries
--   | G.null symmetries = emptySymmetries
--   | (G.head symmetries).size <= 1 = emptySymmetries
--   | otherwise = Symmetries permGroup benesNetwork charactersReal charactersImag symmetries
--   where
--     permutations = (.permutation) <$> symmetries
--     permGroup = PermutationGroup permutations
--     benesNetwork = mkBatchedBenesNetwork $ G.map permutationToBenesNetwork permutations
--     characters = G.convert $ (.character) <$> symmetries
--     charactersReal = G.map realPart characters
--     charactersImag = G.map imagPart characters

-- mkSymmetriesFromRepresentation :: (HasCallStack) => PermutationGroup -> Representation -> Symmetries
-- mkSymmetriesFromRepresentation group (Representation r) =
--   either error id $ mkSymmetries . G.toList $ uncurry Symmetry <$> gs
--   where
--     gs = let (PermutationGroup g) = group in (\(x, φ) -> (g ! x.index, φ)) <$> r

-- instance FromJSON Symmetries where
--   parseJSON =
--     either (fail . toString) (pure . compileGroupRepresentation)
--       <=< fmap groupRepresentationFromGenerators . parseJSON

-- instance ToJSON Symmetries where
--   toJSON s = toJSON s.symmOriginal
-- instance ToJSON Symmetries where
--   toJSON s =
--     object
--       [ "symmOriginal" .= s.symmOriginal
--       , "symmNetwork" .= s.symmNetwork
--       ]

-- instance IsList Symmetries where
--   type Item Symmetries = Symmetry
--   toList = toSymmetryList
--   fromList x = either error id (mkSymmetries x)

-- instance HasField "numberBits" Symmetries (Maybe Int) where
--   getField s
--     | nullSymmetries s = Nothing
--     | otherwise = Just . G.length . unPermutation . G.head $ s.symmGroup.permutations

emptyRepresentation :: Representation a
emptyRepresentation = Representation G.empty

nullRepresentation :: Representation a -> Bool
nullRepresentation = G.null . unRepresentation

-- areSymmetriesReal :: Symmetries -> Bool
-- areSymmetriesReal = G.all (== 0) . symmCharactersImag

permuteBasisState :: Permutation -> BasisState t -> BasisState t
permuteBasisState p (BasisState n bits)
  | p.length == n = BasisState n (permuteBits p bits)
  | otherwise = error $ "permuteBasisState: size mismatch: " <> show p.length <> " != " <> show n

isRepresentativeSlow
  :: Representation Permutation
  -- ^ Symmetry group
  -> Maybe Int
  -- ^ Spin inversion
  -> BasisState t
  -- ^ Basis state
  -> Maybe Int
  -- ^ Norm if the state is a representative, nothing otherwise
isRepresentativeSlow (B.toList . unRepresentation -> group) spinInversion basisState1 = case spinInversion of
  Nothing -> go 0 basisState1 group
  Just _ ->
    let count1 = go 0 basisState1 group
        count2 = go 0 basisState2 group
     in (+) <$> count1 <*> count2
  where
    basisState2 = invertBasisState basisState1
    go !count _ [] = Just count
    go !count bits ((RepElement p φ) : rest) = case compare bits' bits of
      LT -> Nothing
      EQ
        | φ /= 0 -> Nothing
        | otherwise -> go (count + 1) bits rest
      GT -> go count bits rest
      where
        bits' = permuteBasisState p bits

stateInfoSlow
  :: Representation Permutation
  -- ^ Symmetry group
  -> Maybe Int
  -- ^ Spin inversion
  -> BasisState t
  -- ^ Basis state
  -> (BasisState t, Int)
  -- ^ Norm if the state is a representative, nothing otherwise
stateInfoSlow (B.toList . unRepresentation -> group) spinInversion basisState1 = case spinInversion of
  Nothing -> go (basisState1, -1) 0 basisState1 group
  Just _ -> error "not implemented"
  -- let count1 = go 0 basisState1 group
  --     count2 = go 0 basisState2 group
  --  in (+) <$> count1 <*> count2
  where
    basisState2 = invertBasisState basisState1
    go (!rep, !index) _ _ [] = (rep, index)
    go (!rep, !index) !k !original ((RepElement p _) : rest) =
      case compare permuted rep of
        LT -> go (permuted, k) (k + 1) original rest
        _ -> go (rep, index) (k + 1) original rest
      where
        permuted = permuteBasisState p original

-- toSymmetryList :: (HasCallStack) => Symmetries -> [Symmetry]
-- toSymmetryList (Symmetries (PermutationGroup gs) _ λsRe λsIm) =
--   Data.List.zipWith3 toSymmetry (G.toList gs) (G.toList λsRe) (G.toList λsIm)
--   where
--     toSymmetry :: (HasCallStack) => Permutation -> Double -> Double -> Symmetry
--     toSymmetry g λRe λIm
--       | r ≈ 1 && s' ≈ fromIntegral (round s' :: Int) = mkSymmetry g (round s')
--       | otherwise = error $ "failed to reconstruct Symmetry from " <> show (toList g) <> " and λ = " <> show λRe <> " + " <> show λIm <> "j"
--       where
--         -- φ = sector / periodicity, but care needs to be taken because φ is approximate
--         s' = φ * fromIntegral n
--         n = getPeriodicity g
--         -- polar returns the phase in (-π, π], we add 2π to get a positive value
--         (r, _φ) = polar (λRe :+ λIm)
--         φ = _φ / (2 * pi) + (if _φ < 0 then 1 else 0)

-- newCpermutation_group :: Symmetries -> IO (Ptr Cpermutation_group)
-- newCpermutation_group x =
--   let (DenseMatrix numberShifts numberMasks masks) = bbnMasks (symmNetwork x)
--    in case x.numberBits of
--         Just nBits ->
--           S.unsafeWith masks $ \masksPtr ->
--             S.unsafeWith (x.symmNetwork.bbnShifts) $ \shiftsPtr ->
--               S.unsafeWith (symmCharactersReal x) $ \eigvalsRealPtr ->
--                 S.unsafeWith (symmCharactersImag x) $ \eigvalsImagPtr ->
--                   new
--                     $ Cpermutation_group
--                       { cpermutation_group_refcount = 0
--                       , cpermutation_group_number_bits = fromIntegral nBits
--                       , cpermutation_group_number_shifts = fromIntegral numberShifts
--                       , cpermutation_group_number_masks = fromIntegral numberMasks
--                       , cpermutation_group_masks = masksPtr
--                       , cpermutation_group_shifts = shiftsPtr
--                       , cpermutation_group_eigvals_re = (castPtr :: Ptr Double -> Ptr CDouble) eigvalsRealPtr
--                       , cpermutation_group_eigvals_im = (castPtr :: Ptr Double -> Ptr CDouble) eigvalsImagPtr
--                       , cpermutation_group_haskell_payload = nullPtr
--                       }
--         Nothing -> pure nullPtr

newtype AbelianPermutationGroup = AbelianPermutationGroup {unAbelianPermutationGroup :: PermutationGroup}
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

-- destroyCpermutation_group :: Ptr Cpermutation_group -> IO ()
-- destroyCpermutation_group p = when (p /= nullPtr) $ free p
abelianSubgroup :: PermutationGroup -> AbelianPermutationGroup
abelianSubgroup (PermutationGroup g t) = AbelianPermutationGroup (PermutationGroup g' t')
  where
    -- abelianSubset returns candidates; we consider the first 100 and select the largest
    indices = Data.List.maximumBy (comparing G.length) . take 100 $ abelianSubset t
    g' = G.map (g G.!) . G.convert $ indices
    t' = shrinkMultiplicationTable t indices

abelianRepresentations :: AbelianPermutationGroup -> B.Vector (Representation Permutation)
abelianRepresentations (AbelianPermutationGroup (PermutationGroup g t)) =
  G.fromList $ Representation . fmap repToSymmetry . unRepresentation <$> abelianGroupRepresentations t
  where
    repToSymmetry x = RepElement (g ! x.element.index) x.phase

groupRepresentations :: PermutationGroup -> B.Vector (Representation Permutation)
groupRepresentations group0 = G.fromList $ unpackRep <$> abelianReps
  where
    (abelianGroup, abelianTable) = abelianization group0.table
    abelianReps = abelianGroupRepresentations abelianTable
    unpackCoset :: RepElement GroupElement -> [(Int, RepElement Permutation)]
    unpackCoset x = f <$> indices
      where
        (Coset (IntSet.toList -> indices)) = abelianGroup ! x.element.index
        f i = (i, RepElement (group0.permutations ! i) x.phase)
    unpackRep :: Representation GroupElement -> Representation Permutation
    unpackRep = either error id . fromGenerators . fmap snd . sortOn fst . concatMap unpackCoset . unRepresentation

data AbelianSubsetHistory = AbelianSubsetHistory
  { included :: !(U.Vector Bool)
  , mask :: !(U.Vector Bool)
  }
  deriving stock (Show)

commMatrixRow :: MultiplicationTable -> Int -> U.Vector Bool
commMatrixRow (MultiplicationTable matrix) i = U.generate matrix.dmCols $ \k ->
  indexDenseMatrix matrix (i, k) == indexDenseMatrix matrix (k, i)

addToAbelianSubsetHistory :: MultiplicationTable -> AbelianSubsetHistory -> Int -> AbelianSubsetHistory
addToAbelianSubsetHistory t h i = AbelianSubsetHistory included' mask'
  where
    included' = G.modify (\v -> UM.write v i True) h.included
    mask' = G.zipWith (&&) h.mask (commMatrixRow t i)

mergeAbelianSubsetHistories :: [AbelianSubsetHistory] -> [AbelianSubsetHistory]
mergeAbelianSubsetHistories = fmap go . Data.List.groupBy (\a b -> a.mask == b.mask)
  where
    go :: [AbelianSubsetHistory] -> AbelianSubsetHistory
    go hs = runST $ case hs of
      (h : rest) -> do
        mh <- G.thaw h.included
        forM_ rest $ \h' -> G.iforM_ h'.included $ \i v ->
          UM.modify mh (|| v) i
        AbelianSubsetHistory <$> G.unsafeFreeze mh <*> pure h.mask
      [] -> error "should never happen by construction"

abelianSubset :: MultiplicationTable -> [U.Vector Int]
abelianSubset t = (\h -> sortVectorBy compare . G.map fst . G.filter snd . G.indexed $ h.included) <$> go h0
  where
    !n = dmRows . unMultiplicationTable $ t
    h0 = AbelianSubsetHistory (G.replicate n False) (G.replicate n True)

    getChoices h =
      G.toList . G.filter (not . (h.included G.!)) . G.map fst . G.filter snd . G.indexed $ h.mask

    go :: AbelianSubsetHistory -> [AbelianSubsetHistory]
    go h = case getChoices h of
      [] -> [h]
      choices ->
        let hs' =
              Data.List.sortOn (negate . G.length . G.filter id . (.included))
                . mergeAbelianSubsetHistories
                . Data.List.sortOn (.mask)
                $ addToAbelianSubsetHistory t h
                  <$> choices
         in concatMap go hs'

data GroupElement = GroupElement {index :: !Int, table :: !MultiplicationTable}
  deriving stock (Show)

instance Semigroup GroupElement where
  -- NOTE: a.table and b.table should be the same
  a <> b = GroupElement index' a.table
    where
      index' = unsafeIndexDenseMatrix (unMultiplicationTable a.table) (a.index, b.index)

instance Eq GroupElement where
  a == b = a.index == b.index

instance Ord GroupElement where
  compare = comparing (.index)

getGroupElements :: MultiplicationTable -> [GroupElement]
getGroupElements t@(MultiplicationTable (DenseMatrix n _ _)) = [GroupElement i t | i <- [0 .. n - 1]]

primeFactors :: Int -> [Int]
primeFactors = go 2
  where
    go :: Int -> Int -> [Int]
    go !d !n
      | n <= 1 = []
      | d * d > n = [n]
      | n `mod` d == 0 = d : go d (n `div` d)
      | otherwise = go (d + 1) n

combineFactors :: [Int] -> [(Int, Int)]
combineFactors xs = (\l -> (head l, length l)) <$> Data.List.NonEmpty.group xs

groupElementOrder :: GroupElement -> Int
groupElementOrder x0 = go (x0 <> x0)
  where
    go !x
      | x.index == x0.index = 1
      | otherwise = 1 + go (x <> x0)

groupElementCycle :: GroupElement -> [GroupElement]
groupElementCycle x0 = x0 : go (x0 <> x0)
  where
    go !x
      | x.index == x0.index = []
      | otherwise = x : go (x <> x0)

instance HasIdentity GroupElement where
  isIdentity x = (x <> x).index == x.index

selectPrimeElements :: [GroupElement] -> [(GroupElement, Int, Int)]
selectPrimeElements = mapMaybe getInfo
  where
    getInfo x =
      let !order = groupElementOrder x
       in case combineFactors $ primeFactors order of
            [(p, n)] -> Just (x, p, n)
            _ -> Nothing

groupGenerators :: [GroupElement] -> [GroupElement]
groupGenerators group0 = go factors0 [] (selectPrimeElements group0) []
  where
    factors0 = combineFactors . primeFactors . length $ group0

    go [] gens _ _ = gens
    go ((p, n) : otherFactors) gens candidates subgroup =
      go factors' gens' candidates' subgroup'
      where
        (g, _, n') =
          Data.List.maximumBy (comparing third) $
            filter ((== p) . second) candidates
        factors'
          | n == n' = otherFactors
          | otherwise = (p, n - n') : otherFactors
        gens' = g : gens
        subgroup'
          | null subgroup = groupElementCycle g
          | otherwise = [a <> b | a <- subgroup, b <- groupElementCycle g]
        candidates' =
          let seen = IntSet.fromList $ (.index) <$> subgroup'
           in flip filter candidates $ \(cg, cp, cn) ->
                not (IntSet.member cg.index seen)
                  && cn
                    <= Data.List.maximum [fn | (fp, fn) <- factors', fp == cp]

        second (_, b, _) = b
        third (_, _, c) = c

abelianGroupRepresentationsFromGenerators :: B.Vector GroupElement -> [Representation GroupElement]
abelianGroupRepresentationsFromGenerators generators =
  fmap (either error id . fromGenerators)
    $ mapM (uncurry phases)
      . G.toList
    $ G.zip generators orders
  where
    orders = G.map groupElementOrder generators
    phases :: GroupElement -> Int -> [RepElement GroupElement]
    phases g n = [RepElement g (i % n) | i <- [0 .. n - 1]]

abelianGroupRepresentations :: MultiplicationTable -> [Representation GroupElement]
abelianGroupRepresentations t = abelianGroupRepresentationsFromGenerators generators
  where
    generators = G.fromList . groupGenerators . getGroupElements $ t

newtype Coset = Coset {unCoset :: IntSet.IntSet}
  deriving stock (Eq, Ord, Show)

getTimes :: MultiplicationTable -> Int -> Int -> Int
getTimes (MultiplicationTable matrix) = curry (unsafeIndexDenseMatrix matrix)

getInvert :: MultiplicationTable -> Int -> Int
getInvert t@(MultiplicationTable matrix) =
  \i ->
    case G.find (\j -> i `times` j == identityIndex) (U.generate n id) of
      Just j -> j
      Nothing -> error "group element has no inverse"
  where
    !n = matrix.dmRows
    !times = getTimes t
    !identityIndex =
      case G.find (\j -> isIdentity (GroupElement j t)) (U.generate n id) of
        Just k -> k
        Nothing -> error "the group has no or multiple identities"

commutatorSubgroup :: MultiplicationTable -> [Int]
commutatorSubgroup t =
  IntSet.toList $
    IntSet.fromList
      [ g `times` h `times` (inverses ! g) `times` (inverses ! h)
      | g <- [0 .. n - 1]
      , h <- [0 .. n - 1]
      ]
  where
    !n = t.unMultiplicationTable.dmRows
    !times = getTimes t
    !invert = getInvert t
    !inverses = U.generate n invert

abelianization :: MultiplicationTable -> (B.Vector Coset, MultiplicationTable)
abelianization t = (cosets, cosetMultiplicationTable t cosets)
  where
    n = t.unMultiplicationTable.dmRows
    times = getTimes t
    coset g = Coset $ IntSet.fromList [g `times` h | h <- commutatorSubgroup t]
    !cosets = go [] IntSet.empty [0 .. n - 1]
    go !cs !seen (g : gs)
      | g `IntSet.member` seen = go cs seen gs
      | otherwise =
          let !c' = coset g
              seen' = unCoset c' `IntSet.union` seen
           in go (c' : cs) seen' gs
    go cs _ [] = G.fromList cs

cosetMultiplicationTable :: MultiplicationTable -> B.Vector Coset -> MultiplicationTable
cosetMultiplicationTable t cosets =
  MultiplicationTable . generateDenseMatrix n n $ \i j ->
    let x = Data.List.head . IntSet.elems . unCoset $ cosets ! i
        y = Data.List.head . IntSet.elems . unCoset $ cosets ! j
        !z = getTimes t x y
     in case G.findIndex (\(Coset coset) -> z `IntSet.member` coset) cosets of
          Just k -> k
          Nothing -> error "G/[G, G] is not closed under <>"
  where
    n = G.length cosets
