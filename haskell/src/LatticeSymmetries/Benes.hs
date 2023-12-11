{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE TemplateHaskell #-}

module LatticeSymmetries.Benes
  ( Permutation
  , unPermutation
  , mkPermutation
  , permuteVector
  , permuteBits
  , identityPermutation
  , isIdentityPermutation
  , Mapping (..)
  , permutationFromMappings
  , minimalSupport
  -- randomPermutation,
  , BenesNetwork (..)
  , permutationToBenesNetwork
  , benesPermuteBits

    -- * Hm...
  , BatchedBenesNetwork (..)
  , emptyBatchedBenesNetwork
  , mkBatchedBenesNetwork
  )
where

import Control.Monad.ST
import Data.Aeson (FromJSON (..), ToJSON (..), object, (.=))
import Data.Bits
import Data.List qualified as L
import Data.Primitive.Ptr qualified as P
import Data.Set qualified as Set
import Data.Vector qualified as B
import Data.Vector.Generic ((!))
import Data.Vector.Generic qualified as G
import Data.Vector.Generic.Mutable qualified as GM
import Data.Vector.Storable qualified as S
import Data.Vector.Storable.Mutable qualified as SM
import Data.Vector.Unboxed qualified as U
import Data.Vector.Unboxed.Deriving (derivingUnbox)
import Data.Vector.Unboxed.Mutable qualified as UM
import GHC.Exts (IsList (..))
import LatticeSymmetries.BitString
import LatticeSymmetries.Dense
import LatticeSymmetries.Utils
import System.IO.Unsafe (unsafePerformIO)

-- import System.Random

import Control.Arrow qualified
import Data.Validity (Validity (validate), check, prettyValidate)
import GHC.Records (HasField (getField))
import Prelude hiding (cycle)

-- auto const n = _info.source.size();
-- _info.inverse_source.resize(n);
-- _info.inverse_target.resize(n);
-- for (auto i = size_t{0}; i < n; ++i) {
--     _info.inverse_source[_info.source[i]] = static_cast<uint16_t>(i);
--     _info.inverse_target[_info.target[i]] = static_cast<uint16_t>(i);
-- }
-- _cxt.types.resize(n);

-- | A permutation of numbers @[0 .. N - 1]@.
newtype Permutation = Permutation {unPermutation :: U.Vector Int}
  deriving stock (Show, Eq, Ord, Generic)
  deriving anyclass (NFData)

instance ToJSON Permutation where
  toJSON (Permutation p) = toJSON p

instance FromJSON Permutation where
  parseJSON = (eitherToParser . mkPermutation) <=< parseJSON

instance IsList Permutation where
  type Item Permutation = Int
  toList (Permutation p) = G.toList p
  fromList p = either error id $ mkPermutation (G.fromList p)

-- | Rearrange elements of the input vector according to the given permutation.
permuteVector
  :: (HasCallStack, G.Vector v a)
  => Permutation
  -- ^ Specifies how to order elements
  -> v a
  -- ^ Input vector
  -> v a
  -- ^ Rearranged vector
permuteVector (Permutation p) xs
  | G.length p == G.length xs = G.generate (G.length p) (G.unsafeIndex xs . G.unsafeIndex p)
  | otherwise = error $ "length mismatch: " <> show (G.length p) <> " != " <> show (G.length xs)

permuteBits :: Permutation -> Integer -> Integer
permuteBits (Permutation p) x = go 0 zeroBits
  where
    go !i !y
      | i < G.length p =
          let y' = if testBit x (p ! i) then setBit y i else y
           in go (i + 1) y'
      | otherwise = y

instance HasField "length" Permutation Int where
  -- \| Get the length of the permutation. If we are given a permutation of numbers @[0 .. N-1]@, then
  -- this function will return @N@.
  getField (Permutation p) = G.length p

-- | Generate the identity permutation of given length.
identityPermutation
  :: (HasCallStack)
  => Int
  -- ^ Length of the permutation
  -> Permutation
identityPermutation size
  | size >= 0 = Permutation $ G.generate size id
  | otherwise = error $ "invalid size: " <> show size

-- | Create a permutation from vector.
mkPermutation :: U.Vector Int -> Either Text Permutation
mkPermutation p = Control.Arrow.left fromString $ prettyValidate (Permutation p)

instance Validity Permutation where
  validate (Permutation p) =
    mconcat
      [ check (not (G.null p)) "p is not empty"
      , check
          (Set.toAscList (Set.fromList (G.toList p)) == [0 .. G.length p - 1])
          ("p is a permutation of [0 .. " <> show (G.length p - 1) <> "]")
      ]

instance Semigroup Permutation where
  (<>) x (Permutation ys) = Permutation $ permuteVector x ys

data Mapping a = Mapping !a !a
  deriving stock (Show, Eq, Ord)

derivingUnbox
  "Mapping"
  [t|forall a. (UM.Unbox a) => Mapping a -> (a, a)|]
  [|\(Mapping a b) -> (a, b)|]
  [|\(!a, !b) -> Mapping a b|]

-- | Create a permutation from a list of index mappings, e.g. [1->2, 3->5, 2->1].
-- Assumes that the list of mappings is well-formed.
permutationFromMappings :: (HasCallStack) => Maybe Int -> [Mapping Int] -> Permutation
permutationFromMappings maybeN mappings = runST $ do
  p <- UM.generate n id
  forM_ mappings $ \(Mapping x y) ->
    UM.write p x y
  either error id . mkPermutation <$> U.unsafeFreeze p
  where
    !n = fromMaybe ((+ 1) . L.maximum $ (\(Mapping x y) -> max x y) <$> mappings) maybeN

minimalSupport :: Permutation -> Maybe Int
minimalSupport = G.findIndex (uncurry (/=)) . G.indexed . unPermutation

isIdentityPermutation :: Permutation -> Bool
isIdentityPermutation = G.all (uncurry (==)) . G.indexed . unPermutation

{-
randomShuffle :: (RandomGen g, G.Vector v a) => v a -> g -> (v a, g)
randomShuffle xs g₀ = runST $ do
  v <- G.thaw xs
  let n = G.length xs
  buffer <- GM.new n
  let go !i !g
        | i < n = do
            let (j, g') = uniformR (i, n - 1) g
            vᵢ <- GM.read v i
            vⱼ <- GM.read v j
            GM.write v j vᵢ
            GM.write buffer i vⱼ
            go (i + 1) g'
        | otherwise = pure g
  g' <- go 0 g₀
  xs' <- G.unsafeFreeze buffer
  pure (xs', g')

-- | Generate a random permutation
randomPermutation :: RandomGen g => Int -> g -> (Permutation, g)
randomPermutation n g = (mkPermutation p, g')
  where
    (p, g') = randomShuffle (G.enumFromN 0 n) g
-}

newtype Index = Index Int
  deriving stock (Show, Eq, Ord)

data Focus = Source | Target
  deriving stock (Show, Eq, Ord)

otherLocation :: Focus -> Focus
otherLocation Source = Target
otherLocation Target = Source

newtype Value = Value Int
  deriving stock (Show, Eq, Ord)

data Point = Point !Index !Value
  deriving stock (Show, Eq, Ord)

data Edge = Edge !Focus !Point !Point
  deriving stock (Show, Eq, Ord)

data InvertiblePermutation = InvertiblePermutation
  {ipPermutation :: !(U.Vector Int), ipInverse :: !(U.Vector Int)}

ipValue :: InvertiblePermutation -> Index -> Value
ipValue p (Index i) = Value (ipPermutation p ! i)

ipIndex :: InvertiblePermutation -> Value -> Index
ipIndex p (Value i) = Index (ipInverse p ! i)

mkInvertiblePermutation :: Permutation -> InvertiblePermutation
mkInvertiblePermutation (Permutation p) = InvertiblePermutation p inverse
  where
    inverse = runST $ do
      v <- GM.new (G.length p)
      G.iforM_ p $ \i pᵢ ->
        GM.write v pᵢ i
      G.unsafeFreeze v

isSmaller :: Int -> Index -> Bool
isSmaller δ (Index i) = i `mod` (2 * δ) < δ

getNeighbor :: Int -> Index -> Index
getNeighbor δ i@(Index i')
  | isSmaller δ i = Index (i' + δ)
  | otherwise = Index (i' - δ)

getCycle :: Int -> InvertiblePermutation -> InvertiblePermutation -> Edge -> [Edge]
getCycle !δ !src !tgt edge₀@(Edge loc₀ _ b₀) = edge₀ : go loc₀ b₀
  where
    getIndex Source x = ipIndex src x
    getIndex Target x = ipIndex tgt x
    getValue Source i = ipValue src i
    getValue Target i = ipValue tgt i

    go loc (Point _ x₀)
      | edge == edge₀ = [edge₀]
      | x₀ == xᵢ = edge : go loc' pⱼ
      | x₀ == xⱼ = edge : go loc' pᵢ
      | otherwise = error $ "this should not have happened: " <> show loc <> ", " <> show edge
      where
        loc' = otherLocation loc
        i = getIndex loc' x₀
        xᵢ = getValue loc' i
        pᵢ = Point i xᵢ
        j = getNeighbor δ i
        xⱼ = getValue loc' j
        pⱼ = Point j xⱼ
        edge = if isSmaller δ i then Edge loc' pᵢ pⱼ else Edge loc' pⱼ pᵢ

data Swap = Swap !Focus !Index !Index
  deriving stock (Show, Eq)

solveCycle :: [Edge] -> [Swap]
solveCycle = go []
  where
    shouldSwap (Edge _ (Point _ a) (Point _ b)) (Edge _ (Point _ c) (Point _ d))
      | a == c || b == d = False
      | a == d || b == c = True
      | otherwise = error "should not have happened"
    go acc [Edge _ (Point i _) (Point j _)]
      | i < j = acc
      | otherwise = error "unsolvable cycle"
    go acc (e₀ : e₁ : es)
      | shouldSwap e₀ e₁ =
          let (Edge loc₁ pᵢ@(Point i _) pⱼ@(Point j _)) = e₁
           in go (Swap loc₁ i j : acc) (Edge loc₁ pⱼ pᵢ : es)
      | otherwise = go acc (e₁ : es)
    go _ _ = error "should not have happened"

getMask :: [Swap] -> Integer
getMask = go zeroBits
  where
    go !acc [] = acc
    go !acc ((Swap _ (Index i) (Index j)) : others) = go (setBit acc (min i j)) others

applySwaps :: InvertiblePermutation -> [Swap] -> InvertiblePermutation
applySwaps perm swaps = runST $ do
  p' <- G.thaw $ ipPermutation perm
  inverse' <- G.thaw $ ipInverse perm
  forM_ swaps $ \s -> applySwap s p' inverse'
  InvertiblePermutation <$> G.unsafeFreeze p' <*> G.unsafeFreeze inverse'
  where
    applySwap (Swap _ (Index a) (Index b)) p inverse = do
      pa <- GM.read p a
      pb <- GM.read p b
      GM.write p a pb
      GM.write p b pa
      GM.write inverse pb a
      GM.write inverse pa b

data SolverState = SolverState
  { ssSource :: !InvertiblePermutation
  , ssTarget :: !InvertiblePermutation
  , ssDelta :: !Int
  }

initialEdges :: SolverState -> [Edge]
initialEdges s = mkEdge <$> go 0
  where
    n = G.length . ipPermutation $ ssSource s
    δ = ssDelta s
    go !i
      | i + δ <= n = [(Index (i + j), Index (i + j + δ)) | j <- [0 .. δ - 1]] <> go (i + 2 * δ)
      | otherwise = []
    mkEdge (!i, !j) =
      let !pᵢ = Point i (ipValue (ssSource s) i)
          !pⱼ = Point j (ipValue (ssSource s) j)
       in Edge Source pᵢ pⱼ

stageCycles :: SolverState -> [[Edge]]
stageCycles s = cycles Set.empty (initialEdges s)
  where
    cycles _ [] = []
    cycles !visited (!e : es)
      | Set.notMember e visited =
          let cycle = getCycle (ssDelta s) (ssSource s) (ssTarget s) e
           in cycle : cycles (visited `Set.union` Set.fromList cycle) es
      | otherwise = cycles visited es

solveStage :: SolverState -> (Integer, Integer, SolverState)
solveStage !s = (srcMask, tgtMask, s')
  where
    swaps = mconcat $ solveCycle <$> stageCycles s
    (srcSwaps, tgtSwaps) = L.partition (\(Swap loc _ _) -> loc == Source) swaps
    !srcMask = getMask srcSwaps
    !tgtMask = getMask tgtSwaps
    !s' =
      SolverState
        (applySwaps (ssSource s) srcSwaps)
        (applySwaps (ssTarget s) tgtSwaps)
        (2 * ssDelta s)

solve :: InvertiblePermutation -> InvertiblePermutation -> ([Integer], [Integer], [Int])
solve src tgt = unpack $ go s₀
  where
    n = G.length (ipPermutation src)
    unpack = foldr (\(a, b, c) (a', b', c') -> (a : a', b : b', c : c')) ([], [], [])
    s₀ = SolverState src tgt 1
    go !s
      | 2 * ssDelta s <= n =
          let (a, b, s') = solveStage s
           in (a, b, ssDelta s) : go s'
      | otherwise = []

data BenesNetwork = BenesNetwork {bnMasks :: !(B.Vector Integer), bnShifts :: !(U.Vector Int)}
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

instance ToJSON BenesNetwork where
  toJSON x = object ["masks" .= x.bnMasks, "shifts" .= x.bnShifts]

mkBenesNetwork :: [Integer] -> [Integer] -> [Int] -> BenesNetwork
mkBenesNetwork srcMasks tgtMasks δs
  | null δs = BenesNetwork G.empty G.empty
  | isOkay =
      BenesNetwork
        (G.fromList $ L.init srcMasks <> reverse tgtMasks)
        (G.fromList $ δs <> drop 1 (reverse δs))
  | otherwise = error $ "invalid backward masks: " <> show srcMasks <> ", " <> show tgtMasks <> ", " <> show δs
  where
    isOkay = case srcMasks of
      [] -> True
      _ -> L.last srcMasks == zeroBits

extendToPowerOfTwo :: Permutation -> Permutation
extendToPowerOfTwo (Permutation p) = Permutation (G.generate n f)
  where
    !n = ((2 :: Int) ^) $ (ceiling :: Double -> Int) $ logBase 2 (fromIntegral (G.length p))
    f !i
      | i < G.length p = p G.! i
      | otherwise = i

permutationToBenesNetwork :: Permutation -> BenesNetwork
permutationToBenesNetwork p = mkBenesNetwork srcMasks tgtMasks shifts
  where
    p' = extendToPowerOfTwo p
    (srcMasks, tgtMasks, shifts) =
      solve
        (mkInvertiblePermutation (identityPermutation p'.length))
        (mkInvertiblePermutation p')

bitPermuteStep :: Integer -> Integer -> Int -> Integer
bitPermuteStep x θ δ = (x `xor` y) `xor` (y `shiftL` δ)
  where
    y = ((x `shiftR` δ) `xor` x) .&. θ

benesPermuteBits :: BenesNetwork -> Integer -> Integer
benesPermuteBits (BenesNetwork masks δs) = go 0
  where
    n = G.length δs
    go !i !x
      | i < n = go (i + 1) (bitPermuteStep x (G.unsafeIndex masks i) (G.unsafeIndex δs i))
      | otherwise = x

data BatchedBenesNetwork = BatchedBenesNetwork
  { bbnMasks :: {-# UNPACK #-} !(DenseMatrix S.Vector Word64)
  , bbnShifts :: {-# UNPACK #-} !(S.Vector Word64)
  }
  deriving stock (Show, Eq)

emptyBatchedBenesNetwork :: BatchedBenesNetwork
emptyBatchedBenesNetwork = mkBatchedBenesNetwork G.empty

mkBatchedBenesNetwork :: B.Vector BenesNetwork -> BatchedBenesNetwork
mkBatchedBenesNetwork networks
  | G.null networks = BatchedBenesNetwork (DenseMatrix 0 0 G.empty) G.empty
  | sameShifts = unsafePerformIO $ do
      masks <- GM.new (numberShifts * numberMasks * numberWords)
      SM.unsafeWith masks $ \masksPtr ->
        loopM 0 (< numberShifts) (+ 1) $ \i ->
          writeManyBitStrings
            numberWords
            (P.advancePtr masksPtr (i * numberMasks * numberWords))
            (getBitStrings i)
      masks' <- DenseMatrix numberShifts (numberMasks * numberWords) <$> G.unsafeFreeze masks
      pure $ BatchedBenesNetwork masks' (G.convert (G.map fromIntegral shifts))
  | otherwise = error "networks have different shifts"
  where
    getBitStrings i = G.toList $ G.map (\n -> BitString (n.bnMasks ! i)) networks
    sameShifts = G.all ((== shifts) . bnShifts) networks
    shifts = bnShifts (G.head networks)
    numberShifts = G.length shifts
    numberBits = 2 * G.maximum shifts
    numberWords = (numberBits + 63) `div` 64
    numberMasks = G.length networks
