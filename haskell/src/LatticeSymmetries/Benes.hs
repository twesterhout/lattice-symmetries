{-# LANGUAGE OverloadedRecordDot #-}

module LatticeSymmetries.Benes
  ( BenesNetwork (..)
  , permutationToBenesNetwork
  , benesPermuteBits
  , BatchedBenesNetwork (..)
  , emptyBatchedBenesNetwork
  , mkBatchedBenesNetwork
  )
where

import Control.Monad.ST (runST)
import Data.Aeson (ToJSON (..), object, (.=))
import Data.Bits (Bits (setBit, shiftL, shiftR, zeroBits, (.&.)))
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
import LatticeSymmetries.BitString
import LatticeSymmetries.Dense
import LatticeSymmetries.Permutation
import LatticeSymmetries.Utils
import System.IO.Unsafe (unsafePerformIO)
import Prelude hiding (cycle)

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
mkInvertiblePermutation (unPermutation -> p) = InvertiblePermutation p inverse
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

data BenesNetwork = BenesNetwork {bnMasks :: !(B.Vector BitString), bnShifts :: !(U.Vector Int)}
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

instance ToJSON BenesNetwork where
  toJSON x = object ["masks" .= x.bnMasks, "shifts" .= x.bnShifts]

mkBenesNetwork :: [Integer] -> [Integer] -> [Int] -> BenesNetwork
mkBenesNetwork srcMasks tgtMasks δs
  | null δs = BenesNetwork G.empty G.empty
  | isOkay =
      BenesNetwork
        (G.fromList . fmap BitString $ L.init srcMasks <> reverse tgtMasks)
        (G.fromList $ δs <> drop 1 (reverse δs))
  | otherwise = error $ "invalid backward masks: " <> show srcMasks <> ", " <> show tgtMasks <> ", " <> show δs
  where
    isOkay = case srcMasks of
      [] -> True
      _ -> L.last srcMasks == zeroBits

extendToPowerOfTwo :: Permutation -> Permutation
extendToPowerOfTwo (unPermutation -> p) = either error id $ mkPermutation (G.generate n f)
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

bitPermuteStep :: BitString -> BitString -> Int -> BitString
bitPermuteStep x θ δ = (x `xor` y) `xor` (y `shiftL` δ)
  where
    y = ((x `shiftR` δ) `xor` x) .&. θ

benesPermuteBits :: BenesNetwork -> BitString -> BitString
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

instance ToJSON BatchedBenesNetwork where
  toJSON x = object ["bbnMasks" .= denseMatrixToList x.bbnMasks, "bbnShifts" .= G.toList x.bbnShifts]

emptyBatchedBenesNetwork :: BatchedBenesNetwork
emptyBatchedBenesNetwork = mkBatchedBenesNetwork G.empty

mkBatchedBenesNetwork :: HasCallStack => B.Vector BenesNetwork -> BatchedBenesNetwork
mkBatchedBenesNetwork networks
  | G.null networks || numberShifts * numberShifts == 0 = BatchedBenesNetwork (DenseMatrix 0 0 G.empty) G.empty
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
    getBitStrings i = G.toList $ G.map (\n -> n.bnMasks ! i) networks
    sameShifts = G.all ((== shifts) . bnShifts) networks
    shifts = bnShifts (G.head networks)
    numberShifts = G.length shifts
    numberBits = 2 * G.maximum shifts
    numberWords = (numberBits + 63) `div` 64
    numberMasks = G.length networks
