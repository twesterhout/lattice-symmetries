module Main (main) where

import Data.Complex
import qualified Data.Vector.Storable as S
import Data.Yaml.Aeson
import Foreign.Storable
import LatticeSymmetries
import Test.Hspec

anySpinEDException :: Selector SpinEDException
anySpinEDException = const True

anyLatticeSymmetriesException :: Selector LatticeSymmetriesException
anyLatticeSymmetriesException = const True

parseLines :: FromJSON a => [Text] -> Either ParseException a
parseLines xs = decodeEither' . encodeUtf8 . unlines $ xs

checkRight s y = case s of
  Left _ -> s `shouldSatisfy` isRight
  Right x -> x `shouldBe` y

checkJust s y = case s of
  Just x -> x `shouldBe` y
  Nothing -> s `shouldSatisfy` isJust

main :: IO ()
main = hspec $ do
  describe "SymmetrySpec" $ do
    it "parses JSON specifications" $ do
      let s₁ =
            parseLines @SymmetrySpec
              [ "permutation: [3, 1, 2, 4]",
                "sector: 8"
              ]
          s₂ =
            parseLines @SymmetrySpec
              [ "permutation: [1, 1, 2, 4]",
                "sector: -5"
              ]
          s₃ =
            parseLines @SymmetrySpec
              [ "permutation: []",
                "sector: 5"
              ]
          s₄ =
            parseLines @SymmetrySpec
              [ "permutation: [1]",
                "sector: null"
              ]
      case s₁ of
        Left _ -> s₁ `shouldSatisfy` isRight
        Right x -> x `shouldBe` SymmetrySpec (fromList [3, 1, 2, 4]) 8
      case s₂ of
        Left _ -> s₂ `shouldSatisfy` isRight
        Right x -> x `shouldBe` SymmetrySpec (fromList [1, 1, 2, 4]) (-5)
      s₃ `shouldSatisfy` isLeft
      s₄ `shouldSatisfy` isLeft
  describe "denseMatrixFromList" $ do
    it "handles empty lists" $ do
      denseMatrixFromList [] `shouldBe` Just (DenseMatrix (0, 0) S.empty :: DenseMatrix Double)
      denseMatrixFromList [[], []] `shouldBe` Just (DenseMatrix (2, 0) S.empty :: DenseMatrix Double)
    it "handles square matrices" $ do
      denseMatrixFromList [[1, 2], [3, 4]]
        `shouldBe` Just (DenseMatrix (2, 2) (fromList [1, 2, 3, 4]) :: DenseMatrix Double)
    it "handles rectangular matrices" $ do
      denseMatrixFromList [[1, 2, 3], [4, 5, 6]]
        `shouldBe` Just (DenseMatrix (2, 3) (fromList [1, 2, 3, 4, 5, 6]) :: DenseMatrix Double)
      denseMatrixFromList [[1, 2], [3, 4], [5, 6]]
        `shouldBe` Just (DenseMatrix (3, 2) (fromList [1, 2, 3, 4, 5, 6]) :: DenseMatrix Double)
    it "general lists" $ do
      (denseMatrixFromList [[1, 2], [3], [5, 6]] :: Maybe (DenseMatrix Int))
        `shouldBe` Nothing
  describe "DenseMatrixSpec" $ do
    it "parses JSON specifications" $ do
      parseLines @DenseMatrixSpec
        [ "[[1, 2,  1],",
          " [4, 2,  1],",
          " [1, -2, 3]]"
        ]
        `checkRight` DenseMatrixSpec [[1, 2, 1], [4, 2, 1], [1, -2, 3]]
      parseLines @DenseMatrixSpec
        [ "[[1, 2],",
          " [4, 2],",
          " [1, [-2, 3]]]"
        ]
        `checkRight` DenseMatrixSpec [[1, 2], [4, 2], [1, (-2) :+ 3]]
      parseLines @DenseMatrixSpec ["[]"] `checkRight` DenseMatrixSpec []
      parseLines @DenseMatrixSpec ["[[], []]"] `checkRight` DenseMatrixSpec [[], []]
  describe "denseToCSRSquare" $ do
    it "works" $ do
      denseToCSRSquare (fromList [[1, 0, 2], [0, 0, -3]]) `shouldBe` Nothing
      print $ denseToCSRSquare (fromList [[1, 0, 2], [0, 4, 0], [0, 0, -3]])
      True `shouldBe` True
    it ".." $ do
      let m =
            fromList
              [ [1, 0, 0.2, 0],
                [0, 1, (-2) :+ 8, 0],
                [0, -2.1, 1, 0],
                [0 :+ 0.1, 0, 0, 3]
              ]
          expr = denseToCSRSquare m
       in case expr of
            Just m' -> sparseToDense m' `shouldBe` m
            Nothing -> expr `shouldSatisfy` isJust
      let m =
            fromList
              [[3]]
          expr = denseToCSRSquare m
       in case expr of
            Just m' -> sparseToDense m' `shouldBe` m
            Nothing -> expr `shouldSatisfy` isJust
      let m =
            fromList
              []
          expr = denseToCSRSquare m
       in case expr of
            Just m' -> sparseToDense m' `shouldBe` m
            Nothing -> expr `shouldSatisfy` isJust
  describe "Storable instances" $ do
    it "has correct sizeOf" $ do
      sizeOf (undefined :: Csparse_matrix) `shouldBe` trueCsparse_matrixSizeOf
      sizeOf (undefined :: Cbit_index) `shouldBe` trueCbit_indexSizeOf
      sizeOf (undefined :: Cterm) `shouldBe` trueCtermSizeOf
      sizeOf (undefined :: Coutput_buffer) `shouldBe` trueCoutput_bufferSizeOf
      sizeOf (undefined :: Csparse_operator) `shouldBe` trueCsparse_operatorSizeOf
    it "has correct alignment" $ do
      alignment (undefined :: Csparse_matrix) `shouldBe` trueCsparse_matrixAlignment
      alignment (undefined :: Cbit_index) `shouldBe` trueCbit_indexAlignment
      alignment (undefined :: Cterm) `shouldBe` trueCtermAlignment
      alignment (undefined :: Coutput_buffer) `shouldBe` trueCoutput_bufferAlignment
      alignment (undefined :: Csparse_operator) `shouldBe` trueCsparse_operatorAlignment
  describe "applyOperatorTerm'" $ do
    it ".." $ do
      let spec =
            InteractionSpec
              [ [1, 0, 0, 0],
                [0, -1, 2, 0],
                [0, 2, -1, 0],
                [0, 0, 0, 1]
              ]
              [[0, 1]]
          (Right term) = toOperatorTerm spec
      print =<< term `applyOperatorTerm'` [0x1]
      True `shouldBe` True

{-
describe "BasisSpec" $ do
  it "parses JSON specifications" $ do
    let s₁ =
          parseLines @BasisSpec $
            [ "number_spins: 4",
              "hamming_weight: 2",
              "symmetries:",
              "  - permutation: [1, 2, 3, 0]",
              "    sector: 0",
              "  - permutation: [3, 2, 1, 0]",
              "    sector: 0"
            ]
        expected₁ = BasisSpec 4 (Just 2) [SymmetrySpec [1, 2, 3, 0] False 0, SymmetrySpec [3, 2, 1, 0] False 0]
    case s₁ of
      Left _ -> s₁ `shouldSatisfy` isRight
      Right x -> x `shouldBe` expected₁

    let s₂ =
          parseLines @BasisSpec $
            [ "number_spins: 100",
              "hamming_weight: null",
              "symmetries: []"
            ]
        expected₂ = BasisSpec 100 Nothing []
    case s₂ of
      Left _ -> s₂ `shouldSatisfy` isRight
      Right x -> x `shouldBe` expected₂
  it "creates SymmetryGroup" $ do
    s₁ <- toSymmetry $ SymmetrySpec [3, 2, 1, 0] False 0
    s₂ <- toSymmetry $ SymmetrySpec [1, 2, 3, 0] False 0
    s₂' <- toSymmetry $ SymmetrySpec [1, 2, 3, 0] False 1
    s₃ <- toSymmetry $ SymmetrySpec [0, 1, 2, 3] True 0
    g <- mkGroup [s₁, s₂, s₃]
    getGroupSize g `shouldBe` 16
    mkGroup [s₁, s₂', s₃] `shouldThrow` anyLatticeSymmetriesException
  it "creates SpinBasis" $ do
    let s₁ = BasisSpec 4 (Just 2) [SymmetrySpec [1, 2, 3, 0] False 0, SymmetrySpec [3, 2, 1, 0] False 0]
    _ <- toBasis s₁
    let s₂ = BasisSpec (-2) Nothing [SymmetrySpec [1, 2, 3, 0] False 0]
    toBasis s₂ `shouldThrow` anySpinEDException

describe "InteractionSpec" $ do
  it "parses JSON specifications" $ do
    let (s₁ :: Either ParseException InteractionSpec) =
          decodeEither' . encodeUtf8 . unlines $
            [ "matrix: [[1, 2,  1],",
              "         [4, 2,  1],",
              "         [1, -2, 3]]",
              "sites: [[0, 0, 0], [1, 1, 1]]"
            ]
    s₁ `shouldSatisfy` isRight
  it "creates Interaction from InteractionSpec" $ do
    let matrix =
          [ [1.0 :+ 0.0, 0.0 :+ 0.0],
            [0.0 :+ 0.0, (-1.0) :+ 0.0]
          ]
        sites = [[0], [2]]
        spec = InteractionSpec matrix sites
    x <- toInteraction spec
    isRealInteraction x `shouldBe` True
-}
