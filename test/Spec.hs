module Main (main) where

import Data.Complex
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Unboxed as U
import Data.Yaml.Aeson
import Foreign.Storable
import LatticeSymmetries
import LatticeSymmetries.IO
import LatticeSymmetries.Parser
import LatticeSymmetries.Sparse
import LatticeSymmetries.Types
import Test.Hspec
import Text.Parsec (parse)

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
      denseMatrixFromList []
        `shouldBe` Right (DenseMatrix S.empty :: DenseMatrix S.Vector 0 0 Double)
      denseMatrixFromList [[], []]
        `shouldBe` Right (DenseMatrix S.empty :: DenseMatrix S.Vector 2 0 Double)
    it "handles square matrices" $ do
      denseMatrixFromList [[1, 2], [3, 4]]
        `shouldBe` Right (DenseMatrix (fromList [1, 2, 3, 4]) :: DenseMatrix S.Vector 2 2 Double)
    it "handles rectangular matrices" $ do
      denseMatrixFromList [[1, 2, 3], [4, 5, 6]]
        `shouldBe` Right (DenseMatrix (fromList [1, 2, 3, 4, 5, 6]) :: DenseMatrix S.Vector 2 3 Double)
      denseMatrixFromList [[1, 2], [3, 4], [5, 6]]
        `shouldBe` Right (DenseMatrix (fromList [1, 2, 3, 4, 5, 6]) :: DenseMatrix S.Vector 3 2 Double)
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
  -- describe "denseToCSRSquare" $ do
  --   it "works" $ do
  --     denseToSparse (fromList [[1, 0, 2], [0, 0, -3]]) `shouldBe` Nothing
  --     print $ denseToSparse (fromList [[1, 0, 2], [0, 4, 0], [0, 0, -3]])
  --     True `shouldBe` True
  --   it ".." $ do
  --     let m =
  --           fromList
  --             [ [1, 0, 0.2, 0],
  --               [0, 1, (-2) :+ 8, 0],
  --               [0, -2.1, 1, 0],
  --               [0 :+ 0.1, 0, 0, 3]
  --             ]
  --         expr = denseToSparse m
  --      in case expr of
  --           Just m' -> sparseToDense m' `shouldBe` m
  --           Nothing -> expr `shouldSatisfy` isJust
  --     let m =
  --           fromList
  --             [[3]]
  --         expr = denseToSparse m
  --      in case expr of
  --           Just m' -> sparseToDense m' `shouldBe` m
  --           Nothing -> expr `shouldSatisfy` isJust
  --     let m =
  --           fromList
  --             []
  --         expr = denseToSparse m
  --      in case expr of
  --           Just m' -> sparseToDense m' `shouldBe` m
  --           Nothing -> expr `shouldSatisfy` isJust
  -- describe "Storable instances" $ do
  --   it "has correct sizeOf" $ do
  --     sizeOf (undefined :: Csparse_matrix) `shouldBe` trueCsparse_matrixSizeOf
  --     sizeOf (undefined :: Cbit_index) `shouldBe` trueCbit_indexSizeOf
  --     sizeOf (undefined :: Cterm) `shouldBe` trueCtermSizeOf
  --     sizeOf (undefined :: Coutput_buffer) `shouldBe` trueCoutput_bufferSizeOf
  --     sizeOf (undefined :: Csparse_operator) `shouldBe` trueCsparse_operatorSizeOf
  --   it "has correct alignment" $ do
  --     alignment (undefined :: Csparse_matrix) `shouldBe` trueCsparse_matrixAlignment
  --     alignment (undefined :: Cbit_index) `shouldBe` trueCbit_indexAlignment
  --     alignment (undefined :: Cterm) `shouldBe` trueCtermAlignment
  --     alignment (undefined :: Coutput_buffer) `shouldBe` trueCoutput_bufferAlignment
  --     alignment (undefined :: Csparse_operator) `shouldBe` trueCsparse_operatorAlignment
  -- describe "applyOperatorTerm'" $ do
  --   it ".." $ do
  --     let spec =
  --           InteractionSpec
  --             [ [1, 0, 0, 0],
  --               [0, -1, 2, 0],
  --               [0, 2, -1, 0],
  --               [0, 0, 0, 1]
  --             ]
  --             [[0, 1]]
  --         (Right term) = toOperatorTerm spec
  --     print =<< term `applyOperatorTerm'` [0x1]
  --     True `shouldBe` True
  describe "binarySearch" $ do
    it ".." $ do
      binarySearch (S.fromList [1, 3, 4, 5, 7, 9]) (4 :: Int) `shouldBe` (Just 2)
      binarySearch (S.fromList [1, 3, 4, 5, 7, 9]) (2 :: Int) `shouldBe` Nothing
      binarySearch (S.fromList []) (2 :: Int) `shouldBe` Nothing
      binarySearch (S.fromList [1]) (1 :: Int) `shouldBe` Just 0
  describe "cooToCSR" $ do
    it ".." $ do
      let (m :: CSR S.Vector Int 3 4 Int) = cooToCsr @U.Vector @S.Vector . fromList $ [(2, 2, -3), (1, 1, 1), (2, 3, 8), (1, 1 :: Int, 3 :: Int)]
      print m
      m `csrIndex` (1, 1) `shouldBe` 4
  describe "denseToCSR" $ do
    it ".." $ do
      let (m :: CSR S.Vector Int 3 4 Int) = denseToCsr @S.Vector @S.Vector (fromList [[1, 0, 2, 0], [0, 4, 0, 0], [0, 0, -3, 8]])
      print m
      m `csrIndex` (0, 0) `shouldBe` 1
      m `csrIndex` (0, 1) `shouldBe` 0
      m `csrIndex` (0, 2) `shouldBe` 2
      m `csrIndex` (2, 3) `shouldBe` 8
  describe "Num CSR" $ do
    it "(+)" $ do
      let (a :: CSR S.Vector Int 3 4 Int) =
            denseToCsr @S.Vector @S.Vector (fromList [[1, 0, 2, 0], [0, 4, 0, 0], [0, 0, -3, 8]])
          (b :: CSR S.Vector Int 3 4 Int) =
            denseToCsr @S.Vector @S.Vector (fromList [[3, -8, 0, 0], [1, 5, 0, 0], [0, 0, -3, 0]])
          (c :: CSR S.Vector Int 3 4 Int) =
            denseToCsr @S.Vector @S.Vector (fromList [[4, -8, 2, 0], [1, 9, 0, 0], [0, 0, -6, 8]])
      print a
      print b
      print c
      a + b `shouldBe` c
  describe "pPrimitiveOperator" $ do
    it ".." $ do
      parse pPrimitiveOperator "" ("c†↓₁₀" :: Text)
        `shouldBe` Right (SpinfulFermionicOperator FermionicCreationOperator 'c' SpinDown 10)
      parse pPrimitiveOperator "" ("f₀" :: Text)
        `shouldBe` Right (SpinlessFermionicOperator FermionicAnnihilationOperator 'f' 0)
      parse pPrimitiveOperator "" ("n↑₃₈" :: Text)
        `shouldBe` Right (SpinfulFermionicOperator FermionicNumberCountingOperator 'n' SpinUp 38)
      parse pPrimitiveOperator "" ("Sᶻ₃₈" :: Text)
        `shouldBe` Right (SpinOperator SpinZOperator 'S' 38)
      parse pOperatorString "" ("n₃₈ f↓₁₅" :: Text)
        `shouldBe` Right
          ( SpinlessFermionicOperator FermionicNumberCountingOperator 'n' 38
              :| [SpinfulFermionicOperator FermionicAnnihilationOperator 'f' SpinDown 15]
          )
      parse pOperatorString "" ("n↑₃₈f↓₁₅" :: Text)
        `shouldBe` Right
          ( SpinfulFermionicOperator FermionicNumberCountingOperator 'n' SpinUp 38
              :| [SpinfulFermionicOperator FermionicAnnihilationOperator 'f' SpinDown 15]
          )
  describe "sortByWithParity" $ do
    it ".." $ do
      sortByWithParity compare [3 :: Int, 1, 4, 8, 0] `shouldBe` (Odd, [0, 1, 3, 4, 8])
      sortByWithParity compare [3 :: Int, 1] `shouldBe` (Odd, [1, 3])
      sortByWithParity compare [1 :: Int, 3] `shouldBe` (Even, [1, 3])
      sortByWithParity compare ([] :: [Int]) `shouldBe` (Even, [])
  -- describe "csrMatMul" $ do
  --   it ".." $
  --     do
  --       let a =
  --             denseToCsr @Float . fromList $
  --               [ [0.0, 0.0, 0.0, 0.87383316],
  --                 [0.0, 0.0, 0.0, 0.13542224],
  --                 [0.24292645, 0.0, 0.0, 0.0]
  --               ]
  --           b =
  --             denseToCsr @Float . fromList $
  --               [ [0.0, 0.0, 0.0, 0.67190817, 0.0],
  --                 [0.80930562, 0.0, 0.0, 0.94956449, 0.51363394],
  --                 [0.63021572, 0.94316889, 0.0, 0.89156391, 0.0],
  --                 [0.0, 0.0, 0.64573977, 0.05604065, 0.28140917]
  --               ]
  --           c =
  --             denseToCsr @Float . fromList $
  --               [ [0.0, 0.0, 0.56426882, 0.048970178, 0.24590467],
  --                 [0.0, 0.0, 0.08744753, 0.00758915, 0.03810906],
  --                 [0.0, 0.0, 0.0, 0.16322428, 0.0]
  --               ]
  --       csrMatMul a b `shouldBe` c
  describe "combineNeighbors" $ do
    it ".." $ do
      combineNeighbors (==) (+) (S.fromList [1 :: Int, 2, 3]) `shouldBe` (S.fromList [1, 2, 3])
      combineNeighbors (==) (+) (S.fromList [1 :: Int, 1, 3]) `shouldBe` (S.fromList [2, 3])
      combineNeighbors (==) (const) (S.fromList [1 :: Int, 1, 2, 1, 1, 1, 3, 3])
        `shouldBe` (S.fromList [1, 2, 1, 3])
      combineNeighbors (==) (+) (S.fromList ([] :: [Int])) `shouldBe` (S.fromList [])

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
