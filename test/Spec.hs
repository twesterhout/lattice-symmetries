{-# LANGUAGE BinaryLiterals #-}
{-# LANGUAGE OverloadedLists #-}

module Main (main) where

import Data.Bits
import Data.Complex
import Data.Ratio ((%))
import Data.Type.Equality
import qualified Data.Vector as V
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Unboxed as U
import Data.Yaml.Aeson
import Foreign.Storable
import GHC.Exts (IsList (..))
import LatticeSymmetries
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import qualified LatticeSymmetries.CSR as CSR
import LatticeSymmetries.ComplexRational
import qualified LatticeSymmetries.Dense as Dense
import LatticeSymmetries.Generator
import LatticeSymmetries.IO
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Operator
import LatticeSymmetries.Parser
import LatticeSymmetries.Sparse
import LatticeSymmetries.Types
import Test.Hspec
import Text.Parsec (parse)
import Text.PrettyPrint.ANSI.Leijen (hPutDoc, hardline, pretty, putDoc)
import Prelude hiding (Product, Sum, toList)

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

type SpinPolynomial =
  Sum (Scaled ComplexRational (Product (Generator Int SpinGeneratorType)))

type FermionicPolynomial =
  Sum (Scaled ComplexRational (Product (Generator Int FermionGeneratorType)))

type Dense r c = DenseMatrix V.Vector r c ComplexRational

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
      let (m :: CSR S.Vector Int 3 4 Int) = denseToCsr (fromList [[1, 0, 2, 0], [0, 4, 0, 0], [0, 0, -3, 8]])
      print m
      m `csrIndex` (0, 0) `shouldBe` 1
      m `csrIndex` (0, 1) `shouldBe` 0
      m `csrIndex` (0, 2) `shouldBe` 2
      m `csrIndex` (2, 3) `shouldBe` 8
  describe "Num CSR" $ do
    it "(+)" $ do
      let (a :: CSR S.Vector Int 3 4 Int) =
            denseToCsr (fromList [[1, 0, 2, 0], [0, 4, 0, 0], [0, 0, -3, 8]])
          (b :: CSR S.Vector Int 3 4 Int) =
            denseToCsr (fromList [[3, -8, 0, 0], [1, 5, 0, 0], [0, 0, -3, 0]])
          (c :: CSR S.Vector Int 3 4 Int) =
            denseToCsr (fromList [[4, -8, 2, 0], [1, 9, 0, 0], [0, 0, -6, 8]])
      print a
      print b
      print c
      a + b `shouldBe` c
    it "(+)" $ do
      let a₁' =
            fromList @(SomeCSR S.Vector Int Double) $
              [ (0, 0, 0.05052698254161314),
                (0, 1, 0.34388186562432865),
                (0, 2, 0.3736125370797053),
                (0, 4, 0.37422323470438656),
                (1, 1, 0.08956717735368347),
                (1, 3, 0.20028803087704727),
                (2, 0, 0.010520730797154432),
                (2, 1, 0.7913535778080222),
                (2, 3, 0.6557777513743904),
                (2, 4, 0.9609587020099947)
              ]
          b₁' =
            fromList @(SomeCSR S.Vector Int Double) $
              [ (2, 3, 0.3942046813633675),
                (2, 4, 0.0008746925379253812)
              ]
          (c₁ :: [(Int, Int, Double)]) =
            [ (0, 0, 0.05052698254161314),
              (0, 1, 0.34388186562432865),
              (0, 2, 0.3736125370797053),
              (0, 4, 0.37422323470438656),
              (1, 1, 0.08956717735368347),
              (1, 3, 0.20028803087704727),
              (2, 0, 0.010520730797154432),
              (2, 1, 0.7913535778080222),
              (2, 3, 1.049982432737758),
              (2, 4, 0.96183339454792)
            ]
          a₂ =
            fromList @(CSR S.Vector Int 5 2 Double) $
              [(1, 1, 0.9913900269818051), (2, 0, 0.701927418412059)]
          b₂ =
            fromList @(CSR S.Vector Int 5 2 Double) $
              [(2, 0, 0.5563116167541139), (2, 1, 0.5507459794849707), (3, 1, 0.761606773358582)]
          (c₂ :: [(Int, Int, Double)]) =
            [ (1, 1, 0.9913900269818051),
              (2, 0, 1.2582390351661727),
              (2, 1, 0.5507459794849707),
              (3, 1, 0.761606773358582)
            ]
      withSomeCsr a₁' $ \a₁ -> withSomeCsr b₁' $ \b₁ ->
        case sameShape a₁ b₁ of
          Just Refl -> toList (a₁ + b₁) `shouldBe` c₁
          x -> x `shouldSatisfy` isJust
      toList (a₂ + b₂) `shouldBe` c₂
  -- describe "pPrimitiveOperator" $ do
  --   it ".." $ do
  --     parse pPrimitiveOperator "" ("c†↓₁₀" :: Text)
  --       `shouldBe` Right (SpinfulFermionicOperator FermionicCreationOperator 'c' SpinDown 10)
  --     parse pPrimitiveOperator "" ("f₀" :: Text)
  --       `shouldBe` Right (SpinlessFermionicOperator FermionicAnnihilationOperator 'f' 0)
  --     parse pPrimitiveOperator "" ("n↑₃₈" :: Text)
  --       `shouldBe` Right (SpinfulFermionicOperator FermionicNumberCountingOperator 'n' SpinUp 38)
  --     parse pPrimitiveOperator "" ("Sᶻ₃₈" :: Text)
  --       `shouldBe` Right (SpinOperator SpinZOperator 'S' 38)
  --     parse pOperatorString "" ("n₃₈ f↓₁₅" :: Text)
  --       `shouldBe` Right
  --         ( SpinlessFermionicOperator FermionicNumberCountingOperator 'n' 38
  --             :| [SpinfulFermionicOperator FermionicAnnihilationOperator 'f' SpinDown 15]
  --         )
  --     parse pOperatorString "" ("n↑₃₈f↓₁₅" :: Text)
  --       `shouldBe` Right
  --         ( SpinfulFermionicOperator FermionicNumberCountingOperator 'n' SpinUp 38
  --             :| [SpinfulFermionicOperator FermionicAnnihilationOperator 'f' SpinDown 15]
  --         )
  -- describe "sortByWithParity" $ do
  --   it ".." $ do
  --     sortByWithParity compare [3 :: Int, 1, 4, 8, 0] `shouldBe` (Odd, [0, 1, 3, 4, 8])
  --     sortByWithParity compare [3 :: Int, 1] `shouldBe` (Odd, [1, 3])
  --     sortByWithParity compare [1 :: Int, 3] `shouldBe` (Even, [1, 3])
  --     sortByWithParity compare ([] :: [Int]) `shouldBe` (Even, [])
  describe "csrMatMul" $ do
    it ".." $ do
      let (a :: StorableDenseMatrix 3 4 Float) =
            fromList $
              [ [0.0, 0.0, 0.0, 0.87383316],
                [0.0, 0.0, 0.0, 0.13542224],
                [0.24292645, 0.0, 0.0, 0.0]
              ]
          (b :: StorableDenseMatrix 4 5 Float) =
            fromList $
              [ [0.0, 0.0, 0.0, 0.67190817, 0.0],
                [0.80930562, 0.0, 0.0, 0.94956449, 0.51363394],
                [0.63021572, 0.94316889, 0.0, 0.89156391, 0.0],
                [0.0, 0.0, 0.64573977, 0.05604065, 0.28140917]
              ]
          (c :: StorableDenseMatrix 3 5 Float) =
            fromList $
              [ [0.0, 0.0, 0.56426882, 0.048970178, 0.24590467],
                [0.0, 0.0, 0.08744753, 0.00758915, 0.03810906],
                [0.0, 0.0, 0.0, 0.16322428, 0.0]
              ]
      csrMatMul (denseToCsr @S.Vector @Int a) (denseToCsr b) `shouldBe` (denseToCsr c)
    it "multiplies matrices" $ do
      let (a :: CSR S.Vector Int 5 2 Double) = fromList $ [(0, 0, 0.8213116953652885), (0, 1, 0.2721720739350296), (1, 1, 0.0004725096934642403), (2, 0, 0.03350833939255393), (4, 0, 0.2598645975559749)]
          (b :: CSR S.Vector Int 2 6 Double) = fromList $ [(0, 2, 0.7175422685574337), (0, 3, 0.2603509882468078), (0, 5, 0.20233694825026882), (1, 0, 0.9695204042199307), (1, 1, 0.25913885822625227), (1, 3, 0.3241512404858632)]
          (c :: CSR S.Vector Int 5 6 Double) = fromList $ [(0, 1, 0.07053036048059469), (0, 0, 0.2638763791388668), (0, 5, 0.16618170200246693), (0, 3, 0.3020542269386639), (0, 2, 0.589325857085161), (1, 3, 0.00015316460327802846), (1, 1, 0.00012244562246515968), (1, 0, 0.0004581077890052858), (2, 5, 0.006779975133623629), (2, 3, 0.008723929275360856), (2, 2, 0.02404364986332557), (4, 5, 0.05258020962776023), (4, 3, 0.06765600478405706), (4, 2, 0.1864638328480788)]
      csrMatMul a b `shouldBe` c

  describe "combineNeighbors" $ do
    it ".." $ do
      combineNeighbors (==) (+) (S.fromList [1 :: Int, 2, 3]) `shouldBe` (S.fromList [1, 2, 3])
      combineNeighbors (==) (+) (S.fromList [1 :: Int, 1, 3]) `shouldBe` (S.fromList [2, 3])
      combineNeighbors (==) (const) (S.fromList [1 :: Int, 1, 2, 1, 1, 1, 3, 3])
        `shouldBe` (S.fromList [1, 2, 1, 3])
      combineNeighbors (==) (+) (S.fromList ([] :: [Int])) `shouldBe` (S.fromList [])
  describe "pSpinOperator" $ do
    it "parses simple operators" $ do
      parse pSpinOperator "" ("σᶻ₁₀" :: Text)
        `shouldBe` Right (fromList [Scaled 1 (Generator 10 SpinZ)])
      parse pSpinOperator "" ("S⁺₀" :: Text)
        `shouldBe` Right (fromList [Scaled (fromRational (1 % 2)) (Generator 0 SpinPlus)])
    it "parses composite operators" $ do
      parse pSpinOperator "" ("σˣ₅" :: Text)
        `shouldBe` Right
          ( fromList
              [ Scaled 1 (Generator 5 SpinPlus),
                Scaled 1 (Generator 5 SpinMinus)
              ]
          )
      parse pSpinOperator "" ("σʸ₈₉₄₃" :: Text)
        `shouldBe` Right
          ( fromList
              [ Scaled (ComplexRational 0 (-1)) (Generator 8943 SpinPlus),
                Scaled (ComplexRational 0 1) (Generator 8943 SpinMinus)
              ]
          )
  describe "pOperatorString" $ do
    it ".." $ do
      parse (pOperatorString pSpinOperator) "" ("σᶻ₁₀σ⁺₁₁" :: Text)
        `shouldBe` Right [Scaled 1 [Generator 10 SpinZ, Generator 11 SpinPlus]]
      parse (pOperatorString pSpinOperator) "" ("σᶻ₁₀ σ⁺₁₁" :: Text)
        `shouldBe` Right [Scaled 1 [Generator 10 SpinZ, Generator 11 SpinPlus]]
      parse (pOperatorString pSpinOperator) "" ("σᶻ₁₀ Sˣ₅ σ⁺₁₁" :: Text)
        `shouldBe` Right
          ( [ Scaled
                (fromRational (1 % 2))
                [Generator 10 SpinZ, Generator 5 SpinPlus, Generator 11 SpinPlus],
              Scaled
                (fromRational (1 % 2))
                [Generator 10 SpinZ, Generator 5 SpinMinus, Generator 11 SpinPlus]
            ]
          )
  describe "sumToCanonical" $ do
    it ".." $ do
      -- let (Right [x₀]) =
      --       toList <$> parse (pOperatorString pSpinOperator) "" ("σᶻ₁₀ σ⁺₅ σ⁺₁₁" :: Text)
      -- productToCanonical x₀

      let (Right x₀) = parse (pOperatorString pSpinOperator) "" ("σᶻ₁₀ Sˣ₅ σ⁺₁₁" :: Text)
          y₀ :: SpinPolynomial
          y₀ =
            [ Scaled
                (fromRational (1 % 2))
                [Generator 5 SpinPlus, Generator 10 SpinZ, Generator 11 SpinPlus],
              Scaled
                (fromRational (1 % 2))
                [Generator 5 SpinMinus, Generator 10 SpinZ, Generator 11 SpinPlus]
            ]
          x₁ :: SpinPolynomial
          x₁ =
            [ Scaled
                5
                [Generator 5 SpinMinus, Generator 10 SpinZ, Generator 11 SpinPlus],
              Scaled
                3
                [Generator 11 SpinMinus, Generator 5 SpinPlus, Generator 10 SpinZ],
              Scaled
                (-1)
                [Generator 11 SpinMinus, Generator 5 SpinZ]
            ]
          y₁ :: SpinPolynomial
          y₁ =
            [ Scaled
                (-1)
                [Generator 5 SpinZ, Generator 11 SpinMinus],
              Scaled
                3
                [Generator 5 SpinPlus, Generator 10 SpinZ, Generator 11 SpinMinus],
              Scaled
                5
                [Generator 5 SpinMinus, Generator 10 SpinZ, Generator 11 SpinPlus]
            ]
      simplify x₀ `shouldBe` y₀
      simplify x₁ `shouldBe` y₁

  describe "simplify (spin)" $ do
    let parseString ::
          Text ->
          Polynomial ComplexRational (Generator Int SpinGeneratorType)
        -- Sum (Scaled ComplexRational (Product (Generator Int SpinGeneratorType)))
        parseString s = case parse (pOperatorString pSpinOperator) "" s of
          Left e -> error (show e)
          Right x -> x
    it "computes spin commutators" $ do
      simplify (parseString "σ⁺₁₀ σᶻ₁₀" - parseString "σᶻ₁₀ σ⁺₁₀")
        `shouldBe` [Scaled (-2) [Generator 10 SpinPlus]]
      simplify (parseString "σ⁺₁₀ σ⁻₁₀" - parseString "σ⁻₁₀ σ⁺₁₀")
        `shouldBe` [Scaled 1 [Generator 10 SpinZ]]
      simplify (parseString "σ⁻₁₀ σᶻ₁₀" - parseString "σᶻ₁₀ σ⁻₁₀")
        `shouldBe` [Scaled 2 [Generator 10 SpinMinus]]
    it "simplifies products of primitive operators" $ do
      simplify (parseString "σᶻ₁₀ σᶻ₁₀") `shouldBe` [Scaled 1 [Generator 10 SpinIdentity]]
      simplify (parseString "σ⁺₁₀ σ⁺₁₀") `shouldBe` []
      simplify (parseString "σ⁻₁₀ σ⁻₁₀") `shouldBe` []
      simplify (parseString "σ⁺₁₀ σᶻ₁₀") `shouldBe` [Scaled (-1) [Generator 10 SpinPlus]]
      simplify (parseString "σ⁻₁₀ σᶻ₁₀") `shouldBe` [Scaled 1 [Generator 10 SpinMinus]]
      simplify (parseString "σᶻ₁₀ σ⁺₁₀") `shouldBe` [Scaled 1 [Generator 10 SpinPlus]]
      simplify (parseString "σᶻ₁₀ σ⁻₁₀") `shouldBe` [Scaled (-1) [Generator 10 SpinMinus]]
      simplify (parseString "σ⁺₁₀ σ⁻₁₀")
        `shouldBe` [ Scaled (fromRational (1 % 2)) [Generator 10 SpinIdentity],
                     Scaled (fromRational (1 % 2)) [Generator 10 SpinZ]
                   ]
      simplify (parseString "σ⁻₁₀ σ⁺₁₀")
        `shouldBe` [ Scaled (fromRational (1 % 2)) [Generator 10 SpinIdentity],
                     Scaled (fromRational (-1 % 2)) [Generator 10 SpinZ]
                   ]
    it "simplifies products with identities" $ do
      forM_ ([SpinIdentity, SpinPlus, SpinMinus, SpinZ] :: [SpinGeneratorType]) $ \g -> do
        simplify [Scaled 1 [Generator 5 SpinIdentity, Generator 5 g]]
          `shouldBe` [Scaled (1 :: Rational) [Generator (5 :: Int) g]]
        simplify [Scaled 1 [Generator 5 g, Generator 5 SpinIdentity]]
          `shouldBe` [Scaled (1 :: Rational) [Generator (5 :: Int) g]]

  describe "CSR.csrMatrixFromDense" $ do
    it "converts dense matrices to sparse form" $ do
      CSR.csrMatrixFromDense ([[1, 2], [3, 4]] :: Dense.DenseMatrix S.Vector (Complex Double))
        `shouldBe` (CSR.CsrMatrix [0, 1, 2] [1, 0] [2, 3] [1, 4])
  describe "Num CSR" $ do
    it "(+)" $
      do
        let a = CSR.csrMatrixFromDense ([[1, 2], [0, 2]] :: Dense.DenseMatrix S.Vector (Complex Double))
            b = CSR.csrMatrixFromDense ([[0, 0], [1, 3]] :: Dense.DenseMatrix S.Vector (Complex Double))
            c = CSR.csrMatrixFromDense ([[1, 2], [1, 5]] :: Dense.DenseMatrix S.Vector (Complex Double))
        (a + b) `shouldBe` c

  -- describe "builds operators" $ do
  --   it ".." $ do
  --     print $ mkSpinOperator "σ⁺₁ σ⁻₂" [[0, 1], [1, 2], [2, 0]]
  --     True `shouldBe` True

  describe "NonbranchingTerm" $ do
    it ".." $ do
      let a = nonbranchingRepresentation (Generator (3 :: Int) FermionCreate)
          b = nonbranchingRepresentation (Generator (4 :: Int) FermionAnnihilate)
          c = nonbranchingRepresentation (Generator (3 :: Int) FermionAnnihilate)
          d = nonbranchingRepresentation (Generator (4 :: Int) FermionIdentity)
      print a
      print b
      case a <> a of
        (NonbranchingTerm v _ _ _ _ _) -> v `shouldBe` 0
      case b <> b of
        (NonbranchingTerm v _ _ _ _ _) -> v `shouldBe` 0
      (a <> c) `shouldBe` (nonbranchingRepresentation (Generator (3 :: Int) FermionCount))
      (d <> d) `shouldBe` d
      (a <> d) `shouldBe` a
      (d <> a) `shouldBe` a
      (b <> d) `shouldBe` b
      (d <> b) `shouldBe` b
      print (a <> b)
      print (b <> a)

  describe "binomialCoefficient" $ do
    it "computes binomial(n, k)" $ do
      binomialCoefficient 10 5 `shouldBe` 252
      binomialCoefficient 10 0 `shouldBe` 1
      binomialCoefficient 10 24 `shouldBe` 0
      binomialCoefficient 4 1 `shouldBe` 4
  describe "stateIndex" $ do
    it "computes indices of representatives" $ do
      isStateIndexIdentity (LatticeSymmetries.Basis.SpinBasis 10 Nothing) `shouldBe` True
      isStateIndexIdentity (LatticeSymmetries.Basis.SpinBasis 10 (Just 8)) `shouldBe` False
      stateIndex (LatticeSymmetries.Basis.SpinBasis 10 Nothing) "|0000000101⟩" `shouldBe` Just 5
      stateIndex (LatticeSymmetries.Basis.SpinBasis 10 (Just 2)) "|0000000101⟩" `shouldBe` Just 1

  -- describe "BitString" $ do
  --   it ".." $ do
  --     let x =
  describe "applyOperator" $ do
    it ".." $ do
      let basis1 = LatticeSymmetries.Basis.SpinBasis 10 Nothing
          operator1 = operatorFromString basis1 "σ⁺₀" [[5]]
      -- hPutDoc stdout (pretty (opTerms operator1) <> hardline)
      applyOperator operator1 "|0000000000⟩" `shouldBe` [(1, "|0000100000⟩")]
      applyOperator operator1 "|0000100000⟩" `shouldBe` []
      applyOperator operator1 "|0111011111⟩" `shouldBe` [(1, "|0111111111⟩")]
      let basis2 = LatticeSymmetries.Basis.SpinBasis 2 Nothing
          operator2 =
            operatorFromString basis2 "σˣ₀ σˣ₁" [[0, 1]]
              + operatorFromString basis2 "σʸ₀ σʸ₁" [[0, 1]]
              + operatorFromString basis2 "σᶻ₀ σᶻ₁" [[0, 1]]
      -- hPutDoc stdout (pretty (opTerms operator2) <> hardline)
      -- print $ getNonbranchingTerms operator2
      applyOperator operator2 "|00⟩" `shouldBe` [(1, "|00⟩")]
      applyOperator operator2 "|01⟩" `shouldBe` [(2, "|10⟩"), (-1, "|01⟩")]
      applyOperator operator2 "|10⟩" `shouldBe` [(2, "|01⟩"), (-1, "|10⟩")]
      applyOperator operator2 "|11⟩" `shouldBe` [(1, "|11⟩")]

      let basis3 = LatticeSymmetries.Basis.SpinfulFermionicBasis 2 SpinfulNoOccupation
          operator3 =
            operatorFromString basis3 "c†↑₀ c↑₁" [[0, 1], [3, 2]]
      hPutDoc stdout (pretty (opTerms operator3) <> hardline)

-- let basis3 = LatticeSymmetries.Basis.SpinfulFermionicBasis 2 SpinfulNoOccupation
--     operator3 =
