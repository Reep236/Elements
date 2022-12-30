{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoStarIsType #-}
{-# LANGUAGE LambdaCase #-}

module Compounds 
    ( Ionic (..) 
    , Covalent2 
    , mkPolarCov2 
    , mkMaybeIonic
    , ValidIonic
    , Molecule (..)
    , Molecules (..)
    , Molecular (..)
    , KnownMolecule
    , molecValI
    , shatterMolecule
    , shatterSimple 
    ) 
    where 

import Elements 
import Slaters.Unsafe (percentIonic')
import GHC.TypeLits
import Data.Proxy
import Data.Type.Bool (If (..), Not (..), type (&&))
import Data.Type.Equality (type (==))
import Data.Kind (Constraint)

type family ValidIonic_ (e1Ct :: Nat) (e2Ct :: Nat) (bcs :: (Nat, Nat)) :: Constraint where 
    ValidIonic_ c1 c2 '(c1, c2) = () :: Constraint 
    ValidIonic_ _  _   _        = TypeError (Text "No natural ionic form")

type ValidIonic (e1 :: Element) (e1Ct :: Nat) (e2 :: Element) (e2Ct :: Nat) (c1 :: Nat) (c2 :: Nat) = 
    ( Species (Pos c1) ~ IonRep e1 
    , Species (Neg c2) ~ IonRep e2 
    , KnownNat e1Ct 
    , KnownNat e2Ct
    , KnownNat c1 
    , KnownNat c2
    , ValidIonic_ e1Ct e2Ct (BinCoeffs c1 c2)
    , '(e1Ct, e2Ct) ~ BinCoeffs c1 c2
    )

data Ionic (e1Ct :: Nat) (e2Ct :: Nat) =
    forall a b. 
    ( KnownNat a, KnownNat b
    , '(e1Ct, e2Ct) ~ BinCoeffs a b)
    => Ionic (Species (Pos a)) (Species (Neg b)) 

type family BinCoeffs (s1Charge :: Nat) (s2Charge :: Nat) :: (Nat, Nat) where 
    BinCoeffs x x = '(1, 1) 
    BinCoeffs x y = '(y, x)

instance Eq (Ionic a b) where 
    (==) (Ionic a b) (Ionic c d) = (specValue a == specValue c) && (specValue b == specValue d)

instance (KnownNat e1Ct, KnownNat e2Ct) => Show (Ionic e1Ct e2Ct) where 
    show (Ionic a b) =
        let as = case specValue a of
                  Right a       -> show a 
                  Left (a, b) -> "(" ++ show a ++ show b ++ ")" 
            bs = case specValue b of
                  Right a       -> show a 
                  Left  (a, b) -> "(" ++ show a ++ show b ++ ")" 
            an = natVal $ Proxy @e1Ct 
            bn = natVal $ Proxy @e2Ct 
        in as ++ (if an == 1 then "" else show an) ++ bs ++ (if bn == 1 then "" else show bn)

data Covalent2 (e1Ct :: Nat) (e2Ct :: Nat) = UnsafeMkCov2 Element Element deriving Eq

instance (KnownNat a, KnownNat b) => Show (Covalent2 a b) where 
    show (UnsafeMkCov2 e1 e2) = 
        let e1Ct = natVal $ Proxy @a 
            e2Ct = natVal $ Proxy @b 
        in  show e1 ++ (if e1Ct == 1 then "" else show e1Ct) ++ show e2 ++ (if e2Ct == 1 then "" else show e2Ct)

mkPolarCov2 :: 
    forall e1 e1Ct e2 e2Ct.
      ( KnownElem e1, KnownElem e2
      , ValidPolar2 e1 e1Ct e2 e2Ct Neutral) 
    => Covalent2 e1Ct e2Ct 
mkPolarCov2 = UnsafeMkCov2 (elemValI @e1) (elemValI @e2)

-- Position, Bonds to parent, Formal Charge, Element, Children|Count
data Molecule = Central Nat Charge Element [Molecule] | Terminating Nat Charge Element Nat deriving Eq
data Molecules = Molecules Molecule Int deriving Eq

class Molecular a where 
    asMolecule :: a -> Molecule

instance Show Molecule where 
    show (Terminating _ fc e n) = show e ++ (if n == 1 then "" else show n) ++ '|':show fc
    show (Central _ fc e bs) = '(':show e ++ '|':show fc ++ ")" 
                                     ++ concatMap (\x -> '(':show x ++ ")") bs

instance Show Molecules where 
    show (Molecules m n) = '(':show m ++ ')':show n 

instance Molecular Molecule where 
    asMolecule = id 

shatterMolecule :: Molecule -> [Atoms]
shatterMolecule (Terminating _ _ e n) = [Atoms e $ fromIntegral n]
shatterMolecule (Central _ _ e ms) = Atoms e 1 : concatMap shatterMolecule ms

shatterSimple :: Molecule -> [Atoms]
shatterSimple = foldr insert [] . shatterMolecule 
    where

        insert :: Atoms -> [Atoms] -> [Atoms]
        insert a [] = [a]
        insert a@(Atoms e n) (a2@(Atoms e2 n2):as) = if e == e2 then Atoms e (n + n2) : as 
                                                                else a2 : insert a as 


newtype MolecSing (m :: Molecule) = MolecSing Molecule 

class KnownMolecPart m where 
    molecSing :: MolecSing m 

type KnownMolecule (m :: Molecule) = (KnownMolecPart m, FC0 (FormalCharge m))

type family FC0 (s :: Charge) :: Constraint where 
    FC0 Neutral = () :: Constraint 
    FC0 chg = TypeError (Text "Nonzero Formal Charge: " :<>: ShowType chg)

molecValI :: forall m. KnownMolecule m => Molecule 
molecValI = molecPartValI @m

class KnownMolecules ms where
    molecPartValIs :: [Molecule]

instance TypeError (Text "Cannot make empty molecule") => KnownMolecules '[] where 
    molecPartValIs = undefined 

instance KnownMolecPart m => KnownMolecules '[m] where 
    molecPartValIs = [molecPartValI @m]

instance (KnownMolecPart m, KnownMolecules (m_ ': ms)) => KnownMolecules (m ': m_ ': ms) where 
    molecPartValIs = molecPartValI @m : molecPartValIs @(m_ ': ms)

molecPartValI :: forall m. KnownMolecPart m => Molecule 
molecPartValI  = case molecSing @m of 
               MolecSing m -> m 

instance ( KnownNat b
         , KnownSpec fc 
         , KnownNat n
         , KnownElem e 
         , ValidMolecule (Terminating b fc e n)
         ) => KnownMolecPart (Terminating b fc e n) where 
    molecSing = MolecSing $ Terminating (fromIntegral . natVal $ Proxy @b) 
                                        (chargeValI @fc)
                                        (elemValI @e) 
                                        (fromIntegral . natVal $ Proxy @n)

instance ( KnownNat b
         , KnownSpec fc 
         , KnownElem e 
         , KnownMolecules ms
         , ValidMolecule (Central b fc e ms)
         ) => KnownMolecPart (Central b fc e ms) where 
    molecSing = MolecSing $ Central (fromIntegral . natVal $ Proxy @b) (chargeValI @fc) (elemValI @e) (molecPartValIs @ms)

type UnstableState (e :: Element) (n :: Nat) = ShowType e :<>: Text " unstable as " 
                                                          :<>: ShowType e :<>: ShowType n

type family IsDiatomic (e :: Element) :: Bool where 
    IsDiatomic H  = True 
    IsDiatomic N  = True 
    IsDiatomic O  = True 
    IsDiatomic F  = True 
    IsDiatomic Cl = True
    IsDiatomic Br = True 
    IsDiatomic I  = True 
    IsDiatomic At = True 
    IsDiatomic Ts = True
    IsDiatomic e  = False 

isDiatomic :: Element -> Bool 
isDiatomic H  = True 
isDiatomic N  = True 
isDiatomic O  = True 
isDiatomic F  = True 
isDiatomic Cl = True
isDiatomic Br = True 
isDiatomic I  = True 
isDiatomic At = True 
isDiatomic Ts = True
isDiatomic e  = False 

type family BondsToParent (ms :: [Molecule]) :: Nat where 
    BondsToParent '[]                         = 0 
    BondsToParent (Terminating n1 _ _ n2 ': ms) = (n1 * n2) + BondsToParent ms
    BondsToParent (Central n _ _ _ ': ms)        = n + BondsToParent ms

bondsToParent :: forall a. Integral a => [Molecule] -> a 
bondsToParent = foldr (\case 
                         Terminating n1 _ _ n2 -> (+) (fromIntegral $ n1*n2) 
                         Central n _ _ _ -> (+) (fromIntegral n)) 0 

type family AreValidMolecules (y1 :: Bool) (xs :: [Molecule]) :: Bool where 
    AreValidMolecules y '[]       = y 
    AreValidMolecules y (x ': xs) = AreValidMolecules (IsValidMolecule x && y) xs 

type family IsValidMolecule (m :: Molecule) :: Bool where
    IsValidMolecule (Terminating 0 Neutral  S 8) = True  
    IsValidMolecule (Terminating 0 Neutral  S n) = TypeError (UnstableState S n)
    IsValidMolecule (Terminating 0 Neutral  P 4) = True  
    IsValidMolecule (Terminating 0 Neutral  P n) = TypeError (UnstableState P n)
    IsValidMolecule (Terminating 0 Neutral  e 2) = If (IsDiatomic e) True (TypeError (UnstableState e 2))
    IsValidMolecule (Terminating 0 Neutral  e 1) = If (IsDiatomic e) (TypeError (UnstableState e 1)) True
    IsValidMolecule (Terminating n fc e _)       = If (IsOctetExpandable (ValenceE e + NaturalCharge fc - 3) n)
                                                      True 
                                                      (TypeError (InvalidOxidationState e n))
    IsValidMolecule (Central n fc e bs)     
      = AreValidMolecules (If (IsOctetExpandable (ValenceE e + NaturalCharge fc - 3) (n + BondsToParent bs))
                              True 
                              (TypeError (InvalidOxidationState e (n + BondsToParent bs)))
                          ) bs

isValidMolecule :: Molecule -> Bool 
isValidMolecule = \case 
    Terminating 0 Neutral S 8 -> True 
    Terminating 0 Neutral S n -> False 
    Terminating 0 Neutral P 4 -> True 
    Terminating 0 Neutral P n -> False 
    Terminating 0 Neutral e 2 -> isDiatomic e 
    Terminating 0 Neutral e 1 -> not $ isDiatomic e 
    Terminating n fc      e _ -> isOctetExpandable (valenceE e + charge fc) n
    Central n fc e bs -> foldr ((&&) . isValidMolecule) (isOctetExpandable (valenceE e + charge fc) $ n + bondsToParent bs) bs 

type family ValidMolecule (m :: Molecule) :: Constraint where 
    ValidMolecule m = 
        ( If (IsValidMolecule m) 
          (() :: Constraint)
          (TypeError (Text "Invalid Molecule: " :<>: ShowType m)))

type family FormalCharges (y1 :: Charge) (ms :: [Molecule]) where 
    FormalCharges y '[] = y 
    FormalCharges y (x ': xs) = FormalCharges (AddChg y (FormalCharge x)) xs 

type family FormalCharge (m :: Molecule) :: Charge where 
    FormalCharge (Terminating _ c _ n) = MulChg n c
    FormalCharge (Central _ c _ bs)    = FormalCharges c bs

formalCharge :: Molecule -> Charge 
formalCharge = \case  
    Terminating _ fc _ n -> fc * Pos n
    Central _ c _ bs     -> foldr ((+) . formalCharge) c bs
     

mkMaybeIonic :: 
    forall e1 e1Ct e2 e2Ct c1 c2. 
        ( KnownElem e1, KnownElem e2
        , ValidPolar2 e1 e1Ct e2 e2Ct Neutral
        , ValidIonic  e1 e1Ct e2 e2Ct c1 c2)
        => Either (Covalent2 e1Ct e2Ct) (Ionic e1Ct e2Ct)
mkMaybeIonic = let e1 = elemValI @e1 
                   e2 = elemValI @e2 
               in  if percentIonic' e1 e2 < 50 
                     then Left  $ UnsafeMkCov2 e1 e2 
                     else Right $ Ionic (mkIonRep @e1) (mkIonRep @e2)

rawCombineMolecules :: Molecule -> Molecule -> Maybe Molecule 
rawCombineMolecules m1 m2 = case m1 of 
    Terminating b c e n 
      | n == 1    -> let res = Central b c e [m2] 
                      in if isValidMolecule res then Just res else Nothing
      | otherwise -> let res = Central b c e [m2, Terminating b c e $ n-1]
                      in if isValidMolecule res then Just res else Nothing 
    Central b c e ms -> let res = Central b c e 
                         in if isValidMolecule (res $ m2:ms) then Just . res $ m2:ms else Nothing

