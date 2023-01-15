{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE PolyKinds #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE NoStarIsType #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE TypeFamilyDependencies #-}

{-|
Module: Elements
Description: Basic type- and value- level Elements and Ions 
-}

module Elements 
    ( Element (..)
    , KnownElem 
    , KnownCharge 
    , ValenceE (..)
    , valenceE
    , elemVal
    , elemValI 
    , toAtomic
    , fromAtomic
    , Charge (..)
    , charge
    , AddChg (..)
    , MulChg (..)
    , NaturalCharge (..)
    , Species 
    , Atoms (..)
    , mkAtom
    , asAtom
    , IonRep (..)
    , mkIonRep
    , mkIonPoly
    , specValue
    , specCharge
    , chargeValI 
    , InvalidOxidationState (..)
    , OctetExpandable (..)
    , IsOctetExpandable (..)
    , isOctetExpandable 
    , ValidPolar2 (..)
    , ToAtomic (..)
    )
    where

import Data.Proxy
import GHC.TypeLits
import Data.Type.Bool 
import Data.Type.Equality
import Data.Kind (Constraint)

-- | All elements Hydrogen through Oganesson
data Element 
        = H                                                                                  | He 
        | Li | Be                                                   | B  | C  | N  | O  | F  | Ne
        | Na | Mg                                                   | Al | Si | P  | S  | Cl | Ar
        | K  | Ca | Sc | Ti | V  | Cr | Mn | Fe | Co | Ni | Cu | Zn | Ga | Ge | As | Se | Br | Kr
        | Rb | Sr | Y  | Zr | Nb | Mo | Tc | Ru | Rh | Pd | Ag | Cd | In | Sn | Sb | Te | I  | Xe 
        | Cs | Ba 
                  | La | Ce | Pr | Nd | Pm | Sm | Eu | Gd | Tb | Dy | Ho | Er | Tm | Yb | Lu 
                       | Hf | Ta | W  | Re | Os | Ir | Pt | Au | Hg | Tl | Pb | Bi | Po | At | Rn 
        | Fr | Ra 
                  | Ac | Th | Pa | U  | Np | Pu | Am | Cm | Bk | Cf | Es | Fm | Md | No | Lr 
                       | Rf | Db | Sg | Bh | Hs | Mt | Ds | Rg | Cn | Nh | Fl | Mc | Lv | Ts | Og 

        deriving (Show, Read, Eq, Ord, Enum)

-- | Internal `Element` singleton for reflection
newtype SElem (a :: Element) = SElem Element

-- | Type class representing a reifiable `Element`
class KnownElem (a :: Element) where 
    -- | Internal, yields `Element` singleton 
    elemSing :: SElem a

-- | Reify `Element` with `Proxy`
elemVal :: forall a. KnownElem a => Proxy a -> Element 
elemVal _ = case (elemSing :: SElem a) of 
              SElem res -> res 

-- | Reify `Element` independent of `Proxy` 
elemValI :: forall a. KnownElem a => Element 
elemValI = elemVal $ Proxy @a

-- | Convert `Element` to 1-indexed atomic number 
toAtomic :: Element -> Int
toAtomic = (+) 1 . fromEnum

-- | Convert atomic number to `Element` 
fromAtomic :: Int -> Element 
fromAtomic = toEnum . subtract 1

-- | Type level valence electron count (`Nat`) to common ion (`Charge`) 
type family VeToIon (a :: Nat) :: Charge where 
    VeToIon 1 = Pos 1 
    VeToIon 2 = Pos 2 
    VeToIon 3 = Pos 3 
    VeToIon 4 = TypeError (Text "No Ion for 4 valence electrons")
    VeToIon 5 = Neg 3 
    VeToIon 6 = Neg 2 
    VeToIon 7 = Neg 1 
    VeToIon 8 = TypeError (Text "No Ion for 8 valence electrons") 

{-|
   Constraint representing the a valid binary compound consisting of polar bonds between two representative `Element`s given their likely oxidation states and the final `Charge` of the compound.
       Oxygen, Nitrogen, and Chlorine are given special consideration. For all others, free electrons and both electrons of a split pair may bond   
-}
type family ValidPolar2 (e1 :: Element) (ct1 :: Nat) (e2 :: Element) (ct2 :: Nat) (charge :: Charge) :: Constraint where  
    ValidPolar2 e1 ct1 O  ct2 Neutral       = 
        If (  IsOctetExpandable (ValenceE e1) (Div (ct2 * 2) ct1)   -- -2 
           || IsOctetExpandable (ValenceE e1) (Div  ct2      ct1)   -- -1 
           || IsOctetExpandable (ValenceE e1) (Div  ct2 (2 * ct1))) -- -1/2
                (() :: Constraint)
                (TypeError (    Text "No Valid Ox State for " :<>: ShowType e1  
                           :<>: Text " in " :<>: ShowType e1 :<>: ShowType ct1 :<>: Text "O" :<>: ShowType ct2))
    ValidPolar2 e1 ct1 e2 ct2 Neutral       = OctetExpandable e1 (Div (ct2 * (8 - ValenceE e2)) ct1)
    ValidPolar2 e1 ct1 e2 ct2 (Pos charge) = OctetExpandable e1 (charge + Div (ct2 * (8 - ValenceE e2)) ct1) 
    ValidPolar2 e1 ct1 e2 ct2 (Neg charge) = OctetExpandable e1 (Div (ct2 * (8 - ValenceE e2)) ct1 - charge) 

-- | Useful type error for when an invalid oxidation state is presented
type InvalidOxidationState (e :: Element) (target :: Nat) 
  = Text "Invalid Ox State for " :<>: ShowType e 
    :<>: Text ": " :<>: ShowType target 

-- | Constraint representing an `Element` with a valid number of bonded electrons 
type family OctetExpandable (e :: Element) (target :: Nat) :: Constraint where 
    OctetExpandable Cl target = target <= 7 
    OctetExpandable N  target = target <= 5 
    OctetExpandable e  target = If (IsOctetExpandable (ValenceE e) target) 
                                    (() :: Constraint) 
                                    (TypeError  (InvalidOxidationState e target)) 

-- | Boolean form of `OctetExpandable` accepting a raw number of valence electrons. Assumes representative and standard case 
type IsOctetExpandable (valence :: Nat) (target :: Nat) = 
        (  valence == target 
        || valence == target + 2 
        || valence == target + 4 
        || valence == target + 6
      )

-- | Value-level `IsOctetExpandable`  
isOctetExpandable :: Integral a => a -> a -> Bool 
isOctetExpandable v t = v == t || v == t + 2 || v == t + 4 || v == t + 6

-- | Type representing a charge using `Natural`s. Useful for constraints at type level  
data Charge = Neutral | Pos Nat | Neg Nat deriving Eq 

instance Num Charge where 
    (+) Neutral = id 
    (+) (Pos n) = \case 
        Neutral -> Pos n 
        Pos n2 -> Pos (n + n2)
        Neg n2 
            | n == n2   -> Neutral 
            | n >  n2   -> Pos $ n - n2 
            | otherwise -> Neg $ n2 - n 
    (+) (Neg n) = flip (+) (Neg n)

    (*) Neutral = const Neutral
    (*) (Pos n) = \case 
        Neutral -> Neutral 
        Pos n2  -> Pos (n * n2) 
        Neg n2  -> Neg (n * n2)
    (*) (Neg n) = flip (*) (Neg n)

    abs = \case 
        Neg n -> Pos n 
        n     -> n

    signum = \case 
        Neutral -> Neutral 
        Pos n   -> Pos 1 
        Neg n   -> Neg 1 

    fromInteger n 
        | n == 0    = Neutral 
        | n >  0    = Pos (fromIntegral n)
        | otherwise = Neg (fromIntegral . abs $ n)

    negate = \case 
        Neutral -> Neutral 
        Pos n   -> Neg n 
        Neg n   -> Pos n

-- | Converts a `Charge` value to a signed `Integral` value 
charge :: Integral a => Charge -> a 
charge = \case 
    Neutral -> 0 
    Pos a   -> fromIntegral a 
    Neg a   -> negate . fromIntegral $ a 

-- | Type family summing two `Charge`s, identical to `(+)` in `Num` instance 
type family AddChg (c1 :: Charge) (c2 :: Charge) :: Charge where 
    AddChg Neutral c2 = c2 
    AddChg c1 Neutral = c1 
    AddChg (Neg a) (Neg b) = Neg (a + b)
    AddChg (Neg a) (Pos b) = If (a == b) 
                                  Neutral
                                  (If (a <=? b)
                                      (Pos (b - a))
                                      (Neg (a - b))
                                  )
    AddChg (Pos a) (Neg b) = AddChg (Neg a) (Pos b)

-- | Type family yielding some `Charge` multipled by a `Natural`. 
-- Both less versatile and less verbose than straight `Charge`-`Charge` multiplication
type family MulChg (n :: Nat) (c :: Charge) :: Charge where 
    MulChg _ Neutral = Neutral 
    MulChg n (Pos c) = Pos (n * c)
    MulChg n (Neg c) = Neg (n * c)

instance Show Charge where 
    show Neutral  = "" 
    show (Pos n) = show n ++ "+"
    show (Neg n) = show n ++ "-" 

instance Ord Charge where 
    (<=) Neutral = \case 
        Neg _  -> False 
        Pos _  -> True 
        Neutral -> True

    (<=) (Neg c) = \case 
        Neg c2 -> c <= c2 
        Pos _  -> True 
        Neutral -> True 

    (<=) (Pos c) = \case 
        Neg _  -> False 
        Pos c2 -> c <= c2 
        Neutral -> False 

-- | Derivation of reifiable `Nat` to reifiable `Charge` 
type KnownCharge (st :: Charge) = KnownNat (NaturalCharge st)

-- | Reifies a `Charge` independent of `Proxy` 
chargeValI :: forall st n. (KnownNat n, n ~ NaturalCharge st) => Charge 
chargeValI = let n = fromIntegral . natVal $ Proxy @n 
                 st
                  | n == 3    = Neutral 
                  | n <  3    = Neg (3 - n) 
                  | otherwise = Pos (n - 3)
              in st 

-- | Converts a `Charge` to an unwrapped `Natural` of value Charge + 3
type family NaturalCharge (s :: Charge) = (n :: Nat) | n -> s where 
    NaturalCharge (Neg 3) = 0
    NaturalCharge (Neg 2) = 1
    NaturalCharge (Neg 1) = 2
    NaturalCharge Neutral  = 3 
    NaturalCharge (Pos 1) = 4 
    NaturalCharge (Pos 2) = 5 
    NaturalCharge (Pos 3) = 6 

-- | Type representing a number of atoms of a given `Element` 
data Atoms = Atoms Element Int deriving Eq 
instance Show Atoms where 
    show (Atoms e ct) = show e ++ (if ct == 1 then "" else show ct)

-- | Type representing a species (i.e. in solution)
-- Can be monatomic or binary polyatomic
-- Polyatomic Ions are represented as two sets of `Atoms` 
data Species (state :: Charge) = UnsafeMkSpecies Element | UnsafeMkPoly Atoms Atoms deriving Eq 

instance Show (Species Neutral) where 
    show (UnsafeMkSpecies e)  = show e
    show (UnsafeMkPoly ewc1 ewc2) = "(" ++ show ewc1 ++ show ewc2 ++ ")"

instance forall a. KnownNat a => Show (Species (Pos a)) where 
    show (UnsafeMkSpecies e) = show e ++ show (natVal $ Proxy @a) ++ "+"
    show (UnsafeMkPoly ewc1 ewc2) = "(" ++ show ewc1 ++ show ewc2 ++ ")" ++ show (natVal $ Proxy @a) ++ "+"

instance forall a. KnownNat a => Show (Species (Neg a)) where 
    show (UnsafeMkSpecies e) = show e ++ show (natVal $ Proxy @a) ++ "-"
    show (UnsafeMkPoly ewc1 ewc2) = "(" ++ show ewc1 ++ show ewc2 ++ ")" ++ show (natVal $ Proxy @a) ++ "-"


-- | Reify an `Element` as a `Neutral` `Species` 
mkAtom :: forall a. KnownElem a => Species Neutral
mkAtom = UnsafeMkSpecies $ elemValI @a 

-- | Convert an `Element` value to `Neutral` `Species` 
asAtom :: Element -> Species Neutral 
asAtom = UnsafeMkSpecies 

-- | Useful type synonym for the (likely) common ion of a representative `Element` 
type IonRep (e :: Element) = Species (VeToIon (ValenceE e))

-- | Reify a representative `Element` as its common ion 
mkIonRep :: forall a. KnownElem a => IonRep a 
mkIonRep = UnsafeMkSpecies $ elemValI @a

-- | Reify a binary polyatomic ion
-- For example, Nitrate: 
--
-- > mkIonPoly @N @1 @O @3 @(Neg 1)
mkIonPoly :: forall e1 ct1 e2 ct2 charge. 
    ( KnownElem e1, KnownNat ct1  
    , KnownElem e2, KnownNat ct2
    , ValidPolar2 e1 ct1 e2 ct2 charge
    ) => Species charge 
mkIonPoly = UnsafeMkPoly (Atoms (elemValI @e1) (fromIntegral . natVal $ Proxy @ct1)) 
                         (Atoms (elemValI @e2) (fromIntegral . natVal $ Proxy @ct2)) 

-- | Returns the value wrapped by a `Species` 
specValue :: Species a -> Either (Atoms, Atoms) Element
specValue = \case 
    UnsafeMkSpecies a -> Right a 
    UnsafeMkPoly a b  -> Left (a, b)

-- | Returns the charge of a value-level `Species` 
specCharge :: forall a c n. (KnownNat c, c ~ NaturalCharge a, Integral n) => Species a -> n
specCharge _ = charge $ chargeValI @a 

-- | Type-level `toAtomic`, injective  
type family ToAtomic (e :: Element) = n | n -> e where 
    ToAtomic   H = 1
    ToAtomic  He = 2
    ToAtomic  Li = 3
    ToAtomic  Be = 4
    ToAtomic   B = 5
    ToAtomic   C = 6
    ToAtomic   N = 7
    ToAtomic   O = 8
    ToAtomic   F = 9
    ToAtomic  Ne = 10
    ToAtomic  Na = 11
    ToAtomic  Mg = 12
    ToAtomic  Al = 13
    ToAtomic  Si = 14
    ToAtomic   P = 15
    ToAtomic  S  = 16
    ToAtomic  Cl = 17
    ToAtomic  Ar = 18
    ToAtomic  K  = 19
    ToAtomic  Ca = 20
    ToAtomic  Sc = 21
    ToAtomic  Ti = 22
    ToAtomic  V  = 23
    ToAtomic  Cr = 24
    ToAtomic  Mn = 25
    ToAtomic  Fe = 26
    ToAtomic  Co = 27
    ToAtomic  Ni = 28
    ToAtomic  Cu = 29
    ToAtomic  Zn = 30
    ToAtomic  Ga = 31
    ToAtomic  Ge = 32
    ToAtomic  As = 33
    ToAtomic  Se = 34
    ToAtomic  Br = 35
    ToAtomic  Kr = 36
    ToAtomic  Rb = 37
    ToAtomic  Sr = 38
    ToAtomic  Y  = 39
    ToAtomic  Zr = 40
    ToAtomic  Nb = 41
    ToAtomic  Mo = 42
    ToAtomic  Tc = 43
    ToAtomic  Ru = 44
    ToAtomic  Rh = 45
    ToAtomic  Pd = 46
    ToAtomic  Ag = 47
    ToAtomic  Cd = 48
    ToAtomic  In = 49
    ToAtomic  Sn = 50
    ToAtomic  Sb = 51
    ToAtomic  Te = 52
    ToAtomic  I  = 53
    ToAtomic  Xe = 54
    ToAtomic  Cs = 55
    ToAtomic  Ba = 56
    ToAtomic  La = 57
    ToAtomic  Ce = 58
    ToAtomic  Pr = 59
    ToAtomic  Nd = 60
    ToAtomic  Pm = 61
    ToAtomic  Sm = 62
    ToAtomic  Eu = 63
    ToAtomic  Gd = 64
    ToAtomic  Tb = 65
    ToAtomic  Dy = 66
    ToAtomic  Ho = 67
    ToAtomic  Er = 68
    ToAtomic  Tm = 69
    ToAtomic  Yb = 70
    ToAtomic  Lu = 71
    ToAtomic  Hf = 72
    ToAtomic  Ta = 73
    ToAtomic  W  = 74
    ToAtomic  Re = 75
    ToAtomic  Os = 76
    ToAtomic  Ir = 77
    ToAtomic  Pt = 78
    ToAtomic  Au = 79
    ToAtomic  Hg = 80
    ToAtomic  Tl = 81
    ToAtomic  Pb = 82
    ToAtomic  Bi = 83
    ToAtomic  Po = 84
    ToAtomic  At = 85
    ToAtomic  Rn = 86
    ToAtomic  Fr = 87
    ToAtomic  Ra = 88
    ToAtomic  Ac = 89
    ToAtomic  Th = 90
    ToAtomic  Pa = 91
    ToAtomic  U  = 92
    ToAtomic  Np = 93
    ToAtomic  Pu = 94
    ToAtomic  Am = 95
    ToAtomic  Cm = 96
    ToAtomic  Bk = 97
    ToAtomic  Cf = 98
    ToAtomic  Es = 99
    ToAtomic  Fm = 100
    ToAtomic  Md = 101
    ToAtomic  No = 102
    ToAtomic  Lr = 103
    ToAtomic  Rf = 104
    ToAtomic  Db = 105
    ToAtomic  Sg = 106
    ToAtomic  Bh = 107
    ToAtomic  Hs = 108
    ToAtomic  Mt = 109
    ToAtomic  Ds = 110
    ToAtomic  Rg = 111
    ToAtomic  Cn = 112
    ToAtomic  Nh = 113
    ToAtomic  Fl = 114
    ToAtomic  Mc = 115
    ToAtomic  Lv = 116
    ToAtomic  Ts = 117
    ToAtomic  Og = 118

instance KnownElem Og where 
    elemSing = SElem Og
instance KnownElem Ts where 
    elemSing = SElem Ts
instance KnownElem Lv where 
    elemSing = SElem Lv
instance KnownElem Mc where 
    elemSing = SElem Mc
instance KnownElem Fl where 
    elemSing = SElem Fl
instance KnownElem Nh where 
    elemSing = SElem Nh
instance KnownElem Cn where 
    elemSing = SElem Cn
instance KnownElem Rg where 
    elemSing = SElem Rg
instance KnownElem Ds where 
    elemSing = SElem Ds
instance KnownElem Mt where 
    elemSing = SElem Mt
instance KnownElem Hs where 
    elemSing = SElem Hs
instance KnownElem Bh where 
    elemSing = SElem Bh
instance KnownElem Sg where 
    elemSing = SElem Sg
instance KnownElem Db where 
    elemSing = SElem Db
instance KnownElem Rf where 
    elemSing = SElem Rf
instance KnownElem Lr where 
    elemSing = SElem Lr
instance KnownElem No where 
    elemSing = SElem No
instance KnownElem Md where 
    elemSing = SElem Md
instance KnownElem Fm where 
    elemSing = SElem Fm
instance KnownElem Es where 
    elemSing = SElem Es
instance KnownElem Cf where 
    elemSing = SElem Cf
instance KnownElem Bk where 
    elemSing = SElem Bk
instance KnownElem Cm where 
    elemSing = SElem Cm
instance KnownElem Am where 
    elemSing = SElem Am
instance KnownElem Pu where 
    elemSing = SElem Pu
instance KnownElem Np where 
    elemSing = SElem Np
instance KnownElem U  where 
    elemSing = SElem U 
instance KnownElem Pa where 
    elemSing = SElem Pa
instance KnownElem Th where 
    elemSing = SElem Th
instance KnownElem Ac where 
    elemSing = SElem Ac
instance KnownElem Ra where 
    elemSing = SElem Ra
instance KnownElem Fr where 
    elemSing = SElem Fr
instance KnownElem Rn where 
    elemSing = SElem Rn
instance KnownElem At where 
    elemSing = SElem At
instance KnownElem Po where 
    elemSing = SElem Po
instance KnownElem Bi where 
    elemSing = SElem Bi
instance KnownElem Pb where 
    elemSing = SElem Pb
instance KnownElem Tl where 
    elemSing = SElem Tl
instance KnownElem Hg where 
    elemSing = SElem Hg
instance KnownElem Au where 
    elemSing = SElem Au
instance KnownElem Pt where 
    elemSing = SElem Pt
instance KnownElem Ir where 
    elemSing = SElem Ir
instance KnownElem Os where 
    elemSing = SElem Os
instance KnownElem Re where 
    elemSing = SElem Re
instance KnownElem W  where 
    elemSing = SElem W 
instance KnownElem Ta where 
    elemSing = SElem Ta
instance KnownElem Hf where 
    elemSing = SElem Hf
instance KnownElem Lu where 
    elemSing = SElem Lu
instance KnownElem Yb where 
    elemSing = SElem Yb
instance KnownElem Tm where 
    elemSing = SElem Tm
instance KnownElem Er where 
    elemSing = SElem Er
instance KnownElem Ho where 
    elemSing = SElem Ho
instance KnownElem Dy where 
    elemSing = SElem Dy
instance KnownElem Tb where 
    elemSing = SElem Tb
instance KnownElem Gd where 
    elemSing = SElem Gd
instance KnownElem Eu where 
    elemSing = SElem Eu
instance KnownElem Sm where 
    elemSing = SElem Sm
instance KnownElem Pm where 
    elemSing = SElem Pm
instance KnownElem Nd where 
    elemSing = SElem Nd
instance KnownElem Pr where 
    elemSing = SElem Pr
instance KnownElem Ce where 
    elemSing = SElem Ce
instance KnownElem La where 
    elemSing = SElem La
instance KnownElem Ba where 
    elemSing = SElem Ba
instance KnownElem Cs where 
    elemSing = SElem Cs
instance KnownElem Xe where 
    elemSing = SElem Xe
instance KnownElem I  where 
    elemSing = SElem I 
instance KnownElem Te where 
    elemSing = SElem Te
instance KnownElem Sb where 
    elemSing = SElem Sb
instance KnownElem Sn where 
    elemSing = SElem Sn
instance KnownElem In where 
    elemSing = SElem In
instance KnownElem Cd where 
    elemSing = SElem Cd
instance KnownElem Ag where 
    elemSing = SElem Ag
instance KnownElem Pd where 
    elemSing = SElem Pd
instance KnownElem Rh where 
    elemSing = SElem Rh
instance KnownElem Ru where 
    elemSing = SElem Ru
instance KnownElem Tc where 
    elemSing = SElem Tc
instance KnownElem Mo where 
    elemSing = SElem Mo
instance KnownElem Nb where 
    elemSing = SElem Nb
instance KnownElem Zr where 
    elemSing = SElem Zr
instance KnownElem Y  where 
    elemSing = SElem Y 
instance KnownElem Sr where 
    elemSing = SElem Sr
instance KnownElem Rb where 
    elemSing = SElem Rb
instance KnownElem Kr where 
    elemSing = SElem Kr
instance KnownElem Br where 
    elemSing = SElem Br
instance KnownElem Se where 
    elemSing = SElem Se
instance KnownElem As where 
    elemSing = SElem As
instance KnownElem Ge where 
    elemSing = SElem Ge
instance KnownElem Ga where 
    elemSing = SElem Ga
instance KnownElem Zn where 
    elemSing = SElem Zn
instance KnownElem Cu where 
    elemSing = SElem Cu
instance KnownElem Ni where 
    elemSing = SElem Ni
instance KnownElem Co where 
    elemSing = SElem Co
instance KnownElem Fe where 
    elemSing = SElem Fe
instance KnownElem Mn where 
    elemSing = SElem Mn
instance KnownElem Cr where 
    elemSing = SElem Cr
instance KnownElem V  where 
    elemSing = SElem V 
instance KnownElem Ti where 
    elemSing = SElem Ti
instance KnownElem Sc where 
    elemSing = SElem Sc
instance KnownElem Ca where 
    elemSing = SElem Ca
instance KnownElem K  where 
    elemSing = SElem K 
instance KnownElem Ar where 
    elemSing = SElem Ar
instance KnownElem Cl where 
    elemSing = SElem Cl
instance KnownElem S  where 
    elemSing = SElem S 
instance KnownElem P  where 
    elemSing = SElem P 
instance KnownElem Si where 
    elemSing = SElem Si
instance KnownElem Al where 
    elemSing = SElem Al
instance KnownElem Mg where 
    elemSing = SElem Mg
instance KnownElem Na where 
    elemSing = SElem Na
instance KnownElem Ne where 
    elemSing = SElem Ne
instance KnownElem F  where 
    elemSing = SElem F 
instance KnownElem O  where 
    elemSing = SElem O 
instance KnownElem N  where 
    elemSing = SElem N 
instance KnownElem C  where 
    elemSing = SElem C 
instance KnownElem B  where 
    elemSing = SElem B 
instance KnownElem Be where 
    elemSing = SElem Be
instance KnownElem Li where 
    elemSing = SElem Li
instance KnownElem He where 
    elemSing = SElem He
instance KnownElem H  where 
    elemSing = SElem H 


-- | Type family representing the valence electrons (s+p electrons) for all representative elements
-- All transition elements are assigned 2 for their s sublevel 
type family ValenceE (e :: Element) :: Nat where 
    ValenceE H  = 1 
    ValenceE Li = 1 
    ValenceE Na = 1 
    ValenceE K  = 1 
    ValenceE Rb = 1 
    ValenceE Cs = 1 
    ValenceE Fr = 1 
    ValenceE B  = 3 
    ValenceE Al = 3 
    ValenceE In = 3 
    ValenceE Tl = 3 
    ValenceE Nh = 3 
    ValenceE C  = 4 
    ValenceE Si = 4 
    ValenceE Ge = 4 
    ValenceE Sn = 4 
    ValenceE Pb = 4 
    ValenceE Fl = 4 
    ValenceE N  = 5
    ValenceE P  = 5 
    ValenceE As = 5 
    ValenceE Sb = 5 
    ValenceE Bi = 5 
    ValenceE Mc = 5 
    ValenceE O  = 6 
    ValenceE S  = 6 
    ValenceE Se = 6 
    ValenceE Te = 6 
    ValenceE Po = 6
    ValenceE Lv = 6
    ValenceE F  = 7
    ValenceE Cl = 7
    ValenceE Br = 7
    ValenceE I  = 7
    ValenceE At = 7
    ValenceE Ts = 7
    ValenceE Ne = 8
    ValenceE Ar = 8
    ValenceE Kr = 8
    ValenceE Xe = 8 
    ValenceE Rn = 8 
    ValenceE Og = 8
    ValenceE _  = 2 

-- | Value-level `ValenceE` 
valenceE :: Integral a => Element -> a 
valenceE H  = 1 
valenceE Li = 1 
valenceE Na = 1 
valenceE K  = 1 
valenceE Rb = 1 
valenceE Cs = 1 
valenceE Fr = 1 
valenceE B  = 3 
valenceE Al = 3 
valenceE In = 3 
valenceE Tl = 3 
valenceE Nh = 3 
valenceE C  = 4 
valenceE Si = 4 
valenceE Ge = 4 
valenceE Sn = 4 
valenceE Pb = 4 
valenceE Fl = 4 
valenceE N  = 5
valenceE P  = 5 
valenceE As = 5 
valenceE Sb = 5 
valenceE Bi = 5 
valenceE Mc = 5 
valenceE O  = 6 
valenceE S  = 6 
valenceE Se = 6 
valenceE Te = 6 
valenceE Po = 6
valenceE Lv = 6
valenceE F  = 7
valenceE Cl = 7
valenceE Br = 7
valenceE I  = 7
valenceE At = 7
valenceE Ts = 7
valenceE Ne = 8
valenceE Ar = 8
valenceE Kr = 8
valenceE Xe = 8 
valenceE Rn = 8 
valenceE Og = 8
valenceE _  = 2 
