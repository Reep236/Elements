{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE LambdaCase #-}

{-|
Module: ElectronConf
Description: Quantum numbers and electron configuration
-}

module ElectronConf 
    ( QL (..)
    , KnownQL
    , qlVal 
    , qlValI 
    , Sublevel (..)
    , EConf (..)
    , GetN
    , peakN 
    , lastSubL 
    , anumToEConf 
    , anumToSublevel
    , ANum2PQN 
    , ANumLastSL
    , ANumToEConf
    , addElectrons 
    , addElectronsSL 
    , dropElectrons 
    , dropElectronsSL 
    , formIon
    , formIonSL
    , getN
    , getMS
    , getML
    , getL 
    , foldrEC
    , lengthEC
    , slValI 
    )

    where

import GHC.TypeLits
import Data.Type.Bool 
import Data.Type.Equality 
import Data.Proxy


-- | Type representing states of the azimuthal or angular momentum quantum number, l 
data QL = SL | PL | DL | FL deriving (Eq, Ord, Enum)

-- | Internal, singleton `QL` for reflection
newtype QlSing (a :: QL) = QlSing QL 

-- | Type class representing a reifiable azimuthal quantum number 
class KnownQL (a :: QL) where 
    -- | Internal, yields singleton QL 
    qlSing :: QlSing a  

instance KnownQL SL where 
    qlSing = QlSing SL 

instance KnownQL PL where 
    qlSing = QlSing PL 

instance KnownQL DL where 
    qlSing = QlSing DL 

instance KnownQL FL where 
    qlSing = QlSing FL 

-- | Reification of `QL` using `Proxy`
qlVal :: forall q. KnownQL q => Proxy q -> QL
qlVal _ = case (qlSing :: QlSing q) of 
            QlSing q -> q

-- | Reification of `QL` independent of `Proxy`
qlValI :: forall q. KnownQL q => QL
qlValI = qlVal $ Proxy @q

instance Show QL where 
    show SL = "s"
    show PL = "p"
    show DL = "d"
    show FL = "f"

-- | Readable definition of the quantum spin number, ms 
data QMS = Up | Down deriving (Show, Eq, Ord, Enum) 

-- | Type representing a given sublevel
-- For example, 1s2 = `SubL 1 SL 2` 
data Sublevel = SubL Nat QL Nat deriving Eq  

infixr 5 :<-:

-- | Definition of electron configuration equivalent to a singly linked list 
data EConf = Nucleus | Sublevel :<-: EConf  deriving Eq 

-- | Reification of a `Sublevel` independent of `Proxy`
slValI :: forall ec n l e. 
    ( KnownQL l
    , KnownNat n
    , KnownNat e
    , SubL n l e ~ ec) 
    => Sublevel  
slValI = SubL (fromIntegral . natVal $ Proxy @n) (qlValI @l) (fromIntegral . natVal $ Proxy @e)

-- | Extract Principal Quantum Number (n) from `Sublevel` 
getN :: Integral a => Sublevel -> a 
getN = \case 
    SubL n _ _        -> fromIntegral n 

-- | Extract Angular Momentum Quantum Number (l) from `Sublevel` 
getL :: Sublevel -> QL 
getL = \case 
    SubL _ l _        -> l 

-- | Extract Magnetic Quantum Number (ml) from `Sublevel` 
getML :: Sublevel -> Int
getML = \case 
    SubL _ l e -> 
        case l of
          SL -> 1 
          PL -> fromIntegral $ (e - 1) `mod` 3 - 1
          DL -> fromIntegral $ (e - 1) `mod` 5 - 2 
          FL -> fromIntegral $ (e - 1) `mod` 7 - 3 

-- | Extract Spin Quantum Number (ms) from `Sublevel` 
getMS :: Sublevel -> QMS 
getMS = \case  
    SubL _ l e -> 
        case l of 
            SL -> toEnum . fromIntegral $  e - 1 
            PL -> toEnum . fromIntegral $ (e - 1) `div` 3 
            DL -> toEnum . fromIntegral $ (e - 1) `div` 5 
            FL -> toEnum . fromIntegral $ (e - 1) `div` 7

-- | foldr definition as EConf cannot have a `Foldable` instance 
foldrEC :: (Sublevel -> b -> b) -> b -> EConf -> b 
foldrEC f init = \case 
    ec :<-: ecs -> foldrEC f (f ec init) ecs 
    Nucleus -> init

lengthEC :: EConf -> Int 
lengthEC = \case 
    Nucleus   -> 0 
    _ :<-: ec -> 1 + lengthEC ec 

instance Show Sublevel where 
    show (SubL a b c) = show a ++ show b ++ show c

instance Show EConf where 
    show Nucleus            = "Nucleus"
    show (ec1 :<-: Nucleus) = show ec1 
    show (ec1 :<-: ec2)     = show ec2 ++ "|" ++ show ec1


-- | Type level conversion of an atomic number to electron configuration 
-- Far too slow for elements past Sodium 

{-# WARNING ANumToEConf "Greatly lengthens type-checking, avoid if possible " #-}
type family ANumToEConf (n :: Nat) :: EConf  where 
    ANumToEConf 0 = TypeError (Text "No Atomic # 0")
    ANumToEConf 1 = SubL 1 SL 1 :<-: Nucleus 
    ANumToEConf 2 = SubL 1 SL 2 :<-: Nucleus 
    ANumToEConf 3 = SubL 2 SL 1 :<-: (SubL 1 SL 2 :<-: Nucleus) 
    ANumToEConf 4 = SubL 2 SL 2 :<-: (SubL 1 SL 2 :<-: Nucleus)
    ANumToEConf a = AN2ECHelper (a - 4) (ANumToEConf 4) 

type family AN2ECHelper (z :: Nat) (st :: EConf) :: EConf where 
    AN2ECHelper 0 st = st
    AN2ECHelper z (SubL n SL e :<-: rem) = 
        If (e == 1) 
           (AN2ECHelper (z - 1) (SubL n SL 2 :<-: rem))
           (If (n <=? 3) 
               (AN2ECHelper (z - 1) (SubL n PL 1 :<-: (SubL n SL 2 :<-: rem)))
               (If (n <=? 5)
                   (AN2ECHelper (z - 1) (SubL (n - 1) DL 1 :<-: (SubL n SL 2 :<-: rem)))
                   (AN2ECHelper (z - 1) (SubL (n - 2) FL 1 :<-: (SubL n SL 2 :<-: rem)))))

    AN2ECHelper z (SubL n PL e :<-: rem) = 
        If (e <=? 5) 
           (AN2ECHelper (z - 1) (SubL n PL (e + 1) :<-: rem)) 
           (AN2ECHelper (z - 1) (SubL (n + 1) SL 1 :<-: (SubL n PL 6 :<-: rem)))

    AN2ECHelper z (SubL n DL e :<-: rem) = 
        If (e <=? 9)
           (AN2ECHelper (z - 1) (SubL n DL (e + 1) :<-: rem))
           (AN2ECHelper (z - 1) (SubL (n + 1) PL 1 :<-: (SubL n DL 10 :<-: rem)))

    AN2ECHelper z (SubL n FL e :<-: rem) = 
        If (e <=? 13)
           (AN2ECHelper (z - 1) (SubL n FL (e + 1) :<-: rem))
           (AN2ECHelper (z - 1) (SubL (n + 1) DL 1 :<-: (SubL n FL 14 :<-: rem)))

-- | Type-level equivalent of `getN` 
type family GetN (ec :: Sublevel) :: Nat where 
    GetN (SubL n _ _)  = n 

-- | Shortcut type-level conversion of atomic number to PQN 
type family ANum2PQN (z :: Nat) :: Nat where 
    ANum2PQN 1 = 1 
    ANum2PQN 2 = 1
    ANum2PQN 3 = 2
    ANum2PQN 4 = 2
    ANum2PQN z = GetN (ANumLastSL z (SubL 2 SL 2))

-- | Type-level conversion of atomic number to final sublevel, equivalent to 
-- @ case ANumToEConf z of 
--      sl :<-: ec -> sl @ 
-- Without extra expense 

{-# WARNING ANumLastSL "Greatly lengthens type-checking, avoid if possible" #-}
type family ANumLastSL (z :: Nat) (ec :: Sublevel) :: Sublevel where 
    ANumLastSL 0 st = st
    ANumLastSL z (SubL n SL e) = 
        If (e == 1) 
           (ANumLastSL (z - 1) (SubL n SL 2))
           (If (n <=? 3) 
               (ANumLastSL (z - 1) (SubL n PL 1))
               (If (n <=? 5)
                   (ANumLastSL (z - 1) (SubL (n - 1) DL 1))
                   (ANumLastSL (z - 1) (SubL (n - 2) FL 1))))

    ANumLastSL z (SubL n PL e) = 
        If (e <=? 5) 
           (ANumLastSL (z - 1) (SubL n PL (e + 1))) 
           (ANumLastSL (z - 1) (SubL (n + 1) SL 1))

    ANumLastSL z (SubL n DL e) = 
        If (e <=? 9)
           (ANumLastSL (z - 1) (SubL n DL (e + 1)))
           (ANumLastSL (z - 1) (SubL (n + 1) PL 1))

    ANumLastSL z (SubL n FL e) = 
        If (e <=? 13)
           (ANumLastSL (z - 1) (SubL n FL (e + 1)))
           (ANumLastSL (z - 1) (SubL (n + 1) DL 1))

-- | Value-level equivalent of `ANumToEConf` 
anumToEConf :: Int -> EConf 
anumToEConf = \case 
    1 -> SubL 1 SL 1 :<-: Nucleus
    2 -> SubL 1 SL 2 :<-: Nucleus 
    3 -> SubL 2 SL 1 :<-: SubL 1 SL 2 :<-: Nucleus  
    4 -> SubL 2 SL 2 :<-: SubL 1 SL 2 :<-: Nucleus 
    a -> addElectrons (a - 4) (anumToEConf 4) 

-- | Adding electrons to electron configuration
-- Useful for anion formation or stepping through elements 
addElectrons :: Int -> EConf -> EConf 
addElectrons 0 = id 
addElectrons n = \case 
    Nucleus -> addElectrons (n-1) (SubL 1 SL 1 :<-: Nucleus)
    a       -> foldr (\n acc@(b:<-:bs) -> let nx = addElectronsSL n b in nx :<-: (if getL nx == getL b then bs else acc)) a 
             $ replicate n 1

-- | Dropping electrons from an electron configuration 
-- Useful for cation formation or stepping through elements 
-- Yields no warning regarding being left with just the nucleus 
dropElectrons :: Int -> EConf -> EConf 
dropElectrons 0 = id 
dropElectrons n = \case 
    Nucleus   -> Nucleus
    SubL _ _ 1 :<-: ec -> dropElectrons (n - 1) ec 
    SubL p l a :<-: ec -> dropElectrons (n - 1) $ SubL p l (a-1) :<-: ec 

-- | Formation of an ion electron configuration from `Int` charge and electron configuration 
formIon :: Int -> EConf -> EConf
formIon n 
    | n == 0 = id 
    | n >  0 = dropElectrons n
    | n <  0 = addElectrons  $ negate n 

-- | Formation of an ion yielding just the final `Sublevel`
-- Any `Left` value indicates just `Nucleus` 
formIonSL :: Int -> Sublevel -> Either EConf Sublevel 
formIonSL n 
    | n == 0 = Right 
    | n >  0 = dropElectronsSL n
    | n <  0 = Right . addElectronsSL (negate n)  

-- | `anumToEConf` yielding just the final `Sublevel`
-- Value equivalent of `ANumLastSL`
anumToSublevel :: Int -> Sublevel 
anumToSublevel = \case 
    1 -> SubL 1 SL 1
    2 -> SubL 1 SL 2 
    3 -> SubL 2 SL 1   
    4 -> SubL 2 SL 2  
    a -> addElectronsSL (a - 4) (anumToSublevel 4) 

-- | `addElectrons` yielding just the final `Sublevel` 
addElectronsSL :: Int -> Sublevel -> Sublevel
addElectronsSL 0 a = a 
addElectronsSL z (SubL n SL e) 
    | e == 1    = addElectronsSL (z - 1)  $ SubL n SL 2 
    | n <= 3    = addElectronsSL (z - 1)  $ SubL n PL 1  
    | n <= 5    = addElectronsSL (z - 1)  $ SubL (n - 1) DL 1  
    | otherwise = addElectronsSL (z - 1)  $ SubL (n - 2) FL 1  
addElectronsSL z (SubL n PL e) = 
    if e <= 5 then addElectronsSL (z - 1) $ SubL n PL (e + 1) 
              else addElectronsSL (z - 1) $ SubL (n + 1) SL 1 
addElectronsSL z (SubL n DL e) = 
    if e <= 9 then addElectronsSL (z - 1) $ SubL n DL (e + 1)
              else addElectronsSL (z - 1) $ SubL (n + 1) PL 1 
addElectronsSL z (SubL n FL e) = 
    if e <= 13 then addElectronsSL (z - 1) $ SubL n FL (e + 1)
               else addElectronsSL (z - 1) $ SubL (n + 1) DL 1

-- | `dropElectrons` yielding just the final `Sublevel` 
-- `Left` indicates just `Nucleus` 
dropElectronsSL :: Int -> Sublevel -> Either EConf Sublevel  
dropElectronsSL 0 = Right 
dropElectronsSL n = \case 
    SubL 1 SL 1     -> Left Nucleus 
    SubL p SL 1     -> dropElectronsSL (n-1) $ SubL (p - 1) PL 6 
    SubL p PL 1 
        | p <= 3    -> dropElectronsSL (n-1) $ SubL p     SL 2 
        | p <= 5    -> dropElectronsSL (n-1) $ SubL (p-1) DL 10 
        | otherwise -> dropElectronsSL (n-1) $ SubL (p-2) FL 14
    SubL p DL 1     -> dropElectronsSL (n-1) $ SubL (p+1) SL 2 
    SubL p FL 1     -> dropElectronsSL (n-1) $ SubL (p+2) SL 2 
    SubL p l  e | e > 1 -> dropElectronsSL (n - 1) (SubL p l $ e - 1)

-- | Unwraps an `EConf` to `Nucleus` or the final `Sublevel` 
lastSubL :: EConf -> Either EConf Sublevel 
lastSubL = \case 
    (ec :<-: ecs) -> Right ec
    Nucleus  -> Left Nucleus

-- | Should a `DL` or `FL` sublevel lower the PQN, returns the final `SL` sublevel
peakN :: EConf -> Either EConf Sublevel 
peakN = \case 
    ec@(SubL _ l _) :<-: ecs -> if l == PL || l == SL then Right ec else peakN ecs     
    Nucleus                  -> Left Nucleus 
