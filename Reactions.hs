{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE NoStarIsType #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE AllowAmbiguousTypes #-}

module Reactions
    ( Hydrocarbon
    , Molecular
    , MolReaction (..)
    , asMolecule
    , ValidHydrocarbon 
    , mkHydrocarbon 
    , CO2
    , H2O
    , O2 
    , combustHydrocarbon 
    , deltaHBondEnergy
    ) where 

import Elements
import Compounds
import Bonds
import GHC.TypeLits
import Data.Type.Bool (If (..))
import Data.Type.Equality (type (==))
import Data.Proxy 

data MolReaction = MolReaction {inRxn :: [Molecules], outRxn :: [Molecules]} deriving Eq 

instance Show MolReaction where 
    show (MolReaction (i:is) (o:os)) =  show i ++ concatMap ((++) " + " . show) is ++ " -> " 
                                     ++ show o ++ concatMap ((++) " + " . show) os  

cancelCommonElements :: Eq a => [a] -> [a] -> ([a], [a]) 
cancelCommonElements as bs = (as', bs')
    where 
        dropFirst pxs y [] = pxs 
        dropFirst pxs y (x:xs) = if x == y then pxs ++ xs else dropFirst (x:pxs) y xs 

        as' = foldr (dropFirst []) as bs
        bs' = foldr (dropFirst []) bs as 

deltaHBondEnergy :: Floating a => MolReaction -> a
deltaHBondEnergy (MolReaction i o) =
    let bdsInRaw  = concatMap (\(Molecules m n) -> concat . replicate n $ shatterToBonds m) i
        bdsOutRaw = concatMap (\(Molecules m n) -> concat . replicate n $ shatterToBonds m) o 
        (bdsIn, bdsOut) = cancelCommonElements bdsInRaw bdsOutRaw
        breaking  = sum . map bondEnergy $ bdsIn
        forming   = sum . map bondEnergy $ bdsOut
     in breaking - forming

newtype Hydrocarbon = Hydrocarbon Molecule deriving (Eq, Show)

instance Molecular Hydrocarbon where 
    asMolecule (Hydrocarbon m) = m 

type family ValidHydrocarbon (cCt :: Nat) (hCt :: Nat) :: Nat where 
    ValidHydrocarbon cCt hCt = 
        If (hCt == 2*cCt + 2)
            1 
           (If (hCt == 2*cCt)
               2 
               (If (hCt == 2*cCt - 2)
                   3 
                   (TypeError (Text "Invalid Hydrocarbon: C" :<>: ShowType cCt 
                          :<>: Text "H" :<>: ShowType hCt
            ))))

type family Replicate4H (e :: Molecule) (n :: Nat) :: Molecule where 
    Replicate4H (Central bc Neutral e '[]) 0 = Central bc Neutral e '[Terminating 1 Neutral H (4-bc)] 
    Replicate4H (Central 3  Neutral e '[]) n = Central 3 Neutral e '[Replicate4H (Central 1 Neutral e '[]) (n - 1)]
    Replicate4H (Central bc Neutral e '[]) n = 
        Central bc Neutral e '[Terminating 1 Neutral H (3 - bc), Replicate4H (Central 1 Neutral e '[]) (n - 1)]

mkHydrocarbon :: 
    forall cCt hCt bondClass m.
    ( KnownNat cCt, KnownNat hCt
    , bondClass ~ ValidHydrocarbon cCt hCt
    , m ~ If (cCt == 1)
             (Central 0 Neutral C '[Terminating 1 Neutral H hCt])
             (Central 0 Neutral C 
                (  Terminating 1 Neutral H (4 - bondClass)
                ': '[Replicate4H (Central bondClass Neutral C '[]) (cCt - 2)]))
    , KnownMolecule m) => 
    Hydrocarbon

mkHydrocarbon = Hydrocarbon $ molecValI @m

type O2  = Terminating 0 Neutral O 2 
type CO2 = Central 0 Neutral C '[Terminating 2 Neutral O 2]
type H2O = Central 0 Neutral O '[Terminating 1 Neutral H 2]


combustHydrocarbon :: Hydrocarbon -> MolReaction 
combustHydrocarbon hc = let (nc_, nh_) = case shatterSimple . asMolecule $ hc of 
                                [Atoms C ncRaw, Atoms H nhRaw] -> (ncRaw, nhRaw)
                                [Atoms H nhRaw, Atoms C ncRaw] -> (ncRaw, nhRaw)
                                other -> error $ "Impossible: Non-Hydrocarbon " ++ show other
                             
                            (nc, nh, n) = if nh_ `mod` 4 == 0 then (nc_, nh_, 1) else (nc_*2, nh_*2, 2)
                         in MolReaction [Molecules (asMolecule hc) n, Molecules (molecValI @O2) $ nc + (nh `div` 4)]
                                        [Molecules (molecValI @CO2) nc, Molecules (molecValI @H2O) (nh `div` 2)]

