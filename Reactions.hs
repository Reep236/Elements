{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE NoStarIsType #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE AllowAmbiguousTypes #-}

{-|
Module: Reactions 
Description: Molecular (for the time) reactions and rough thermochemistry
-}

module Reactions
    ( Hydrocarbon
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

-- | Data type representing some reaction with reactants and products. 
-- Termed "inRxn" and "outRxn" to leave "reactants" and "products" as future getters when reactions are validated on creation.
data MolReaction = MolReaction {inRxn :: [Molecules], outRxn :: [Molecules]} deriving Eq 

instance Show MolReaction where 
    show (MolReaction (i:is) (o:os)) =  show i ++ concatMap ((++) " + " . show) is ++ " -> " 
                                     ++ show o ++ concatMap ((++) " + " . show) os  

-- | Helper function to remove the proper *number* of duplicates between each of two lists 
-- Seems to have high time complexity for what is aiming to be accomplished 
cancelCommonElements :: Eq a => [a] -> [a] -> ([a], [a]) 
cancelCommonElements as bs = (as', bs')
    where 
        dropFirst pxs y [] = pxs 
        dropFirst pxs y (x:xs) = if x == y then pxs ++ xs else dropFirst (x:pxs) y xs 

        as' = foldr (dropFirst []) as bs
        bs' = foldr (dropFirst []) bs as 

-- | Uses `bondEnergy` to compute approximate ethalpy change for a gaseous reaction 
-- Relatively inaccurate at present 
deltaHBondEnergy :: Floating a => MolReaction -> a
deltaHBondEnergy (MolReaction i o) =
    let bdsInRaw  = concatMap (\(Molecules m n) -> concat . replicate n $ shatterToBonds m) i
        bdsOutRaw = concatMap (\(Molecules m n) -> concat . replicate n $ shatterToBonds m) o 
        (bdsIn, bdsOut) = cancelCommonElements bdsInRaw bdsOutRaw
        breaking  = sum . map bondEnergy $ bdsIn
        forming   = sum . map bondEnergy $ bdsOut
     in breaking - forming

-- | newtype wrapper for a Hydrocarbon `Molecule`  
newtype Hydrocarbon = Hydrocarbon Molecule deriving (Eq, Show)

instance Molecular Hydrocarbon where 
    asMolecule (Hydrocarbon m) = m 

-- | Type family yielding the order for single- (-ane), double- (-ene), or triple- (-yne) bonded hydrocarbon 
-- Emits a `TypeError` for none of the above 
type family ValidHydrocarbon (cCt :: Nat) (hCt :: Nat) :: Nat where 
    ValidHydrocarbon cCt hCt = 
        If (hCt == 2*cCt + 2)
            1 
           (If (hCt == 2*cCt)
               2 
               (If (hCt == 2*cCt - 2)
                   3 
                   (TypeError (Text "Invalid simple Hydrocarbon: C" :<>: ShowType cCt 
                          :<>: Text "H" :<>: ShowType hCt
            ))))

-- | Fills a Carbon (or any other `Element` accepting 4 covalent bonds) chain with Hydrogens to length @n@
type family Replicate4H (e :: Molecule) (n :: Nat) :: Molecule where 
    Replicate4H (Central bc Neutral e '[]) 0 = Central bc Neutral e '[Terminating 1 Neutral H (4-bc)] 
    Replicate4H (Central 3  Neutral e '[]) n = Central 3 Neutral e '[Replicate4H (Central 1 Neutral e '[]) (n - 1)]
    Replicate4H (Central bc Neutral e '[]) n = 
        Central bc Neutral e '[Terminating 1 Neutral H (3 - bc), Replicate4H (Central 1 Neutral e '[]) (n - 1)]

-- | Reifies a `Hydrocarbon` given a number of `C`s and `H`s 
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

-- | Helpful type synonym for the `Molecule` representing Oxygen 
type O2  = Terminating 0 Neutral O 2 
-- | Helpful type synonym for the `Molecule` representing Carbon Dioxide 
type CO2 = Central 0 Neutral C '[Terminating 2 Neutral O 2]
-- | Helpful type synonym for the `Molecule` representing Water 
type H2O = Central 0 Neutral O '[Terminating 1 Neutral H 2]


-- | Yields the `MolReaction` representing the complete combustion of some `Hydrocarbon` 
combustHydrocarbon :: Hydrocarbon -> MolReaction 
combustHydrocarbon hc = let (nc_, nh_) = case shatterSimple . asMolecule $ hc of 
                                [Atoms C ncRaw, Atoms H nhRaw] -> (ncRaw, nhRaw)
                                [Atoms H nhRaw, Atoms C ncRaw] -> (ncRaw, nhRaw)
                                other -> error $ "Impossible: Non-Hydrocarbon " ++ show other
                             
                            (nc, nh, n) = if nh_ `mod` 4 == 0 then (nc_, nh_, 1) else (nc_*2, nh_*2, 2)
                         in MolReaction [Molecules (asMolecule hc) n, Molecules (molecValI @O2) $ nc + (nh `div` 4)]
                                        [Molecules (molecValI @CO2) nc, Molecules (molecValI @H2O) (nh `div` 2)]
