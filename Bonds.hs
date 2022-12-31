{-# LANGUAGE DataKinds #-}
{-# LANGUAGE LambdaCase #-}

module Bonds where 

import Elements
import Compounds 
import Slaters.Unsafe
import Slaters
import GHC.TypeLits

-- | Purely value-level type (`Int` not `Nat`) representing a bond between two `Element`s and the order of said bond 
data Bond = Bond Int Element Element deriving Eq

instance Show Bond where 
    show (Bond n e1 e2) 
        = let bondType
                | n == 1 = "-"
                | n == 2 = "="
                | n == 3 = "-=-"
           in show e1 ++ bondType ++ show e2 

-- | Yields the length of some `Bond` approximated by the sum of covalent radii 
bondLength :: Floating a => Bond -> a 
bondLength (Bond _ e1 e2) = rCov' e1 + rCov' e2 

-- | Very rough approximation of "reflexive" single bond energy (A-A)
-- Regression on covalent radius, valence electrons, and effective nuclear charge mirroring Coulomb's law 
reflexBondE :: Floating a => Element -> a 
reflexBondE e = (*) 2.4226377319
              . (*) (zEff  e ** (-2.22550756547)) 
              . (/ (rCov' e **   (-1.2185578373))) 
              $ fromIntegral (valenceE e) ** 1.46107382193

-- | Helper function, eV ~= kJ/96.485
kJ2eV :: Floating a => a -> a 
kJ2eV = (/ 96.485)

-- | Helper function, kJ ~= eV*96.485
eV2kJ :: Floating a => a -> a 
eV2kJ = (*) 96.485

-- | Helper function, kCal ~= kJ/4.184
kJ2kCal :: Floating a => a -> a 
kJ2kCal = (/ 4.184)

-- | Helper function, kJ ~= kCal*4.184
kCal2kJ :: Floating a => a -> a 
kCal2kJ = (*) 4.184

-- | Single bond energy based on Pauling electronegativity difference for dissociation energy 
--
-- > Ed(AB) = 0.5(Ed(AA) + Ed(BB)) + (XA - XB)^2
-- Either very close or not at all 
singBondE :: Floating a => Element -> Element -> a 
singBondE e1 e2 = eV2kJ 
                $ (kJ2eV (reflexBondE e1) + kJ2eV (reflexBondE e2))/2
                + (eNeg' e1 - eNeg' e2)^^2

-- | Higher order bond energies, small dataset used for a regression between single bond energy and bond length 
bondEnergy :: Floating a => Bond -> a 
bondEnergy = \case 
    Bond 1 e1 e2     -> singBondE e1 e2 
    b@(Bond 2 e1 e2) -> (/) (singBondE e1 e2 * 1857.55285277)
                      $ bondLength b ** 1.34627103448
    b@(Bond 3 e1 e2) -> (+) 387.249296239
                      . (/) (singBondE e1 e2 * 5.47446611285e7)
                      $ bondLength b ** 3.41631365146
    _ -> error "Variable bond order not implemented"

-- | Sum of `bondEnergy`/ies within some `Molecular` compound 
totalBondEnergy :: (Molecular m, Floating a) => m -> a 
totalBondEnergy = sum . map bondEnergy . shatterToBonds . asMolecule 

-- | Break apart some `Molecule` into component `Bond`s
shatterToBonds :: Molecule -> [Bond] 
shatterToBonds a@(Terminating 0 _ e n) 
  | e == H    = [Bond 1 H H]
  | e == S    = replicate 8 $ Bond 1 S S   
  | e == P    = replicate 6 $ Bond 1 P P
  | n == 2    = [Bond (8 - valenceE e) e e] 
  | n == 1    = [] 
  | otherwise = error "Impossible! Invalid molecule"
shatterToBonds (Terminating {}) = []
shatterToBonds (Central _ _ ep ms) 
  = flip concatMap ms $ \case 
        Terminating b _ e n  -> replicate (fromIntegral n) $ Bond (fromIntegral b) ep e
        c@(Central b _ e _)  -> Bond (fromIntegral b) ep e : shatterToBonds c
