{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}

{-|
Module: Slaters.Unsafe 
Description: Value-level interface for Slaters-derived quantities 
-}
module Slaters.Unsafe where 

import Elements 
import ElectronConf 
import Slaters
import GHC.TypeLits 

-- | Yields n* from n (Slater's effective principal quantum number)
effPQN :: (Integral i, Fractional a) => i -> a   
effPQN = \case  
           1 -> 1 
           2 -> 2 
           3 -> 3 
           4 -> 3.7
           5 -> 4.0 
           6 -> 4.2
           7 -> 4.3

-- | Orbital Exponent 
--
-- > (Z-S)/n* 
orbitalExp :: Fractional a => Element -> a 
orbitalExp e = (zEff e)/(effPQN . getN . anumToSublevel . toAtomic $ e)

-- | `orbitalExp` using a given `Species` with associated `zEffIon` 
-- `Maybe` accounts for the potential of an invalid anion (`Left` `Nucleus`)
orbitalExpIon :: (Fractional a, KnownCharge st) => Species st -> Maybe a 
orbitalExpIon s = case specValue s of 
                    Right e  -> case formIonSL (specCharge s) . anumToSublevel . toAtomic $ e of 
                                  Right sl -> (/ (effPQN . getN $ sl)) <$> zEffIon s  
                                  Left  _  -> Nothing 
                    Left  _  -> Nothing  

-- | Returns max/absolute atomic radius in picometers 
-- (53 pm / a0) * (n / `orbitalExp`)
rMax :: Fractional a => Element -> a 
rMax e = (* 53) $ (fromIntegral . getN . anumToSublevel . toAtomic $ e) / orbitalExp e

-- | `rCov` given just the `Element` 
rCov' :: Floating a => Element -> a 
rCov' = rCov . rMax

-- | `rMax` using a given `Species` with associated `orbitalExpIon` 
rMaxIon :: (Fractional a, KnownCharge st) => Species st -> Maybe a  
rMaxIon s = case specValue s of 
              Right e -> case formIonSL (specCharge s) . anumToSublevel . toAtomic $ e of 
                           Right sl -> (* 53) . (/) (fromIntegral . getN $ sl) <$> orbitalExpIon s 
                           Left  _  -> Nothing 
              Left  _ -> Nothing  

-- | `rCov` using a given `Species` with associated `rMaxIon`
rCovIon :: (Floating a, KnownCharge st) => Species st -> Maybe a 
rCovIon = fmap rCov . rMaxIon

-- | `eNeg` given just the `Element`
eNeg' :: Floating a => Element -> a 
eNeg' e = eNeg (zEff e) (rCov' e)

-- | `eNeg` using a given `Species` and associated `zEffIon` and `rCovIon`
eNegIon :: (Floating a, KnownCharge st) => Species st -> Maybe a 
eNegIon s = eNeg <$> zEffIon s <*> rCovIon s

-- | `eNegAR` given just the `Element` 
eNegAR' :: Floating a => Element -> a 
eNegAR' e = eNegAR (zEff e) (rCov' e)

-- | `eNegAR` using a given `Species` and associated `zEffIon` and `rCovIon`
eNegARIon :: (Floating a, KnownCharge st) => Species st -> Maybe a 
eNegARIon s = eNegAR <$> zEffIon s <*> rCovIon s

-- | `percentIonic` given just two `Element`s 
percentIonic' :: Floating a => Element -> Element -> a 
percentIonic' e1 e2 = percentIonic (eNeg' e1) (eNeg' e2)
