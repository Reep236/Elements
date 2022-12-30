{-# LANGUAGE LambdaCase #-}

module Slaters.Unsafe where 

import Elements 
import ElectronConf 
import Slaters
import GHC.TypeLits 


effPQN :: forall i a. (Integral i, Fractional a) => i -> a   
effPQN = \case  
           1 -> 1 
           2 -> 2 
           3 -> 3 
           4 -> 3.7
           5 -> 4.0 
           6 -> 4.2
           7 -> 4.3


orbitalExp :: forall a. Fractional a => Element -> a 
orbitalExp e = (zEff e)/(effPQN . getN . anumToSublevel . toAtomic $ e)

orbitalExpIon :: forall a st. (Fractional a, KnownCharge st) => Species st -> Maybe a 
orbitalExpIon s = case specValue s of 
                    Right e  -> case formIonSL (specCharge s) . anumToSublevel . toAtomic $ e of 
                                  Right sl -> (/ (effPQN . getN $ sl)) <$> zEffIon s  
                                  Left  _  -> Nothing 
                    Left  _  -> Nothing  

rMax :: forall a. Fractional a => Element -> a 
rMax e = (* 53) $ (fromIntegral . getN . anumToSublevel . toAtomic $ e) / orbitalExp e

rCov' :: forall a. Floating a => Element -> a 
rCov' = rCov . rMax

rMaxIon :: forall a st. (Fractional a, KnownCharge st) => Species st -> Maybe a  
rMaxIon s = case specValue s of 
              Right e -> case formIonSL (specCharge s) . anumToSublevel . toAtomic $ e of 
                           Right sl -> (* 53) . (/) (fromIntegral . getN $ sl) <$> orbitalExpIon s 
                           Left  _  -> Nothing 
              Left  _ -> Nothing  

rCovIon :: forall a st. (Floating a, KnownCharge st) => Species st -> Maybe a 
rCovIon = fmap rCov . rMaxIon

eNeg' :: forall a. Floating a => Element -> a 
eNeg' e = eNeg (zEff e) (rCov' e)

eNegIon :: forall a st. (Floating a, KnownCharge st) => Species st -> Maybe a 
eNegIon s = eNeg <$> zEffIon s <*> rCovIon s

eNegAR' :: forall a. Floating a => Element -> a 
eNegAR' e = eNegAR (zEff e) (rCov' e)

eNegARIon :: forall a st. (Floating a, KnownCharge st) => Species st -> Maybe a 
eNegARIon s = eNegAR <$> zEffIon s <*> rCovIon s

percentIonic' :: forall a. Floating a => Element -> Element -> a 
percentIonic' e1 e2 = percentIonic (eNeg' e1) (eNeg' e2)
