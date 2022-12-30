module Slaters where 

import Elements
import ElectronConf 
import GHC.TypeLits

screenConst :: forall a. Fractional a => EConf -> a 
screenConst ec = 
    case peakN ec of 
      Left Nucleus -> 0 
      Right (SubL 1 SL e) -> 0.3 * fromIntegral (e - 1) 
      Right (SubL n SL e) -> foldrEC (evalEff (n, SL) (0.35, 0, 0.85, 1)) 0 ec
      Right (SubL n PL e) -> foldrEC (evalEff (n, PL) (0.35, 0, 0.85, 1)) 0 ec
      Right (SubL n DL e) -> foldrEC (evalEff (n, DL) (0.35, 1, 1   , 1)) 0 ec
      Right (SubL n FL e) -> foldrEC (evalEff (n, FL) (0.35, 1, 1   , 1)) 0 ec
    where 
        evalEff :: (Natural, QL) -> (a, a, a, a) -> Sublevel -> a -> a 
        evalEff (n, l) (sameSL, lowerL, nSub1, nSub2) ec a 
          = case ec of
              SubL n' l' e | n' == n && l' == l -> a + sameSL * fromIntegral (e - 1) 
              SubL n' SL e | n' == n && PL == l -> a + sameSL * fromIntegral e 
              SubL n' l' e | n' == n && l' < l  -> a + lowerL * fromIntegral e 
              SubL n' l' e | n' == n - 1        -> a + nSub1  * fromIntegral e 
              SubL n' l' e                      -> a + nSub2  * fromIntegral e 

zEff :: forall a. Fractional a => Element -> a 
zEff e = let atm = toAtomic e 
             ec  = anumToEConf atm 
         in fromIntegral atm - screenConst ec

zEffIon :: forall a st. (KnownCharge st, Fractional a) => Species st -> Maybe a 
zEffIon st = case specValue st of
               Right e -> let atm = toAtomic e
                              ec  = formIon (specCharge st) . anumToEConf $ atm 
                           in Just $ fromIntegral atm - screenConst ec 
               Left _ -> Nothing

rCov :: forall a. Floating a => a -> a
rCov = (*) 570.114890425 . subtract 1.23657508201 . (** 0.0774677838022) 

eNegAR :: Fractional a => a -> a -> a 
eNegAR zeff rcov = 3590*zeff/(rcov ^^ 2) + 0.744

eNeg :: Floating a => a -> a -> a 
eNeg zeff = (+) 1.26455649347 . (*) 1.42106566995 . log . eNegAR zeff

percentIonic :: Floating a => a -> a -> a 
percentIonic en1 en2 = (*) 100 . (1 -) . exp . (/ 4) . negate . (^ 2) $ en2 - en1 
