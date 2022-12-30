{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}

module Slaters.Safe where 

import Elements 
import ElectronConf 
import Slaters 
import GHC.TypeLits
import Data.Proxy 


type ValidPQN (n :: Nat) = (1 <= n, n <= 7, KnownNat n)

effPQN :: forall n a. (ValidPQN n, Fractional a) => a 
effPQN = case natVal $ Proxy @n of 
           1 -> 1 
           2 -> 2 
           3 -> 3 
           4 -> 3.7
           5 -> 4.0 
           6 -> 4.2
           7 -> 4.3 

orbitalExp :: forall e a n. (Fractional a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
orbitalExp = (zEff $ elemValI @e) / (effPQN @n)

rMax :: forall e a n. (Fractional a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
rMax = (* 53) $ (fromIntegral . natVal $ Proxy @n)/(orbitalExp @e)

rCov' :: forall e a n. (Floating a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
rCov' = rCov $ rMax @e  

eNegAR' :: forall e a n. (Floating a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
eNegAR' = eNegAR (zEff $ elemValI @e) (rCov' @e)

-- Way too slow for elements past Na
eNeg' :: forall e a n. (Floating a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
eNeg' = eNeg (zEff $ elemValI @e) (rCov' @e)

percentIonic' :: 
    forall e1 e2 a n1 n2. 
        ( Floating a, KnownElem e1, KnownElem e2
        , ValidPQN n1, n1 ~ ANum2PQN (ToAtomic e1)
        , ValidPQN n2, n2 ~ ANum2PQN (ToAtomic e2))
        => a 
percentIonic' = percentIonic (eNeg' @e1) (eNeg' @e2) 
