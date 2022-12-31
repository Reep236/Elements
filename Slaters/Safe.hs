{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}

{-|
Module: Slaters.Safe
Description: Type-level interface for Slaters-derived values 
-}

module Slaters.Safe 
    {-# WARNING "There is little benefit to a type-level interface given all elements through Oganesson are considered 'Valid.' Subject to removal unless further functionality neccessitates a type-level interface. All functions are extremely slow past Sodium, since they calculate the full electron configuration" #-} where 

import Elements 
import ElectronConf 
import Slaters 
import GHC.TypeLits
import Data.Proxy 

-- | Constraint representing valid @n@ for which there is an @n*@
type ValidPQN (n :: Nat) = (1 <= n, n <= 7, KnownNat n)

-- | Reifies some @n@ as its associated @n*@ 
effPQN :: forall n a. (ValidPQN n, Fractional a) => a 
effPQN = case natVal $ Proxy @n of 
           1 -> 1 
           2 -> 2 
           3 -> 3 
           4 -> 3.7
           5 -> 4.0 
           6 -> 4.2
           7 -> 4.3 

-- | Yields an orbital exponent, @(Z-S)/n*@, given some type-level `Element` satisfying `ValidPQN` 
orbitalExp :: forall e a n. (Fractional a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
orbitalExp = (zEff $ elemValI @e) / (effPQN @n)

-- | Reifies some `Element` satisfying `ValidPQN` as its max atomic radius in pm
-- (53 pm/a0)*(n / `orbitalExp`)
rMax :: forall e a n. (Fractional a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
rMax = (* 53) $ (fromIntegral . natVal $ Proxy @n)/(orbitalExp @e)

-- | Approximate covalent radius for some type-level `Element`
rCov' :: forall e a n. (Floating a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
rCov' = rCov $ rMax @e  

-- | Allred-Rochow electronegativity given some type-level `Element` 
eNegAR' :: forall e a n. (Floating a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
eNegAR' = eNegAR (zEff $ elemValI @e) (rCov' @e)

-- | Approximate Pauling electronegativity given some type-level `Element` 
eNeg' :: forall e a n. (Floating a, KnownElem e, ValidPQN n, n ~ ANum2PQN (ToAtomic e)) => a 
eNeg' = eNeg (zEff $ elemValI @e) (rCov' @e)

-- | Percent Ionic character using Pauling electronegativity for two type-level `Element`s 
percentIonic' :: 
    forall e1 e2 a n1 n2. 
        ( Floating a, KnownElem e1, KnownElem e2
        , ValidPQN n1, n1 ~ ANum2PQN (ToAtomic e1)
        , ValidPQN n2, n2 ~ ANum2PQN (ToAtomic e2))
        => a 
percentIonic' = percentIonic (eNeg' @e1) (eNeg' @e2) 
