{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE AllowAmbiguousTypes #-}

{-|
Module: Orbitals
Description: Better atomic orbitals using Energy Wave Theory
-}

module Orbitals 
    ( orbitalRadii
    , orbitalRadiiIon
    , atomicRadius
    , atomicRadiusIon
    , iEnergy
    , ampFacE
    ) where 

import ElectronConf
import Elements
import Slaters.Unsafe 
import GHC.TypeLits
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV 
import Control.Monad ((>=>))
import Data.Foldable (foldrM)
import System.IO.Unsafe (unsafePerformIO)
import Data.Fixed (Fixed (..), Nano)


rho :: Fractional a => a 
rho = 3.859764640e22 

ke  :: Fractional a => a 
ke  = 10

al  :: Fractional a => a 
al  = 9.215405708e-19

c   :: Fractional a => a 
c   = 299792458 

oe  :: Fractional a => a 
oe  = 2.138743820 

glamb :: Fractional a => a 
glamb = 0.9873318320

lambl :: Fractional a => a 
lambl = 2.854096501e-17

a0pm :: Fractional a => a 
a0pm = 52.917724899

a0  :: Fractional a => a 
a0  = 5.2917724899e-11 -- In m  

nA  :: Fractional a => a 
nA  = 6.02214e23 


theta :: Fractional a => QL -> a
theta = \case 
    SL -> 0.5  -- (In paper) cos 60 
    PL -> 0.75 -- (In paper) (cos 60 + cos 0) / 2
    DL -> 2/5 
    FL -> 2/3  

thetas :: Fractional a => EConf -> V.Vector ((Nat, Nat), a)
thetas (SubL 1 SL n :<-: Nucleus) = V.singleton ((1, n), 1)
thetas (SubL 3 SL n :<-: rem)     = V.cons ((3, n), 2/3) $ thetas rem
thetas rem = flip V.unfoldr rem $ \case 
                SubL pqn l e :<-: ec -> Just (((pqn, e), theta l), ec)
                Nucleus              -> Nothing 

thetasE :: Fractional a => Element -> V.Vector ((Nat, Nat), a)
thetasE = thetas . anumToEConf . toAtomic

t1Z  :: (Integral i, Fractional a) => i -> a -> a 
t1Z  z rx = fromIntegral z / (rx ^^ 2)

t1Z' :: (Integral i, Fractional a) => i -> a -> a 
t1Z' z rx = (*) (negate 2) $ fromIntegral z / (rx ^^ 3)

t2  :: Fractional a => Nat -> a -> a 
t2  n rx = (fromIntegral n ^^ 2) / (rx ^^ 3) 

t2' :: Fractional a => Nat -> a -> a 
t2' n rx = negate (3 * fromIntegral n ^^ 2) / (rx ^^ 4)

repTerm  :: Fractional a => a -> ((Nat, Nat), a) -> a -> a 
repTerm  rx ((_, e), t) ry = fromIntegral e / ((rx + t*ry)^^2)

repTerm' :: Fractional a => a -> ((Nat, Nat), a) -> a -> a 
repTerm' rx ((_, e), t) ry = (*) (negate 2) $ fromIntegral e / ((rx + t*ry)^^3)

rptDry   :: Fractional a => a -> ((Nat, Nat), a) -> a -> a 
rptDry   rx ((_, e), t) ry = (*) (negate 2) . (*) t $ fromIntegral e / ((rx + t*ry)^^3)

rOrbSingExprZEC :: (Integral i, Fractional a) => i -> EConf -> Int -> V.Vector a -> a
rOrbSingExprZEC z ec n rs = 
    let ts = thetas ec 
        (ps, rls) = V.splitAt n rs 
        (rx, ls) = (V.head rls, V.tail rls)
        rem = ps V.++ ls 
        (pts, tlts) = V.splitAt n ts 
        (((pqn, nx1), tx), lts) = (V.head tlts, V.tail tlts) 
        (SubL _ l _ :<-: _) = ec
        tem
          | pqn < 4   = pts V.++ V.map (tx <$) lts 
          | l == PL   = V.map (0.475 <$) (pts V.++ lts)
          | otherwise = pts V.++ V.map (0.75 <$) lts 
     in negate (t1Z z rx) 
        + t2 pqn rx 
        + repTerm rx ((pqn, nx1 - 1), tx) rx 
        + (sum . V.map (uncurry . flip $ repTerm rx) $ V.zip rem tem)

rOrbSingExprZEC' :: (Fractional a, Integral i) => i -> EConf -> Int -> Int -> V.Vector a -> a 
rOrbSingExprZEC' z ec n dn rs
  | n == dn =  
    let ts = thetas ec
        (ps, rls) = V.splitAt n rs 
        (rx, ls) = (V.head rls, V.tail rls)
        rem = ps V.++ ls 
        (pts, tls) = V.splitAt n ts 
        (((pqn, nx1), tx), lts) = (V.head tls, V.tail tls) 
        (SubL _ l _ :<-: _) = ec
        tem
          | pqn < 4   = pts V.++ V.map (tx <$) lts 
          | l == PL   = V.map (0.475 <$) (pts V.++ lts)
          | otherwise = pts V.++ V.map (0.75 <$) lts 
     in negate (t1Z' z rx) 
        + t2' pqn rx 
        + repTerm' rx ((pqn, nx1 - 1), tx) rx 
        + (sum . V.map (uncurry . flip $ repTerm' rx) $ V.zip rem tem)
  | otherwise = rptDry (rs V.! n) (ts V.! dn) (rs V.! dn)
        
        where 
            ts = thetas ec

rOrbExprZEC :: (Integral i, Fractional a) => i -> EConf -> V.Vector a -> V.Vector a
rOrbExprZEC z ec rs = let ts = thetas ec 
                       in V.generate (V.length ts) $ flip (rOrbSingExprZEC z ec) rs

rOrbExpr :: Fractional a => Element -> V.Vector a -> V.Vector a 
rOrbExpr e = let z = toAtomic e in rOrbExprZEC z (anumToEConf z)

rOrbJCZEC :: (Integral i, Fractional a) => i -> EConf -> V.Vector a -> V.Vector (V.Vector a)
rOrbJCZEC z ec rs = V.generate (V.length rs) (\x -> V.generate (V.length rs) $ flip (rOrbSingExprZEC' z ec x) rs) 

rOrbJC :: Fractional a => Element -> V.Vector a -> V.Vector (V.Vector a)
rOrbJC e = let z = toAtomic e in rOrbJCZEC z (anumToEConf z)

gauss :: (Fractional a, Eq a) 
      => V.Vector (V.Vector a) 
      -> V.Vector a 
      -> IO (V.Vector a)
gauss a b = do
    aM <- V.mapM V.thaw a >>= V.thaw 
    mapM_ (MV.modifyM aM (`MV.grow` 1)) ([0 .. MV.length aM - 1] :: [Int])
    let lelem = V.length b 
    V.imapM_ (\idx val -> MV.read aM idx >>= \vec -> MV.write vec lelem val) b 
    solveColumnsM aM
    bktkM aM
    where       
        solveRowM :: Fractional a
                  => MV.IOVector a 
                  -> MV.IOVector a 
                  -> Int 
                  -> IO () 
        solveRowM pivot row idx = MV.read row idx -- Row val to zero  
                              >>= \h -> MV.read pivot idx  -- Corresponding pivot value 
                              >>= \c -> MV.imapM_ (\i v -> MV.modify row (subtract $ h*v/c) i) pivot 

        solveColumnM :: Fractional a
                     => MV.IOVector (MV.IOVector a)
                     -> Int  
                     -> IO () 
        solveColumnM vec n = MV.read vec n 
                         >>= \pivot -> mapM_ (MV.read vec >=> flip (solveRowM pivot) n) 
                                             ([n + 1 .. MV.length vec - 1] :: [Int])

        solveColumnsM :: (Fractional a, Eq a)
                      => MV.IOVector (MV.IOVector a)
                      -> IO () 
        solveColumnsM mat = mapM_ (\n -> swapPivot n mat >> solveColumnM mat n) ([0..MV.length mat - 1] :: [Int])

        swapPivot :: (Fractional a, Eq a) 
                  => Int 
                  -> MV.IOVector (MV.IOVector a) 
                  -> IO ()
        swapPivot n vec = let l = MV.length vec 
                           in MV.read vec n 
                          >>= \e -> MV.read e n 
                          >>= \x -> if x == 0 then MV.swap vec n (l - 1) >> swapPivot n vec 
                                              else mapM_ (MV.modify e (/ x)) ([0..l] :: [Int])

        bktkM :: forall a. Fractional a 
              => MV.IOVector (MV.IOVector a) 
              -> IO (V.Vector a)
        bktkM mat = MV.foldrM f V.empty mat  
              where 
                  l = MV.length mat + 1 -- Rows are 1 longer than mat 
                  f :: MV.IOVector a -> V.Vector a -> IO (V.Vector a)
                  f row vec = 
                    let n = V.length vec -- Number of x vals accumulated 
                        yM = V.ifoldM' (\yv i xv -> (yv -) . (*) xv <$> MV.read row (l - 1 - n + i)) 0 vec 
                     in yM -- Subtract x*c for prior rows    
                    >>= \y -> MV.read row (l - 2 - n) -- read next x coefficient 
                    >>= \c -> (`V.cons` vec) . (/ c) . (+) y -- prepend (y + yM)/c
                    <$> MV.read row (l - 1) -- final y value 

newtonRaphson :: 
    (Ord a, Fractional a) 
      => (V.Vector a -> V.Vector (V.Vector a)) 
      -> (V.Vector a -> V.Vector a) 
      -> Int 
      -> a 
      -> V.Vector a 
      -> V.Vector a
newtonRaphson j f mxI eps st = unsafePerformIO $ -- Unfortunately, a STVector of STVector is not a thing 
    foldrM (\_ v -> (\n -> if maximum (V.toList n) < eps then v else V.zipWith (+) v n) 
                <$> gauss (j v) (V.map negate $ f v)) st ([1 .. mxI] :: [Int]) 

nrSoln :: forall a. (Fractional a, Ord a) => Element -> Int -> a -> V.Vector a -> V.Vector a
nrSoln e mxI eps = V.reverse . newtonRaphson (rOrbJC e) (rOrbExpr e) mxI eps . V.reverse

nrSolnZEC :: forall a i. (Integral i, Fractional a, Ord a) 
          => i -> EConf -> Int -> a -> V.Vector a -> V.Vector a 
nrSolnZEC z ec mxI eps = V.reverse . newtonRaphson (rOrbJCZEC z ec) (rOrbExprZEC z ec) mxI eps . V.reverse 


approxVals :: (Integral i, Fractional a) => i -> EConf -> [a]
approxVals z = flip foldrEC []
           $ \case 
                SubL n PL _ -> (:) (2 * fromIntegral n / fromIntegral z)
                SubL n _  _ -> (:) (fromIntegral n / fromIntegral z)

-- | Lowest precision that doesn't cause a div by zero error for Hydrogen 
type E8 = 100000000 

orbitalRadiiZECBohr :: forall i. Integral i => i -> EConf -> [Fixed E8]
orbitalRadiiZECBohr z ec 
  | z >= 101  = V.toList . V.map realToFrac
              . nrSolnZEC @Nano z ec 1000 1e-8 . V.fromList 
              $ approxVals z ec 
  | otherwise = V.toList 
              . nrSolnZEC z ec 1000 1e-8 . V.fromList 
              $ approxVals z ec 

-- | Yields orbital radii (s1 - highest `Sublevel`) for any `Element` 
orbitalRadii :: Element -> [Fixed E8]
orbitalRadii = map (* a0pm) . orbitalRadiiBohr

orbitalRadiiBohr :: Element -> [Fixed E8]
orbitalRadiiBohr e = let z = toAtomic e in orbitalRadiiZECBohr z (anumToEConf z)

-- | Yields obirtal radii (s1 - highest `Sublevel`) for any ion 
orbitalRadiiIon :: forall st. KnownCharge st => Species st -> [Fixed E8]
orbitalRadiiIon = map (* a0pm) . orbitalRadiiIonBohr

orbitalRadiiIonBohr :: forall st. KnownCharge st => Species st -> [Fixed E8]
orbitalRadiiIonBohr ion = 
    let (z, ec) = case specValue ion of 
                        Right e -> ( toAtomic e
                                   , formIon (charge $ chargeValI @st) $ anumToEConf (toAtomic e)) 
                        Left  _ -> error "Polyatomic ion has no orbital radii" 
     in orbitalRadiiZECBohr z ec

atomicRadiusZECBohr :: forall a i. (Fractional a, Integral i) => i -> EConf -> a 
atomicRadiusZECBohr z = realToFrac . maximum . orbitalRadiiZECBohr z

atomicRadiusZEC :: forall a i. (Fractional a, Integral i) => i -> EConf -> a 
atomicRadiusZEC z = (*) a0pm . atomicRadiusZECBohr z 

-- | Better estimate of atomic radius than Slater's rules can provide 
-- Based on EWT paper; often between current calculated and empirical values in most cases 
atomicRadius :: Fractional a => Element -> a 
atomicRadius = (*) a0pm . atomicRadiusBohr

atomicRadiusBohr :: Fractional a => Element -> a
atomicRadiusBohr = realToFrac . maximum . orbitalRadiiBohr

-- | `atomicRadius` for ions
atomicRadiusIon :: forall st a. (KnownCharge st, Fractional a) => Species st -> a 
atomicRadiusIon = (*) a0pm . atomicRadiusIonBohr

atomicRadiusIonBohr :: forall st a. (KnownCharge st, Fractional a) => Species st -> a 
atomicRadiusIonBohr = realToFrac . maximum . orbitalRadiiIonBohr

amplitudeFactorRX :: forall a. Fractional a => a -> Sublevel -> a 
amplitudeFactorRX rx = \case 
    SubL 1 SL _     -> 1 
    SubL _ SL _     -> 1 + (1/rx)
    SubL _ PL e
        | e < 4     -> 0.5*(1   + (1/rx))
        | otherwise -> 0.5*(0.5 + (1/rx))
    SubL _ DL e 
        | e < 6     -> 0.25*(1  + (1/rx))
        | otherwise -> 0.25*(0.25 + (1/rx))
    SubL _ FL e
        | e < 8     -> 0.125*(1 + (1/rx))
        | otherwise -> 0.125*(0.125 + (1/rx))

amplitudeFactorZEC :: forall a i. (Fractional a, Integral i) => i -> EConf -> a 
amplitudeFactorZEC z = \case 
    ec@(sl :<-: _) -> amplitudeFactorRX (realToFrac . last $ orbitalRadiiZECBohr z ec) sl

-- | Amplitude factor for a given `Element` (dimensionless)
ampFacE :: forall a. Fractional a => Element -> a 
ampFacE e = let z = toAtomic e 
             in case anumToEConf z of 
                  sl :<-: ec -> amplitudeFactorRX (realToFrac . last $ orbitalRadiiBohr e) sl 

constantSkeleton :: Floating a => a
constantSkeleton = glamb*(pi*rho*(ke^^7)*(al^^6)*(c^^2)*oe)/(3*(lambl^^2))

iEnergyAmp :: Floating a => a -> a -> a -- J
iEnergyAmp r0 delta = constantSkeleton*(2*delta)/(a0*r0)

iEnergyZEC :: forall a i. (Floating a, Integral i) => i -> EConf -> a 
iEnergyZEC z = \case 
        ec@(sl :<-: _) -> let r0 = realToFrac . last $ orbitalRadiiZECBohr z ec
                           in iEnergyAmp r0 (amplitudeFactorRX r0 sl)

-- | First ionization energy for a given `Element` (kJ/mol)
iEnergy :: Floating a => Element -> a 
iEnergy e = let z = toAtomic e in (/ 1000) . (*) nA $ iEnergyZEC z (anumToEConf z)
