{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE AllowAmbiguousTypes #-}

{-|
Module: Orbitals
Description: Better atomic orbitals using Electron Wave Theory
-}

module Orbitals 
    ( orbitalRadii
    , orbitalRadiiIon
    , atomicRadius
    , atomicRadiusIon
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

theta :: Fractional a => QL -> a
theta = \case 
    SL -> 0.5 -- (In paper) cos 60 
    PL -> 0.75 -- (In paper) (cos 60 + cos 0) / 2
    DL -> (1 + 1/2 + 0 - realToFrac (sqrt 3) / 2)/4 -- (Unofficial) (cos 0 + cos 60 + 2*cos 90 + cos 150) / 5
    FL -> (realToFrac (sqrt 2) + 2 + realToFrac (sqrt 3))/6 -- (Unofficial) (2*cos 90 + 2*cos 0 + 2*cos 30) / 6 

thetas :: Fractional a => EConf -> V.Vector (Nat, Nat, a)
thetas (SubL 1 SL n :<-: Nucleus) = V.singleton (1, n, 1)
thetas (SubL 3 SL n :<-: rem)     = V.cons (3, n, 1/3) $ thetas rem
thetas rem = flip V.unfoldr rem $ \case 
                SubL pqn l e :<-: ec -> Just ((pqn, e, theta l), ec)
                Nucleus              -> Nothing 

thetasE :: Fractional a => Element -> V.Vector (Nat, Nat, a)
thetasE = thetas . anumToEConf . toAtomic

t1Z :: (Integral i, Fractional a) => i -> a -> a 
t1Z z rx = fromIntegral z / (rx ^^ 2)

t1Z' :: (Integral i, Fractional a) => i -> a -> a 
t1Z' z rx = (*) (negate 2) $ fromIntegral z / (rx ^^ 3)

t2 :: Fractional a => Nat -> a -> a 
t2 n rx = (fromIntegral n ^^ 2) / (rx ^^ 3) 

t2' :: Fractional a => Nat -> a -> a 
t2' n rx = negate (3 * fromIntegral n ^^ 2) / (rx ^^ 4)

repTerm :: Fractional a => a -> (Nat, Nat, a) -> a -> a 
repTerm rx (_, n, t) ry = fromIntegral n / ((rx + t*ry)^^2)

repTerm' :: Fractional a => a -> (Nat, Nat, a) -> a -> a 
repTerm' rx (_, n, t) ry = (*) (negate 2) $ fromIntegral n / ((rx + t*ry)^^3)

rptDry :: Fractional a => a -> (Nat, Nat, a) -> a -> a 
rptDry rx (_, n, t) ry = (*) (negate 2) . (*) t $ fromIntegral n / ((rx + t*ry)^^3)

rOrbSingExprZEC :: (Integral i, Fractional a) => i -> EConf -> Int -> V.Vector a -> a
rOrbSingExprZEC z ec n rs = 
    let ts = thetas ec 
        (ps, rls) = V.splitAt n rs 
        (rx, ls) = (V.head rls, V.tail rls)
        rem = ps V.++ ls 
        (pts, tlts) = V.splitAt n ts 
        ((pqn, nx1, tx), lts) = (V.head tlts, V.tail tlts) 
        tem = pts V.++ lts 
     in negate (t1Z z rx) 
        + t2 pqn rx 
        + repTerm rx (pqn, nx1 - 1, tx) rx 
        + (sum . V.map (uncurry . flip $ repTerm rx) $ V.zip rem tem)

rOrbSingExprZEC' :: (Fractional a, Integral i) => i -> EConf -> Int -> Int -> V.Vector a -> a 
rOrbSingExprZEC' z ec n dn rs
  | n == dn =  
    let ts = thetas ec 
        (ps, rls) = V.splitAt n rs 
        (rx, ls) = (V.head rls, V.tail rls)
        rem = ps V.++ ls 
        (pts, tlts) = V.splitAt n ts 
        ((pqn, nx1, tx), lts) = (V.head tlts, V.tail tlts) 
        tem = pts V.++ lts 
     in negate (t1Z' z rx) 
        + t2' pqn rx 
        + repTerm' rx (pqn, nx1 - 1, tx) rx 
        + (sum . V.map (uncurry . flip $ repTerm' rx) $ V.zip rem tem)
  | otherwise = rptDry (rs V.! n) (ts V.! dn) (rs V.! dn)
        
        where 
            ts = thetas ec

rOrbExprZEC :: (Integral i, Fractional a) => i -> EConf -> V.Vector a -> V.Vector a
rOrbExprZEC z ec rs = let ts = thetas ec 
                       in V.generate (length ts) (\n -> rOrbSingExprZEC z ec n rs)

rOrbExpr :: Fractional a => Element -> V.Vector a -> V.Vector a 
rOrbExpr e = let z = toAtomic e in rOrbExprZEC z (anumToEConf z)

rOrbJCZEC :: (Integral i, Fractional a) => i -> EConf -> V.Vector a -> V.Vector (V.Vector a)
rOrbJCZEC z ec rs = V.generate (length rs) (\x -> V.generate (length rs) 
                                           (\y -> rOrbSingExprZEC' z ec x y rs))

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
        mulRowM :: Fractional a => a -> MV.IOVector a -> IO () 
        mulRowM c vec = mapM_ (MV.modify vec $ (*) c) ([0..MV.length vec - 1] :: [Int])

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

-- | Yields orbital radii (s1 - highest `Sublevel`) for any `Element` 
orbitalRadii :: Element -> [Fixed E8]
orbitalRadii e 
  | e >= Md    = V.toList . V.map ((*) 52.9177 . realToFrac) 
               . nrSoln @Nano e 1000 1e-3 . V.fromList 
               . approxVals (toAtomic e) 
               . anumToEConf 
               . toAtomic 
               $ e

  | otherwise  = V.toList . V.map (* 52.9177)
               . nrSoln e 1000 1e-3 . V.fromList 
               . approxVals (toAtomic e)
               . anumToEConf 
               . toAtomic 
               $ e

-- | Yields obirtal radii (s1 - highest `Sublevel`) for any ion 
orbitalRadiiIon :: forall st. KnownCharge st => Species st -> [Fixed E8]
orbitalRadiiIon ion = 
    let (z, ec) = case specValue ion of 
                        Right e -> ( toAtomic e
                                   , formIon (charge $ chargeValI @st) $ anumToEConf (toAtomic e)) 
                        Left  _ -> error "Polyatomic ion has no orbital radii" 
     in if z >= 101 then V.toList . V.map ((*) 52.9177 . realToFrac) 
                      . nrSolnZEC @Nano z ec 1000 1e-3 . V.fromList 
                      $ approxVals z ec

                    else V.toList . V.map (* 52.9177)
                       . nrSolnZEC z ec 1000 1e-3 . V.fromList 
                       $ approxVals z ec

-- | Better estimate of atomic radius than Slater's rules can provide 
-- Based on EWT paper; often between current calculated and empirical values in most cases 
atomicRadius :: Fractional a => Element -> a
atomicRadius = realToFrac . maximum . orbitalRadii 

-- | `atomicRadius` for ions
atomicRadiusIon :: forall st a. (KnownCharge st, Fractional a) => Species st -> a 
atomicRadiusIon = realToFrac . maximum . orbitalRadiiIon

