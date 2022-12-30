{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE LambdaCase #-}

import Elements
import ElectronConf
import Compounds
import Reactions 
import Bonds 
import Slaters.Unsafe
import Slaters 
import GHC.TypeLits
import Data.Proxy
import Data.Type.Bool 
import Data.Type.Equality


readBond :: String -> (Bond, Double)
readBond = \case
    x1:x2:x3:x4:x5:xs 
     | x2 `elem` ['–','=','≡'] -> ( Bond (pbt x2) (read [x1]) (rdNxt x3 x4)
                                  , if x5 `elem` ['1'..'9'] then read $ x5:xs else read xs)
     | otherwise               -> (Bond (pbt x3) (read [x1, x2]) (rdNxt x4 x5)
                                  , if head xs `elem` ['1'..'9'] then read xs else read $ drop 1 xs)

    where 
        pbt :: Char -> Int 
        pbt '–' = 1 
        pbt '=' = 2 
        pbt '≡' = 3 

        rdNxt :: Char -> Char -> Element
        rdNxt x1 x2 = read $ x1 : [x2 | x2 `elem` ['a'..'z']]

bondEnergy' :: Floating a => Bond -> (Bond, a)
bondEnergy' b = (b, kJ2kCal $ bondEnergy b)

main :: IO()
main = do
    print . bondEnergy' $ Bond 1 C C  
    print . bondEnergy' $ Bond 2 C C  
    print . bondEnergy' $ Bond 3 C C  
    print . bondEnergy' $ Bond 1 C H 
    print . bondEnergy' $ Bond 2 C O
    print . bondEnergy' $ Bond 1 O H 
