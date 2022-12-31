# Elements: Chemistry DSL in Haskell 

A wide variety of from-scratch Chemistry taking advantage of type- and value- level computation. Currently depends just on base and has no requirement for Stack or Cabal. GHC and GHCi should handle everything fine themselves. Haddock (WIP) is available through Github pages. 

# Elements # 

Contains the definitions of elements H-Og along with basic type- and value- functions dealing with atomic structure and ions

# Electron Configuration # 

Electron configuration and quantum numbers. Independent of `Elements,` all calculations done by atomic number 

Type-level computation is reasonably effective up to sodium. 

# Slaters # 

Calculation of screening constant, effective nuclear charge, covalent radius (from atomic radius), and electronegativity using Slater's rules. 

Has two child modules which apply Elements on the type- and value- level 

## Slaters.Safe ## 

    "Safe" computation, reification from type-level values. Enables compile-time type errors. Likely to be removed due to negligible benefit. 

## Slaters.Unsafe ##

    "Unsafe" computation, value-level. Additionally includes computation for ionic radii. 

# Compounds # 

Molecules, Ionic compounds, and associated functions. Currently only accepts a pre-prepared Lewis-like structure for construction. 

# Bonds # 

Bond types, experimental work on bond energy calculation for thermodynamics 

# Reactions # 

Framework for molecular reactions at present 
