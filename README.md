# DisorderedPointClusters.jl

We consider point configurations in the 3-torus $T^3\subset\mathbb{R}^3$. Placing $2n$ parallel planes in the torus, we separate the torus into $2n$ disjoint regions. On each plane, $p$ points are placed. At these points, we intend to place catenoidal necks to connect alternating regions and in doing so create 2 separate, connected labyrinths, as is the case in 3-periodic bicontinuous minimal surfaces. 

## Installation

```
julia> ]
(@v1.8) pkg> add DisorderedPointClusters
```

## Usage

The method `generateGridLayers` takes as input the fineness of the grid (50), the amount of points that are placed per layer (7), the amount of parallel planes (6) and the radius of the catenoidal necks (3). In addition, there are a few options that may be changed, such as `maxIter` determining the amount of Monte Carlo steps, `monteCarloStartPoints` determining the amount of repetitions, `MD_Method` which chooses between the three discussed options of energy associated to the catenoids or `NGrid2` allowing for non-square grids.

```
julia> generateGridLayers(50, 7, 4, 3; MD_Method="2D-3", maxIter = 35000, monteCarloStartPoints = 2)
```

This input results in the following image: 
<p align="center">
  <img src="https://github.com/matthiashimmelmann/DisorderedPointClusters.jl/blob/main/pictures/MolecularConfiguration1.678112526442e9.png" width="1200", height="500">
</p>
