# DisorderedPointClusters.jl

We consider point configurations in the 3-torus <img src="https://render.githubusercontent.com/render/math?math=T^3\subset\mathbb{R}^3">. 

## Installation

```
julia> ]
(@v1.8) pkg> add DisorderedPointClusters
```

## Usage

There are two main methods in this package: `plot_implicit_surface` and `plot_implicit_curve`. Let us first consider an example of the former:

```
julia> generateGridLayers(50,7, 6, 3; MD_Method="2D-3", maxIter = 35000, monteCarloStartPoints = 2)
```

The result of this can be seen in the following image: 
<p align="center">
  <img src="https://github.com/matthiashimmelmann/DisorderedPointClusters.jl/blob/main/pictures/MolecularConfiguration1.686438587406e9.png" width="900", height="450">
</p>
