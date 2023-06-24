# DisorderedPointClusters.jl

We consider point configurations in the 3-torus <img src="https://render.githubusercontent.com/render/math?math=T^3\subset\mathbb{R}^3">. 

## Installation

```
julia> ]
(@v1.8) pkg> add Implicit3DPlotting
```

## Usage

There are two main methods in this package: `plot_implicit_surface` and `plot_implicit_curve`. Let us first consider an example of the former:

```
julia> generateGridLayers(50,7, 6, 3; MD_Method="2D-3", maxIter = 35000, monteCarloStartPoints = 2)
```

The result of this can be seen in the following image: 
<p align="center">
  <img src="https://user-images.githubusercontent.com/65544132/114864346-2b0ec700-9df1-11eb-8ad4-4ef2d4e1c9f3.png" width="600", height="600">
</p>
