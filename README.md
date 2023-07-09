# DisorderedPointClusters.jl

We consider point configurations in the 3-torus $T^3\subset\mathbb{R}^3$. Placing $2n$ parallel planes in the torus, we separate the torus into $2n$ disjoint regions. On each plane, $m$ points $p$ are placed. At these points, we intend to place catenoidal necks (with different software) to connect alternating regions and in doing so create 2 separate, connected labyrinths, as is the case in 3-periodic bicontinuous minimal surfaces. 

As energy functional, we choose the Riesz energy, which is monotonic and purely dependend on the particles' diestance: $q_{ij} = \frac{1}{\lvert\lvert p_i - p_j \lvert\lvert}$ for $\textit{comparable}$ points $p_i$ and $p_j$. The question then becomes, which points we want to compare. There are three different options:

1. `3D`: The energy functional is considered between all points.
2. `2D`: All points are projected into the same plane by setting $z=0$. Then, we consider distances between all points in that plane, forgetting about the third dimension.
3. `2D-3`: Since our ultimate goal is to place catenoidal necks where the points lie, for any point only adjacent layers should play a role in determining the energy functional. To model this, for each point a copy is placed in the layer below and above, shadowing its behavior. Then, interactions occuring between points inside a layer are considered.

In each of these cases, the last step of the optimization procedure is to sufficiently separate the points so the catenoidal necks don't intersect. This is done by a heavy penalization when points are too close in the Manhattan distance.

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
