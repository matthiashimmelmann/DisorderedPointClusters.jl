module DisorderedPointClusters
#TODO Circular Necks instead of squares

import GLMakie: scatter!, Axis3, Figure, Point3f0, Point2f0, Figure, cam3d!, Axis, xlims!, ylims!, save
import LinearAlgebra: norm, det, pinv, svd
import HomotopyContinuation: Variable, Expression, evaluate, differentiate
import Combinatorics: multiexponents

export generateGridLayers

#=
Creates a square of size `NeckSize` based on the Manhattan-distance around the center point `neckCenter`
in the grid defined by `NGrid` and `NGrid2`.
=#
function createStar(neckCenter, neckSize, NGrid, NGrid2)
    directions = collect(Set(vcat([[[part[1],part[2]], [-part[1],part[2]], [part[1],-part[2]], [-part[1],-part[2]]] for part in vcat([collect(multiexponents(2, exp)) for exp in 0:neckSize]...)]...)))
    #directions = collect(Set(vcat([[[i,j],[i,-j],[-i,j],[-i,-j]] for i in 0:neckSize for j in 0:neckSize]...)))
    star = [[mod(((dir .+ neckCenter)[1]-1),NGrid)+1, mod(((dir .+ neckCenter)[2]-1),NGrid2)+1] for dir in directions]
    return star
end

#=
In the arrray `layerColors` the calculated positions `neckConfig` of the catenoidal necks are coloured
according to their size `neckSize`. 
=#
function color_Neck_Positions(NeckArray, neckConfig, xvarz, xs, totalGrid, layerColours, neckSize; MD_Method)   
    xs_eval = [[] for _ in 1:size(layerColours)[1]]
    for layer in 1:size(NeckArray)[1]
        xs_eval[layer] = [vcat(evaluate(xs[layer][pos][1:2], xvarz=>neckConfig), (layer)/size(layerColours)[1]) for pos in 1:NeckArray[layer]]
        for pos in 1:length(xs_eval[layer])
            neckCenter = filter(t->isapprox(xs_eval[layer][pos][1:2], totalGrid[layer,t[1],t[2],1:2]; atol=1e-3), [(i,j) for i in 1:size(layerColours)[2] for j in 1:size(layerColours)[3]])[1]
            indices = createStar(neckCenter, neckSize, size(layerColours)[2], size(layerColours)[3])
            for index in indices
                layerColours[layer, index[1], index[2]] = 1 .- layerColours[layer, index[1], index[2]]
            end
        end
    end
    return(layerColours)
end

#=
`createListOfRelevantDistances` generates a vector of lists of torus distances, 
i.e. for two molecules all 9 possible distance pairs, in the 2-torus,
for all molecules that match what `MD_Method` dictates. 
No duplicates are calculated because of the if-statements.
=#
function createListOfRelevantDistances(xs, NGrid, NGrid2; MD_Method, distance="euclidean")
    list_of_relevant_molecule_interactions = []
    for layer1 in 1:length(xs), layer2 in layer1:length(xs), pos_torus1 in 1:length(xs[layer1]), pos_torus2 in 1:length(xs[layer2])
        if MD_Method == "2D-3" && layer1==layer2 && pos_torus1<pos_torus2 && all(t->(layer1!=t[1]||pos_torus1!=t[4]||t[3]!=pos_torus2) && (layer1!=t[1]||pos_torus1!=t[3]||t[4]!=pos_torus2), list_of_relevant_molecule_interactions)
            push!(list_of_relevant_molecule_interactions, (layer1, layer1, pos_torus1, pos_torus2))
        elseif (MD_Method == "2D" || MD_Method == "3D") && (layer1!=layer2 || pos_torus1<pos_torus2)
            push!(list_of_relevant_molecule_interactions, (layer1, layer2, pos_torus1, pos_torus2))
        end
    end
    
    multiPeriods = MD_Method == "3D" ? [[mod(i-1,3)-1, NGrid2/NGrid*(mod(div(i-1,3),3)-1),0] for i in 1:9] : [[mod(i-1,3)-1, NGrid2/NGrid*(mod(div(i-1,3),3)-1)] for i in 1:9]
    if distance=="euclidean"
        return [[sum((period + xs[layer1][pos1] - xs[layer2][pos2]).^2) for period in multiPeriods] for (layer1, layer2, pos1, pos2) in list_of_relevant_molecule_interactions]
    elseif distance=="sameplane"
        return [[sum((period + xs[layer1][pos1] - xs[layer2][pos2]).^2) for period in multiPeriods] for (layer1, layer2, pos1, pos2) in filter(t->t[1]==t[2],list_of_relevant_molecule_interactions)]
    elseif distance=="manhattan"
        return [[(period + xs[layer1][pos1] - xs[layer2][pos2]) for period in multiPeriods] for (layer1, layer2, pos1, pos2) in list_of_relevant_molecule_interactions]
    else
        throw(error("Distance type needs to be specified!"))
    end
  end

#=
Does what the method's name suggests: Calculate an approximation of the `energyFunction`'s gradient at `point`.
=#
function calculateGradient(energyFunction, point; t = 1e-5)
    gradient = []
    for i in 1:length(point)
        push!(gradient, (energyFunction(point+[i==j ? t : 0 for j in 1:length(point)])-energyFunction(point))/t)
    end
    return gradient
end


#=
The method `createNeckMatrix` generates a matrix with necks given by `NeckArray`.
Depending on the parameter MD_Method, the looks of the matrix may differ.
=#
function createNeckMatrix(NeckArray; MD_Method)
    xs = [[] for _ in NeckArray]
    for layer in 1:size(NeckArray)[1]
        xs[layer] = [[Variable(:x,layer,neck,xy) for xy in 1:2] for neck in 1:NeckArray[layer]]
        if MD_Method == "2D-3"
            xs[layer] = vcat(xs[layer], [[Variable(:x,(layer<length(NeckArray) ? layer+1 : 1),neck,xy) for xy in 1:2] for neck in 1:NeckArray[(layer<length(NeckArray) ? layer+1 : 1)]], [[Variable(:x,(layer>1 ? layer-1 : length(NeckArray)),neck,xy) for xy in 1:2] for neck in 1:NeckArray[(layer>1 ? layer-1 : length(NeckArray))]])
        elseif MD_Method == "3D"
            xs[layer] = [[entry[1],entry[2],layer/length(NeckArray)]  for entry in xs[layer]]
        end
    end
    return xs
end

#=
The method takeMonteCarloStep calculates one random step on the grid 
in direction 0/1/-1 for each entry of the solution vector.
=#
function takeMonteCarloStep(sol, NGrid, NGrid2)
    sol = sol + [rand([-1,0,1])/NGrid for _ in sol]
    scaling = [i%2==1 ? 1 : NGrid2/NGrid for i in 1:length(sol)]
    return sol - scaling .* floor.(sol ./ scaling)
end


#=
`monteCarlo` is the framework for a Monte Carlo Optimization.
It takes the neck point matrix `xs`, the array of initial guesses `initialPoints`,
the list of variables `xvarz` and the rectangular grid sizes `NGrid` and `NGrid2`
as input. It outputs the least energy configuration that was found along random walks
from the starting positions.
=#
function monteCarlo(xs, initialPoints, xvarz, NGrid, NGrid2, NeckSize; MD_Method, maxIter)
  distanceList = createListOfRelevantDistances(xs, NGrid, NGrid2; MD_Method = MD_Method, distance="euclidean")
  #distancesSameplane = createListOfRelevantDistances(xs, NGrid, NGrid2; MD_Method = MD_Method, distance="sameplane")
  manhattanDistance = createListOfRelevantDistances(xs, NGrid, NGrid2; MD_Method = MD_Method, distance="manhattan")
  #PrioritizeSamePlane?
  energyFunction = sol -> sum(1 ./ map(t->minimum(evaluate(t, xvarz=>sol)), distanceList) )  
  energyFunctionManhattan = sol -> sum(1 ./ map(t->minimum(evaluate(t, xvarz=>sol)), distanceList)) + sum(minimum(map(t->sum(abs.(t)), evaluate.(t,xvarz=>sol)))<=(2*NeckSize+2)/NGrid ? 10000/minimum(map(t->sum(abs.(t)), evaluate.(t,xvarz=>sol))) : 0 for t in manhattanDistance)
  #energyFunctionSquare = sol -> energyFunction(sol) + sum(minimum(map(t->maximum(abs.(t)), evaluate.(t,xvarz=>sol)))<=(2*NeckSize+1)/NGrid ? 100/minimum(map(t->maximum(abs.(t)), evaluate.(t,xvarz=>sol))) : 0 for t in manhattanDistance)

  outputList = []
  for initialPoint in initialPoints
    push!(outputList, initialPoint)
    prevsol = Base.copy(initialPoint)
    saveEnergy = energyFunction(prevsol)
    for iter in 1:maxIter
        cursol = takeMonteCarloStep(prevsol, NGrid, NGrid2)

        energy = energyFunction(cursol)
        iter%50==0 && println(iter," ", energy)
        if energy <= saveEnergy
            outputList[end] = cursol
            saveEnergy = energy
            prevsol = cursol
        end

        if iter%1000==0
            gradient = calculateGradient(energyFunction, prevsol)
            avg = sum([norm(grd) for grd in gradient]) / length(gradient)
            gradient = round.(NGrid .* round.(gradient))
            gradient = [grd > avg ? 1 : (grd < -avg ? -1 : 0) for grd in gradient] ./ NGrid
            scaling = [i%2==1 ? 1 : NGrid2/NGrid for i in 1:length(prevsol)]
            cursol = prevsol - gradient
            cursol = cursol - scaling .* floor.(cursol ./ scaling)
            if energyFunction(cursol) <= saveEnergy
                println("HMC improved!")
                outputList[end] = cursol
                saveEnergy = energyFunction(cursol)
                prevsol = cursol
            end
        end
    end

    saveEnergy = energyFunctionManhattan(prevsol)
    for iter in maxIter+1:maxIter+maxIter/10
        cursol = takeMonteCarloStep(prevsol, NGrid, NGrid2)

        energy = energyFunctionManhattan(cursol)
        iter%50==0 && println(iter," ", energy, " (Manhattan distance)")
        if energy <= saveEnergy
            outputList[end] = cursol
            saveEnergy = energy
            prevsol = cursol
        end
    end
    println("Final Energy: ", saveEnergy)
  end
  return argmin(sol->energyFunction(sol), outputList)
end

#=
There are three ways to generate point clouds: "2D-3", "2D" and "3D". 
They can be chosen with the parameter "MD_Method".
The input is composed of the number of grid points in each layer, "NGrid",
the number of Necks in each layer, summarized in the array "NeckArray",
and the size of the catenoidal necks in the output, "NeckSize".
WARNING don't choose NeckSize too large to avoid overlap.
The maxIter parameter denotes how many MD-steps are supposed to be performed.
=#
function generateGridLayers(NGrid::Int, NeckArray::Vector, NeckSize::Int; NGrid2 = NGrid, MD_Method = "2D-3", maxIter = 15000, monteCarloStartPoints = 3)
    length(NeckArray)%2==0 || throw(error("There must be an even number of layers!"))
    MD_Method=="2D-3" || MD_Method=="2D" || MD_Method=="3D" || throw(error("Please choose a permissible MD Method!"))
    (MD_Method!="2D-3" || length(NeckArray)>=4) || throw(error("When using the method 2D-3 you need to use at least 4 layers!"))
    NGrid2 >= NGrid || throw(error("For some reason, NGrid2 needs to be at least as large as `NGrid`. Please abide by the rules of the code god."))
    sum([NeckArray[i] for i in 1:2:length(NeckArray)])==sum([NeckArray[i] for i in 2:2:length(NeckArray)]) || throw(error("Please choose an amount of necks per layer that can lead to a balanced minimal surface!"))

    pointGrid = [[Point3f0(a,b,c) for a in 1/(NGrid*2):1/NGrid:1-1/(NGrid*2) for b in 1/(NGrid*2):1/NGrid:NGrid2/NGrid-1/(NGrid*2)]  for c in 0:1/length(NeckArray):1-1/length(NeckArray)]
    xvarz = [Variable(:x, layer, neck, pos) for layer in 1:length(NeckArray) for neck in 1:NeckArray[layer] for pos in 1:2]
    gridPositions = [[(i,j) for i in 1:NGrid for j in 1:NGrid2] for _ in 1:length(NeckArray)]
    initialPoints = []
    for _ in 1:monteCarloStartPoints
      initialPoint = []
      for layer in 1:length(NeckArray), _ in 1:NeckArray[layer]
        randPos = rand(gridPositions[layer])
        gridPositions[layer] = filter!(t->t!=(randPos), gridPositions[layer])
        gridPositions[layer!=length(NeckArray) ? layer+1 : 1] = filter!(t->t!=(randPos), gridPositions[layer!=length(NeckArray) ? layer+1 : 1])
        gridPositions[layer!=1 ? layer-1 : length(NeckArray)] = filter!(t->t!=(randPos), gridPositions[layer!=1 ? layer-1 : length(NeckArray)])
        append!(initialPoint, Vector{Float64}(pointGrid[layer][(randPos[1]-1)*NGrid+(randPos[2]-1)%NGrid2+1][1:2]))
      end
      push!(initialPoints, initialPoint)
    end
    #TODO points closer than necksize => push apart
    xs = createNeckMatrix(NeckArray; MD_Method = MD_Method)
    neckConfig = monteCarlo(xs, initialPoints, xvarz, NGrid, NGrid2, NeckSize; MD_Method = MD_Method, maxIter = maxIter)

    totalGrid = Array{Any,4}(undef, length(NeckArray), NGrid, NGrid2, 3)
    layerColours = Array{Any,3}(undef, length(NeckArray), NGrid, NGrid2)
    for layer in 1:length(NeckArray), pos1 in 0:NGrid-1, pos2 in 0:NGrid2-1
        totalGrid[layer, pos1+1, pos2+1,:] = [1/(NGrid*2)+pos1/NGrid, 1/(NGrid*2)+pos2/NGrid, (layer-1)/length(NeckArray)]
        layerColours[layer, pos1+1, pos2+1] = mod(layer, 2)
    end
    layerColours = color_Neck_Positions(NeckArray, neckConfig, xvarz, xs, totalGrid, layerColours, NeckSize; MD_Method = MD_Method)
    open("pomelocoordinates.txt", "w") do io
        for k in [(layer,pos1,pos2) for pos1 in 1:NGrid for pos2 in 1:NGrid2 for layer in 1:length(NeckArray)]
            write(io, string(totalGrid[k[1],k[2],k[3],1], " ", totalGrid[k[1],k[2],k[3],2], " ", totalGrid[k[1],k[2],k[3],3], " ", layerColours[k[1],k[2],k[3]] == 1 ? "1" : "2", "\n"))
        end
    end;
    #=
    open("pomelocoordinates2.txt", "w") do io
        for i in 1:length([(layer,pos1,pos2) for pos1 in 1:NGrid for pos2 in 1:NGrid2 for layer in 1:length(NeckArray)])
            k = [(layer,pos1,pos2) for pos1 in 1:NGrid for pos2 in 1:NGrid2 for layer in 1:length(NeckArray)][i]
            write(io, "$(i): ", string(totalGrid[k[1],k[2],k[3],1], " ", totalGrid[k[1],k[2],k[3],2], " ", totalGrid[k[1],k[2],k[3],3], " ", layerColours[k[1],k[2],k[3]] == 1 ? "c(1, 0, 0, 1)" : "c(0, 0, 1, 1)", "\n"))
        end
    end;=#

    scene=Figure(resolution = (650, 650), fontsize=22)
    if MD_Method == "2D"
        ax = Axis(scene[1,1])
        xlims!(ax, (0,1))
        ylims!(ax, (0,1))
        xs_eval = [[evaluate(xs[layer][pos], xvarz=>neckConfig) for pos in 1:NeckArray[layer]] for layer in 1:length(NeckArray)]
        foreach(k -> scatter!(ax, Point2f0(xs_eval[k[1]][k[2]][1:2]), color = mod(k[1],2)==0 ? :red : :blue, grid=true, markersize=17), [(i,j) for i in 1:size(xs_eval)[1] for j in 1:length(xs_eval[i])])
    elseif MD_Method == "2D-3"
        scene=Figure(resolution = (1800, 900), fontsize = 25, ticksize=18)
        colors = repeat([:red, :green, :blue, :purple, :orange, :pink, :saddlebrown, :grey30], 20)
        xs_eval = [[evaluate(xs[layer][pos], xvarz=>neckConfig) for pos in 1:NeckArray[layer]] for layer in 1:length(NeckArray)]
        ax = vcat([Axis(scene[div(layer-1,2)+1, mod(layer-1,2)+1], title = "Layer = $(layer)", titlecolor=colors[layer]) for layer in 1:length(NeckArray)], Axis(scene[:,3:4], title = "All Layers together"))
        for layer in 1:length(NeckArray)
            layer1 = (layer==1) ? length(NeckArray) : layer-1
            layer2 = (layer==length(NeckArray)) ? 1 : layer+1
            xlims!(ax[layer], (0,1))
            ylims!(ax[layer], (0,1))
            foreach(k -> scatter!(ax[layer], Point2f0(xs_eval[layer][k][1:2]), color = :black, markersize=17), 1:NeckArray[layer])
            foreach(k -> scatter!(ax[layer], Point2f0(xs_eval[layer1][k][1:2]), color = :red, markersize=17), 1:NeckArray[layer1])
            foreach(k -> scatter!(ax[layer], Point2f0(xs_eval[layer2][k][1:2]), color = :green, markersize=17), 1:NeckArray[layer2])
        end
        xlims!(ax[end], (0,1))
        ylims!(ax[end], (0,1))
        foreach(k -> scatter!(ax[end], Point2f0(xs_eval[k[1]][k[2]][1:2]), color = colors[k[1]], markersize=17), [(i,j) for i in 1:size(xs_eval)[1] for j in 1:length(xs_eval[i])])
    else
        ax = Axis3(scene[1,1])
        xlims!(ax, (0,1))
        ylims!(ax, (0,1))
        foreach(k -> mod(layerColours[k[1],k[2],k[3]],2) != mod(k[1],2) ? scatter!(ax, Point3f0(totalGrid[k[1],k[2],k[3],:]), color=layerColours[k[1],k[2],k[3]] == 1 ? :red : :blue, grid=true) : nothing, [(layer,pos1,pos2) for pos1 in 1:NGrid for pos2 in 1:NGrid2 for layer in 1:length(NeckArray)])
    end
    save("MolecularConfiguration$(Base.time()).png", scene)
    display(scene)
end

function generateGridLayers(NGrid::Int, NNecks::Int, NLayers::Int, NeckSize::Int; NGrid2 = NGrid, MD_Method = "2D-3", maxIter = 25000, monteCarloStartPoints = 2)
    generateGridLayers(NGrid, [NNecks for _ in 1:NLayers], NeckSize; NGrid2 = NGrid2, MD_Method = MD_Method, maxIter = maxIter, monteCarloStartPoints = monteCarloStartPoints)
end

generateGridLayers(60, 6, 8, 3; MD_Method="2D-3", maxIter = 25000, monteCarloStartPoints = 2)

end