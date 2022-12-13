module DisorderedPointClusters
#Interactions are taken per 3 consecutive 2D-layers

import GLMakie: scatter!, Axis3, Figure, Point3f0, Scene, cam3d!
import LinearAlgebra: norm, det, pinv, svd
import HomotopyContinuation: Variable, Expression, evaluate, differentiate
import Combinatorics: multiexponents
#import Zygote: gradient

export generateGridLayers

function createStar(neckCenter, neckSize, NGrid, NGrid2)
    directions = collect(Set(vcat([[[part[1],part[2]], [-part[1],part[2]], [part[1],-part[2]], [-part[1],-part[2]]] for part in vcat([collect(multiexponents(2, exp)) for exp in 0:neckSize]...)]...)))
    star = [[mod(((dir .+ neckCenter)[1]-1),NGrid)+1, mod(((dir .+ neckCenter)[2]-1),NGrid2)+1] for dir in directions]
    return star
end

function color_Neck_Positions(neckConfig, xvarz, xs, totalGrid, layerColours, neckSize; MD_Method)
    xs_eval = MD_Method=="2D-3" ? Array{Float64,3}(undef, size(layerColours)[1], Int(size(xs)[2]/3), 3) : Array{Float64,3}(undef, size(layerColours)[1], size(xs)[2], 3)
    xs_eval[:,:,1:2] = MD_Method=="2D-3" ? evaluate.(xs, xvarz=>neckConfig)[:,1:Int(size(xs)[2]/3),1:2] : evaluate.(xs[:,:,1:2], xvarz=>neckConfig)
    for layer in 1:size(layerColours)[1]
        xs_eval[layer,:,3] .= (layer)/size(layerColours)[1]
        for pos in 1:size(xs_eval)[2]
            neckCenter = filter(t->isapprox(xs_eval[layer,pos,1:2], totalGrid[layer,t[1],t[2],1:2]; atol=1e-3), [(i,j) for i in 1:size(layerColours)[2] for j in 1:size(layerColours)[3]])[1]
            indices = createStar(neckCenter, neckSize, size(layerColours)[2], size(layerColours)[3])
            for index in indices
                layerColours[layer, index[1], index[2]] = 1 .- layerColours[layer, index[1], index[2]]
            end
        end
    end
    return(layerColours)
end

function create_list_of_relevant_distances(xs, NGrid, NGrid2; MD_Method)
  list_of_relevant_molecule_interactions = []
  for layer1 in 1:size(xs)[1], layer2 in layer1:size(xs)[1], pos_torus1 in 1:size(xs)[2], pos_torus2 in 1:size(xs)[2]
      if MD_Method == "2D-3" && pos_torus1<pos_torus2 && all(t->(layer1!=t[1]||pos_torus1!=t[4]||t[3]!=pos_torus2) && (layer1!=t[1]||pos_torus1!=t[3]||t[4]!=pos_torus2), list_of_relevant_molecule_interactions)
          push!(list_of_relevant_molecule_interactions, (layer1, layer1, pos_torus1, pos_torus2))
      elseif (MD_Method == "2D" || MD_Method == "3D") && (layer1!=layer2 || pos_torus1<pos_torus2)
          push!(list_of_relevant_molecule_interactions, (layer1, layer2, pos_torus1, pos_torus2))
      end
  end

  multiPeriods = MD_Method == "3D" ? [[mod(i-1,3)-1, NGrid2/NGrid*(mod(div(i-1,3),3)-1),0] for i in 1:9] : [[mod(i-1,3)-1, NGrid2/NGrid*(mod(div(i-1,3),3)-1)] for i in 1:9]
  return [[sum((period + xs[layer1,pos1,:] - xs[layer2,pos2,:]).^2) for period in multiPeriods] for (layer1, layer2, pos1, pos2) in list_of_relevant_molecule_interactions]
end

function create_Periods(NLayers, NNecks; MD_Method)
    xs = MD_Method == "2D-3" ? Array{Expression,3}(undef, NLayers, NNecks*3, 2) : 
           Array{Expression,3}(undef, NLayers, NNecks, MD_Method == "2D" ? 2 : 3)

    for layer in 1:NLayers
        xs[layer,1:NNecks,1:2] = [Variable(:x,layer,neck,xy) for xy in 1:2 for neck in 1:NNecks]
        if MD_Method == "2D-3"
          xs[(layer<NLayers ? layer+1 : 1),NNecks+1:2*NNecks,1:2] = [Variable(:x,layer,neck,xy) for xy in 1:2 for neck in 1:NNecks]
          xs[(layer>1 ? layer-1 : NLayers),2*NNecks+1:3*NNecks,1:2] = [Variable(:x,layer,neck,xy) for xy in 1:2 for neck in 1:NNecks]
        elseif MD_Method == "3D"
          xs[layer,:,3] .= layer/NLayers
        end
    end
    #=
    periodic_xs = MD_Method == "3D" ? Array{Expression,3}(undef, NLayers, size(xs)[2]*9, 3) : Array{Expression,3}(undef, NLayers, size(xs)[2]*9, 2)
    #NOTE left to right 1-2-3, bottom to top 1-4-7
    for torus in 1:9
        periodic_xs[:, (torus-1)*size(xs)[2]+1:torus*size(xs)[2], :] = xs
        periodic_xs[:, (torus-1)*size(xs)[2]+1:torus*size(xs)[2], 1] .+= (torus-1)%3-1
        periodic_xs[:, (torus-1)*size(xs)[2]+1:torus*size(xs)[2], 2] .+= Int(floor((torus-1)/3))-1
    end
    =#

    return xs#, periodic_xs
end

function takeMonteCarloStep(sol, NGrid, NGrid2)
  sol = sol + [rand([-1,0,1])/NGrid for _ in sol]
  scaling = [i%2==1 ? 1 : NGrid2/NGrid for i in 1:length(sol)]
  return sol - scaling .* floor.(sol ./ scaling)
end

function monteCarlo(xs, initialPoints, xvarz, NGrid, NGrid2; MD_Method, maxIter)
  distanceList = create_list_of_relevant_distances(xs, NGrid, NGrid2; MD_Method = MD_Method)
  energyfunction = sol->sum(1 ./ (map(t->minimum(evaluate(t, xvarz=>sol)), distanceList).^2))
  outputList = []
  for initialPoint in initialPoints
    push!(outputList, initialPoint)
    prevsol = Base.copy(initialPoint)
    saveEnergy = energyfunction(prevsol)
    for iter in 1:maxIter
      cursol = takeMonteCarloStep(prevsol, NGrid, NGrid2)
      energy = energyfunction(cursol)
      iter%100==0 && println(iter," ", energy)
      if energy < saveEnergy
        outputList[end] = cursol
        saveEnergy = energy
        prevsol = cursol
      end
    end
  end

  return argmin(sol->energyfunction(sol), outputList)
end


#=
function gradientDescent(xs, initialPoint, variables, NGrid; energy, MD_Method, maxIter)
    cursol = Base.copy(initialPoint)

    distances = create_list_of_relevant_distances(xs; MD_Method = MD_Method)
    Q=0;
    if energy == "riesz"
        Q = sum([1/d^2 for d in distances])
    elseif energy == "lennardjones"
        σ1 = MD_Method=="2D" ? sqrt(2)/(sqrt(size(xs)[1]*size(xs)[2]/9)) : sqrt(2)/(sqrt(size(xs)[2]/9))
        Q = sum(((σ1/d)^6-(σ1/d)^3 for d in distances))
    end
    ∇Q = differentiate(Q, variables)
    HessQ = differentiate(∇Q, variables)
    isNoMin = 0
    α = 0.5
    prevsol = cursol
    solutionarray = []

    for iter in 1:maxIter
        println(iter, "\tMinima found: ", isNoMin,"\t stepsize: ", round(α, digits=3), "\t\t gradient norm: ", round(norm(evaluate(∇Q, variables=>cursol)), digits=1), "\t\t energy: ", round(evaluate(Q, variables=>cursol), digits=1))
        stepdirection = pinv(evaluate.(HessQ, variables=>cursol))*evaluate(∇Q, variables=>cursol)
        cursol = prevsol - α*stepdirection
        α = iter%80 == 0 ? 0.5 : α
        cursol = iter%530 == 0 ? mod.(Int.(round.(NGrid.*cursol)), NGrid) ./ NGrid .+ 1/(2*NGrid) : cursol
        cursol = iter%190 == 0 ? cursol + 1e-2*(rand(Float64, length(cursol)) .- 0.5) : cursol

        cursol = cursol - floor.(cursol)

        if norm(evaluate(∇Q, variables => cursol)) > norm(evaluate(∇Q, variables => prevsol))
            α = α/3
        else
            α = 1.1*α
            α>0.7 ? α = 0.7 : nothing
        end


        if norm(evaluate(∇Q, variables=>cursol)) <= 1e-2
            isNoMin += 1;
            push!(solutionarray, cursol)
            cursol = cursol - 2e-1*(rand(Float64, length(cursol)) .- 0.5)
            cursol = cursol - floor.(cursol)
        end

        prevsol = cursol
        iter%5==0 && push!(solutionarray, cursol)
    end
    cursol = solutionarray[argmin([evaluate(Q, variables=>sol) for sol in solutionarray])]
    #Round only in the end to grid
    cursol = mod.(Int.(round.(NGrid.*cursol)), NGrid) ./ NGrid .+ 1/(2*NGrid)
    cursol = cursol - floor.(cursol)

    return cursol
end
=#

#=
There are three ways to generate point clouds: "2D-3", "2D" and "3D". 
They can be chosen with the parameter "MD_Method".
The input is composed of the number of grid points in each layer, "NGrid",
the number of Necks in each layer, "NNecks", the number of layers, "NLayers",
and the size of the catenoidal necks in the output, "NeckSize".
WARNING don't choose NeckSize too large to avoid overlap.
The maxIter parameter denotes how many MD-steps are supposed to be performed.
=#
function generateGridLayers(NGrid::Int, NNecks::Int, NLayers::Int, NeckSize::Int; NGrid2 = NGrid, MD_Method = "2D-3", maxIter = 25000, monteCarloStartPoints = 7)
    NLayers%2==0 || throw(error("There must be an even number of layers!"))
    MD_Method=="2D-3" || MD_Method=="2D" || MD_Method=="3D" || throw(error("Please choose a permissible MD Method!"))
    (MD_Method!="2D-3" || NLayers>=4) || throw(error("When using the method 2D-3 you need to use at least 4 layers!"))

    pointGrid = [[Point3f0(a,b,c) for a in 1/(NGrid*2):1/NGrid:1-1/(NGrid*2) for b in 1/(NGrid*2):1/NGrid:NGrid2/NGrid-1/(NGrid*2)]  for c in 0:1/NLayers:1-1/NLayers]
    xvarz = [Variable(:x, layer, neck, pos) for layer in 1:NLayers for neck in 1:NNecks for pos in 1:2]
    gridPositions = [[(i,j) for i in 1:NGrid for j in 1:NGrid2] for _ in 1:NLayers]
    initialPoints = []
    for _ in 1:monteCarloStartPoints
      initialPoint = []
      for layer in 1:NLayers, _ in 1:NNecks
        randPos = rand(gridPositions[layer])
        gridPositions[layer] = filter!(t->t!=(randPos), gridPositions[layer])
        gridPositions[layer!=NLayers ? layer+1 : 1] = filter!(t->t!=(randPos), gridPositions[layer!=NLayers ? layer+1 : 1])
        gridPositions[layer!=1 ? layer-1 : NLayers] = filter!(t->t!=(randPos), gridPositions[layer!=1 ? layer-1 : NLayers])
        append!(initialPoint, Vector{Float64}(pointGrid[layer][(randPos[1]-1)*NGrid+(randPos[2]-1)%NGrid2+1][1:2]))
      end
      push!(initialPoints, initialPoint)
    end

    xs#=, periodic_xs=# = create_Periods(NLayers, NNecks; MD_Method = MD_Method)
    #neckConfig = gradientDescent(periodic_xs, initialPoint, xvarz, NGrid; energy = energy, MD_Method = MD_Method, maxIter = maxIter)
    neckConfig = monteCarlo(xs, initialPoints, xvarz, NGrid, NGrid2; MD_Method = MD_Method, maxIter = maxIter)

    totalGrid = Array{Any,4}(undef, NLayers, NGrid, NGrid2, 3)
    layerColours = Array{Any,3}(undef, NLayers, NGrid, NGrid2)
    for layer in 1:NLayers, pos1 in 0:NGrid-1, pos2 in 0:NGrid2-1
        totalGrid[layer, pos1+1, pos2+1,:] = [1/(NGrid*2)+pos1/NGrid, 1/(NGrid*2)+pos2/NGrid, layer/NLayers]
        layerColours[layer, pos1+1, pos2+1] = mod(layer, 2)
    end

    scene2 = Scene()
    cam3d!(scene2)
    layerColours = color_Neck_Positions(neckConfig, xvarz, xs, totalGrid, layerColours, NeckSize; MD_Method = MD_Method)
    open("pomelocoordinates.txt", "w") do io
        for k in [(layer,pos1,pos2) for pos1 in 1:NGrid for pos2 in 1:NGrid2 for layer in 1:NLayers]
            write(io, string(totalGrid[k[1],k[2],k[3],1], " ", totalGrid[k[1],k[2],k[3],2], " ", totalGrid[k[1],k[2],k[3],3], " ", layerColours[k[1],k[2],k[3]] == 1 ? "1" : "2", "\n"))
        end
    end;
    foreach(k -> scatter!(scene2, Point3f0(totalGrid[k[1],k[2],k[3],:]), color=layerColours[k[1],k[2],k[3]] == 1 ? :red : :blue), [(layer,pos1,pos2) for pos1 in 1:NGrid for pos2 in 1:NGrid2 for layer in 1:NLayers ])
    display(scene2)
end

end
