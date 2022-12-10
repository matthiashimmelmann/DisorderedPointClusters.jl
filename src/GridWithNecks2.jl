module GridWithNecks2
#Interactions are only considered in 3D with vertices pinned to planes

import GLMakie: scatter!, Axis3, Figure, Point3f0, Scene, cam3d!
import LinearAlgebra: norm, det, pinv, svd
import HomotopyContinuation: Variable, Expression, evaluate, differentiate
import Random: shuffle
import Combinatorics: multiexponents

export generateGridLayers

function createStar(neckCenter, neckSize, NGrid)
    directions = collect(Set(vcat([[[part[1],part[2]], [-part[1],part[2]], [part[1],-part[2]], [-part[1],-part[2]]] for part in vcat([collect(multiexponents(2, exp)) for exp in 0:neckSize]...)]...)))
    star = [[mod(((dir .+ neckCenter)[1]-1),NGrid)+1, mod(((dir .+ neckCenter)[2]-1),NGrid)+1] for dir in directions]
    return star
end

function colorNeckPositions(neckConfig, xvarz, xs, totalGrid, layerColours, neckSize)
    xs_eval = Array{Any,3}(undef, size(layerColours)[1], size(xs)[2], 3)
    xs_eval[:,:,1:3] = evaluate.(xs, xvarz=>neckConfig)
    for layer in 1:size(layerColours)[1]
        xs_eval[layer,:,3] .= layer/size(layerColours)[1]
        for pos in 1:size(xs_eval)[2] 
            neckCenter = filter(t->isapprox(xs_eval[layer,pos,1:2], totalGrid[layer,t[1],t[2],1:2]; atol=1e-3), [(i,j) for i in 1:size(layerColours)[2] for j in 1:size(layerColours)[3]])[1]
            display(xs_eval[layer, pos, :])
            display(totalGrid[layer, neckCenter[1], neckCenter[2], :])
            indices = createStar(neckCenter, neckSize, size(layerColours)[2])
            for index in indices
                layerColours[layer, index[1], index[2]] = 1 .- layerColours[layer, index[1], index[2]]
            end
        end
    end
    return(layerColours)
end


function backToTorus!(cursol)
    if any(t->t<0||t>1,cursol)
        #display("Left the Torus")
        cursol = cursol - floor.(cursol)
    end
    return cursol
end

function gradientDescent(periodic_xs, initialPoint, variables, NGrid; energy)
    cursol = Base.copy(initialPoint)
    list_of_relevant_molecule_interactions = []

    for layer1 in 1:size(periodic_xs)[1], layer2 in 1:size(periodic_xs)[1], pos_torus in Int(4*size(periodic_xs)[2]/9)+1:Int(5*size(periodic_xs)[2]/9), pos_general in 1:size(periodic_xs)[2]
      if (pos_torus!=pos_general || layer1!=layer2) && all(t->t[3]!=pos_general||t[4]!=pos_torus||t[1]!=layer2||t[2]!=layer1, list_of_relevant_molecule_interactions)
        push!(list_of_relevant_molecule_interactions, (layer1, layer2, pos_torus, pos_general))
      end
    end

    distances = [sum((periodic_xs[layer1,pos1,:] - periodic_xs[layer2,pos2,:]).^2) for (layer1, layer2, pos1, pos2) in list_of_relevant_molecule_interactions]
    Q=0;
    if energy == "riesz"
        Q = sum([1/d^2 for d in distances])
    elseif energy == "lennardjones"
        σ1 = 1/((sqrt(size(periodic_xs)[2]/9)))
        Q = sum(((σ1/d)^6-(σ1/d)^3 for d in distances))
    end
    ∇Q = differentiate(Q, variables)
    HessQ = differentiate(∇Q, variables)
    isNoMin = 0

    for iter in 1:1000000
        println(iter, " ", isNoMin," ",norm(evaluate(∇Q, variables=>cursol)))
        cursol = cursol - 0.4*pinv(evaluate(HessQ, variables=>cursol))*evaluate(∇Q, variables=>cursol)
        cursol = iter%250 == 0 ? backToTorus!(mod.(Int.(round.(NGrid.*cursol)), NGrid) ./ NGrid .+ 1/(2*NGrid)) : cursol
        cursol = iter%75 == 0 ? cursol - 5e-2*rand(Float64, length(cursol)) : cursol
        backToTorus!(cursol)

        if norm(evaluate(∇Q, variables=>cursol))<=1e-6 && (all(t-> t>1e-8, svd(evaluate(HessQ, variables=>cursol)).S))
            display("Minimum Found!")
            break;
        elseif norm(evaluate(∇Q, variables=>cursol))<=1e-6
            isNoMin += 1;
            if isNoMin == 10
                display("Found no Minimum for 3 times.")
                break;
            end
            cursol = cursol - 5e-2*rand(Float64, length(cursol))
            backToTorus!(cursol)
        end
    end

    #Round only in the end to grid
    return backToTorus!(mod.(Int.(round.(NGrid.*cursol)), NGrid) ./ NGrid .+ 1/(2*NGrid))
end

#=NGrid is the edge length of the basic square grid
NNecks is the number of necks placed per layer
NLayers is the number of layers (needs to be even)=#
function generateGridLayers(NGrid, NNecks, NLayers, neckSize)
    NLayers%2==0 || throw(error("There must be an even number of layers!"))

    pointGrid = [[Point3f0(a,b,c) for a in 1/(NGrid*2):1/NGrid:1-1/(NGrid*2) for b in 1/(NGrid*2):1/NGrid:1-1/(NGrid*2)]  for c in 0:1/NLayers:1-1/NLayers]
    xvarz = [Variable(:x, layer, neck, pos) for layer in 1:NLayers for neck in 1:NNecks for pos in 1:2]
    gridPositions = [[(i,j) for i in 1:NGrid for j in 1:NGrid] for _ in 1:NLayers]
    initialPoint = []
    for layer in 1:NLayers
        for _ in 1:NNecks
            randPos = rand(gridPositions[layer])
            gridPositions[layer] = filter!(t->t!=(randPos), gridPositions[layer])
            gridPositions[layer!=NLayers ? layer+1 : 1] = filter!(t->t!=(randPos), gridPositions[layer!=NLayers ? layer+1 : 1])
            gridPositions[layer!=1 ? layer-1 : NLayers] = filter!(t->t!=(randPos), gridPositions[layer!=1 ? layer-1 : NLayers])
            append!(initialPoint, Vector{Float64}(pointGrid[layer][(randPos[1]-1)*NGrid+(randPos[2]-1)%NGrid+1][1:2]))
        end
    end

    xs = Array{Expression,3}(undef, NLayers, NNecks, 3)
    for layer in 1:NLayers
        xs[layer,1:NNecks,1:2] = [Variable(:x,layer,neck,xy) for xy in 1:2 for neck in 1:NNecks]
        xs[layer,:,3] .= layer/NLayers
    end
    
    periodic_xs = Array{Expression,3}(undef, NLayers, size(xs)[2]*9, 3)
    #NOTE left to right 1-2-3, bottom to top 1-4-7
    for torus in 1:9
        periodic_xs[:, (torus-1)*size(xs)[2]+1:torus*size(xs)[2], :] = xs
        periodic_xs[:, (torus-1)*size(xs)[2]+1:torus*size(xs)[2], 1] .+= (torus-1)%3-1
        periodic_xs[:, (torus-1)*size(xs)[2]+1:torus*size(xs)[2], 2] .+= Int(floor((torus-1)/3))-1
    end
    
    neckConfig = gradientDescent(periodic_xs, initialPoint, xvarz, NGrid; energy = "riesz")#"lennardjones"

    totalGrid = Array{Any,4}(undef, NLayers, NGrid, NGrid, 3)
    layerColours = Array{Any,3}(undef, NLayers, NGrid, NGrid)

    for layer in 1:NLayers, pos1 in 0:NGrid-1, pos2 in 0:NGrid-1
        totalGrid[layer, pos1+1, pos2+1,:] = [1/(NGrid*2)+pos1/NGrid, 1/(NGrid*2)+pos2/NGrid, layer/NLayers]
        layerColours[layer, pos1+1, pos2+1] = mod(layer, 2)
    end

    scene2 = Scene()
    cam3d!(scene2)
    layerColours = colorNeckPositions(neckConfig, xvarz, xs, totalGrid, layerColours, neckSize)
    open("pomelocoordinates.txt", "w") do io
        for k in [(layer,pos1,pos2) for pos1 in 1:NGrid for pos2 in 1:NGrid for layer in 1:NLayers]
            write(io, string(totalGrid[k[1],k[2],k[3],1], " ", totalGrid[k[1],k[2],k[3],2], " ", totalGrid[k[1],k[2],k[3],3], " ", layerColours[k[1],k[2],k[3]] == 1 ? "1" : "2", "\n"))
        end
    end;
    foreach(k -> scatter!(scene2, Point3f0(totalGrid[k[1],k[2],k[3],:]), color=layerColours[k[1],k[2],k[3]] == 1 ? :red : :blue), [(layer,pos1,pos2) for pos1 in 1:NGrid for pos2 in 1:NGrid for layer in 1:NLayers ])
    display(scene2)
end

generateGridLayers(50, 4, 4, 4)

end