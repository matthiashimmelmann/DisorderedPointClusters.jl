module DisorderedPointClusters

import GLMakie: scatter!, Axis3, Figure, Point3f
import LinearAlgebra: norm, det
import HomotopyContinuation: Variable, Expression, evaluate, differentiate

export generateDisorderedPointClusters
    function insertNecks(layers, curLayer, NGrid, Necks)
      pos = rand(1:NGrid^2)
      layerAbove = layers[curLayer==length(layers) ? 1 : curLayer+1][pos] == layers[curLayer][pos]
      layerBelow = layers[curLayer==1 ? length(layers) : curLayer-1][pos] == layers[curLayer][pos]
      #TODO at the moment only orthogonal adjacency.
      anyNeckAdjacent = any(t-> t==(curLayer%2==1 ? 2 : 1), layers[curLayer][vcat(pos, pos%NGrid==0 ? Int((pos/NGrid-1)*NGrid)+1 : pos+1, (pos-1)%NGrid==0 ? Int(((pos-1)/NGrid+1)*NGrid) : (pos-1), pos-NGrid<=0 ? pos-NGrid+NGrid^2 : pos-NGrid, pos+NGrid>NGrid^2 ? pos+NGrid-NGrid^2 : pos+NGrid)])
      if anyNeckAdjacent || layerAbove || layerBelow
        println("Insertion did not work, trying again...")
        insertNecks(layers, curLayer, NGrid, Necks)
      else
        layers[curLayer][pos] = 3-layers[curLayer][pos]
        push!(Necks, [curLayer,pos])
      end
    end

    #=
    function generateMatrix(NLayers, NGrid, NNecks, NeckCables)
      xs=Array{Expression,3}(undef, NLayers, NGrid^2, 3)
      for (i,j) in [(a,b) for a in 1:NLayers for b in 1:NGrid^2]
        xs[i,j,1] = Variable(:x, i, j, 1)
        xs[i,j,2] = Variable(:x, i, j, 2)
        xs[i,:,3] .= (0:1/NLayers:1-1/NLayers)[i]
      end

      for cable in NeckCables
        xs[cable[1]==NLayers ? 1 : cable[1]+1, cable[2], 1:2] = xs[cable[1], cable[2], 1:2]
      end
      return(xs)
    end=#

    function generateMatrix(NLayers, NGrid, NNecks)
      xs=Array{Expression,3}(undef, NLayers, NGrid^2+2*NNecks, 3)
      for (i,j) in [(a,b) for a in 1:NLayers for b in 1:NGrid^2+NNecks]
        xs[i,j,1] = Variable(:x, i, j, 1)
        xs[i,j,2] = Variable(:x, i, j, 2)
        xs[i,:,3] .= (0:1/NLayers:1-1/NLayers)[i]
      end

      for layer in 1:NLayers, i in NGrid^2+NNecks+1:NGrid^2+2*NNecks
        xs[layer==NLayers ? 1 : layer+1, i, 1:2] = xs[layer, i-NNecks, 1:2]
      end
      return(xs)
    end

    function gradientDescent(initialPoint, periodic_xs, xvarz, σ1)
      cursol = Base.copy(initialPoint)
      global ∇Q, HessQ = energyFunction(periodic_xs, xvarz, σ1, cursol)
      while norm(∇Q)>1e-5
        cursol = det(HessQ)==0 ? cursol - 1e-6*∇Q : cursol - HessQ\∇Q
        cursol = cursol - floor.(cursol)
        global ∇Q, HessQ = energyFunction(periodic_xs, xvarz, σ1, cursol)
        display(norm(∇Q))
      end
      return(cursol)
    end

    function energyFunction(periodic_xs, xvarz, σ1, pt)
      eval_periodic_xs = evaluate.(periodic_xs, xvarz=>pt)
      list_of_relevant_molecule_interactions = []
      for layer in 1:size(eval_periodic_xs)[1], i in Int(4*size(periodic_xs)[2]/9)+1:Int(5*size(periodic_xs)[2]/9), j in 1:size(periodic_xs)[2]
        if sqrt(sum((eval_periodic_xs[layer,i,:]-eval_periodic_xs[layer,j,:]).^2)) < 2^(1/6)*sqrt(σ1) && i!=j && all(t->t[1]!=layer||t[2]!=j||t[3]!=i, list_of_relevant_molecule_interactions)
          push!(list_of_relevant_molecule_interactions, (layer,i,j))
        end
      end
      distances = [sum((periodic_xs[layer,i,:] - periodic_xs[layer,j,:]).^2) for (layer,i,j) in list_of_relevant_molecule_interactions]
      lennardJones = 4*sum((σ1/d)^6-(σ1/d)^3 for d in distances) #TODO ADD DIFFERENT NUMBER FOR NECKS -> REPULSION
      ∇Q = differentiate(lennardJones, xvarz)
      HessQ = differentiate(∇Q, xvarz)
      return evaluate(∇Q, xvarz=>pt), evaluate(HessQ, xvarz=>pt)
    end

    function simulateRepulsion(unitGrid, randomGrid, NLayers, NGrid, NNecks)
      initialPoint = Vector{Float64}([])
      xvarz = Vector{Variable}([])
      for i in 1:NLayers, j in 1:NGrid^2+NNecks
        append!(xvarz, [Variable(:x, i, j, 1), Variable(:x, i, j, 2)])
        j<=NGrid^2 ? append!(initialPoint, unitGrid[(i-1)*NGrid^2+j][1:2]) : append!(initialPoint, randomGrid[(i-1)*NNecks+j-NGrid^2][1:2])
      end


      xs = generateMatrix(NLayers, NGrid, NNecks)
      display(evaluate.(xs,xvarz=>initialPoint))
      periodic_xs = Array{Expression,3}(undef, size(xs)[1], size(xs)[2]*9, 3)
      #left to right, bottom to top 1-9
      for i in 1:9
        periodic_xs[:,(i-1)*size(xs)[2]+1:i*size(xs)[2], :] = xs
        periodic_xs[:,(i-1)*size(xs)[2]+1:i*size(xs)[2], 1] .+= (i-1)%3-1
        periodic_xs[:,(i-1)*size(xs)[2]+1:i*size(xs)[2], 2] .+= Int(floor((i-1)/3))-1
      end

      #TODO improve performance by only looking at close neighbors
      σ1 = 1/((1+NNecks/NGrid)*NGrid^2)
      result = gradientDescent(initialPoint, periodic_xs, xvarz, σ1)

      return(evaluate.(xs, xvarz=>result))
    end
    #=
    function simulateRepulsion(unitGrid, NLayers, NGrid, NeckCables, NNecks)
      initialPoint = Vector{Float64}([])
      xvarz = Vector{Variable}([])
      for i in 1:NLayers, j in 1:NGrid^2
        if !any(t->(t[1]==i-1 || (i==1&&t[1]==NLayers)) && t[2]==j, NeckCables)
          append!(initialPoint, unitGrid[(i-1)*NGrid^2+j][1:2])
          append!(xvarz, [Variable(:x, i, j, 1), Variable(:x, i, j, 2)])
        end
      end

      xs = generateMatrix(NLayers, NGrid, NNecks, NeckCables)
      periodic_xs = Array{Expression,3}(undef, size(xs)[1], size(xs)[2]*9, 3)
      #left to right, bottom to top 1-9
      for i in 1:9
        periodic_xs[:,(i-1)*NGrid^2+1:i*NGrid^2, :] = xs
        periodic_xs[:,(i-1)*NGrid^2+1:i*NGrid^2, 1] .+= (i-1)%3-1
        periodic_xs[:,(i-1)*NGrid^2+1:i*NGrid^2, 2] .+= Int(floor((i-1)/3))-1
      end

      #TODO improve performance by only looking at close neighbors
      σ1 = 1/(NGrid^2)
      result = gradientDescent(initialPoint, periodic_xs, xvarz, σ1)

      return(evaluate.(xs, xvarz=>result))
    end=#

    function generateDisorderedPointClusters(NGrid::Int, NNecks::Int, NLayers::Int)
      NLayers%2==0 || throw(error("There must be an even number of layers!"))

      layers = [(i%2==1 ? ones(NGrid^2+2*NNecks) : 2*ones(NGrid^2+2*NNecks)) for i in 1:NLayers]
      for layer in 1:NLayers
        layers[layer][NGrid^2+1:NGrid^2+NNecks] = 3 .- layers[layer][NGrid^2+1:NGrid^2+NNecks]
      end
      unitGrid = [Point3f(a,b,c) for c in 0:1/NLayers:1-1/NLayers for a in 1/(2*NGrid):1/NGrid:1-1/(2*NGrid) for b in 1/(2*NGrid):1/NGrid:1-1/(2*NGrid) ]
      #TODO add NLayer many random points (Necks)
      randomGrid=[]
      for i in 1:NLayers
        for j in 1:NNecks
          push!(randomGrid, Point3f(rand(),rand(),(i-1)/NLayers))
        end
      end
      #=NeckCables=[]
      foreach(curLayer->
              foreach(curNeck->insertNecks(layers, curLayer, NGrid, NeckCables), 1:NNecks),
          1:NLayers)
      correctedGrid = simulateRepulsion(unitGrid, NLayers, NGrid, NeckCables, NNecks)=#

      correctedGrid = simulateRepulsion(unitGrid, randomGrid, NLayers, NGrid, NNecks)
      fig = Figure()
      ax = Axis3(fig[1,1])

      foreach(k -> scatter!(ax, Point3f(correctedGrid[k[1],k[2],:]), color=layers[k[1]][k[2]] == 1.0 ? :red : :blue), [(i,j) for i in 1:size(correctedGrid)[1] for j in 1:size(correctedGrid)[2]])
      display(fig)

      foreach(k->println(correctedGrid[k[1],k[2],1], " ", correctedGrid[k[1],k[2],2], " ", correctedGrid[k[1],k[2],3], " ", layers[k[1]][k[2]] == 1.0 ? 1 : 2), [(i,j) for i in 1:size(correctedGrid)[1] for j in 1:size(correctedGrid)[2]])
    end


    generateDisorderedPointClusters(3, 2, 2)

end
