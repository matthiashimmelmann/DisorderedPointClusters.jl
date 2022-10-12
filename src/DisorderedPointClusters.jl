module DisorderedPointClusters

import GLMakie: scatter!, Axis3, Figure, Point3f
import LinearAlgebra: norm, det, pinv
import HomotopyContinuation: Variable, Expression, evaluate, differentiate
import Random: shuffle

export generateDisorderedPointClusters

#TODO disordered surface with >=3 labyrinths
#TODO calculate one optimum. Transform the points upward in some way (e.g. rotation/reflection). Reflection is easy.
#For rotation, multiply by e^ipi/2 = cos(pi/2)+isin(pi/2) with standard multiplication in C.
#TODO 3D optimization on the Torus...
#TODO put a repulsive grid around the necks
#TODO For Layers: for each neck: look for nearby points and color them differently=> list 1:gridlength for each layer, delete in top/bottom layer and nearby (<2^(1/6)*sigma_1) and continue insertion.

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

    function generateMatrix(NGrid, NNecks, layerValue, totalLayers, numFixed)
      xs=Array{Expression,2}(undef, NGrid+2*NNecks, 3)
      for j in 1:NGrid+2*NNecks-numFixed
        xs[j,1] = Variable(:x, j, 1)
        xs[j,2] = Variable(:x, j, 2)
      end
      xs[:,3] .= layerValue
      return(xs)
    end

    function gradientDescent(initialPoint, periodic_xs, xvarz, σ1)
      cursol = Base.copy(initialPoint)
      eval_periodic_xs = evaluate.(periodic_xs, xvarz=>cursol)
      list_of_relevant_molecule_interactions = []
      for i in Int(4*size(periodic_xs)[1]/9)+1:Int(5*size(periodic_xs)[1]/9), j in 1:size(periodic_xs)[1]
        if i!=j && #=&& sqrt(sum((eval_periodic_xs[layer,i,:]-eval_periodic_xs[layer,j,:]).^2)) < 2^(1/6)*σ1 &&=# all(t->t[1]!=j||t[2]!=i, list_of_relevant_molecule_interactions)
          push!(list_of_relevant_molecule_interactions, (i,j))
        end
      end
      distances = [sum((periodic_xs[i,:] - periodic_xs[j,:]).^2) for (i,j) in list_of_relevant_molecule_interactions]
      #boundarydistances = vcat([[(periodic_xs[i,1]-1.)^2, (periodic_xs[i,1])^2, (periodic_xs[i,2]-1.)^2, (periodic_xs[i,2])^2] for i in Int(4*size(periodic_xs)[1]/9)+1:Int(5*size(periodic_xs)[1]/9)]...)
      lennardJones = sum(((σ1/d)^6-(σ1/d)^3 for d in distances)) #TODO ADD DIFFERENT NUMBER FOR NECKS -> REPULSION
      #Boundaries = sum(1/d^2 for d in boundarydistances)
      ∇Q = differentiate(lennardJones, xvarz)
      HessQ = differentiate(∇Q, xvarz)
      while norm(evaluate(∇Q, xvarz=>cursol))>1e-4
        cursol = cursol - 0.5*pinv(evaluate(HessQ, xvarz=>cursol))*evaluate(∇Q, xvarz=>cursol)
        if any(t->t<0||t>1,cursol)
          display("Left Torus")
          cursol = cursol - floor.(cursol)
        end
        display(norm(evaluate(∇Q, xvarz=>cursol)))
      end
      return(cursol)
    end

    function simulateRepulsion(unitGrid, NGrid, NNecks, fixedPoints, curLayer, NLayers)
      initialPoint = Vector{Float64}([])
      xvarz = Vector{Variable}([])
      previousLayer = []
      layerList = (1:Int(length(unitGrid)/NLayers))[1:NGrid+2*NNecks-length(fixedPoints)]#shuffle(1:Int(length(unitGrid)/NLayers))[1:NGrid+2*NNecks-length(fixedPoints)]

      for j in 1:NGrid+2*NNecks-length(fixedPoints)
        append!(xvarz, [Variable(:x, j, 1), Variable(:x, j, 2)])
        append!(initialPoint, unitGrid[Int((curLayer-1)*length(unitGrid)/NLayers)+layerList[j]][1:2])
      end

      xs = generateMatrix(NGrid, NNecks, (curLayer-1)/NLayers, NLayers, length(fixedPoints))
      for i in 1:length(fixedPoints)
        xs[i+NGrid+2*NNecks-length(fixedPoints), 1:2] = fixedPoints[i][1:2]
      end
      eval_xs = evaluate(xs, xvarz=>initialPoint)

      periodic_xs = Array{Expression,2}(undef, size(xs)[1]*9, 3)
      #left to right, bottom to top 1-9
      for i in 1:9
        periodic_xs[(i-1)*size(xs)[1]+1:i*size(xs)[1], :] = xs
        periodic_xs[(i-1)*size(xs)[1]+1:i*size(xs)[1], 1] .+= (i-1)%3-1
        periodic_xs[(i-1)*size(xs)[1]+1:i*size(xs)[1], 2] .+= Int(floor((i-1)/3))-1
      end

      σ1 = 1/((sqrt(2*NNecks+NGrid)))
      result = gradientDescent(initialPoint, periodic_xs, xvarz, σ1)

      return(evaluate.(xs, xvarz=>result))
    end

    function calculateAdjacents(pos, grid)
      σ = 2^(1/6)/((sqrt(size(grid)[1])))
      indices = findall(t->
        any(q->norm(t-grid[pos,:]+q)<σ, [[a%3-1, floor((a-1)/3)-1, 0.] for a in 1:9]),
        [grid[i,:] for i in 1:size(grid)[1]])
      return indices
    end

    function generateDisorderedPointClusters(NGrid::Int, NNecks::Int, NLayers::Int)
      #TODO calcuate min only once, rest follows.
      NLayers%2==0 || throw(error("There must be an even number of layers!"))
      maxNumberPoints = NGrid+2*NNecks
      nextGridNumber = ceil(sqrt(maxNumberPoints)) #Calculate the next highest perfect square's square root.
      gridArray = []
      unitGrid = [Point3f(a,b,c) for c in 0:1/NLayers:(NLayers-1)/NLayers for a in 1/(2*nextGridNumber):1/nextGridNumber:1-1/(2*nextGridNumber) for b in 1/(2*nextGridNumber):1/nextGridNumber:1-1/(2*nextGridNumber)]#[vcat([Int(nextGridNumber^2*(layer-1)+1):Int(nextGridNumber^2*(layer-1)+totalNumberPoints) for layer in 1:NLayers]...)]

      correctedGrid = simulateRepulsion(unitGrid, NGrid, NNecks, [], 1, NLayers)
      for z in 0:1/NLayers:(NLayers-1)/NLayers
        helperGrid = Array{Float64,2}(undef, size(correctedGrid)[1], size(correctedGrid)[2])
        helperGrid .= correctedGrid
        helperGrid[:,3] .= z
        push!(gridArray, helperGrid)
      end

      fig = Figure()
      ax = Axis3(fig[1,1])

      layers = [[] for _ in 1:length(gridArray)]
      numbersInGrid = [[i for i in 1:size(gridArray[layer])[1]] for layer in 1:length(gridArray)]
      for i in 1:length(gridArray)
        layers[i] = [i%2 for _ in 1:size(gridArray[i])[1]]
        for j in 1:NNecks
          randPos = rand(numbersInGrid[i])
          layers[i][randPos] = 1-layers[i][randPos]
          numbersInGrid[i] = filter(t-> !(t in calculateAdjacents(randPos, gridArray[i])), numbersInGrid[i])
          if i!=length(gridArray)
            numbersInGrid[i+1] = filter(t-> t!=randPos, numbersInGrid[i+1])
          end
        end
        println(numbersInGrid[i])
        #layers[i][shuffle(1:size(gridArray[i])[1])[1:NNecks]] .= 1-i%2
      end
      println(layers)

      foreach(k -> scatter!(ax, Point3f(gridArray[k[2]][k[1],:]), color=layers[k[2]][k[1]] == 1 ? :red : :blue), [(i,j) for i in 1:size(gridArray[1])[1] for j in 1:length(gridArray)])
      display(fig)

      foreach(k->println(gridArray[k[2]][k[1],1], " ", gridArray[k[2]][k[1],2], " ", gridArray[k[2]][k[1],3], " ", layers[k[2]][k[1]] == 1 ? "1" : "2"), [(i,j) for i in 1:size(gridArray[1])[1] for j in 1:length(gridArray)])
    end


    generateDisorderedPointClusters(16, 6, 8)
end
