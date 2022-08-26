module DisorderedPointClusters

import GLMakie: scatter!, Axis3, Figure, Point3f
import LinearAlgebra: norm
import HomotopyContinuation: Variable, Expression, evaluate, differentiate
export generateDisorderedPointClusters
    #TODO add analytic differentiation

    function insertNecks(layers, curLayer, NGrid, Necks)
      pos = rand(1:NGrid^2)
      layerAbove = layers[curLayer==length(layers) ? 1 : curLayer+1][pos] == layers[curLayer][pos]
      layerBelow = layers[curLayer==1 ? length(layers) : curLayer-1][pos] == layers[curLayer][pos]
      #At the moment only orthogonal adjacency asked for
      anyNeckAdjacent = any(t-> t==(curLayer%2==1 ? 2 : 1), layers[curLayer][vcat(pos, pos%NGrid==0 ? Int((pos/NGrid-1)*NGrid)+1 : pos+1, (pos-1)%NGrid==0 ? Int(((pos-1)/NGrid+1)*NGrid) : (pos-1), pos-NGrid<=0 ? pos-NGrid+NGrid^2 : pos-NGrid, pos+NGrid>NGrid^2 ? pos+NGrid-NGrid^2 : pos+NGrid)])
      if anyNeckAdjacent || layerAbove || layerBelow
        println("Insertion did not work, trying again...")
        insertNecks(layers, curLayer, NGrid, Necks)
      else
        layers[curLayer][pos] = 3-layers[curLayer][pos]
        push!(Necks, [curLayer,pos])
      end
    end

    function generateMatrix(NLayers, NGrid, NNecks, NeckCables)
      xs=Array{Expression,3}(undef, NLayers, NGrid^2, 3)
      for (i,j) in [(a,b) for a in 1:NLayers for b in 1:NGrid^2]
        xs[i,j,1] = Variable(:x, i, j, 1)
        xs[i,j,2] = Variable(:x, i, j, 2)
        xs[i,:,3] .= (0:1/NLayers:1-1/NLayers)[i]
      end
      #TODO add NeckCableEquations

      for cable in NeckCables
        xs[cable[1]==NLayers ? 1 : cable[1]+1, cable[2], 1:2] = xs[cable[1], cable[2], 1:2]
      end
      return(xs)
    end

    function gradientdescent(gradient, HessQ, xvarz, curpoint)
      cursol = Base.copy(curpoint)
      while norm(evaluate(gradient, xvarz=>cursol))>1e-5
        Jac = evaluate(HessQ, xvarz=>cursol)
        cursol = cursol - Jac\evaluate(gradient, xvarz=>cursol)
        #cursol = cursol - floor.(cursol)
        display(norm(evaluate(gradient, xvarz=>cursol)))
      end
      return(cursol)
    end

    function simulateRepulsion(unitGrid, NLayers, NGrid, NeckCables, NNecks)
      initialPoint = []
      xvarz = []
      for i in 1:NLayers
        for j in 1:NGrid^2
          if any(t->(t[1]==i-1||i==1&&t[1]==NLayers) && t[2]==j,NeckCables) || i==1 && (NLayers,j) in NeckCables
            continue
          else
            append!(initialPoint, unitGrid[(i-1)*NGrid^2+j][1:2])
            append!(xvarz, [Variable(:x, i, j, 1), Variable(:x, i, j, 2)])
          end
        end
      end
      xvarz = Vector{Variable}(xvarz)
      xs = generateMatrix(NLayers, NGrid, NNecks, NeckCables)

      dist = vcat([(sum((xs[layer,i,:] - xs[layer,j,:]).^2)) for layer in 1:NLayers for i in 1:NGrid^2 for j in i+1:NGrid^2]...)
      distBoundary = vcat([sum([(xs[layer,i,1]-0)^2,(xs[layer,i,1]-1)^2,(xs[layer,i,2]-0)^2,(xs[layer,i,2]-1)^2]) for layer in 1:NLayers for i in 1:NGrid^2]...)
      NeckEquations = [sum((xs[cable[1],cable[2],:]-xs[cable[1],cable[3],:]).^2) for cable in filter(t->t[2]!=t[3], [vcat(cable,j) for cable in NeckCables for j in 1:NGrid^2])]
      σ1 = 1/(1*NGrid^2)
      σ2 = 1/(2*NGrid^2)
      #TODO add a higher/lower factor for necks.
      lennardJones = 4*sum((σ1/d)^6-(σ1/d)^3 for d in dist) + 4*sum((σ2/d)^6-(σ2/d)^3 for d in distBoundary) + 0.02*sum((1.8*σ1/d)^6-(1.8*σ1/d)^3 for d in NeckEquations)
      ∇Q = differentiate(lennardJones, xvarz)
      HessQ = differentiate(differentiate(lennardJones, xvarz), xvarz)
      result = gradientdescent(∇Q, HessQ, xvarz, initialPoint)

      return(evaluate.(xs, xvarz=>result))
    end

    function generateDisorderedPointClusters(NGrid::Int, NNecks::Int, NLayers::Int)
      NLayers%2==0 || throw(error("There must be an even number of layers!"))

      layers = [(i%2==1 ? ones(NGrid^2) : 2*ones(NGrid^2)) for i in 1:NLayers]
      unitGrid = [Point3f(a,b,c) for c in 0:1/NLayers:1-1/NLayers for a in 1/(2*NGrid):1/NGrid:1-1/(2*NGrid) for b in 1/(2*NGrid):1/NGrid:1-1/(2*NGrid) ]
      #TODO add NLayer many random points (Necks)
      #=randomGrid=[]
      for i in 1:NLayers
        for j in 1:NNecks
          push!(randomGrid, Point3f(rand(),rand(),(i-1)/NLayers))
        end
      end=#
      NeckCables=[]
      foreach(curLayer->
              foreach(curNeck->insertNecks(layers, curLayer, NGrid, NeckCables), 1:NNecks),
          1:NLayers)
      correctedGrid = simulateRepulsion(unitGrid, NLayers, NGrid, NeckCables, NNecks)

      #correctedGrid = simulateRepulsion(unitGrid, randomGrid, NLayers, NGrid, NNecks)
      fig = Figure()
      ax = Axis3(fig[1,1])

      foreach(k -> scatter!(ax, Point3f(correctedGrid[k[1],k[2],:]), color=layers[k[1]][k[2]] == 1.0 ? :red : :blue), [(i,j) for i in 1:NLayers for j in 1:NGrid^2])
      display(fig)

      foreach(k->println(correctedGrid[k[1],k[2],1], " ", correctedGrid[k[1],k[2],2], " ", correctedGrid[k[1],k[2],3], " ", layers[k[1]][k[2]] == 1.0 ? 1 : 2), [(i,j) for i in 1:NLayers for j in 1:NGrid^2])
    end


    generateDisorderedPointClusters(6, 10, 6)

end
