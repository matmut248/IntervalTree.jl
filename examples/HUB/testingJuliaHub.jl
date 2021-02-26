using Pkg
Pkg.resolve()
Pkg.instantiate()
Pkg.activate(".")
Pkg.add("LinearAlgebraicRepresentation")
Pkg.add("IntervalTrees")
Pkg.add("BenchmarkTools")
Pkg.add("OrderedCollections")
Pkg.add("Test")

println("ciao")

using IntervalTrees
using BenchmarkTools
using OrderedCollections
using Base.Threads
using Test
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

function coordintervals(coord,bboxes)
    boxdict = OrderedDict{Array{Float64,1},Array{Int64,1}}()
    for (h,box) in enumerate(bboxes)
        key = box[coord,:]
        if haskey(boxdict,key) == false
            boxdict[key] = [h]
        else
            push!(boxdict[key], h)
        end
    end
    return boxdict
end
function boxcovering(bboxes, index, tree)
    covers = [[zero(eltype(Int64))] for k=1:length(bboxes)]		#zero(eltype(Int64)) serve per rendere covers type stable
    @threads for (i,boundingbox) in collect(enumerate(bboxes))
        extent = bboxes[i][index,:]
        iterator = IntervalTrees.intersect(tree, tuple(extent...))
        addIntersection!(covers, i, iterator)
    end
    return covers
end
function addIntersection!(covers::Array{Array{Int64,1},1}, i::Int64, iterator)
    splice!(covers[i],1)		#splice serve a togliere gli zeri iniziali all'interno di covers
    @threads for x in collect(iterator)
        append!(covers[i],x.value)
    end
end
function createIntervalTree(boxdict::AbstractDict{Array{Float64,1},Array{Int64,1}})
    tree = IntervalTrees.IntervalMap{Float64,Array}()
    for (key, boxset) in boxdict
        tree[tuple(key...)] = boxset
    end
    return tree
end
function removeSelfIntersection!(covers::Array{Array{Int64,1},1})
	@threads for k=1:length(covers)
        covers[k] = setdiff(covers[k],[k])	#toglie le intersezioni con se stesso 
    end
end
function boundingbox(vertices::Lar.Points)
    firstDim = vertices[1,:]
    secondDim = vertices[2,:]
    if (size(vertices,1)==3)
        thirdDim = vertices[3,:]
         minimum = Threads.@spawn hcat([min(firstDim...), min(secondDim...), min(thirdDim...)])
         maximum = Threads.@spawn hcat([max(firstDim...), max(secondDim...), max(thirdDim...)])
    else
         minimum = Threads.@spawn hcat([min(firstDim...), min(secondDim...)])
         maximum = Threads.@spawn hcat([max(firstDim...), max(secondDim...)])
    end
    return fetch(minimum),fetch(maximum)
end
function createModel(numFaces)
    numPoints = div(numFaces,3)
    V=(zeros(Float64,(3,numPoints)))
    FV=[[zero(eltype(Int64))]]
    pop!(FV)

    for i = 1:numPoints
        x,y,z = rand(Float64, (1, 3))
        V[:,i] = [x,y,z]
    end

    range = 1:numFaces 
    for i = 1:numFaces
        f=Int(rand((1:numPoints)))
        s=Int(rand((1:numPoints)))
        t=Int(rand((1:numPoints)))
        if f!=s && s!=t && t!=f && !in([f,s,t],FV)
            push!(FV,[f,s,t])
        else
            i = i-1
        end
    end
    return (V,FV)
end
function spaceindex(model::Lar.LAR)::Array{Array{Int,1},1}
    V,CV = model[1:2]
    dim = size(V,1)
    
    cellpoints = [ V[:,CV[k]]::Lar.Points for k=1:length(CV) ]		    #calcola le celle
    bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]    #calcola i boundingbox delle celle
    
    xboxdict = Threads.@spawn coordintervals(1,bboxes)
    yboxdict = Threads.@spawn coordintervals(2,bboxes)

    # xs,ys sono di tipo IntervalTree
    xs = Threads.@spawn createIntervalTree(fetch(xboxdict))
    ys = Threads.@spawn createIntervalTree(fetch(yboxdict))
    
    xcovers = Threads.@spawn boxcovering(bboxes, 1, fetch(xs))                        #lista delle intersezioni dei bb sulla coordinata x
    ycovers = Threads.@spawn boxcovering(bboxes, 2, fetch(ys))                        #lista delle intersezioni dei bb sulla coordinata x
    covers = [intersect(pair...) for pair in zip(fetch(xcovers),fetch(ycovers))]      #lista delle intersezioni dei bb su entrambe le coordinate

    if dim == 3
		zboxdict = Threads.@spawn coordintervals(3,bboxes)
		zs = Threads.@spawn createIntervalTree(fetch(zboxdict))
		zcovers = Threads.@spawn boxcovering(bboxes, 3, fetch(zs))
		covers = [intersect(pair...) for pair in zip(fetch(zcovers),covers)]
    end
    
    removeSelfIntersection!(covers)       #rimozione delle intersezioni con se stesso
    return covers
end

@testset "BigData test" begin
	"""
	funzione che ricava tutte le coppie di vertici opposti del bounding box
	a partire dai punti che lo definiscono
	"""
	function buildPoints(bb)
		x1,y1,z1 = bb[1:3]
		x2,y2,z2 = bb[4:6]
		p1 = [x1 x2; y1 y2; z1 z2]
		p2 = [x2 x1; y1 y2; z1 z2]
		p3 = [x1 x2; y2 y1; z1 z2]
		p4 = [x1 x2; y1 y2; z2 z1]
		return p1,p2,p3,p4
	end

	"""
	funzione che calcola se c'è intersezione tra due boundingbox
	c'è intersezione se esiste almeno un punto nel boundingbox le cui dimensioni sono comprese
	all'interno del secondo boundingbox, o viceversa 
	"""
	function hasIntersection(bb1,bb2)
		#punti del boundingbox bb1
		x11,y11,z11 = bb1[1:3]
		x12,y12,z12 = bb1[4:6]
		#punti del boundingbox bb2
		x21,y21,z21 = bb2[1:3]
		x22,y22,z22 = bb2[4:6]

		#si controlla re il primo vertice di bb1 è contenuto in bb2
		#oppure se bb1 contiene bb2 
		first = ((x11>=x21 && x11<=x22) || (x11<=x21 && x12>=x22)) && ((y11>=y21 && y11<=y22) || (y11<=y21 && y12>=y22)) && ((z11>=z21 && z11<=z22) || (z11<=z21 && z12>=z22))
		#si controlla se il secondo vertice di bb1 è contenuto in bb2 
		second = (x12>=x21 && x12<=x22) && (y12>=y21 && y12<=y22) && (z12>=z21 && z12<=z22)

		return first || second
	end

	#creazione del modello complesso
	(V,FV) = createModel(1000)
	cover = spaceindex((V,FV))
	cells = [V[:,FV[k]] for k = 1:length(FV)]
	bb = [hcat(boundingbox(c)...) for c in cells]

	#testing delle intersezioni
	for j = 1:length(bb)
		for i = 1:length(cover[j])
			p1,p2,p3,p4 = buildPoints(bb[j])
			p5,p6,p7,p8 = buildPoints(bb[cover[j][i]])
			@test (hasIntersection(p1,bb[cover[j][i]]) || hasIntersection(p5,bb[j])
				|| hasIntersection(p2,bb[cover[j][i]]) || hasIntersection(p6,bb[j]) 
				|| hasIntersection(p3,bb[cover[j][i]]) || hasIntersection(p7,bb[j]) 
				|| hasIntersection(p4,bb[cover[j][i]]) || hasIntersection(p8,bb[j]) ) == true
		end
	end

	#testing delle non intersezioni
	for j = 1:length(bb)
		noIntersection = setdiff((h for h = 1:length(FV)),cover[j])		#tutti i bb che non intersecano il j-esimo bb
		for i = 1:length(noIntersection)
			if bb[j] != bb[noIntersection[i]]
				p1,p2,p3,p4 = buildPoints(bb[j])
				p5,p6,p7,p8 = buildPoints(bb[noIntersection[i]])
				@test (hasIntersection(p1,bb[noIntersection[i]]) || hasIntersection(p5,bb[j])
					|| hasIntersection(p2,bb[noIntersection[i]]) || hasIntersection(p6,bb[j]) 
					|| hasIntersection(p3,bb[noIntersection[i]]) || hasIntersection(p7,bb[j]) 
					|| hasIntersection(p4,bb[noIntersection[i]]) || hasIntersection(p8,bb[j]) ) == false
			end		
		end
	end
end