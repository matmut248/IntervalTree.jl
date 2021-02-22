using Test
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using DataStructures,IntervalTrees
include("../examples/tetGrid.jl")



#Questi sono i test del package LinearAlgebraicRepresentation.
#Il codice rivisitato continua a soddisfare tutti i test.
@testset "spaceindex tests" begin

	# 2x2x2 cuboidal grid for 1-, 2-, and 3-dim tests
	V,(VV,EV,FV,CV) = Lar.cuboidGrid([2,2,2],true)
	W,_ = Lar.apply(Lar.r(1,1,pi/6),(V,[VV,EV,FV,CV]))

	function test_bboxes(bboxes)
		# initialize accumulator
		accumulator = BitArray{1}()
		for k=1:size(bboxes[1],1)
			push!(accumulator, true)
		end
		# testing data 
		for h=1:length(bboxes)
			accumulator = (bboxes[h][:,1] .< bboxes[h][:,2]) .& accumulator
		end
		return (&)(accumulator...)
	end
	
	@testset "boundingbox Tests" begin
			
		@testset "Edge tests" begin # 
			cellpoints = [ W[:,EV[k]]::Lar.Points for k=1:length(EV) ]
			bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
			@test true == test_bboxes(bboxes)
		end
		@testset "Face tests" begin # 
			cellpoints = [ W[:,FV[k]]::Lar.Points for k=1:length(FV) ]
			bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
			@test true == test_bboxes(bboxes)
		end
		@testset "Cell tests" begin # 
			cellpoints = [ W[:,CV[k]]::Lar.Points for k=1:length(CV) ]
			bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
			@test true == test_bboxes(bboxes)
		end
	end

	@testset "coordintervals Tests" begin
        # 2x2x2 cuboidal grid for 1-, 2-, and 3-dim tests
        V,(VV,EV,FV,CV) = Lar.cuboidGrid([2,2,2],true)
        W,_ = Lar.apply(Lar.r(1,1,pi/6),(V,[VV,EV,FV,CV]))
			
		cellpoints = [ W[:,EV[k]]::Lar.Points for k=1:length(EV) ]
		bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
		@testset "Edge tests" begin # 
			@test typeof(coordintervals(1,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(coordintervals(2,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(coordintervals(3,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
		end
		cellpoints = [ W[:,FV[k]]::Lar.Points for k=1:length(FV) ]
		bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
		@testset "Face tests" begin # 
			@test typeof(coordintervals(1,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(coordintervals(2,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(coordintervals(3,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
		end
		cellpoints = [ W[:,CV[k]]::Lar.Points for k=1:length(CV) ]
		bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
		@testset "Cell tests" begin # 
			@test typeof(coordintervals(1,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(coordintervals(2,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(coordintervals(3,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
		end
	end
	
	@testset "boxcovering Tests" begin
		V,(VV,EV,FV,CV) = Lar.cuboidGrid([2,2,2],true)
		W,_ = Lar.apply(Lar.r(1,1,pi/6),(V,[VV,EV,FV,CV]))
		cellpoints = [ W[:,EV[k]]::Lar.Points for k=1:length(EV) ]
		bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
		dict = coordintervals(1,bboxes)
		@test typeof(dict) == OrderedDict{Array{Float64,1},Array{Int64,1}}
		@test length(coordintervals(1,bboxes)) == 54
		@test length(coordintervals(2,bboxes)) == 54
		@test length(coordintervals(3,bboxes)) == 54

		V,(VV,EV,FV) = Lar.cuboidGrid([2,1],true)
		cellpoints = [ V[:,EV[k]]::Lar.Points for k=1:length(EV) ]
		bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
		@test bboxes == [[0.0 0.0; 0.0 1.0],
        [1.0 1.0; 0.0 1.0],
        [2.0 2.0; 0.0 1.0],
        [0.0 1.0; 0.0 0.0],
        [0.0 1.0; 1.0 1.0],
        [1.0 2.0; 0.0 0.0],
        [1.0 2.0; 1.0 1.0]]
        xboxdict = Dict(
         [0.0, 0.0] => [1],
         [1.0, 1.0] => [2],
         [2.0, 2.0] => [3],
         [0.0, 1.0] => [4, 5],
         [1.0, 2.0] => [6, 7])
        @test xboxdict == coordintervals(1,bboxes)
		xs = IntervalTrees.IntervalMap{Float64, Array}()
		for (key,boxset) in xboxdict
			xs[tuple(key...)] = boxset
		end
       @test typeof(xs) ==
		IntervalTrees.IntervalBTree{Float64,
		IntervalValue{Float64,Array},64}
	end
	
	@testset "Refactoring spaceindex tests" begin
		V,(VV,EV,FV) = Lar.cuboidGrid([2,1],true)
		EV = [[1, 2], [3, 4], [5, 6], [1, 3], [2, 4], [3, 5], [4, 6]]
		cellpoints = [ V[:,EV[k]]::Lar.Points for k=1:length(EV) ]
		bboxes = [hcat(boundingbox(cell)...) for cell in cellpoints]
		xboxdict = coordintervals(1,bboxes)
		yboxdict = coordintervals(2,bboxes)
		xs = IntervalTrees.IntervalMap{Float64, Array}()
		for (key,boxset) in xboxdict
			xs[tuple(key...)] = boxset
		end
		ys = IntervalTrees.IntervalMap{Float64, Array}()
		for (key,boxset) in yboxdict
			ys[tuple(key...)] = boxset
		end
		xcovers = boxcovering(bboxes, 1, xs)
		ycovers = boxcovering(bboxes, 2, ys)
		covers = [intersect(pair...) for pair in zip(xcovers,ycovers)]
		
		@test covers == Array{Int64,1}[[1, 4, 5], [4, 5, 2, 6, 7], [6, 7, 3], 
			[1, 4, 2, 6], [1, 5, 2, 7], [4, 2, 6, 3], [5, 2, 7, 3]]
	end
end


#test delle nuove funzioni
@testset "test nuove funzioni" begin
	
	@testset "createIntervalTree test" begin
		dict = OrderedDict([0.0, 1.0] => [1, 3],[1.0, 1.0] => [2],[0.0, 0.0] => [4],[0.0, 2.0] => [5])
		t = createIntervalTree(dict)
		@test typeof(t) == IntervalTrees.IntervalBTree{Float64,IntervalValue{Float64,Array},64}
		@test t.root.entries[1].first == 0.0
		@test t.root.entries[1].last == 0.0
		@test t.root.entries[1].value == [4]

		@test t.root.entries[2].first == 0.0
		@test t.root.entries[2].last == 1.0
		@test t.root.entries[2].value == [1,3]

		@test t.root.entries[3].first == 0.0
		@test t.root.entries[3].last == 2.0
		@test t.root.entries[3].value == [5]

		@test t.root.entries[4].first == 1.0
		@test t.root.entries[4].last == 1.0
		@test t.root.entries[4].value == [2]
	end

	@testset "removeSelfIntersection! test" begin
		covers = [[4, 1, 3, 5, 2], [1, 3, 5, 2], [4, 1, 3, 5, 2], [4, 1, 3, 5], [4, 1, 3, 5, 2]]
		removeSelfIntersection!(covers)
		@test typeof(covers) == Array{Array{Int64,1},1}
		@test covers[1] == [4, 3, 5, 2]
		@test covers[2] == [1, 3, 5]
		@test covers[3] == [4, 1, 5, 2]
		@test covers[4] == [1, 3, 5]
		@test covers[5] == [4, 1, 3, 2]

	end
	@testset "addIntersection! test" begin
		bb = [[0.0 1.0; 0.0 0.0], [1.0 1.0; 0.0 1.0], [0.0 1.0; 1.0 1.0], [0.0 0.0; 0.0 1.0], [0.0 2.0; 0.0 1.0]];
		dict = OrderedDict([0.0, 1.0] => [1, 3],[1.0, 1.0] => [2],[0.0, 0.0] => [4],[0.0, 2.0] => [5])
		t = createIntervalTree(dict)
		c = boxcovering(bb,1,t)
		@test c == [[4, 1, 3, 5, 2], [1, 3, 5, 2], [4, 1, 3, 5, 2], [4, 1, 3, 5], [4, 1, 3, 5, 2]];
	end

end

#test di un caso complesso: griglia tridimensionale con 2 tetraedri randomici su ogni punto intero della griglia
@testset "tetgrid test" begin
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
	(V,FV,EV) = tetgrid.randomTetGrid()
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