using Test
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
include("intervaltree.jl")
using DataStructures,IntervalTrees


"""
Questi sono i test del package LinearAlgebraicRepresentation.
Il codice rivisitato continua a soddisfare tutti i test.
"""

@testset "Old spaceindex tests" begin

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

#=
"""
@testset "New spaceindex tests" begin

    @testset "spaceindex tests" begin
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
"""
=#