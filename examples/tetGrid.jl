#=
1. Create a regular $N^3$ grid T of random tetrahedra:
	(a) Generate four vertices for each tetrahedron randomly in a half-open fixed-size cube.
	(b) Let these cubes form a regular $N^3$ grid.
2. Create a regular $(N − 1)^3$ grid C of such cubes.
3. Align T and C such that the grid nodes of C are at the centers of the grid cells of T .
4. Measure time for T ∪ C.
=#
module tetgrid

    using ViewerGL; GL = ViewerGL
    using LinearAlgebraicRepresentation
    Lar = LinearAlgebraicRepresentation;

    function t1t2(i,j,k)
        v1 = rand(3,4).+ [2i,2j,2k]
        v2 = rand(3,4).+ [2i,2j,2k]
        ev = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
        fv = [[2,3,4],[1,3,4],[1,2,4],[1,2,3]]
        t1 = [v1,fv,ev] 
        t2 = [v2,fv,ev]
        return t1,t2
    end

    function randomTetGrid()
        
        N = 3  # N ≥ 2
        tets = []
        for i=0:N-1, j=0:N-1, k=0:N-1
            push!(tets, t1t2(i,j,k)...)
        end
        V,FV,EV = Lar.struct2lar(Lar.Struct(tets))
        #GL.VIEW([ GL.GLGrid(V,FV, GL.COLORS[1],1), GL.GLFrame2 ]);

        W,(_,EW,FW,_) = Lar.cuboidGrid([N-1,N-1,N-1],true);
        U = (W .* 2) .+ 0.5
        #GL.VIEW([ GL.GLGrid(V,FV,GL.COLORS[1],0.5), GL.GLGrid(U,FW,GL.COLORS[1],0.5), GL.GLFrame2 ]);

        V,FV,EV = Lar.struct2lar(Lar.Struct([Lar.Struct(tets), Lar.Struct([(U,FW,EW)]) ]))

        GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],0.5), GL.GLFrame2 ]);
        return (V,FV,EV)
    end
end