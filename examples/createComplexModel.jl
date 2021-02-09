#=
funzione che crea un modello geometrico complesso per il testing.
il modello viene creato a partire da un numero di facce passate come parametro,
e un numero di punti pari al numero di facce /3 generati in modo randomico.
=#

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