
using LinearAlgebra

#= BEN YAHIA Bilail 
23/03/2021 =#
function monLU(A)
   
   #initiliasation 
    AA = copy(A)
    n,m = size(AA)
    
    #début de l'agorithme
    for k =1:n-1
        for i = k+1:n
            AA[i,k] = AA[i,k]/AA[k,k]
        end
        for i = k+1:n
            for j = k+1:n
                AA[i,j] = AA[i,j] - AA[i,k]*AA[k,j]  
            end
        end
    end
    
    #détermination de la matrice L et U 
    Identite = Matrix{Float64}(I, n,n)
    L = tril(AA,-1)+Identite
    U = triu(AA)
    return L,U
end         

#création des tests
A = [1 1 2; 
    -1 -2 3;
    3 -7 4]

b = [8,1,10]

L,U = monLU(A)
L

U

L*U

function descente(L,b)

    n,m = size(L)
    x = zeros(Int, n)
    t = 0
    
    for i = 1:n
        for j = 1:i-1
            t = t+L[i,j]*x[j]
        end
        x[i] = b[i]-t
        t = 0
    end
    return x
end


y = descente(L,b)

L*y

function remontee(U,z)
    n,m = size(U)
    x = zeros(Int, n)
    t = 0
    
    for i = 1:n
        for j = 1:i-1
            t = t+U[n-i,j]*x[n-j]
        end
        x[n-i] = 1/U[i,i]*(z[n-i] - t) 
    end
    return x
end


remontee(U,b)


