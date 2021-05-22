using LinearAlgebra

#= BEN YAHIA Bilail 
06/04/2021 =#

println("\t \t TP 2 \n")
println("-------------------------\n")


function monLU(A) #retourne L et U 
   
   #initiliasation 
    AA = copy(A)
    n,m = size(AA) #pour récupérer les dimensions de la matrice 
    
    #début de l'agorithme
    for k =1:n-1
        for i = k+1:n #parcour toutes les lignes sous le pivot
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
    3 -7 4.]

b = [8;1;10.]

println("La matrice A : ")
println(A)
println("\nle vecteur b : ")
println(b)

L,U = monLU(A)
println("\n La décomposition : ")
println("L : ")
println(L)
println("U : ")
println(U)

println("\nvérificaton L*U : ") # L*U = A 
println(L*U)

println("\n\tRésolution de système : \n-------------------------\n")
function descente(L,b)

    n,m = size(L)
    x = zeros(Int, n) # x = Vector{Float64}(undef,n)
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


z = descente(L,b)
println("résoltion descendante de L*z = b, avec z =")
println(z)


function remonte(U, z)
    x = zeros(1, size(b,1))

    for i in 1:size(z,1)
        sum = 0
        for j in i+1:size(z,1)
            sum += x[j] * U[i,j]
        end
        x[i] = (z[i] - sum) / U[i,i]
    end
    return x
end

x = remonte(U,z) 
println("\nrésoltion ascedante de U*z")
println(x)


function monPivotLu(A::Matrix{Float64}) # retourne perm, UP et LP 

   #initiliasation 
    AA = copy(A)
    n = size(A,1) #pour récupérer les dimensions de la matrice 
    Identite = Matrix{Float64}(I, n,n)
    Q = Matrix{Float64}(I, n,n) #matrice de permutation

    for k = 1:(n-1)
        imax=k #ligne du pivot max
        pivot=abs(AA[k,k])

        for i = k+1:n
            if abs(AA[i,k]) > pivot
                pivot = abs(AA[i,k])
                imax = i
            end
        end
        AA[k,:],AA[imax,:] = AA[imax,:],AA[k,:] #on interverti deux lignes entre elles
        Q[k,:],Q[imax,:] = Q[imax,:],Q[k,:]

        for i = k+1:n                  #on refait l'algo de monLU
            AA[i,k] = AA[i,k]/AA[k,k]
        end
        for i = k+1:n
            for j = k+1:n
                AA[i,j] = AA[i,j] - AA[i,k]*AA[k,j]  
            end
        end
    end
    
    #détermination de la matrice L et U 
    LP = tril(AA,-1)+Identite
    UP = triu(AA)
    return LP,UP,Q
end

println("\nMéthode avec pivot :")
LP,UP,Q = monPivotLu(A)
print("LP = ")
println(LP)
print("UP = ")
println(UP)
print("Q = ")
println(Q)

println("\nvérification : ")
println(Q*A - LP*UP)
# doit vérifier Q*A = LP*UP

#Les systèmes à résoudre 
  #1 
  A1 = [1. 1 2;
        -1 -2 3;
        3 -7 4]
  b1 = [8.;1;10]

  #2
  A2 = [1. -1 2 -1;
        2 1 -2 -2;
        -1 2 -4 1;
        3 0 0 -3]
  b2 = [-1.;2;1;-3]

#=
L2,U2,Q2 = monPivotLu(A2)
y2 = descente(L2,b2)
x2 = remonte(U2,y2)
print(" les solutions du systèmes sont : ")
println(x2)
=#

#=script -------------------------------------------------------
A = [1 1 2; 
    -1 -2 3;
    3 -7 4.]
L,U = monLU(A)
L
U
L*U    #on doit retrouver A
L*U-AA #on doit trouver 0 
b = [8;1;10]
y = descente(L,b)
L*y   #on doit retomber sur b 
--------------------------------------------------------------=#