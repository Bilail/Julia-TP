using LinearAlgebra

#= BEN YAHIA Bilail 
11/05/2021 
Methode de dicrétisation : on a une barre de longueur 1m que l'on découpe en N morceaux (N+1 points)
=#

println("\t \t TP 4")
println("-------------------------\n")

function calcul_sol_exact(N) #retourne un vecteur
    sol_exact = []
    vec_x = []
    for i in range(0, 1, length=N+1) #on peut soit utilisé la longueur ou le pas pour range
        push!(sol_exact,-sin(2*pi*i))
        push!(vec_x,i)
    end
    return vec_x, sol_exact
end

#= Autre méthode
vec_x = range(0, 1, length=N+1)
sol_exact = -sin.vec_x =#
N = 10
println("Soit N = 10, discrétisation : ")
x,sol_exact = calcul_sol_exact(N)
println(" le vecteur x : ", x)
println(" le vecteur sol_excat : " , sol_exact)

function calcul_A(N) #retourne la matrice A à partir de N
    A = zeros(Int,N-1, N-1)
    for i = 1:N-2
        A[i,i] = -2
        A[i+1,i] = 1
        A[i,i+1] = 1
    end
    return A
end 

println("Calcul de la matrice A", calcul_A(N))

function calcul_f(N,X,u,h)
    f = []
    for i = 1:N-1
        push!(f,(u[i-1]-2*u[i]+u[i+1])/h*h)
    end 
end 
h = 1/N
println("Calcul de f ", calcul_f(N,x,sol_exact,h))
