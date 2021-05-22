# AC.Letournel
# TP4 Equation de la chaleur

# programmer la solution exacte u pour f(x)=4*pi^2*sin(2*pi*x) en (N+1) points
function calcul_sol_exacte(N)
	x = range(0,1,length=N+1)
	u = -sin.(2*pi*x)
	return x,u
end

# construire la matrice tri-diagonale A, carrée (N-1)
using LinearAlgebra
function calcul_A(N)
    h = 1/N
    A = 1/h^2 * (diagm(0=>-2*ones(N-1))+diagm(1=>ones(N-2))+diagm(-1=>ones(N-2)))
    return A
end

A = calcul_A(N)

# programmer f(x) sur (N+1) points et tronquer f(0) et f(N+1)
function calcul_f(N)
	x = range(0,1,length=N+1)
	x = x[2:N]
	f = 4*pi^2*sin.(2*pi*x)
	return x,f
end

# résoudre le système linéaire A.v=f où v est de dimension N-1
function calcul_v(N)
	A = calcul_A(N)
	f = calcul_f(N)
	v = A\f
  # ajouter v0 et v15
  push!(v,0)
  pushfirst!(v,0)
	return v
end

N = 10
x10, u10 = calcul_sol_exacte(N)
f10 = calcul_f(N)
v10 = calcul_v(N)
println("la norme (u-v) pour N=10 vaut ", norm(u10-v10))


N = 15
x15, u15 = calcul_sol_exacte(N)
f15 = calcul_f(N)
v15 = calcul_v(N)
println("la norme (u-v) pour N=15 vaut  ", norm(u15-v15))

N = 100
x100, u100 = calcul_sol_exacte(N)
f100 = calcul_f(N)
v100 = calcul_v(N)
println("la norme (u-v) pour N=100 vaut  ", norm(u100-v100))





