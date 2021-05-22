using LinearAlgebra
using Plots

function u(x)
    return -sin(2*pi*x)
end

function f(x)
    return 4*pi^2*sin(2*pi*x)
end

function calcul_sol_exact(N)
    x = [k/N for k in 0:N] 
    sol_exact  = u.(x)
    return x, sol_exact
end

function calcul_A(N)
    d = [1 for k in 1:N-2]
    D = [-2 for k in 1:N-1]
    return Tridiagonal(d,D,d)
end

function calcul_f(N)
    return [f(k/N) for k in 1:N-1]
end

function calcul_v(N)
    A = calcul_A(N)
    v = A\calcul_f(N)
    v = 1/(N^2) * v
    return v
end

#Pour la résolution du système tridiagonal, on pourrait aussi utiliser l'algorithme de Thomas.
#Ce devrait être plus efficace que l'utilisation de la méthode A\.
function calcul_v_Thomas(N)
    a = [1. for k in 1:N-2]
    b = [-2. for k in 1:N-1]
    c = [1. for k in 1:N-2]
    d = calcul_f(N)
    x = [0. for k in 1:N-1]
    for i in 2:N-1
        w = a[i-1]/b[i-1]
        b[i] = b[i] - w*c[i-1]
        d[i] = d[i] - w*d[i-1]
    end
    x[N-1] = d[N-1]/b[N-1]
    for i in N-2:-1:1
        x[i] = (d[i] - c[i]*x[i+1])/b[i]
    end
    x = 1/(N^2) * x
    return x
end

function calcul_erreur_approx(N)
    v = vcat([0],calcul_v_Thomas(N),[0])
    x, sol_exact = calcul_sol_exact(N)
    erreur_approx = max(abs.(v-sol_exact)...)
    return erreur_approx
end

#On trace les valeurs de u(x_k) et des solutions approchées v_k pour N=10.
#On peut rajouter l'argument seriestype = :scatter dans la fonction plot pour n'afficher que les points.
N = 10;
x, sol_exact = calcul_sol_exact(N);
v = vcat([0],calcul_v_Thomas(N),[0]);
p = plot(x, sol_exact, label = "Solution exacte", ylabel = "y", xlabel = "x", color=:red, title = "N = 10")
plot!(p, x, v, color=:green, label="Solution approchée")

#Même chose pour N=100.
N = 100;
x, sol_exact = calcul_sol_exact(N);
v = vcat([0],calcul_v_Thomas(N),[0]);
p = plot(x, sol_exact, label = "Solution exacte", ylabel = "y", xlabel = "x", color=:red, title = "N = 100")
plot!(p, x, v, color=:green, label="Solution approchée")

#On trace la valeur de erreur_approx pour N allant de 3 à 100.
y = [calcul_erreur_approx(N) for N in 3:100];
plot(3:100, y, title = "Erreur en fonction de N")

#Le graphique ci-dessous permet de constater que la diminution de l'erreur est inversement quadratique en fonction de N. 
z = [1/sqrt(t) for t in y];
plot(3:100, z)
