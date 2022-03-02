using SymPy
using NLsolve
using PyPlot
using Parameters

@vars M C
@vars a g y r m N₀ Nₜ c

f(M, C) = M * ( a * C * ( Nₜ / (N₀ + Nₜ) ) - g / (M + (1 - M - C)) + y * (1 - M - C))
h(M, C) = C * ( r * (1 - c * Nₜ) * (1 - M - C) - m - a * M * ( Nₜ / (N₀ + Nₜ) )) 

SymPy.solve(f(M,C), M)
SymPy.solve(f(M,C), C)
SymPy.solve(h(M,C), M)
SymPy.solve(h(M,C), C)

function coral_model!(F, x, par)
    @unpack a, g, y, r, m, N₀, Nₜ, c = par
    F[1] = x[1] * ( r * (1 - c * Nₜ) * (1 - x[2] - x[1]) - m - a * x[2] * ( Nₜ / (N₀ + Nₜ) ))
    F[2] = x[2] * ( a * x[1] * ( Nₜ / (N₀ + Nₜ) ) - g / (x[2] + (1 - x[2] - x[1])) + y * (1 - x[2] - x[1]))
end

@with_kw mutable struct CoralPar
    a = 2.3
    g = 0.4
    y = 0.7
    r = 1.8
    m = 0.15
    N₀ = 0.5
    Nₜ = 0.5
    c = 0.25
end

function macro_isocline(C, par)
    @unpack a, g, y, r, m, N₀, Nₜ, c = par
    return 1 - C + (a * C / y) * (Nₜ/ (N₀+Nₜ)) - (g / (y * (1 - C)))
end

function coral_isocline(C, par)
    @unpack a, g, y, r, m, N₀, Nₜ, c = par
    return ((N₀+Nₜ) * (m - r * (1 - c * Nₜ) * (1-C))) / (r * (N₀ + Nₜ) * ( c * Nₜ -1 ) - a * Nₜ)
end

let 
    Crange = 0.0:0.001:1.0
    macrodata = [macro_isocline(C, CoralPar(Nₜ=0.7)) for C in Crange]
    coraldata = [coral_isocline(C, CoralPar(Nₜ=0.7)) for C in Crange]
    isoclines = figure()
    plot(Crange, macrodata)
    plot(Crange, coraldata)
    xlim(0.0,1.0)
    ylim(0.0,1.0)
    return isoclines
end

function coral_iso_zero(M, par)
    @unpack a, g, y, r, m, N₀, Nₜ, c = par
    (-M*N₀*Nₜ*c*r + M*N₀*r - M*(Nₜ^2)*c*r + M*Nₜ*a + M*Nₜ*r + N₀*Nₜ*c*r + N₀*m - N₀*r + (Nₜ^2)*c*r + Nₜ*m - Nₜ*r)/(r*(N₀*Nₜ*c - N₀ + (Nₜ^2)*c - Nₜ))
end


function coral_equil(par, stable_guess, unstable_guess)
    coralaxis = coral_iso_zero(0.0, par)
    macroaxis = macro_isocline(0.0, par)
    if stable_guess==0.0
        return [coralaxis, macroaxis]
    else
        stableinterior = nlsolve((F, x) -> coral_model!(F, x, par), stable_guess, autodiff = :forward).zero
        unstableinterior = nlsolve((F, x) -> coral_model!(F, x, par), unstable_guess, autodiff = :forward).zero
        return [coralaxis, macroaxis, stableinterior, unstableinterior]
    end
end

coral_equil(CoralPar(), [0.2; 0.4], [0.5; 0.3])[4]

coral_equil(CoralPar(Nₜ=0.7), [0.2; 0.4], [0.5; 0.3])[4]

coral_equil(CoralPar(Nₜ=0.45), 0.0, 0.0)[1]