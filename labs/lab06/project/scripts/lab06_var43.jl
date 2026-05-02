# lab06_var43.jl — Лабораторная работа №6 (SIR-модель эпидемии), вариант 43
# Автор: Митрофанов Тимур Александрович, НПИбд-02-23
#
# Параметры варианта 43: N = 5505, I0 = 45, R0 = 3, S0 = N - I0 - R0 = 5457.
# Рассматриваются два случая критического порога I*:
#   Случай 1: I* = 50  ->  I(0) = 45 ≤ I*   (эпидемия не разгорается)
#   Случай 2: I* = 30  ->  I(0) = 45 > I*   (эпидемия разгорается)

using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))

# === Параметры варианта 43 ===
const N      = 5505
const I0     = 45.0
const R0_init = 3.0
const S0     = N - I0 - R0_init     # = 5457
const α      = 0.01                  # коэффициент заболеваемости
const β      = 0.02                  # коэффициент выздоровления
const tspan  = (0.0, 200.0)

# === Правая часть системы ОДУ с критическим порогом I* ===
# Если I(t) ≤ I* — заражение не происходит (dS/dt = 0, dI/dt = -βI).
# Если I(t) > I* — заражение идёт (dS/dt = -αS, dI/dt = αS - βI).
function sir!(du, u, p, t)
    S, I, R = u[1], u[2], u[3]
    I_star  = p[1]
    if I > I_star
        du[1] = -α * S
        du[2] =  α * S - β * I
    else
        du[1] = 0.0
        du[2] = -β * I
    end
    du[3] = β * I
end

u0 = [S0, I0, R0_init]

# === Случай 1: I(0) ≤ I*  (I* = 50) ===
prob1 = ODEProblem(sir!, u0, tspan, [50.0])
sol1  = solve(prob1, Tsit5(), saveat = 0.5)

t1  = sol1.t
S1  = [u[1] for u in sol1.u]
I1  = [u[2] for u in sol1.u]
R1  = [u[3] for u in sol1.u]

p1 = plot(t1, S1, label = "S(t)", color = :blue,  lw = 1.7)
plot!(p1, t1, I1, label = "I(t)", color = :green, lw = 1.7)
plot!(p1, t1, R1, label = "R(t)", color = :red,   lw = 1.7)
hline!(p1, [50.0], label = "I* = 50", color = :grey, ls = :dot)
title!(p1, "Случай 1: I(0) ≤ I*  (I* = 50)")
xlabel!(p1, "Время t");  ylabel!(p1, "Численность")
savefig(p1, plotsdir(script_name, "case1_I0_le_Istar.png"))

# === Случай 2: I(0) > I*  (I* = 30) ===
prob2 = ODEProblem(sir!, u0, tspan, [30.0])
sol2  = solve(prob2, Tsit5(), saveat = 0.5)

t2  = sol2.t
S2  = [u[1] for u in sol2.u]
I2  = [u[2] for u in sol2.u]
R2  = [u[3] for u in sol2.u]

p2 = plot(t2, S2, label = "S(t)", color = :blue,  lw = 1.7)
plot!(p2, t2, I2, label = "I(t)", color = :green, lw = 1.7)
plot!(p2, t2, R2, label = "R(t)", color = :red,   lw = 1.7)
hline!(p2, [30.0], label = "I* = 30", color = :grey, ls = :dot)
title!(p2, "Случай 2: I(0) > I*  (I* = 30)")
xlabel!(p2, "Время t");  ylabel!(p2, "Численность")
savefig(p2, plotsdir(script_name, "case2_I0_gt_Istar.png"))

# === Сводная картинка ===
p_sum = plot(p1, p2, layout = (1, 2), size = (1300, 500))
savefig(p_sum, plotsdir(script_name, "summary.png"))

# === Информация в консоль ===
println("Параметры варианта 43:  N = $N, I0 = $I0, R0 = $R0_init, S0 = $S0")
println("                       α = $α, β = $β")
println("Случай 1 (I* = 50): max I = ", maximum(I1), ",  R(end) = ", R1[end])
println("Случай 2 (I* = 30): max I = ", maximum(I2), ",  R(end) = ", R2[end])
println("Графики сохранены в ", plotsdir(script_name))
