# lab07_var43.jl — Лабораторная работа №7 (Эффективность рекламы), вариант 43
# Автор: Митрофанов Тимур Александрович, НПИбд-02-23
#
# Модель: dn/dt = (α1(t) + α2(t)·n)·(N - n).
# Параметры варианта 43: N = 3310, n(0) = 22.
# Три случая:
#   1) α1 = 0.211,         α2 = 0.000011    — преобладает платная реклама (Мальтус)
#   2) α1 = 0.0000311,     α2 = 0.21        — преобладает «сарафан» (логистический рост)
#   3) α1 = 0.511·sin(t),  α2 = 0.311·sin(t) — периодические коэффициенты

using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))

# === Параметры варианта 43 ===
const N  = 3310.0
const n0 = 22.0

# === Правая часть ОДУ с переменными коэффициентами через параметр ===
function rhs!(du, u, p, t)
    a1_fn, a2_fn = p
    n = u[1]
    du[1] = (a1_fn(t) + a2_fn(t) * n) * (N - n)
end

# === Хелпер: считаем и сохраняем графики n(t) и dn/dt ===
function run_case(a1_fn, a2_fn, t_end::Float64, label::String, fname::String;
                  saveat::Float64 = 0.01, max_marker::Bool=false)
    prob = ODEProblem(rhs!, [n0], (0.0, t_end), (a1_fn, a2_fn))
    sol  = solve(prob, Tsit5(), saveat = saveat)
    t  = sol.t
    nv = [u[1] for u in sol.u]
    d  = [(a1_fn(tt) + a2_fn(tt)*nn) * (N - nn) for (tt, nn) in zip(t, nv)]

    p1 = plot(t, nv, label = "n(t)", color = :blue, lw = 1.7)
    hline!(p1, [N], label = "N = $(Int(N))", color = :grey, ls = :dot)
    xlabel!(p1, "Время t");  ylabel!(p1, "n(t)")
    title!(p1, "$label — n(t)")

    p2 = plot(t, d, label = "dn/dt", color = :red, lw = 1.7)
    xlabel!(p2, "Время t");  ylabel!(p2, "dn/dt")
    title!(p2, "$label — скорость dn/dt")
    if max_marker
        idx = argmax(d)
        scatter!(p2, [t[idx]], [d[idx]], color = :black, ms = 6, label = "max")
        scatter!(p1, [t[idx]], [nv[idx]], color = :black, ms = 6, label = "n(t*)")
        println("$label: max dn/dt = $(d[idx]) при t = $(t[idx]), n = $(nv[idx])")
    end

    p_pair = plot(p1, p2, layout = (1, 2), size = (1200, 450))
    savefig(p_pair, plotsdir(script_name, fname))
    return t, nv, d
end

# === Случай 1: преобладает α1 (Мальтус) ===
run_case(
    t -> 0.211, t -> 0.000011,
    30.0, "Случай 1: α₁=0.211, α₂=0.000011",
    "case1.png")

# === Случай 2: преобладает α2 (логистическая кривая) ===
# Динамика очень быстрая: интегрируем на маленьком t_end с малым saveat.
run_case(
    t -> 0.0000311, t -> 0.21,
    0.05, "Случай 2: α₁=0.0000311, α₂=0.21",
    "case2_zoom.png";
    saveat = 0.0001, max_marker = true)

# === Случай 3: периодические коэффициенты ===
run_case(
    t -> 0.511*sin(t), t -> 0.311*sin(t),
    30.0, "Случай 3: α₁=0.511·sin(t), α₂=0.311·sin(t)",
    "case3.png";
    saveat = 0.005)

# === Случай 4: только α1 (только платная реклама) ===
run_case(
    t -> 0.211, t -> 0.0,
    30.0, "Только платная реклама (α₂=0)",
    "case4_only_paid.png")

# === Случай 5: только α2 (только «сарафан») ===
run_case(
    t -> 0.0, t -> 0.21,
    0.05, "Только «сарафан» (α₁=0)",
    "case5_only_word.png";
    saveat = 0.0001)

println("Графики сохранены в: ", plotsdir(script_name))
