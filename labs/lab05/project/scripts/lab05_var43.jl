using DrWatson


using DifferentialEquations
using Plots


mkpath(plotsdir("script"))
mkpath(datadir("script"))
 

x0 = 3.0      
y0 = 8.0      
u0 = [x0, y0]
tspan = (0.0, 200.0)
 

a = 0.19      
b = 0.026     
c = 0.18      
d = 0.032     
 

function predator_prey!(du, u, p, t)
    x = u[1]   # хищники
    y = u[2]   # жертвы
    du[1] = -a * x + b * x * y
    du[2] =  c * y - d * x * y
end
 
prob = ODEProblem(predator_prey!, u0, tspan)
sol  = solve(prob, Tsit5(), dtmax = 0.08)
 
x_val    = sol[1, :]
y_val    = sol[2, :]
time_val = sol.t
 

p1 = plot(time_val, x_val, label="Хищники (x)", color=:red, lw=1.5)
plot!(p1, time_val, y_val,  label="Жертвы (y)", color=:blue, lw=1.5)
title!(p1, "Изменение численности хищников и жертв во времени")
xlabel!(p1, "Время t")
ylabel!(p1, "Численность популяции")
 

x_eq = c / d          
y_eq = a / b          
println("Стационарная точка: x_eq = ", x_eq, ", y_eq = ", y_eq)
 

p2 = plot(y_val, x_val, label="Фазовая траектория", color=:green, lw=1.5)
scatter!(p2, [y0],   [x0],   label="Начальная точка",    color=:black, ms=6)
scatter!(p2, [y_eq], [x_eq], label="Стационарная точка", color=:red,   ms=6, marker=:diamond)
title!(p2, "Фазовый портрет: хищники от жертв")
xlabel!(p2, "Жертвы (y)")
ylabel!(p2, "Хищники (x)")
 

solutions = plot(p1, p2, layout=(2,1), size=(900, 800))
savefig(solutions, plotsdir("script", "lab05_var43.png"))
println("Графики сохранены: ", plotsdir("script", "lab05_var43.png"))
