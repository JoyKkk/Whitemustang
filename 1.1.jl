using ControlSystems
using Plots

T = 0.5
K = 1.0

H1 = tf([K], [1.0])
H2 = tf([K], [T, 0.0])
H4 = tf([1.0, 0.0], [T, 1.0])
H5 = tf([K], [T, 1.0])
H6_1 = tf([K], [T^2, 0*T, 1.0])
H6_2 = tf([K], [T^2, 0.5*T, 1.0])
H6_3 = tf([K], [T^2, 1.5*T, 1.0])

systems = Dict(
    "H1 (Пропорц.)"   => H1,
    "H2 (Интегр.)"    => H2,
    "H4 (Реал. Дифф.)" => H4,
    "H5 (Апериод.)"   => H5,
    "H6 (a=0)"        => H6_1,
    "H6 (a=0.5)"      => H6_2,
    "H6 (a=1.5)"      => H6_3,
)

println("=== Устойчивость ===")
for name in sort(collect(keys(systems)))
    sys = systems[name]
    p = poles(sys)
    z = tzeros(sys)
    
    if isempty(p) && all(real.(z) .<= 0)
        status = "Устойчива (нет динамических полюсов)"
    elseif all(real.(p) .< 0)
        status = "Устойчива"
    elseif any(real.(p) .> 0)
        status = "НЕустойчива"
    else
        status = "На границе устойчивости"
    end
    
    println("$name:")
    println("  Нули: ", isempty(z) ? "нет" : z)
    println("  Полюса: ", isempty(p) ? "нет" : p)
    println("  Статус: $status\n")
end

plt = plot(title="Карта нулей и полюсов", xlabel="Re(s)", ylabel="Im(s)", grid=true)
colors = [:blue, :red, :green, :purple, :orange, :magenta, :cyan]

for (i, name) in enumerate(sort(collect(keys(systems))))
    sys = systems[name]
    # Для H7 (системы с задержкой) pzmap строит полюса рациональной части
    pzmap!(sys, label=name, color=colors[i], markersize=6)
end
display(plt)