using ControlSystems
using Plots
using LinearAlgebra

num_s = [1.0]
den_s = [1.0, 3.2361, 5.2361, 5.2361, 3.2361, 1.0]
sys_s = tf(num_s, den_s)

println("- Пункты 1 и 2")
println("H(s) = ", sys_s)
z_s = tzeros(sys_s)
p_s = poles(sys_s)
println("Нули H(s): ", z_s)
println("Полюса H(s): ", p_s)
println("Порядок системы: ", length(p_s))
println("Система устойчива: ", all(real.(p_s) .< 0))

println("\n- Пункт 3")
# Группируем корни: 1 вещественный и 2 пары комплексно-сопряженных
real_p = [p for p in p_s if abs(imag(p)) < 1e-5][1]
comp_p1 = [p for p in p_s if imag(p) > 1e-5][1]
comp_p2 = [p for p in p_s if imag(p) > 1e-5][2]

sys_s1 = tf(zpk(Float64[], [real_p], 1.0))
sys_s2 = tf(zpk(Float64[], [comp_p1, conj(comp_p1)], 1.0))
sys_s3 = tf(zpk(Float64[], [comp_p2, conj(comp_p2)], 1.0))

println("Звено 1: ", sys_s1)
println("Звено 2: ", sys_s2)
println("Звено 3: ", sys_s3)

Ts = 0.1
sys_z = c2d(sys_s, Ts, :tustin)

println("\n- Пункт 4")
println("H(z) = ", sys_z)

println("\n- Пункт 5: Нули и полюса H(z)")
z_z = tzeros(sys_z)
p_z = poles(sys_z)
println("Нули H(z): ", z_z)
println("Полюса H(z): ", p_z)
println("Система устойчива: ", all(abs.(p_z) .< 1))

println("\n- Пункт 6")
num_z = numvec(sys_z)[1]   
den_z = denvec(sys_z)[1]   

# Нормируем на первый коэффициент знаменателя (a0)
b = num_z ./ den_z[1]
a = den_z ./ den_z[1]

println("y[n] = ", b[1], "*x[n] + ", b[2], "*x[n-1] + ", b[3], "*x[n-2] + ", b[4], "*x[n-3] + ", b[5], "*x[n-4] + ", b[6], "*x[n-5]")
println("       - (", a[2], "*y[n-1] + ", a[3], "*y[n-2] + ", a[4], "*y[n-3] + ", a[5], "*y[n-4] + ", a[6], "*y[n-5])")

println("\n- Пункт 7")
sys_z1 = c2d(sys_s1, Ts, :tustin)
sys_z2 = c2d(sys_s2, Ts, :tustin)
sys_z3 = c2d(sys_s3, Ts, :tustin)

# Каскадное соединение
sys_z_cascade = sys_z1 * sys_z2 * sys_z3

println("Звено 1 в Z: ", sys_z1)
println("Звено 2 в Z: ", sys_z2)
println("Звено 3 в Z: ", sys_z3)

function compute_impulse(b::Vector, a::Vector, N::Int)
    x = zeros(N)
    x[1] = 1.0  # Единичный импульс
    y = zeros(N)
    for n in 1:N
        # Сумма по числителю
        for i in 1:length(b)
            if n - i + 1 >= 1
                y[n] += b[i] * x[n - i + 1]
            end
        end
        # Сумма по знаменателю (рекурсивная часть)
        for i in 2:length(a)
            if n - i + 1 >= 1
                y[n] -= a[i] * y[n - i + 1]
            end
        end
    end
    return y
end

t = 0:Ts:10
N = length(t)

y_full = compute_impulse(b, a, N)

num_c = numvec(sys_z_cascade)[1]
den_c = denvec(sys_z_cascade)[1]
b_c = num_c ./ den_c[1]
a_c = den_c ./ den_c[1]
y_cascade = compute_impulse(b_c, a_c, N)

max_error = maximum(abs.(y_full .- y_cascade))
println("Максимальная ошибка между полной моделью и каскадом: ", max_error)

# Построение графиков
plt = plot(t, y_full, label="Полная система H(z)", linewidth=3, linestyle=:solid)
plot!(plt, t, y_cascade, label="Каскад звеньев", linewidth=2, linestyle=:dash, color=:red)
title!(plt, "Сравнение импульсных характеристик")
xlabel!(plt, "Время (с)")
ylabel!(plt, "Амплитуда")
display(plt)

println("\n- Пункт 8")
K = 1.0
T = 0.5
a_vals = [0.0, 0.5, 1.5]

H1 = tf([K], [1.0])
H2 = tf([K], [T, 0.0])
H3 = tf([T, 0.0], [1.0])
H4 = tf([1.0, 0.0], [T, 1.0])
H5 = tf([K], [T, 1.0])
H6a = tf([K], [T^2, a_vals[1]*T, 1.0])
H6b = tf([K], [T^2, a_vals[2]*T, 1.0])
H6c = tf([K], [T^2, a_vals[3]*T, 1.0])
H7 = pade(2.0, 2)

H_cont = Dict("H1"=>H1, "H2"=>H2, "H3"=>H3, "H4"=>H4, "H5"=>H5, "H6a"=>H6a, "H6b"=>H6b, "H6c"=>H6c, "H7"=>H7)

m = 100
H1_z_12 = tf([1.0], [1.0])
H2_z_12 = tf([1.0, 0.0], [m, -m])
H3_z_12 = tf([m, -m], [1.0, 0.0])
H5_z_12 = tf([m], [1.0, m+1])
H6_z_12 = tf([1.0], [m^2, m, 1.0])

H_disc_12 = Dict("H1"=>H1_z_12, "H2"=>H2_z_12, "H3"=>H3_z_12, "H5"=>H5_z_12, "H6a"=>H6_z_12, "H6b"=>H6_z_12, "H6c"=>H6_z_12)

for (name, sys) in H_cont
    try
        sys_z_tustin = c2d(sys, Ts, :tustin)
        if haskey(H_disc_12, name)
            println("\nСравнение для $name:")
            println("  Билинейное: ", sys_z_tustin)
            println("  Из пункта 1.2:       ", H_disc_12[name])
            println("  Полюса Tustin: ", poles(sys_z_tustin))
            println("  Полюса 1.2:    ", poles(H_disc_12[name]))
        else
            println("\n$name (билин): ", sys_z_tustin)
        end
    catch e
        println("\nОшибка при преобразовании $name: ", e)
    end
end