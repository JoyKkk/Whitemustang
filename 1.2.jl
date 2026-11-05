using ControlSystems
using Plots

function main()
    m = 100
    K = 1
    Ts = 1.0
    z = tf("z", Ts)

    H1 = K * z / z
    H2 = (K * z) / (m * (z - 1))
    H3 = K * m * (z - 1) / z
    H5 = K * m / (z + m + 1)
    H6 = K / (m^2 * z^2 + m * z + 1)
    H7 = K / (z^m)

    systems = [H1, H2, H3, H5, H6, H7]
    names = ["H1", "H2", "H3", "H5", "H6", "H7"]

    for i in 1:length(systems)
        println("=== $(names[i]) ===")
        println("Нули:  ", tzero(systems[i]))
        println("Полюса: ", poles(systems[i]))
        println()
    end

    plots_list = []
    for i in 1:length(systems)
        zs = ComplexF64.(tzero(systems[i]))
        ps = ComplexF64.(poles(systems[i]))

        all_pts = vcat(zs, ps)
        r = isempty(all_pts) ? 1.5 : max(1.5, maximum(abs.(all_pts)) * 1.3)

        p = plot(title=names[i], grid=true, legend=false,
                 xlabel="Re(z)", ylabel="Im(z)",
                 xlim=(-r, r), ylim=(-r, r),
                 aspect_ratio=1, margin=8Plots.mm,
                 titlefontsize=11)

        θ = range(0, 2π, length=300)
        plot!(p, cos.(θ), sin.(θ),
              color=:gray, linestyle=:solid, linewidth=1, label="")

        if !isempty(zs)
            scatter!(p, real(zs), imag(zs),
                     markershape=:circle,
                     markercolor=:white, markerstrokecolor=:blue,
                     markersize=8, markerstrokewidth=1.5, label="")
        end
        if !isempty(ps)
            scatter!(p, real(ps), imag(ps),
                     markershape=:xcross,
                     markercolor=:red, markersize=6, label="")
        end

        push!(plots_list, p)
    end

    legend_plot = plot(
        xlim=(0, 1), ylim=(0, 1),
        framestyle=:none,
        xticks=false, yticks=false,
        legend=:inside,
        legendfontsize=12,
        showaxis=false,
    )
    scatter!(legend_plot, [-10.0], [-10.0],
             markershape=:circle,
             markercolor=:white, markerstrokecolor=:blue,
             markersize=8, markerstrokewidth=1.5,
             label="нули")
    scatter!(legend_plot, [-10.0], [-10.0],
             markershape=:xcross,
             markercolor=:red,
             markersize=8,
             label="полюса")

    l = @layout [grid(2, 3); b{0.12h}]
    display(plot(plots_list..., legend_plot, layout=l, size=(1400, 900)))
end

main()