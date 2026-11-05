using Plots
plotlyjs()

w = collect(range(0.0, 100.0, length=2001))

function draw_system(name, curves)
    p1 = plot(
        title="ЛАЧХ: $name",
        xlabel="Частота ω (рад/с)",
        ylabel="L(ω), дБ",
        xlims=(0, 100),
        lw=2,
        grid=true,
        legend=false
    )

    p2 = plot(
        title="ФЧХ: $name",
        xlabel="Частота ω (рад/с)",
        ylabel="Фаза (рад)",
        xlims=(0, 100),
        lw=2,
        grid=true,
        legend=false
    )

    for (i, c) in enumerate(curves)
        A = copy(c.mag)

        A[.!isfinite.(A) .| (A .<= 0)] .= NaN

        L = 20 .* log10.(A)
        plot!(p1, w, L, label=c.label, lw=2)
        plot!(p2, w, c.phase, label=false, lw=2)

        if hasproperty(c, :cutoff) && c.cutoff !== nothing
            wc = c.cutoff

            level = occursin("H<sub>4</sub>", name) ? 20log10(2) - 3.0103 : -3.0103

            hline!(p1, [level],
                linestyle=:dash, color=:gray,
                label=i == 1 ? "Уровень −3 дБ" : false)

            vline!(p1, [wc],
                linestyle=:dot, color=:red, label=false)

            scatter!(p1, [wc], [level],
                color=:red, markersize=5, label=false)

            annotate!(p1, wc, level,
                text("ωс = $(round(wc, digits=2))", 9, :left))
        end
    end

    p = plot(
        p1, p2,
        layout=(1, 2),
        size=(1100, 450),
        legend=:topright
    )

    display(p)
end

draw_system("H<sub>1</sub>", [
    (label="H1", mag=ones(length(w)),
     phase=zeros(length(w)))
])

draw_system("H<sub>2</sub>", [
    (label="H2", mag=2 ./ w,
     phase=fill(-π/2, length(w)))
])

draw_system("H<sub>3</sub>", [
    (label="H3", mag=0.5 .* w,
     phase=fill(π/2, length(w)))
])

draw_system("H<sub>4</sub>", [
    (label="H4",
     mag=w ./ sqrt.(1 .+ (0.5 .* w).^2),
     phase=π/2 .- atan.(0.5 .* w),
     cutoff=2.0)
])

draw_system("H<sub>5</sub>", [
    (label="H5",
     mag=1 ./ sqrt.(1 .+ (0.5 .* w).^2),
     phase=-atan.(0.5 .* w),
     cutoff=2.0)
])

function H6_curve(a)
    x = 1 .- w.^2 ./ 4
    y = a .* w ./ 2

    A = 1 ./ sqrt.(x.^2 .+ y.^2)

    if a == 0
        phase = [
            wi == 2 ? NaN : (wi < 2 ? 0.0 : -π)
            for wi in w
        ]
    else
        phase = -atan.(y, x)
    end

    u = (2 - a^2 + sqrt((a^2 - 2)^2 + 4)) / 2
    wc = 2 * sqrt(u)

    return (label="a = $a", mag=A, phase=phase, cutoff=wc)
end

draw_system("H<sub>6</sub>", [
    H6_curve(0.0),
    H6_curve(0.5),
    H6_curve(1.5)
])

draw_system("H<sub>7</sub>", [
    (label="H7",
     mag=ones(length(w)),
     phase=-2 .* w)
])