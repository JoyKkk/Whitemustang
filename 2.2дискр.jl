using Plots
plotlyjs()

w = collect(range(0.001, π, length=2001))

function get_H(w)
    z = exp(im * w)

    H1 = 1
    H2 = z / (100 * (z - 1))
    H3 = 100 * (z - 1) / z
    H5 = 100 / (z + 101)
    H6 = 1 / (10000 * z^2 + 100 * z + 1)
    H7 = z^(-100)

    return [H1, H2, H3, H5, H6, H7]
end

function unwrap_phase(phi)
    result = copy(phi)
    offset = 0.0

    for i in 2:length(result)
        delta = phi[i] - phi[i-1]

        if delta > π
            offset -= 2π
        elseif delta < -π
            offset += 2π
        end

        result[i] += offset
    end

    return result
end

function draw_system(name, H, w)
    A = abs.(H)
    L = name == "H<sub>7</sub>" ? zeros(length(w)) : 20 .* log10.(abs.(H))
    phase = unwrap_phase(angle.(H))

    p1 = plot(
        w, L,
        title="ЛАЧХ: $name",
        xlabel="Частота ω (рад)",
        ylabel="L(ω), дБ",
        lw=2,
        grid=true,
        legend=false,
        color = "#ff7514"
    )

    p2 = plot(
        w, phase,
        title="ФЧХ: $name",
        xlabel="Частота ω (рад)",
        ylabel="Фаза (рад)",
        lw=2,
        grid=true,
        legend=false,
        color = "#ff7514"
    )

    display(plot(p1, p2, layout=(1, 2), size=(1100, 400)))
end

names = ["H<sub>1</sub>", "H<sub>2</sub>", "H<sub>3</sub>", "H<sub>5</sub>", "H<sub>6</sub>", "H<sub>7</sub>"]

for (i, name) in enumerate(names)
    H = [get_H(wi)[i] for wi in w]
    draw_system(name, H, w)
end


# ===============================
A3 = abs.([get_H(wi)[3] for wi in w])

wc = π / 2
Ac = 200 / sqrt(2)

p = plot(
    w, A3,
    title="АЧХ H<sub>3</sub>, срез",
    xlabel="Частота ω (рад)",
    ylabel="A(ω)",
    lw=2,
    label="H3",
    grid=true
)

hline!(p, [Ac],
    linestyle=:dash,
    label="Уровень −3 дБ"
)

scatter!(p, [wc], [Ac],
    markersize=5,
    label="Частота среза"
)

display(p)