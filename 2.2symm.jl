w2 = collect(range(-2π, 2π, length=4001))

function check_symmetry(name, index)
    H = [get_H(wi)[index] for wi in w2]

    A = index == 6 ? ones(length(w2)) : abs.(H)
    phase = unwrap_phase(angle.(H))

    p1 = plot(
        w2, A,
        title="АЧХ: $name",
        xlabel="ω (рад)",
        ylabel="A(ω)",
        lw=2,
        grid=true,
        color = "#ff4f00",
        legend=false
    )

    p2 = plot(
        w2, phase,
        title="ФЧХ: $name",
        xlabel="ω (рад)",
        ylabel="Фаза (рад)",
        lw=2,
        grid=true,
        color = "#ff4f00",
        legend=false
    )

    display(plot(p1, p2, layout=(1, 2), size=(1100, 400)))
end

for (i, name) in enumerate(names)
    check_symmetry(name, i)
end
