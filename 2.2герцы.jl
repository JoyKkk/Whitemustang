function draw_frequency_hz(fs)
    f = w .* fs ./ (2π)

    p1 = plot(
        title="ЛАЧХ, fs = $fs Гц",
        xlabel="",
        ylabel="L(f), дБ",
        grid=true,
        legend=:outerright,
        size=(1100, 450)
    )

    p2 = plot(
        title="ФЧХ, fs = $fs Гц",
        xlabel="Частота f (Гц)",
        ylabel="Фаза (рад)",
        grid=true,
        legend=:outerright
    )

    for (i, name) in enumerate(names)
        H = [get_H(wi)[i] for wi in w]

        A = abs.(H)
        L = 20 .* log10.(A)
        phase = unwrap_phase(angle.(H))

        plot!(p1, f, L, label=name, lw=2)
        plot!(p2, f, phase, label="", lw=2)
    end

    display(plot(p1, p2, layout=(2, 1), size=(1100, 800)))
end

draw_frequency_hz(100)
draw_frequency_hz(10000)