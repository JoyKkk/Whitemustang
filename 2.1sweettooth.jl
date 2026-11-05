using FFTW
using Plots
plotlyjs()

fs = 1000.0
T = 4.0
t = 0:1/fs:T-1/fs
N = length(t)

f0 = 5.0
ϕ = π/4

x1 = sin.(2π * f0 .* t .+ ϕ)
x2 = cos.(2π * f0 .* t .+ ϕ)
x3 = exp.(im .* (2π .* f0 .* t .+ ϕ))

freqs = fftshift(fftfreq(N, fs))
get_spectrum(x) = fftshift(fft(x)) ./ N

X1 = get_spectrum(x1)
X2 = get_spectrum(x2)
X3 = get_spectrum(x3)

function plot_signal_spectrums(freqs, X, title_name)
    mask = -15.0 .<= freqs .<= 15.0

    f_plot = freqs[mask]
    X_plot = X[mask]

    A = abs.(X_plot)

    Phase = [
        A[i] > 1e-5 ? angle(X_plot[i]) : 0.0
        for i in 1:length(X_plot)
    ]

    p1 = stem(
        f_plot,
        real.(X_plot),
        imag.(X_plot),
        st=:line,
        lw=3,
        label="",
        title="Трёхмерный спектр: $title_name",
        xlabel="Частота (Гц)",
        ylabel="Вещественная (Re)",
        zlabel="Мнимая (Im)"
    )

    scatter!(
        p1,
        f_plot,
        real.(X_plot),
        imag.(X_plot),
        msize=3,
        label=""
    )

    p2 = plot(
        f_plot,
        A,
        lw=2,
        color=:blue,
        legend=false,
        title="Амплитудный спектр",
        xlabel="",
        ylabel="Амплитуда"
    )

    p3 = plot(
        f_plot,
        Phase,
        lw=2,
        color=:red,
        legend=false,
        title="Фазовый спектр",
        xlabel="Частота (Гц)",
        ylabel="Фаза (рад)"
    )

    display(
        plot(
            p1,
            plot(p2, p3, layout=(2,1)),
            layout=(1,2),
            size=(1000,500)
        )
    )
end

plot_signal_spectrums(freqs, X1, "sin(2πft + ϕ)")
plot_signal_spectrums(freqs, X2, "cos(2πft + ϕ)")
plot_signal_spectrums(freqs, X3, "exp(j2πft + ϕ)")