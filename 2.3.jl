using SpecialFunctions
using EngeeDSP
import Plots

fs = 200
T = 2
N = fs*T
t = (0:N-1)/fs

x = sin.(2*pi*t) + 0.5*sin.(2*pi*10*t) + 0.3*sin.(2*pi*30*t)

function spectrum(x, fs, nfft)
    N = length(x)
    X = EngeeDSP.Functions.fft(x, nfft)
    A = abs.(X)/N
    A[2:nfft÷2] .*= 2
    f = (0:nfft÷2)*fs/nfft
    return f, A[1:nfft÷2+1]
end

f, A = spectrum(x, fs, N)

p1 = Plots.plot(f, A, seriestype=:stem,
    xlim=(0,50), ylim=(0,1.1),
    xlabel="Частота, Гц", ylabel="Амплитуда",
    title="Амплитудный спектр сигнала",
    color = "#ff00aa",
    legend=false,
    linewidth=2,
    markersize=4)
display(p1)

x125 = sin.(2*pi*1.25*t) +
       0.5*sin.(2*pi*10*t) +
       0.3*sin.(2*pi*30*t)

f125, A125 = spectrum(x125, fs, N)

p2 = Plots.plot(f125, A125, seriestype=:stem,
    xlim=(0,50), ylim=(0,1.1),
    xlabel="Частота, Гц", ylabel="Амплитуда",
    title="Амплитудный спектр (1.25 Гц у 1-й комп.)",
    color = "#ff00aa",
    legend=false,
    linewidth=2,
    markersize=4)
display(p2)

nfft = 4096
f0, A0 = spectrum(x, fs, nfft)

p3 = Plots.plot(f0, A0, seriestype=:stem,
    xlim=(0,50), ylim=(0,1.1),
    xlabel="Частота, Гц", ylabel="Амплитуда",
    title="Спектр с дополнением нулями (4096)",
    color = "#ff00aa",
    legend=false)
display(p3)

T8 = 8
N8 = fs*T8
t8 = (0:N8-1)/fs

x8 = sin.(2*pi*t8) + 0.5*sin.(2*pi*10*t8) + 0.3*sin.(2*pi*30*t8)

f8, A8 = spectrum(x8, fs, 16384)

p4 = Plots.plot(f8, A8, seriestype=:stem,
    xlim=(0,50), ylim=(0,1.1),
    xlabel="Частота, Гц", ylabel="Амплитуда",
    title="Сигнал 8 секунд и с 16384 нулями",
    color = "#ff00aa",
    legend=false)
display(p4)


N = length(x)
n = 0:N-1

wh = 0.5 .* (1 .- cos.(2*pi*n/(N-1)))

wb = 0.42 .- 0.5*cos.(2*pi*n/(N-1)) .+
      0.08*cos.(4*pi*n/(N-1))

wk = [besseli(0,5*sqrt(1-(2*i/(N-1)-1)^2))/
      besseli(0,5) for i in n]

fh, Ah = spectrum(x .* wh, fs, 4096)
fb, Ab = spectrum(x .* wb, fs, 4096)
fk, Ak = spectrum(x .* wk, fs, 4096)

Ah .*= N/sum(wh)
Ab .*= N/sum(wb)
Ak .*= N/sum(wk)

p5 = Plots.plot(fh, Ah, seriestype=:stem,
    xlim=(0,50), ylim=(0,1.1),
    xlabel="Частота, Гц", ylabel="Амплитуда",
    title="Спектр через окно Ханна",
    label="Ханна",
    color = "#ff00aa")
display(p5)

p7 = Plots.plot(fb, Ab, seriestype=:stem,
    xlim=(0,50), ylim=(0,1.1),
    xlabel="Частота, Гц", ylabel="Амплитуда",
    title="Спектр через окно Блэкмана",
    label="Блэкман",
    color = "#ff00aa")
display(p7)

p8 = Plots.plot(fk, Ak, seriestype=:stem,
    xlim=(0,50), ylim=(0,1.1),
    xlabel="Частота, Гц", ylabel="Амплитуда",
    title="Спектр через окно Кайзера",
    label="Кайзер",
    color = "#ff00aa")
display(p8)

p6 = Plots.plot(fh, Ah, label="Ханна",
    xlim=(0,50), ylim=(0,1.1),
    xlabel="Частота, Гц",
    ylabel="Амплитуда",
    title="Сравнение оконных функций",
    color = "#ff00aa")

Plots.plot!(p6, fb, Ab, label="Блэкман")
Plots.plot!(p6, fk, Ak, label="Кайзер")

display(p6)