using PyPlot
using Distributions

function distribution_plot()
    pred = Normal(5, 2)
    obs = Normal(9, 1)
    opt = Normal(7.5, 0.7)

    x = range(0, 15, length=200)
    plt.fill_between(x, pdf.(pred, x), label="Predicted", alpha=0.7, edgecolor="k")
    plt.fill_between(x, pdf.(obs, x), label="Observed", alpha=0.7, edgecolor="k")
    plt.fill_between(x, pdf.(opt, x), label="Combined", alpha=0.7, edgecolor="k")

    # plt.plot(x, pdf.(obs, x), label="Observed")
    # plt.plot(x, pdf.(opt, x), label="Combined")

    plt.xlim(0, 12)
    plt.ylim(0, 0.6)

    plt.legend()
    plt.grid()

    plt.xlabel("Position")
    plt.ylabel("PDF")
    savetikz("./figures/KalmanDistribution/KalmanDistribution.tex")
end
