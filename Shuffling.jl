"""

Figure 1 of "A Fundamental Tradeoff Between Computation and
Communication in Distributed Computing".

"""
function coded_shuffling_tradeoff_plot(K=10, Q=10, N=2520)
    r = 1:K
    uncoded = 1 .- r./K
    coded = 1 ./ r .* (1 .- r./K)
    plt.plot(r, uncoded, "-o", label="Uncoded")
    plt.plot(r, coded, "-s", label="Coded")
    plt.grid()
    plt.xlim(1, 10)
    plt.ylim(0, 1)
    plt.legend()
    plt.xlabel("Repetition factor")
    plt.ylabel("Communication load")
    # savetikz("./figures/ShufflingTradeoff/ShufflingTradeoff.tex")
end
