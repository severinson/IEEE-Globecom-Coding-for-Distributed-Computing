function Bj(j; K, q, μ)
    μq = floor(Int, μ*q)
    rv = binomial(q-1, j) * binomial(K-q, μq-j)
    rv /= q/K * binomial(K, μq)
    return rv
end

function sq(;K, q, μ)
    μq = floor(Int, μ*q)
    μ_bar = μq / q
    c = 0.0
    for s in μq:-1:1
        if sum(Bj(j, K=K, q=q, μ=μ) for j in s:μq) > 1 - μ_bar
            println("$q <= $s <= $K")
            return max(min(μq, s+1), 1)
        end
    end
    println("$q <= $μq <= $K")
    return 1
end

"""

Computational delay of the unified scheme. Eq. (8) of "A Unified
Coding Framework for Distributed Computing with Straggling Servers".

"""
function computational_delay(q; K, μ, N)
    rv = ExponentialOrder(1.0, K, q)
    # return μ*N * (1 + mean(rv))
    return μ * (1 + mean(rv))
end

"""

Communication load of the unified scheme. Eq. (9) of "A Unified Coding
Framework for Distributed Computing with Straggling Servers".

"""
function communication_load(q; K, μ, N)
    μq = floor(Int, μ*q)
    μ_bar = μq / q
    s = sq(K=K, q=q, μ=μ)
    rv = sum(Bj(j, K=K, q=q, μ=μ)/j for j in s:μq)
    rv += min(
        1 - μ_bar - sum(Bj(j, K=K, q=q, μ=μ) for j in s:μq),
        Bj(s-1, K=K, q=q, μ=μ) / (s-1),
    )
    # rv *= N
    return rv
end

"""

Lower bound on the communication load of the unified scheme. Eq. (11)
of "A Unified Coding Framework for Distributed Computing with
Straggling Servers".

"""
function communication_load_lower(q; μ, N)
    # return N*maximum(q*(1 - min(t*μ, 1)) / (ceil(q/t) * (q-t)) for t in 1:q-1)
    return maximum(q*(1 - min(t*μ, 1)) / (ceil(q/t) * (q-t)) for t in 1:q-1)
end


"""

Figure 1 of "A Unified Coding Framework for Distributed Computing with
Straggling Servers"

"""
function delay_load_tradeoff_plot(K=18, N=180, μ=1/3)
    q = round(Int, 1/μ):round(Int, 1/μ):K-1
    D = computational_delay.(q, K=K, μ=μ, N=N)
    L = communication_load.(q, K=K, μ=μ, N=N)
    Ll = communication_load_lower.(q, μ=μ, N=N)
    println(Ll)
    # println(L)
    plt.plot(D, L)
    plt.plot(D, Ll)
    plt.grid()
    # plt.ylim(20, 120)
    # plt.xlim(60, 160)
    plt.ylabel("Communication load")
    plt.xlabel("Computational latency")
    # savetikz("./figures/UnifiedTradeoff/UnifiedTradeoff.tex")
    return
end
