using Distributions
using PyPlot
using QuadGK
using LambertW
using Test
using PyCall
matplotlib2tikz = pyimport("matplotlib2tikz")

"""

savetikz(path; fig = PyPlot.gcf(), extra::Vector{String})`

"""
function savetikz(path; fig=PyPlot.gcf(), extra=[])
    if extra == []
        matplotlib2tikz.save(
            path, fig,
            figureheight = "\\figureheight",
            figurewidth = "\\figurewidth",
        )
    else
        matplotlib2tikz.save(
            path, fig,
            figureheight="\\figureheight",
            figurewidth="\\figurewidth",
            extra_tikzpicture_parameters=PyCall.pybuiltin("set")(extra),
        )
    end
end

"""

    savetable(xs, ys, path)

Save x and y data as a table that can be included in pgf plots.

"""
function savetable(xs, ys, path)
    open(path, "w") do f
        for (x, y) in zip(xs, ys)
            write(f, "$x $y\n")
        end
    end
end

"""
    ExponentialOrder(scale::Real, total::Int, order::Int)

Random variable representing the order-th largest value out of total
realizations of an exponential random variable with given scale.

"""
function ExponentialOrder(scale::Real, total::Int, order::Int)
    scale > 0 || throw(DomainError(scale, "scale must be positive"))
    total > 0 || throw(DomainError(scale, "total must be positive"))
    order > 0 || throw(DomainError(scale, "order must be positive"))
    order <= total || throw(DomainError(scale, "order must be <= total"))
    var = sum(1/(i^2)  for i=(total-order+1):total) * scale^2
    mean = sum(1/i for i=(total-order+1):total) * scale
    alpha = mean^2 / var # shape parameter
    theta = var / mean # scale
    return Gamma(alpha, theta)
end

function test_exporder()
    # test against values from p74 of A First Course in Order Statistics by Barry Arnold
    @test isapprox(round(mean(ExponentialOrder(1, 10, 10)), digits=6), 2.928968)
    @test isapprox(round(var(ExponentialOrder(1, 10, 10)), digits=6), 1.549768)
    @test isapprox(round(mean(ExponentialOrder(1, 10, 1)), digits=6), 0.1)
    @test isapprox(round(var(ExponentialOrder(1, 10, 1)), digits=6), 0.01)
    @test isapprox(round(mean(ExponentialOrder(1, 6, 4)), digits=6), 0.95)
    @test isapprox(round(var(ExponentialOrder(1, 6, 4)), digits=6), 0.241389)
end

function validate_exporder(n=6, k=4, β=2.0, nsamples=10000000)
    s = zeros(n)
    de = Exponential(β)
    samples = zeros(nsamples)
    for i in 1:nsamples
        s .= rand(de, n)
        sort!(s)
        samples[i] = s[k]
    end
    println(mean(samples), " / ", var(samples))
    plt.hist(samples, 200, density=true, cumulative=true)

    d = ExponentialOrder(β, n, k)
    println(mean(d), " / ", var(d))
    t = range(minimum(samples), maximum(samples), length=100)
    plt.plot(t, cdf.(d, t))
    # plt.ylim(0, 1)
    return
end

function delay_ccdf(n, k, β; l=k, t=range(0, β, length=1000))
    d = ExponentialOrder(β, n, k)
    return t, 1 .- cdf.(d, l.*t.-1/l)
end

"""

Plot the CCDF of the delay as a function of time (Fig. 4a of Lee).

"""
function delay_ccdf_plot(n=10, k=5, β=3.0)
    t, v = delay_ccdf(n, k, β)
    plt.semilogy(t, v, label="MDS \$(10, 5)\$")

    t, v = delay_ccdf(n, n, β)
    plt.semilogy(t, v, label="Uncoded")

    t, v = delay_ccdf(n, 1, β)
    plt.semilogy(t, v, label="Repetition")

    plt.legend()
    plt.xlim(0, 2)
    plt.ylim(1e-3, 1e0)
    plt.xlabel("\$t\$")
    plt.ylabel("\$\\Pr(t \\leq T)\$")
    plt.grid()
    ax = plt.axes()
    # ax.set_xticks(0:0.025:0.5, minor=true)
    plt.grid(true, which="major", ls="-")
    plt.grid(true, which="minor", ls=":")
    savetikz("./figures/delay_ccdf.tex")
    return
end

ET_uncoded(n, β) = (1+β*log(n))/n
ET_mds(n, k, β) = (1+β*log(n/(n-k)))/k

"""

Expected delay computations from Prop. 3 of Lee.

"""
function expected_delays(n=10, k=round(Int, n/2), β=1.0)
    println("Uncoded approx.\t\t", ET_uncoded(n, β))
    d = ExponentialOrder(β/n, n, n)
    println("Uncoded exact\t\t", mean(d)+1/n)

    println("MDS ($n, $k) approx.\t", ET_mds(n, k, β))
    d = ExponentialOrder(β/k, n, k)
    println("MDS ($n, $k) exact\t", mean(d)+1/k)
end

"""optimal expected MDS delay, multiplied by n ((13) of Lee)"""
γ(μ) = -lambertw(-exp(-μ-1), -1)/μ

"""

Plot the normalized optimal expected computing time as a function of μ
(Fig. 5a of Lee).

"""
function γ_plot(n=10)
    μs = 10.0.^range(-1, 1, length=100)
    βs = 1 ./ μs
    γs = γ.(μs)
    plt.semilogx(βs, γs./n)
    plt.xlim(1e-1, 1e1)
    plt.ylim(0, 1.6)
    plt.grid()
    # plt.ylabel("\$\\gamma^*\$")
    plt.ylabel("\$T^*\$")
    plt.xlabel("\$\\beta\$")
    plt.grid(true, which="major", ls="-")
    plt.grid(true, which="minor", ls=":")
    savetikz("./figures/gamma.tex")
    return
end

"""optimal MDS code rate ((12) of Lee)"""
mds_rate_opt(μ) = 1 + 1 / lambertw(-exp(-μ-1), -1)

"""

Return the optimal MDS code rate evaluated exhaustively.

"""
function mds_rate_opt_ex(μ)
    opt_ET = Inf
    opt_k = 0
    for k in range(eps(Float64), 1-eps(Float64), length=100)
        v = ET_mds(1, k, 1/μ)
        if v < opt_ET
            opt_ET = v
            opt_k = k
        end
    end
    return opt_k
end

function opt_α_lee(μ)
    opt_α = 0
    opt_err = Inf
    for α in range(eps(Float64), 1-eps(Float64), length=1000)
        err = abs(1 / (1-α) - log(1 / (1-α)) - μ - 1)
        if err < opt_err
            opt_err = err
            opt_α = α
        end
    end
    return opt_α
end

function opt_α(μ)
    opt_α = 0
    opt_err = Inf
    for α in range(eps(Float64), 1-eps(Float64), length=1000)
        err = abs(α / (1-α) - log(1 / (1-α)) - μ)
        if err < opt_err
            opt_err = err
            opt_α = α
        end
    end
    return opt_α
end

"""

Plot the optimal rate as a function of μ (Fig. 5b of Lee).

"""
function opt_rate_plot()
    μs = 10.0.^range(-4, 4, length=50)
    println(μs)
    βs = 1 ./ μs
    vs = mds_rate_opt.(μs)
    vs_ex = mds_rate_opt_ex.(μs)
    plt.semilogx(βs, vs)
    plt.semilogx(βs, vs_ex, ".")
    plt.xlim(1/1e4, 1/1e-4)
    plt.ylim(0, 1)
    plt.ylabel("\$k^*/n\$")
    plt.xlabel("\$\\beta\$")
    plt.grid(true, which="major", ls="-")
    plt.grid(true, which="minor", ls=":")
    # savetikz("./figures/optrate.tex")
    return
end

function main()

    # MDS
    β = 100
    println(mds_rate_opt(1/β))
    return
    # μ = 1 / β
    # n = 1:1000
    # μs = range()
    # plt.plot(n, exp_uncoded.(n, μ))
    # return
    # println(mean(d) / n)
    # return
    # println(1 + log(n) / μ)
    # return
    n = 10
    k = 5


    n = 10
    k = 2
    d_mds = ExponentialOrder(β, n, k)
    t = range(0, β, length=1000)
    v_mds = 1 .-cdf.(d_mds, t.-1/k)
    plt.semilogy(t, v_mds, label="MDS (10, 2)")

    d_uncoded = ExponentialOrder(β, n, n)
    v_uncoded = 1 .-cdf.(d_uncoded, n.*t)
    plt.semilogy(t, v_uncoded, label="Uncoded")

    # plt.xlim(0, 0.5)
    # plt.ylim(1e-3, 1e0)
    # println(quadgk(x->pdf(d, x), 0, 10))

    d_uncoded = ExponentialOrder(β, n, 1)
    t = range(0, β, length=1000)
    v = 1 .-cdf.(d_uncoded, t)
    # println(quadgk(x->pdf(d_uncoded, x), 0, 10))
    # println(quadgk(, 0, 100))
    plt.plot(t, v, "--", label="Rep.")

    # n = 1
    # d_uncoded = ExponentialK(β, n, 1)
    # t = range(0, 1, length=1000)
    # v = pdf.(d_uncoded, t)
    # # println(quadgk(x->pdf(d_uncoded, x), 0, 10))
    # # println(quadgk(, 0, 100))
    # plt.plot(t, v, "--", label="Rep. 1 worker")

    plt.legend()
    plt.show()
end
