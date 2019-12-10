"""

    matrix_vector_operations(nrows, ncols)

return a tuple (nadditions, nmultiplications) required to multiply a
matrix by a vector

"""
matrix_vector_operations(nrows, ncols) = (nrows-1)*ncols, nrows*ncols

"""

    rs_fft_decoding_operations(n)

return a tuple (nadditions, nmultiplications) required to decode a
Reed-Solomon code of length n using the FFT-based method from "Novel
Polynomial Basis With Fast Fourier Transform and its Application to
Reed-Solomon Erasure Codes".

"""
function rs_fft_decoding_operations(n)
    f = (n, a, b, c) -> a+b*n*log2(c*n)
    return f(n, 2, 8.5, 0.86700826), f(n, 2, 1, 4)
end

"""

    rs_bm_decoding_operations(n)

return a tuple (nadditions, nmultiplications) required to decode a
Reed-Solomon code of length n using the Berlekamp-Massey algorithm
when k symbols are available to the decoder.

"""
rs_bm_decoding_operations(n, k) = n*(n-k-1), n^2*(1-k/n)

function decoding_delay_plot()
    ns =  10.0.^range(1, 4, length=100)
    # ks = ns ./ 2
    plt.loglog(ns, [arithmetic_delay.(rs_fft_decoding_operations(n)...) for n in ns], label="FFT")
    plt.loglog(ns, [arithmetic_delay.(rs_bm_decoding_operations(n, n/2)...) for n in ns], label="BM")
    plt.legend()
    plt.grid()
    plt.show()
end

"""

return the time needed to compute na additions and nm multiplications.

# ta=64, tm=64*log2(64)

"""
function arithmetic_delay(na, nm; ta=1e-9, tm=1e-9)
    return na*ta+nm*tm
end

"""

Return the minimum delay for the MDS-coded scheme.

"""
function min_delay_mds(n, k, β, s, nrows, ncols; enc=true, dec=true)

    # computation delay
    delay_mul = arithmetic_delay(matrix_vector_operations(nrows, ncols)...)
    delay_mul /= k
    delay = delay_mul

    # encoding delay
    rows_per_worker = nrows / k
    # delay_enc = arithmetic_delay(rs_bm_decoding_operations(n, k)...)
    delay_enc = arithmetic_delay(rs_fft_decoding_operations(n)...)
    delay_enc *= rows_per_worker * ncols
    # delay_enc = nrows*ncols*64 / 1e9 # RQ
    if enc
        delay += delay_enc
    end

    # decoding delay
    # delay_dec = arithmetic_delay(rs_bm_decoding_operations(n, k)...)
    delay_dec = arithmetic_delay(rs_fft_decoding_operations(n)...)
    delay_dec *= rows_per_worker
    # delay_dec = nrows*64 / 1e9 # RQ
    if dec
        delay += delay_dec
    end
    return delay
end

"""
    overall_delay_mds(n, k, β, nrows, ncols)

return the overall delay of the MDS-coded scheme.

"""
function mean_delay_mds(n, k, β, s, nrows, ncols; enc=true, dec=true)

    # straggler delay
    d = ExponentialOrder(β, n, k)
    delay = mean(d)
    delay += min_delay_mds(n, k, β, s, nrows, ncols, enc=enc, dec=dec)

    # println("str. ", mean(d))
    # println("mul. ", delay_mul)
    # println("enc. ", delay_enc)
    # println("dec. ", delay_dec)

    return delay
end

"""

return the overall delay of the MDS-coded scheme minimized over all k.

"""
function overall_delay_mds(n, β, s, nrows, ncols; enc=true, dec=true)
    k_opt = 1
    delay_mds = Inf
    for k in 1:n
        delay = overall_delay_mds(n, k, β, s, nrows, ncols, enc=enc, dec=dec)
        if delay < delay_mds
            delay_mds = delay
            k_opt = k
        end
    end
    # println("k ", k_opt)
    return delay_mds
end

"""

return the mean communication delay. bandwidth and latency numbers are
for Google data centers, as suggested in "The Datacenter as a
Computer: Designing Warehouse-Scale Machines, Third Edition".

"""
function communication_delay(nrows, ncols; symbolsize=64)
    nsymbols = nrows*ncols
    nbits = nsymbols*symbolsize
    bandwidth = 40e9
    delay = nbits/bandwidth
    return delay
end

"""

plot the shifted exponential.

In [1], they measure the delay of retrieving files in AWS S3 and
propose modeling it by a shifted exponential s+exp(β). As the file
size goes to 0, β approaches about 5-15 ms and the sum β+s about 25-50
ms.

[1] TOFEC: Achieving optimal throughput-delay trade-off of cloud
storage using erasure codes

"""
function shifted_exponential_plot(β=10e-3, s=1e-4)
    s = 0
    β = 10
    ts = range(0, 10/β, length=200)
    d = Exponential(1/β)
    # plt.plot(ts./1e-3, cdf.(d, ts.-s))
    plt.plot(ts, pdf.(d, ts.-s))
    plt.xlim(0, 0.6)
    plt.ylim(0, β)
    plt.xlabel("Time [s]")
    # plt.xlabel("\$t\$ [s]")
    # plt.ylabel("\$\\Pr(t \\leq T)\$")
    plt.ylabel("PDF")
    plt.grid()
    pltsavetikz(("./exponential.tex"))
end

function task_cdf_plot(n=100, k=25, β=10e-3, s=1e-4, nrows=1e4, ncols=1e3, enc=false, dec=true)
    n = 10
    k = 5
    # overall_delay_mds(n, β, s, nrows, ncols, enc=false)
    # return

    # computation delay
    delay_mul = arithmetic_delay(matrix_vector_operations(nrows, ncols)...)
    delay_uncoded = delay_mul / n

    # MDS-coded
    d = Exponential(β)
    ts = range(0, 10β, length=200)
    delay_mds = min_delay_mds(n, k, β, s, nrows, ncols; enc=false, dec=false)
    vs = pdf.(d, ts.-delay_mds)./1e3
    savetable(ts./1e-3, vs, "./figures/MatrixVectorPDF/data/task_mds.txt")
    plt.plot(ts./1e-3, vs, label="MDS w/o enc. and dec.")

    # with encoding
    delay_mds = min_delay_mds(n, k, β, s, nrows, ncols; enc=false, dec=true)
    vs = pdf.(d, ts.-delay_mds)./1e3
    savetable(ts./1e-3, vs, "./figures/MatrixVectorPDF/data/task_mds_dec.txt")
    plt.plot(ts./1e-3, vs, label="MDS w/o enc.")

    # with encoding and decoding
    ts = range(2500e-3, 3000e-3, length=200)
    delay_mds = min_delay_mds(n, k, β, s, nrows, ncols; enc=true, dec=true)
    vs = pdf.(d, ts.-delay_mds)./1e3
    savetable(ts./1e-3, vs, "./figures/MatrixVectorPDF/data/task_mds_enc_dec.txt")
    plt.plot(ts./1e-3, vs, label="MDS")

    ts = range(0, 10β, length=200)
    vs = pdf.(d, ts.-delay_uncoded)./1e3
    savetable(ts./1e-3, vs, "./figures/MatrixVectorPDF/data/task_uncoded.txt")
    plt.plot(ts./1e-3, vs, label="Uncoded")
    # plt.plot(ts./1e-3, cdf.(d, ts.-s))
    # delay = mean(d) + s

    plt.xlim(0, 20)
    plt.ylim(0, 0.12)
    plt.xlabel("\$t\$ [ms]")
    # plt.ylabel("\$\\Pr(t \\leq T_i)\$")
    plt.grid()
    plt.legend()

    # plt.yscale("log")
    # plt.xscale("log")

    # savetikz("./figures/MatrixVectorCDF/MatrixVectorCDF.tex")
    return
end

function overall_cdf_plot(n=100, k=25, β=10e-3, s=1e-4, nrows=1e4, ncols=1e3, enc=false, dec=true)
    n = 10
    k = 5

    # overall_delay_mds(n, β, s, nrows, ncols, enc=false)
    # return

    # computation delay
    delay_mul = arithmetic_delay(matrix_vector_operations(nrows, ncols)...)
    delay_uncoded = delay_mul / n

    # MDS-coded
    d = ExponentialOrder(β, n, k)
    ts = range(0, 60e-3, length=200)
    delay_mds = min_delay_mds(n, k, β, s, nrows, ncols; enc=false, dec=false)
    vs = pdf.(d, ts.-delay_mds)./1e3
    savetable(ts./1e-3, vs, "./figures/MatrixVectorPDF/data/overall_mds.txt")
    plt.plot(ts./1e-3, vs, label="MDS w/o enc. and dec.")

    # with encoding
    delay_mds = min_delay_mds(n, k, β, s, nrows, ncols; enc=false, dec=true)
    vs = pdf.(d, ts.-delay_mds)./1e3
    savetable(ts./1e-3, vs, "./figures/MatrixVectorPDF/data/overall_mds_dec.txt")
    plt.plot(ts./1e-3, vs, label="MDS w/o enc.")

    # with encoding and decoding
    ts = range(2500e-3, 3000e-3, length=100)
    delay_mds = min_delay_mds(n, k, β, s, nrows, ncols; enc=true, dec=true)
    vs = pdf.(d, ts.-delay_mds)./1e3
    savetable(ts./1e-3, vs, "./figures/MatrixVectorPDF/data/overall_mds_enc_dec.txt")
    plt.plot(ts./1e-3, vs, label="MDS")

    ts = range(0, 10β, length=100)
    d = ExponentialOrder(β, n, n)
    vs = pdf.(d, ts.-delay_uncoded)./1e3
    savetable(ts./1e-3, vs, "./figures/MatrixVectorPDF/data/overall_uncoded.txt")
    plt.plot(ts./1e-3, vs, label="Uncoded")
    # plt.plot(ts./1e-3, cdf.(d, ts.-s))
    # delay = mean(d) + s

    plt.xlim(0, 60)
    plt.ylim(0, 0.16)
    plt.xlabel("\$t\$ [ms]")
    plt.ylabel("\$\\Pr(t \\leq T_i)\$")
    plt.grid()
    plt.legend()

    # plt.yscale("log")
    # plt.xscale("log")

    # savetikz("./figures/MatrixVectorCDF/MatrixVectorCDF.tex")
    return
end

# function overall_cdf_plot(n=100, k=25, β=10e-3, s=1e-4, nrows=1e4, ncols=1e3, enc=false, dec=true)

#     # network latency
#     delay_mds = s
#     delay_uncoded = s

#     # computation delay
#     delay_mul = arithmetic_delay(matrix_vector_operations(nrows, ncols)...)
#     delay_mds += delay_mul / k
#     delay_uncoded += delay_mul / n

#     # encoding delay
#     rows_per_worker = nrows / k
#     delay_enc = arithmetic_delay(rs_fft_decoding_operations(n)...)
#     delay_enc *= rows_per_worker * ncols
#     if enc
#         delay_mds += delay_enc
#     end

#     # decoding delay
#     delay_dec = arithmetic_delay(rs_fft_decoding_operations(n)...)
#     delay_dec *= rows_per_worker
#     if dec
#         delay_mds += delay_dec
#     end

#     # straggler delay
#     ts = range(0, 20e-3, length=100)
#     # ts = range(0, 200e-3, length=100)
#     d = ExponentialOrder(β, n, k)
#     vs = cdf.(d, ts.-delay_mds)
#     plt.plot(ts./1e-3, vs, label="MDS")
#     # savetable(ts, vs, "./figures/MatrixVectorCDF/data/overall_mds.txt")

#     # ts = range(0, 60e-3, length=100)
#     ts = range(0, 200e-3, length=100)
#     d = ExponentialOrder(β, n, n)
#     vs = cdf.(d, ts.-delay_uncoded)
#     plt.plot(ts./1e-3, vs, label="Uncoded")
#     # savetable(ts, vs, "./figures/MatrixVectorCDF/data/overall_uncoded.txt")
#     # plt.plot(ts./1e-3, cdf.(d, ts.-s))
#     # delay = mean(d) + s

#     plt.xlim(0, 60)
#     plt.ylim(0, 1)
#     plt.xlabel("\$t\$ [ms]")
#     plt.ylabel("\$\\Pr(t \\leq T)\$")
#     plt.grid()
#     plt.legend()
# end

"""

β=1e-4 is suggested in the WSC book

"""
function overall_delay_plot(n=100, β=10e-3, s=1e-4, nrows=1e4, ncols=1e3)
    nrowss = round.(Int, 10.0.^range(log10(ncols), log10(ncols)+2, length=100))

    cdelay = communication_delay(nrows, ncols)
    println("Comm. delay: ", cdelay)

    # uncoded delay
    # β = communication_delay(nrows, ncols, n)
    d = ExponentialOrder(β, n, n)
    delays_uncoded = [arithmetic_delay(matrix_vector_operations(nrows, ncols)...)/n
                      for nrows in nrowss]
    delays_uncoded .+= mean(d) + s

    # MDS codes
    delays_mds = [overall_delay_mds(n, β, s, nrows, ncols) for nrows in nrowss]
    plt.plot(nrowss, delays_mds./delays_uncoded, "-", label="MDS enc., dec.")
    savetable(nrowss, delays_mds./delays_uncoded, "./figures/MatrixVectorDelay/data/mds_enc_dec.txt")

    delays_mds .= [overall_delay_mds(n, β, s, nrows, ncols, enc=false) for nrows in nrowss]
    plt.plot(nrowss, delays_mds./delays_uncoded, "-", label="MDS dec.")
    savetable(nrowss, delays_mds./delays_uncoded, "./figures/MatrixVectorDelay/data/mds_dec.txt")

    # delays_mds .= [overall_delay_mds(n, β, nrows, ncols, dec=false) for nrows in nrowss]
    # plt.plot(nrowss, delays_mds./delays_uncoded, "--", label="MDS, enc.")

    delays_mds .= [overall_delay_mds(n, β, s, nrows, ncols, enc=false, dec=false) for nrows in nrowss]
    plt.plot(nrowss, delays_mds./delays_uncoded, "-", label="MDS")
    savetable(nrowss, delays_mds./delays_uncoded, "./figures/MatrixVectorDelay/data/mds.txt")

    plt.yscale("log")
    plt.xscale("log")

    # full repetition delay
    d = ExponentialOrder(β, n, 1)
    delays_rep = [arithmetic_delay(matrix_vector_operations(nrows, ncols)...)
                  for nrows in nrowss]
    delays_rep .+= mean(d) + s

    # plt.semilogx(nrowss, (delays_mds.+cdelay)./delays_uncoded, label="MDS (incl. comm. delay)")
    # plt.semilogx(nrowss, delays_rep./delays_uncoded, label="Repetition")
    # plt.semilogx(nrowss, delays_uncoded./delays_uncoded, label="Uncoded")
    plt.grid()
    plt.legend()
    plt.ylim(1e-2, 1e2)
    plt.xlim(1e3, 1e5)
    plt.grid(true, which="major", ls="-")
    plt.grid(true, which="minor", ls=":")
    plt.ylabel("\$\\mathbb{E}[T^*_\\mathsf{MDS}] / \\mathbb{E}[T_\\mathsf{uncoded}]\$")
    plt.xlabel("\$N_\\mathsf{rows}\$")

    # savetikz("./figures/data/mds_delay.tex")
    plt.show()
    return
end

"""

Plot the fraction of time spent on communication.

"""
function comm_fraction_plot()
    comp_time = 1e-3
    n = round.(Int, 10.0.^range(0, 1, length=100))
    for f in [1, 10, 100]
        comm_time = 1e-9 ./n .+ 100e-6 .*(1 .- 1.0./n)
        comm_time *= f # memory accesses per ms
        plt.plot(n, comm_time ./ (comp_time .+ comm_time), "-o")
    end
    plt.ylim(0, 1)
    plt.xlim(n[1], n[end])
    plt.grid()
    return
end

"""

Plot the fraction of time spent on communication for matrix-vector
multiplication.

"""
function mv_comm_fraction_plot()
    bw = 10e9 # bits per second
    l = 100e-6 # network latency
    n = round.(Int, 10.0.^range(0, 3, length=100)) # number of nodes
    # for comp_delay in [1e-1, 1e-2, 1e-3, 1e-4]
    comp_delay = 1e-3
    m = sqrt.(n.*comp_delay ./ 2e-9) # matrix size
    # comm_delay = l.*n # due to delay
    rvs = ExponentialOrder.(l, n, n)
    comm_delay = mean.(rvs)
    # comm_delay .= m.*64 ./ bw # data transfer
    # comm_delay .+= l
    plt.semilogx(n, comm_delay ./ (comp_delay .+ comm_delay))
                 # label="\$T_\\mathsf{C}=$comp_delay\$")
    # end
    plt.ylim(0, 0.5)
    plt.xlim(n[1], n[end])
    plt.grid()
    # plt.xlabel("\$N_\\mathsf{servers}\$")
    plt.xlabel("Number of servers")
    plt.ylabel("Frac. of time spent on comm.")
    # plt.legend()
    return
end

function latency_plot()
    latency = [1e-9, 100e-9, 2*100e-6, 150e-3]
    labels = ["Arithmetic", "Main mem. ref.", "Intra-WSC RTT", "Inter-WSC RTT"]
    xs = collect(0:length(labels)-1)
    plt.bar(xs, latency)
    plt.yscale("log")
    plt.ylim(1e-10, 1e-0)
    plt.xticks(xs, labels)
    plt.ylabel("Latency [s]")
    plt.grid(true, axis="y")
    # savetikz("./figures/DelayBars/DelayBars.tikz")
end

function straggler_plot()
    # ns = round.(Int, 10.0.^range(0, 3, length=10e))
    ns = [1, 10, 100, 1000]
    β = 100e-6
    t = range(0, 10β, length=100)
    for n in ns
        rv = ExponentialOrder(β, n, n)
        plt.plot(t./1e-3, cdf.(rv, t), label="n=$n")
    end
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid()
    plt.xlabel("time [ms]")
    plt.ylabel("CDF")
    plt.legend()
    savetikz("./figures/Stragglers/Stragglers.tex")
    # rvs = ExponentialOrder.(β, n, n)
    # delays = mean.(rvs)
    # delays ./= delays[1]
    # plt.semilogx(n, delays)
end

"""

Plot the delay from the "Read operations in distributed file system
client" in Achieving Rapid Response Times in Large Online Services.

"""
function google_fs_straggler_plot()
    quantiles = [.5, .9, .99, .999]
    standard_idle = [19, 38, 67, 98]
    standard_busy = [24, 56, 108, 159]
    backup_idle = [16, 28, 38, 51]
    backup_busy = [19, 35, 67, 108]
    plt.plot(standard_idle, quantiles, label="Idle")
    plt.plot(standard_busy, quantiles, label="Busy")
    plt.plot(backup_idle, quantiles, label="Idle, \$2\$ ms rep.")
    plt.plot(backup_busy, quantiles, label="Busy, \$2\$ ms rep.")
    plt.grid()
    plt.ylim(0.5, 1)
    plt.xlim(0, 150)
end

function gd_straggler_plot()
    # % naive, 10: 0.23072916666666668
    # % ignore 1, 10: 0.1453125
    # % naive, 20: 0.2210526315789474
    # % ignore 3, 20: 0.1610047846889952
    # % naive, 30: 0.20703125
    # % ignore 5, 30: 0.18463541666666666

    n = [10, 20, 30]
    naive = [0.23072916666666668, 0.2210526315789474, 0.20703125]
    ignore = [0.1453125, 0.1610047846889952, 0.18463541666666666]
    plt.plot(n, naive)
    plt.plot(n, ignore)
    plt.xlim(10, 30)
    plt.ylim(0.14, 0.24)
end
