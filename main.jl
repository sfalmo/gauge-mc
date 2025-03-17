include("simulation.jl")

using Dates, JLD2, LaTeXStrings
using PythonPlot

rc("text", usetex="True")
rc("xtick", direction="in")
rc("ytick", direction="in")

ϕ_HR(r) = r < 1 ? Inf : 0

function ϕ_LJ(r)
    if r > 2.5
        return 0
    end
    rinv = 1 / r
    rinv3 = rinv * rinv * rinv
    rinv6 = rinv3 * rinv3
    rinv12 = rinv6 * rinv6
    4 * (rinv12 - rinv6)
end

Vext_wall(x; xw=0.5, L=10) = x < xw || x > L - xw ? Inf : 0

function Vext_LJ(x; L=10)
    result = 0
    for d in [x, L - x]
        dinv = 1 / d
        dinv3 = dinv * dinv * dinv
        dinv6 = dinv3 * dinv3
        dinv12 = dinv6 * dinv6
        result += 4 * (dinv12 - dinv6)
    end
    result
end

function Vext_doublewell(x; a=2.5, ϵ_w=2)
    ϵ_w * ((x - 5)^2 - a^2)^2 / a^4
end

ϕ_map = Dict("HR" => ϕ_HR, "LJ" => ϕ_LJ)
Vext_map = Dict("wall" => Vext_wall, "LJ" => Vext_LJ, "double-well" => Vext_doublewell)

function run_sim(ϕ_key, Vext_key, N, ϵ, dϵ)
    L = 10.0
    T = 1.0
    ϕ = ϕ_map[ϕ_key]
    Vext = Vext_map[Vext_key]
    simulate(N, L, T, Vext, ϕ, ϵ, dϵ; equilibration_time=Dates.Second(10), production_time=Dates.Second(60))
end

ρ_x(ps, ρ) = vec(sum(ρ; dims=1)) * (ps[2] - ps[1])

ρ_p(xs, ρ) = vec(sum(ρ; dims=2)) * (xs[2] - xs[1])


function sim_no_shift(ϕ_key, Vext_key, N)
    ϵ(x) = 0
    dϵ(x) = 0

    xs, ps, ρ, stats = run_sim(ϕ_key, Vext_key, N, ϵ, dϵ)
    println(stats["accept"] / stats["move"])
    jldsave("data/$(ϕ_key)_$(Vext_key)_N$(N)_no_shift.jld2"; ϵ=ϵ.(xs), xs, ps, ρ, stats)
end


function sim_shift(ϕ_key, Vext_key, N)
    L = 10
    A = 0.5
    ϵ(x) = A * sin(4π * x / L)
    dϵ(x) = A * cos(4π * x / L) * 4π / L

    xs, ps, ρ, stats = run_sim(ϕ_key, Vext_key, N, ϵ, dϵ)
    println(stats["accept"] / stats["move"])
    jldsave("data/$(ϕ_key)_$(Vext_key)_N$(N)_shift.jld2"; ϵ=ϵ.(xs), xs, ps, ρ, stats)
end


function sim_shift_broken(ϕ_key, Vext_key, N)
    L = 10
    A = 1.5
    ϵ(x) = A * sin(4π * x / L)
    dϵ(x) = A * cos(4π * x / L) * 4π / L

    xs, ps, ρ, stats = run_sim(ϕ_key, Vext_key, N, ϵ, dϵ)
    println(stats["accept"] / stats["move"])
    jldsave("data/$(ϕ_key)_$(Vext_key)_N$(N)_shift_broken.jld2"; ϵ=ϵ.(xs), xs, ps, ρ, stats)
end


function do_plot(ax, axVext, ϕ_key, Vext_key, N, shift_key)
    ϵ, xs, ps, ρ, stats = load("data/$(ϕ_key)_$(Vext_key)_N$(N)_$(shift_key).jld2", "ϵ", "xs", "ps", "ρ", "stats")

    ax[0].plot(xs, ϵ, color="tab:orange")
    ax[1].plot(xs, ρ_x(ps, ρ), color="tab:blue")
    Vext = replace(Vext_map[Vext_key].(xs), Inf => 1000)
    axVext.plot(xs, Vext, color="tab:grey", alpha=0.5)
    ax[2].pcolormesh(xs, ps, ρ, cmap="BuPu", clim=(0, 0.5), rasterized=true)
end


function plot_no_shift_vs_shift(ϕ_key, Vext_key, N)
    fig, ax = subplots(3, 2, sharex=true, sharey="row", figsize=(3.5, 4.5), layout="constrained")
    ax = ax.T
    axVext0 = ax[0][1].twinx()
    axVext1 = ax[1][1].twinx()
    cmap = do_plot(ax[0], axVext0, ϕ_key, Vext_key, N, "no_shift")
    cmap = do_plot(ax[1], axVext1, ϕ_key, Vext_key, N, "shift_broken")
    axVext0.set_yticklabels([])
    axVext0.tick_params(axis="y", colors="tab:grey")
    axVext1.tick_params(axis="y", colors="tab:grey")
    axVext0.set_ylim(-1.5, 10)
    axVext1.set_ylim(-1.5, 10)
    axVext1.set_ylabel(L"$\beta V_\mathrm{ext}$", labelpad=-2, color="tab:grey")
    ax[0][0].set_title(LaTeXString("unshifted"))
    ax[1][0].set_title(LaTeXString("shifted"))
    ax[0][0].set_ylabel(L"$\epsilon / a$", labelpad=-3)
    ax[0][1].set_yticks([0, 0.5, 1, 1.5])
    ax[0][1].set_ylim(0, 1.6)
    ax[0][1].set_ylabel(L"$\rho a$")
    ax[0][2].set_xlim(0, 10)
    ax[0][2].set_ylim(-3, 3)
    ax[0][2].set_ylabel(L"$p / \sqrt{m k_B T}$")
    ax[0][2].set_xlabel(L"$x / a$")
    ax[1][2].set_xlabel(L"$x / a$")
    cbar = fig.colorbar(cmap, orientation="horizontal", ax=[ax[0][2], ax[1][2]], aspect=30)
    cbar.set_label(L"$f a \sqrt{m k_B T}$")
    fig.savefig("figures/$(ϕ_key)_$(Vext_key)_N$(N)_no-shift_vs_shift.pdf", pad_inches=0, bbox_inches="tight")
end


function plot_no_shift_vs_shift_vs_shift_broken(ϕ_key, Vext_key, N)
    fig, ax = subplots(3, 3, sharex=true, sharey="row", figsize=(6, 4.5), layout="constrained")
    ax = ax.T
    axVext0 = ax[0][1].twinx()
    axVext1 = ax[1][1].twinx()
    axVext2 = ax[2][1].twinx()
    cmap = do_plot(ax[0], axVext0, ϕ_key, Vext_key, N, "no_shift")
    cmap = do_plot(ax[1], axVext1, ϕ_key, Vext_key, N, "shift")
    cmap = do_plot(ax[2], axVext2, ϕ_key, Vext_key, N, "shift_broken")
    axVext0.set_yticklabels([])
    axVext0.tick_params(axis="y", colors="tab:grey")
    axVext1.set_yticklabels([])
    axVext1.tick_params(axis="y", colors="tab:grey")
    axVext2.tick_params(axis="y", colors="tab:grey")
    axVext0.set_ylim(0, 10)
    axVext1.set_ylim(0, 10)
    axVext2.set_ylim(0, 10)
    axVext2.set_ylabel(L"$\beta V_\mathrm{ext}$", labelpad=-2, color="tab:grey")
    # ax[0][0].set_title(LaTeXString("unshifted"))
    # ax[1][0].set_title(LaTeXString("shifted"))
    # ax[2][0].set_title(LaTeXString("shifted (broken)"))
    ax[0][0].set_ylabel(L"$\epsilon / a$", labelpad=-3)
    ax[0][1].set_yticks([0, 0.5, 1, 1.5])
    ax[0][1].set_ylim(0, 1.6)
    ax[0][1].set_ylabel(L"$\rho a$")
    ax[0][2].set_xlim(0, 10)
    ax[0][2].set_ylim(-3, 3)
    ax[0][2].set_ylabel(L"$p / \sqrt{m k_B T}$")
    ax[0][2].set_xlabel(L"$x / a$")
    ax[1][2].set_xlabel(L"$x / a$")
    ax[2][2].set_xlabel(L"$x / a$")
    cbar = fig.colorbar(cmap, orientation="horizontal", ax=[ax[0][2], ax[1][2], ax[2][2]], aspect=30, shrink=0.7)
    cbar.set_label(L"$f a \sqrt{m k_B T}$")
    fig.savefig("figures/$(ϕ_key)_$(Vext_key)_N$(N)_no-shift_vs_shift_vs_shift_broken.pdf", pad_inches=0, bbox_inches="tight")
end


ϕ_key = "HR"
Vext_key = "wall"
N = 6

sim_no_shift(ϕ_key, Vext_key, N)
sim_shift(ϕ_key, Vext_key, N)
sim_shift_broken(ϕ_key, Vext_key, N)

plot_no_shift_vs_shift_vs_shift_broken(ϕ_key, Vext_key, N)
