using Dates

mutable struct Particle
    x::Float64
    p::Float64
end

Base.copy(particle::Particle) = Particle(particle.x, particle.p)

struct System
    N::Int
    L::Float64
    β::Float64
    Vext::Function
    ϕ::Function
    particles::Vector{Particle}
    ϵ::Function
    dϵ::Function
    stats::Dict{String,Int}
    function System(N::Int, L::Number, T::Number, Vext::Function, ϕ::Function, ϵ::Function=x -> 0, dϵ::Function=x -> 0; E_max=1000, max_init_tries=10000)
        particles = []
        i = 0
        while length(particles) != N
            push!(particles, Particle(rand() * L, randn() * sqrt(T)))
            E = calc_energy(particles, Vext, ϕ, L, ϵ, dϵ, length(particles))
            if E > E_max
                deleteat!(particles, length(particles))
            end
            i += 1
            if i >= max_init_tries
                error("Could not fit $(N) particles in the box after $(i) tries.")
            end
        end
        new(N, L, 1 / T, Vext, ϕ, particles, ϵ, dϵ, Dict("accept" => 0, "move" => 0))
    end
end

sample_ρ(system::System) = 1

mutable struct Histograms
    dx::Float64
    dp::Float64
    pmax::Float64
    ρ::Matrix{Float64}
    count::Int
    function Histograms(system::System; dx=0.01, dp=0.05, pmax=5.0)
        new(dx, dp, pmax, zeros(Int(system.L / dx), Int(2*pmax / dp)), 0)
    end
end

function bin(histograms::Histograms, system::System, particle::Particle)
    xbins, pbins = size(histograms.ρ)
    x = particle.x
    p = particle.p
    x, p = shift(x, p, system.ϵ, system.dϵ)
    xind = ceil(Int, x / system.L * xbins)
    pind = ceil(Int, (p + histograms.pmax) / (2 * histograms.pmax) * pbins)
    pind = clamp(pind, 1, pbins)
    CartesianIndex((xind, pind))
end

function sample(histograms::Histograms, system::System)
    for particle in system.particles
        histograms.ρ[bin(histograms, system, particle)] += 1
    end
    histograms.count += 1
end

function results(histograms::Histograms, system::System)
    dx = histograms.dx
    dp = histograms.dp
    pmax = histograms.pmax
    xs = dx/2:dx:system.L
    ps = -pmax+dp/2:dp:pmax
    ρ = histograms.ρ / (dx * dp * histograms.count)
    xs, ps, ρ', system.stats
end

function pbc!(system::System, i)
    system.particles[i].x -= floor(system.particles[i].x / system.L) * system.L
end

function dist(xi, xj, L)
    result = xj - xi
    result -= round(result / L) * L
    abs(result)
end

function shift(x, p, ϵ, dϵ)
    x + ϵ(x), p / (1 + dϵ(x))
end

function calc_energy(particles, Vext, ϕ, L, ϵ, dϵ, i)
    xi = particles[i].x
    pi = particles[i].p
    xi, pi = shift(xi, pi, ϵ, dϵ)
    E = Vext(xi)
    E += pi^2 / 2
    for (j, particle_j) in enumerate(particles)
        if i == j
            continue
        end
        xj = particle_j.x
        pj = particle_j.p
        xj, pj = shift(xj, pj, ϵ, dϵ)
        E += ϕ(dist(xi, xj, L))
        if isinf(E)
            break
        end
    end
    E
end

calc_energy(system::System, i) = calc_energy(system.particles, system.Vext, system.ϕ, system.L, system.ϵ, system.dϵ, i)

function trial_move(system::System; Δxmax=0.5, Δpmax=0.5)
    if isempty(system.particles)
        return
    end
    i = rand(1:length(system.particles))
    particle_before = copy(system.particles[i])
    Ebefore = calc_energy(system, i)
    system.particles[i].x += Δxmax * (2 * rand() - 1)
    pbc!(system, i)
    system.particles[i].p += Δpmax * (2 * rand() - 1)
    Eafter = calc_energy(system, i)
    if rand() > exp(-system.β * (Eafter - Ebefore))
        system.particles[i] = particle_before
        return false
    else
        return true
    end
end

function sweep(system::System; transitions=10, stats=true)
    for _ in 1:transitions
        accept = trial_move(system)
        if stats
          system.stats["move"] += 1
          if accept
              system.stats["accept"] += 1
          end
        end
    end
end

function simulate(N::Int, L::Number, T::Number, Vext::Function, ϕ::Function, ϵ::Function=x -> 0, dϵ::Function=x -> 0; equilibration_time=Dates.Second(1), production_time=Dates.Second(2), sweep_transitions=10)
    system = System(N, L, T, Vext, ϕ, ϵ, dϵ)
    histograms = Histograms(system)
    equilibration_start = now()
    while now() - equilibration_start < equilibration_time
        sweep(system; transitions=sweep_transitions, stats=false)
    end
    production_start = now()
    while now() - production_start < production_time
        sweep(system; transitions=sweep_transitions)
        sample(histograms, system)
    end
    results(histograms, system)
end
