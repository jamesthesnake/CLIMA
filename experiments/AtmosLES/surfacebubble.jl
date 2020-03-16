using Distributions
using Random
using StaticArrays
using Test
using DocStringExtensions
using LinearAlgebra

using CLIMA
using CLIMA.Atmos
using CLIMA.ConfigTypes
using CLIMA.GenericCallbacks
using CLIMA.DGmethods.NumericalFluxes
using CLIMA.ODESolvers
using CLIMA.Mesh.Filters
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters
using CLIMA.VariableTemplates

using CLIMA.Parameters
const clima_dir = dirname(pathof(CLIMA))
include(joinpath(clima_dir, "..", "Parameters", "Parameters.jl"))

# -------------------- Surface Driven Bubble ----------------- #
# Rising thermals driven by a prescribed surface heat flux.
# 1) Boundary Conditions:
#       Laterally periodic with no flow penetration through top
#       and bottom wall boundaries.
#       Momentum: Impenetrable(FreeSlip())
#       Energy:   Spatially varying non-zero heat flux up to time t₁
# 2) Domain: 1250m × 1250m × 1000m
# Configuration defaults are in `src/Driver/Configurations.jl`


"""
  Surface Driven Thermal Bubble
"""
function init_surfacebubble!(bl, state, aux, (x, y, z), t)
    FT = eltype(state)
    R_gas::FT = R_d
    c_p::FT = cp_d
    c_v::FT = cv_d
    γ::FT = c_p / c_v
    p0::FT = MSLP

    xc::FT = 1250
    yc::FT = 1250
    zc::FT = 1250
    θ_ref::FT = 300
    Δθ::FT = 0

    #Perturbed state:
    θ = θ_ref + Δθ # potential temperature
    π_exner = FT(1) - grav / (c_p * θ) * z # exner pressure
    ρ = p0 / (R_gas * θ) * (π_exner)^(c_v / R_gas) # density

    q_tot = FT(0)
    ts = LiquidIcePotTempSHumEquil(θ, ρ, q_tot, bl.param_set)
    q_pt = PhasePartition(ts)

    ρu = SVector(FT(0), FT(0), FT(0))
    # energy definitions
    e_kin = FT(0)
    e_pot = gravitational_potential(bl.orientation, aux)
    ρe_tot = ρ * total_energy(e_kin, e_pot, ts)
    state.ρ = ρ
    state.ρu = ρu
    state.ρe = ρe_tot
    state.moisture.ρq_tot = ρ * q_pt.tot
end

function config_surfacebubble(FT, N, resolution, xmax, ymax, zmax)
    # Boundary conditions
    # Heat Flux Peak Magnitude
    F₀ = FT(100)
    # Time [s] at which `heater` turns off
    t₁ = FT(500)
    # Plume wavelength scaling
    x₀ = xmax
    function energyflux(state, aux, t)
        x = aux.coord[1]
        y = aux.coord[2]
        MSEF = F₀ * (cospi(2 * x / x₀))^2 * (cospi(2 * y / x₀))^2
        t < t₁ ? MSEF : zero(MSEF)
    end

    C_smag = FT(0.23)

    ode_solver =
        CLIMA.ExplicitSolverType(solver_method = LSRK144NiegemannDiehlBusch)

    model = AtmosModel{FT}(
        AtmosLESConfigType;
        turbulence = SmagorinskyLilly{FT}(C_smag),
        source = (Gravity(),),
        boundarycondition = (
            AtmosBC(energy = PrescribedEnergyFlux(energyflux)),
            AtmosBC(),
        ),
        moisture = EquilMoist{FT}(),
        init_state = init_surfacebubble!,
        param_set = ParameterSet{FT}(),
    )
    config = CLIMA.AtmosLESConfiguration(
        "SurfaceDrivenBubble",
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        init_surfacebubble!,
        solver_type = ode_solver,
        model = model,
    )

    return config
end

function main()
    CLIMA.init()
    FT = Float64
    # DG polynomial order
    N = 4
    # Domain resolution and size
    Δh = FT(50)
    Δv = FT(50)
    resolution = (Δh, Δh, Δv)
    xmax = FT(2000)
    ymax = FT(2000)
    zmax = FT(2000)
    t0 = FT(0)
    timeend = FT(2000)

    driver_config = config_surfacebubble(FT, N, resolution, xmax, ymax, zmax)
    solver_config =
        CLIMA.setup_solver(t0, timeend, driver_config, init_on_cpu = true)

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do (init = false)
        Filters.apply!(solver_config.Q, 6, solver_config.dg.grid, TMARFilter())
        nothing
    end

    result = CLIMA.invoke!(
        solver_config;
        user_callbacks = (cbtmarfilter,),
        check_euclidean_distance = true,
    )

    @test isapprox(result, FT(1); atol = 1.5e-3)
end

main()
