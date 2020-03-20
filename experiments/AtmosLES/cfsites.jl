using ArgParse
using Distributions
using Random
using StaticArrays
using Test
using Printf
using NCDatasets
using Dierckx
using LinearAlgebra
using DocStringExtensions

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

"""
CMIP6 Test Dataset - cfsites
@Article{gmd-10-359-2017,
AUTHOR = {Webb, M. J. and Andrews, T. and Bodas-Salcedo, A. and Bony, S. and Bretherton, C. S. and Chadwick, R. and Chepfer, H. and Douville, H. and Good, P. and Kay, J. E. and Klein, S. A. and Marchand, R. and Medeiros, B. and Siebesma, A. P. and Skinner, C. B. and Stevens, B. and Tselioudis, G. and Tsushima, Y. and Watanabe, M.},
TITLE = {The Cloud Feedback Model Intercomparison Project (CFMIP) contribution to CMIP6},
JOURNAL = {Geoscientific Model Development},
VOLUME = {10},
YEAR = {2017},
NUMBER = {1},
PAGES = {359--384},
URL = {https://www.geosci-model-dev.net/10/359/2017/},
DOI = {10.5194/gmd-10-359-2017}
}
"""

const seed = MersenneTwister(0)


struct GCMRelaxation{FT} <: Source
    "Relaxation timescale `[s]`"
    τ_relax::FT
end
function atmos_source!(
    s::GCMRelaxation,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    source.ρe -= (state.ρe - aux.ref_state.ρe) / s.τ_relax
    source.moisture.ρq_tot -=
        (state.moisture.ρq_tot - aux.ref_state.ρq_tot) / s.τ_relax
    source.ρu -= (state.ρu - aux.ref_state.ρu) / s.τ_relax
end
# ------------------------ End Boundary Condition --------------------- # 
# We first specify the NetCDF file from which we wish to read our 
# GCM values.
const data = Dataset("./datasets/HadGEM2-A_amip.2004-2008.07.nc", "r");


# Utility function to read and store variables directly from the 
# NetCDF file
function str2var(str::String, var::Any)
    str = Symbol(str)
    @eval(($str) = ($var))
end

# Additional utility function to allow group-id to be read in 
# from a command line argument. If no argument is specified 
# site1 is read by default. 
function parse_groupid()
    exc_handler =
        isinteractive() ? ArgParse.debug_handler : ArgParse.default_handler
    s = ArgParseSettings(exc_handler = exc_handler)
    @add_arg_table! s begin
        "--group-id"
        help = "Set up the GroupID for the CFSite"
        # If no argument is passed choose site1
        default = "site1"
    end
    return parse_args(s)
end

# Call the argument parser and store the command-line 
# value into the groupid variable 
parsed_args = parse_groupid()
groupid = parsed_args["group-id"]

# Define the get_gcm_info function
function get_gcm_info(groupid)

    @printf("--------------------------------------------------\n")
    @info @sprintf("""
       ____ _     ___ __  __    _                                  
      / ___| |   |_ _|  \\/  |  / \\                                 
     | |   | |    | || |\\/| | / _ \\                                
     | |___| |___ | || |  | |/ ___ \\                               
      \\____|_____|___|_| _|_/_/___\\_\\_  __       _     _____ ____  
     | | | | __ _  __| |/ ___| ____|  \\/  |     | |   | ____/ ___| 
     | |_| |/ _` |/ _` | |  _|  _| | |\\/| |_____| |   |  _| \\___ \\ 
     |  _  | (_| | (_| | |_| | |___| |  | |_____| |___| |___ ___) |
     |_| |_|\\__,_|\\__,_|\\____|_____|_|  |_|     |_____|_____|____/ 
     """)


    @printf("\n")
    @printf("CFSite experiment site ID = %s\n", groupid)
    @printf("--------------------------------------------------\n")
    filename = "/Users/asridhar/research/codes/CLIMA/datasets/HadGEM2-A_amip.2004-2008.07.nc"
    req_varnames = ("zg", "ta", "hus", "ua", "va", "pfull")
    # Load NETCDF dataset (HadGEM information)
    # Load the NCDataset (currently we assume all time-stamps are 
    # in the same NCData file). We store this information in `data`. 
    data = NCDataset(filename)
    # To assist the user / inform them of the data processing step
    # we print out some useful information, such as groupnames 
    # and a list of available variables
    @printf("Storing information for group %s ...", groupid)
    for (varname, var) in data.group[groupid]
        for reqvar in req_varnames
            if reqvar == varname
                var = mean(var, dims = 2)
                str2var(varname, var)
            end
        end
        # Store key variables
    end
    @printf("Complete\n")
    @printf("--------------------------------------------------\n")
    @printf("Group data storage complete\n")
    return (zg ./ grav, ta, hus, ua, va, pfull)
end

# Initialise the CFSite experiment :D! 
function init_cfsites!(bl, state, aux, (x, y, z), t, splines)
    FT = eltype(state)
    (spl_temp, spl_pfull, spl_ucomp, spl_vcomp, spl_sphum) = splines

    T = FT(spl_temp(z))
    q_tot = FT(spl_sphum(z))
    u = FT(spl_ucomp(z))
    v = FT(spl_vcomp(z))
    P = FT(spl_pfull(z))

    ρ = air_density(T, P, PhasePartition(q_tot))
    e_int = internal_energy(T, PhasePartition(q_tot))
    e_kin = (u^2 + v^2) / 2
    e_pot = grav * z
    # Assignment of state variables
    state.ρ = ρ
    state.ρu = ρ * SVector(u, v, 0)
    state.ρe = ρ * (e_kin + e_pot + e_int)
    if z <= FT(600)
        state.ρe += rand(seed) * FT(1 / 100) * (state.ρe)
    end
    state.moisture.ρq_tot = ρ * q_tot
    # Assignment of state variables
end

function config_cfsites(FT, N, resolution, xmax, ymax, zmax)

    # Boundary Conditions

    model = AtmosModel{FT}(
        AtmosLESConfigType;
        ref_state = GCMForcedState(),
        turbulence = SmagorinskyLilly{FT}(0.23),
        source = (Gravity(),),
        init_state = init_cfsites!,
        param_set = ParameterSet{FT}(),
    )

    imex_solver = CLIMA.DefaultSolverType()
    exp_solver =
        CLIMA.ExplicitSolverType(solver_method = LSRK144NiegemannDiehlBusch)

    config = CLIMA.AtmosLESConfiguration(
        "CFSites Experiments",
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        init_cfsites!,
        solver_type = exp_solver,
        model = model,
    )
    return config
end

function main()
    CLIMA.init()

    # Working precision
    FT = Float64
    # DG polynomial order
    N = 4
    # Domain resolution and size
    Δh = FT(50)
    Δv = FT(50)
    resolution = (Δh, Δh, Δv)
    # Domain extents
    xmax = FT(2000)
    ymax = FT(2000)
    zmax = FT(2000)
    # Simulation time
    t0 = FT(0)
    timeend = FT(3600 * 6)
    # Courant number
    CFL = FT(0.25)

    # Execute the get_gcm_info function
    (z, ta, hus, ua, va, pfull) = get_gcm_info(groupid)
    # Dropdims for compatibility with the interpolation module
    z = dropdims(z; dims = 2)
    ta = dropdims(ta; dims = 2)
    hus = dropdims(hus; dims = 2)
    ua = dropdims(ua; dims = 2)
    va = dropdims(va; dims = 2)
    pfull = dropdims(pfull; dims = 2)

    splines = (
        spl_temp = Spline1D(z, ta),
        spl_pfull = Spline1D(z, pfull),
        spl_ucomp = Spline1D(z, ua),
        spl_vcomp = Spline1D(z, va),
        spl_sphum = Spline1D(z, hus),
    )

    driver_config = config_cfsites(FT, N, resolution, xmax, ymax, zmax)
    solver_config = CLIMA.setup_solver(
        t0,
        timeend,
        driver_config,
        splines;
        init_on_cpu = true,
        Courant_number = CFL,
    )

    # User defined filter (TMAR positivity preserving filter)
    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do (init = false)
        Filters.apply!(solver_config.Q, 6, solver_config.dg.grid, TMARFilter())
        nothing
    end

    # Invoke solver (calls solve! function for time-integrator)
    result = CLIMA.invoke!(
        solver_config;
        user_callbacks = (cbtmarfilter,),
        check_euclidean_distance = true,
    )

end
main()
