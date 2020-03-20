using CLIMA.PlanetParameters
using DocStringExtensions
using CLIMA.PlanetParameters
using NCDatasets
using Dierckx

export GCMInput, NoCFSites, GCMForcedState

abstract type GCMInput end

vars_state(::GCMInput, FT) = @vars()
vars_aux(::GCMInput, FT) = @vars()
vars_integrals(::GCMInput, FT) = @vars()

function atmos_init_aux!(
    ::GCMInput,
    ::AtmosModel,
    aux::Vars,
    geom::LocalGeometry,
) end
function atmos_nodal_update_aux!(
    ::GCMInput,
    ::AtmosModel,
    state::Vars,
    aux::Vars,
    t::Real,
) end
function preodefun!(::GCMInput, aux::Vars, state::Vars, t::Real) end
function integrate_aux!(::GCMInput, integ::Vars, state::Vars, aux::Vars) end
function flux_radiation!(
    ::GCMInput,
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) end

# --------------------------- Do nothing --------------------------------------- # 
struct NoCFSites <: GCMInput end #

struct CFSource <: Source end

function atmos_source!(
    ::CFSource,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    # Populate with source terms here.
end

# 
# Get initial condition from NCData 
#
function get_ncdata()
    data = Dataset(
        "/home/asridhar/CLIMA/datasets/cfsites_forcing.2010071518.nc",
        "r",
    )
    # Load specific site group via numeric ID in NetCDF file (requires generalisation)
    siteid = data.group["site22"]
    # Allow strings to be read as varnames
    function str2var(str::String, var::Any)
        str = Symbol(str)
        @eval(($str) = ($var))
    end
    # Load all variables
    for (varname, var) in siteid
        str2var(varname, var[:, 1])
    end
    initdata = [height pfull temp ucomp vcomp sphum]
    return initdata
end

# --------------------------- Initialise --------------------------------------- #

"""
    GCMForcedState <: ReferenceState
A reference state prescribed by time-averaged data from GCM (various) forcing information
"""
struct GCMForcedState <: ReferenceState end
vars_aux(m::GCMForcedState, FT) =
    @vars(ρ::FT, p::FT, T::FT, ρu::SVector{3, FT}, ρe::FT, ρq_tot::FT)
# Idea is to pass additional args and store them in here for access through the source function 
function atmos_init_aux!(
    m::GCMForcedState,
    atmos::AtmosModel,
    aux::Vars,
    geom::LocalGeometry,
    splines,
)
    FT = eltype(aux)

    (x, y, z) = aux.coord

    (spl_temp, spl_pfull, spl_ucomp, spl_vcomp, spl_sphum) = splines

    T = FT(spl_temp(z))
    q_tot = FT(spl_sphum(z))
    u = FT(spl_ucomp(z))
    v = FT(spl_vcomp(z))
    p = FT(spl_pfull(z))

    aux.ref_state.T = T
    aux.ref_state.p = p

    q_pt = PhasePartition(q_tot)

    ρ = air_density(T, p, q_pt)

    e_kin = (u^2 + v^2) / 2
    e_pot = gravitational_potential(atmos.orientation, aux)
    e_int = internal_energy(T, q_pt)

    # Assignment of state variables
    aux.ref_state.ρ = ρ
    aux.ref_state.ρu = ρ * SVector(u, v, 0)
    aux.ref_state.ρe = ρ * (e_kin + e_pot + e_int)
    aux.ref_state.ρq_tot = ρ * q_tot

end
