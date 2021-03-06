module CLIMA

include(joinpath("Utilities", "TicToc", "TicToc.jl"))
include(joinpath("Utilities", "ParametersType", "ParametersType.jl"))
include(joinpath("InputOutput", "Writers", "Writers.jl"))
include(joinpath("Common", "ConfigTypes", "ConfigTypes.jl"))
include(joinpath("Common", "PlanetParameters", "UniversalConstants.jl"))
include(joinpath("Common", "Parameters", "Parameters.jl"))
include(joinpath("Common", "PlanetParameters", "PlanetaryConstants.jl"))
include(joinpath("Common", "PlanetParameters", "PlanetParameters.jl"))
include(joinpath("Utilities", "RootSolvers", "RootSolvers.jl"))
include(joinpath("Utilities", "VariableTemplates", "VariableTemplates.jl"))
include(joinpath("Common", "MoistThermodynamics", "MoistThermodynamics.jl"))
include(joinpath(
    "Atmos",
    "Parameterizations",
    "CloudPhysics",
    "MicrophysicsParameters.jl",
))
include(joinpath(
    "Atmos",
    "Parameterizations",
    "CloudPhysics",
    "Microphysics.jl",
))
include(joinpath(
    "Atmos",
    "Parameterizations",
    "SurfaceFluxes",
    "SurfaceFluxes.jl",
))
include(joinpath(
    "Atmos",
    "Parameterizations",
    "TurbulenceConvection",
    "TurbulenceConvection.jl",
))
include(joinpath(
    "Atmos",
    "Parameterizations",
    "SubgridScaleTurbulence",
    "SubgridScaleParameters.jl",
))
include(joinpath("Arrays", "MPIStateArrays.jl"))
include(joinpath("Mesh", "Mesh.jl"))
include(joinpath("DGmethods", "Courant.jl"))
include(joinpath("DGmethods", "SpaceMethods.jl"))
include(joinpath("DGmethods", "DGmethods.jl"))
include(joinpath("Ocean", "ShallowWater", "ShallowWaterModel.jl"))
include(joinpath(
    "Ocean",
    "HydrostaticBoussinesq",
    "HydrostaticBoussinesqModel.jl",
))
include(joinpath("DGmethods_old", "DGBalanceLawDiscretizations.jl"))
include(joinpath("LinearSolvers", "LinearSolvers.jl"))
include(joinpath("LinearSolvers", "GeneralizedConjugateResidualSolver.jl"))
include(joinpath("LinearSolvers", "GeneralizedMinimalResidualSolver.jl"))
include(joinpath("LinearSolvers", "ColumnwiseLUSolver.jl"))
include(joinpath("ODESolvers", "ODESolvers.jl"))
include(joinpath("ODESolvers", "GenericCallbacks.jl"))
include(joinpath("Atmos", "Model", "AtmosModel.jl"))
include(joinpath("InputOutput", "VTK", "VTK.jl"))
include(joinpath("Diagnostics", "Diagnostics.jl"))
include(joinpath("Driver", "Driver.jl"))

end
