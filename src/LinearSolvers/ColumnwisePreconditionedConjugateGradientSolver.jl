module ColumnwisePreconditionedConjugateGradientSolver

export ColumnwisePreconditionedConjugateGradient

using ..LinearSolvers
const LS = LinearSolvers
using ..MPIStateArrays: device, realview

using LinearAlgebra
using LazyArrays
using StaticArrays
using GPUifyLoops

struct ColumnwisePreconditionedConjugateGradient{AT, FT, U, IT} <: LS.AbstractIterativeLinearSolver

  rtol::FT
  atol::FT

  r0::AT
  z0::AT
  p0::AT

  r1::AT
  z1::AT
  p1::AT
  Lp::AT

  dims::U

  max_iter::IT

end

# Define the outer constructor for the
# ColumnwisePreconditionedConjugateGradient
# struct
"""
function ColumnwisePreconditionedConjugateGradient(Q::AT; rtol = eps(eltype(Q)), atol = eps(eltype(Q)), dims = :) where AT

# Description
- Outer constructor for the ColumnwisePreconditionedConjugateGradient struct

# Arguments
- `Q`:(array). The kind of object that linearoperator! acts on.

# Keyword Arguments
- `rtol`: (float). relative tolerance
- `atol`: (float). absolute tolerance
- `dims`: (tuple or : ). the dimensions to compute norms over

# Return
- ColumnwisePreconditionedConjugateGradient struct
"""
function ColumnwisePreconditionedConjugateGradient(Q::AT; rtol = eps(eltype(Q)), atol = eps(eltype(Q)), max_iter = length(Q), dims = :) where AT
    container = []


    # allocate arrays
    r0 = similar(Q)
    z0 = similar(Q)
    p0 = similar(Q)
    r1 = similar(Q)
    z1 = similar(Q)
    p1 = similar(Q)
    Lp = similar(Q)

    # push to container
    push!(container, rtol)
    push!(container, atol)

    push!(container, r0)
    push!(container, z0)
    push!(container, p0)

    push!(container, r1)
    push!(container, z1)
    push!(container, p1)

    push!(container, Lp)

    push!(container, dims)

    push!(container, max_iter)
    # create struct instance
    return ColumnwisePreconditionedConjugateGradient{typeof(Q), eltype(Q), typeof(dims), typeof(max_iter)}(container...)
end


"""
LS.initialize!(linearoperator!, Q, Qrhs, solver::ColumnwisePreconditionedConjugateGradient, args...)

# Description

- This function initializes the iterative solver. It is called as part of the AbstractIterativeLinearSolver routine. SEE CODEREF for documentation on AbstractIterativeLinearSolver

# Arguments

- `linearoperator!`: (function). This applies the predefined linear operator on an array. Applies a linear operator to objecy "y" and overwrites object "z". linearoperator!(z,y)
- `Q`: (array). This is an object that linearoperator! outputs
- `Qrhs`: (array). This is an object that linearoperator! acts on
- `solver`: (struct). This is a scruct for dispatch, in this case for ColumnwisePreconditionedConjugateGradient
- `args...`: (arbitrary). This is optional arguments that can be passed into the function for flexibility.

# Keyword Arguments

- There are no keyword arguments

# Return
- `converged`: (bool). A boolean to say whether or not the iterative solver has converged.
- `threshold`: (float). The value of the residual for the first timestep

# Comment
- This function does nothing for conjugate gradient

"""
function LS.initialize!(linearoperator!, Q, Qrhs,
                        solver::ColumnwisePreconditionedConjugateGradient,
                        args...)

    return false, Inf
end


"""
LS.doiteration!(linearoperator!, Q, Qrhs, solver::ColumnwisePreconditionedConjugateGradient, threshold, applyPC!, args...)

# Description

- This function enacts the iterative solver. It is called as part of the AbstractIterativeLinearSolver routine. SEE CODEREF for documentation on AbstractIterativeLinearSolver

# Arguments

- `linearoperator!`: (function). This applies the predefined linear operator on an array. Applies a linear operator to objecy "y" and overwrites object "z". linearoperator!(z,y)
- `Q`: (array). This is an object that linearoperator outputs
- `Qrhs`: (array). This is an object that linearoperator acts on
- `solver`: (struct). This is a scruct for dispatch, in this case for ColumnwisePreconditionedConjugateGradient
- `threshold`: (float). Either an absolute or relative tolerance
- `applyPC!`: (function). Applies a preconditioner to objecy "y" and overwrites object "z". applyPC!(z,y)
- `args...`: (arbitrary). This is optional arguments that can be passed into the function for flexibility.

# Keyword Arguments

- There are no keyword arguments

# Return
- `converged`: (bool). A boolean to say whether or not the iterative solver has converged.
- `iteration`: (int). Iteration number for the iterative solver
- `threshold`: (float). The value of the residual for the first timestep

# Comment
- This function does conjugate gradient

"""
function LS.doiteration!(linearoperator!, Q, Qrhs,
                         solver::ColumnwisePreconditionedConjugateGradient,
                         threshold, applyPC!, args...)


    # unroll names for convenience

    rtol = solver.rtol
    atol = solver.atol
    residual_norm = typemax(eltype(Q))
    dims = solver.dims
    converged = false

    max_iter = solver.max_iter

    r0 = solver.r0
    z0 = solver.z0
    p0 = solver.p0

    r1 = solver.r1
    z1 = solver.z1
    p1 = solver.p1

    Lp = solver.Lp

    # Smack residual by linear operator
    linearoperator!(r0, Q)
    r0 .= Qrhs - r0
    applyPC!(z0, r0)
    p0 .= z0

    # TODO: FIX THIS
    absolute_residual = maximum( sqrt.( sum(r0 .* r0, dims=dims) ) )
    relative_residual = absolute_residual / maximum( sqrt.(sum(Qrhs .* Qrhs, dims=dims)) )
    # TODO: FIX THIS
    if (absolute_residual <= atol) || (relative_residual <= rtol)
        # wow! what a great guess
        converged = true
        return converged, 1, absolute_residual
    end

    for j in 1:max_iter

        linearoperator!(Lp, p0)

        α = sum(r0 .* z0, dims=dims) ./ sum(p0 .* Lp, dims = dims)

        # Update along preconditioned direction
        @. Q += α * p0

        @. r1 = r0 - α * Lp

        # TODO: FIX THIS
        absolute_residual = maximum( sqrt.( sum(r1 .* r1, dims=dims) ) )
        relative_residual = absolute_residual / maximum( sqrt.(sum(Qrhs .* Qrhs, dims=dims)) )
        # TODO: FIX THIS
        converged = false
        if (absolute_residual <= atol) || (relative_residual <= rtol)
            converged = true
            return converged, j, absolute_residual
        end

        applyPC!(z1, r1)

        β = sum(z1 .* r1, dims=dims) ./ sum(z0 .* r0, dims=dims)

        # Update
        @. p0 = z1 + β * p0
        @. z0 = z1
        @. r0 = r1

    end

    # TODO: FIX THIS
    converged = true
    return converged, max_iter, absolute_residual
end


end #module
