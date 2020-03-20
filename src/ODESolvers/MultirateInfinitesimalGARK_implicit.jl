export MRIGARKImplicit

"""
TODO: Document
"""
mutable struct MRIGARKImplicit{T, RT, AT, LT, Nstages, NΓ, FS, Nstages_sq} <:
               AbstractODESolver
    "time step"
    dt::RT
    "time"
    t::RT
    "rhs function"
    slowrhs!
    "implicit operator, pre-factorized"
    implicitoperator!
    "linear solver"
    linearsolver::LT
    "Storage for solution during the `MRIGARKImplicit` update"
    Qstages::NTuple{Nstages, AT}
    "Storage for RHS during the `MRIGARKImplicit` update"
    Rstages::NTuple{Nstages, AT}
    "Storage for the linear solver rhs vector"
    Qhat::AT
    "RK coefficient matrices for coupling coefficients"
    Γs::NTuple{NΓ, SArray{NTuple{2, Nstages}, RT, 2, Nstages_sq}}
    "RK coefficient matrices for embedded scheme"
    γ̂s::NTuple{NΓ, SArray{NTuple{1, Nstages}, RT, 1, Nstages}}
    "RK coefficient vector C (time scaling)"
    Δc::SArray{NTuple{1, Nstages}, RT, 1, Nstages}
    "fast solver"
    fastsolver::FS

    function MRIGARKImplicit(
        slowrhs!,
        linearsolver::AbstractLinearSolver,
        fastsolver,
        Γs,
        γ̂s,
        Q::AT,
        dt = nothing,
        t0 = 0,
    ) where {AT <: AbstractArray}

        @assert dt != nothing

        NΓ = length(Γs)
        Nstages = size(Γs[1], 1)
        T = eltype(Q)
        RT = real(T)
        LT = typeof(linearsolver)

        # Compute the Δc coefficients
        Δc = sum(Γs[1], dims = 2)[:]

        # Scale in the Δc to the Γ and γ̂, and convert to real type
        Γs = ntuple(k -> RT.(Γs[k] ./ Δc), NΓ)
        γ̂s = ntuple(k -> RT.(γ̂s[k] / Δc[Nstages]), NΓ)

        # Convert to real type
        Δc = RT.(Δc)

        # create storage for the stage values
        Qstages = (Q, ntuple(i -> similar(Q), Nstages - 1)...)
        Rstages = ntuple(i -> similar(Q), Nstages)
        Qhat = similar(Q)

        # create implicit operator. Term comes from the following
        # expression in (2.7) from Sandu (2019):
        #
        # v(0) = y_n;  v' = g(v) + f(y_n),  θ ∈ [0, H]; Y^slow_2 = v(H),
        # y_{n+1} = Y^slow_2 - 1/2 g(y_n) + 1/2 f(y_{n+1}) (Implicit term);
        #
        # where 1/2 is the term in the Γ0 matrix (implicit part, 2,2 entry)
        # TODO: Is this right?
        Γ0, = Γs
        α = dt * Γ0[2, 2]
        # Here we are passing NaN for the time since prefactorization assumes the
        # operator is time independent.  If that is not the case the NaN will
        # surface.
        implicitoperator! = prefactorize(
            EulerOperator(slowrhs!, α),
            linearsolver,
            Q,
            nothing,
            T(NaN),
        )

        FS = typeof(fastsolver)
        new{T, RT, AT, LT, Nstages, NΓ, FS, Nstages^2}(
            RT(dt),
            RT(t0),
            slowrhs!,
            implicitoperator!,
            linearsolver,
            Qstages,
            Rstages,
            Qhat,
            Γs,
            γ̂s,
            Δc,
            fastsolver,
        )
    end
end

"""
    MRIGARKIRK21aSandu(f!, fastsolver, Q; dt, t0 = 0)

The 2rd order, 2 stage implicit scheme from Sandu (2019).
"""
function MRIGARKIRK21aSandu(slowrhs!, linearsolver, fastsolver, Q; dt, t0 = 0)
    T = eltype(Q)
    RT = real(T)
    #! format: off
    Γ0 = [ 1 // 1   0 // 1
          -1 // 2   1 // 2
           0 // 1   0 // 1 ]
    γ̂0 = [-1 // 2   1 // 2]
    #! format: on
    MRIGARKImplicit(slowrhs!, linearsolver, fastsolver, (Γ0,), (γ̂0,), Q, dt, t0)
end

function dostep!(
    Q,
    mrigark::MRIGARKImplicit,
    param,
    timeend::Real,
    adjustfinalstep::Bool,
)
    time, dt = mrigark.t, mrigark.dt
    @assert dt > 0
    if adjustfinalstep && time + dt > timeend
        dt = timeend - time
        @assert dt > 0
    end

    dostep!(Q, mrigark, param, time, dt)

    if dt == mrigark.dt
        mrigark.t += dt
    else
        mrigark.t = timeend
    end
    return mrigark.t
end

function dostep!(
    Q,
    mrigark::MRIGARKImplicit,
    param,
    time::Real,
    dt::AbstractFloat,
)
    fast = mrigark.fastsolver

    implicitoperator!, linearsolver = mrigark.implicitoperator!, mrigark.linearsolver
    Qs = mrigark.Qstages
    Rs = mrigark.Rstages
    Δc = mrigark.Δc
    Qhat = mrigark.Qhat

    Nstages = length(Δc)
    groupsize = 256

    slowrhs! = mrigark.slowrhs!
    Γs = mrigark.Γs
    NΓ = length(Γs)

    rv_Q = realview(Q)
    rv_Qstages = realview.(Qs)
    rv_Rstages = realview.(Rs)
    rv_Qhat = realview(Qhat)

    ts = time
    for s in 1:Nstages
        # Stage dt
        dts = Δc[s] * dt
        stagetime = ts + dts

        p = param isa MRIParam ? param.p : param
        slowrhs!(Rs[s], Qs[s], param, ts, increment = false)

        if param isa MRIParam
            # fraction of the step slower stage increment we are on
            τ = (ts - param.ts) / param.Δts
            event = Event(device(Q))
            event = mri_update_rate!(device(Q), groupsize)(
                realview(Rs[s]),
                τ,
                param.γs,
                param.Rs;
                ndrange = length(realview(Rs[s])),
                dependencies = (event,),
            )
            wait(device(Q), event)
        end

        # fast solver
        γs = ntuple(k -> ntuple(j -> Γs[k][s, j], s), NΓ)
        mriparam = MRIParam(param, γs, realview.(Rs[1:s]), ts, dts)
        solve!(Q, mrigark.fastsolver, mriparam; timeend = stagetime)

        # this is a mess. need to figure out what the hell is actually going on
        # FIXME: DO WE EVEN NEED TO DO THIS?
        Γ0, = Γs
        α = dt * Γ0[s, s]
        linearoperator! = function (LQ, Q)
            slowrhs!(LQ, Q, p, stagetime; increment = false)
            @. LQ = Q - α * LQ
        end
        # slow stage continue
        linearsolve!(implicitoperator!, linearsolver, Qtt, Qhat, p, stagetime)
        slowrhs!(Rstages[s], Qstages[s], p, stagetime, increment = false)
        # update time
        ts += dts
    end
    # Need to compose solution together?!
end
