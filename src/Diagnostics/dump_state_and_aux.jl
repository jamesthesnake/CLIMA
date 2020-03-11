function dump_state_and_aux_init(dgngrp, currtime)
    @assert dgngrp.interpol !== nothing
end

function dump_state_and_aux_collect(dgngrp, currtime)
    DA = CLIMA.array_type()
    dg = Settings.dg
    Q = Settings.Q
    FT = eltype(Q.data)
    bl = dg.balancelaw

    istate = DA(Array{FT}(undef, dgngrp.interpol.Npl, num_state(bl, FT)))
    iaux = DA(Array{FT}(undef, dgngrp.interpol.Npl, num_aux(bl, FT)))

    # interpolate and save
    interpolate_local!(
        dgngrp.interpol,
        Q.data,
        istate,
        project = dgngrp.project,
    )
    interpolate_local!(
        dgngrp.interpol,
        dg.auxstate.data,
        iaux,
        project = dgngrp.project,
    )

    # filename (may also want to take out)
    nprefix = @sprintf(
        "%s_%s-%s-step%04d",
        dgngrp.out_prefix,
        dgngrp.name,
        Settings.starttime,
        dgngrp.step
    )
    statefilename = joinpath(Settings.output_dir, "$(nprefix).nc")
    auxfilename = joinpath(Settings.output_dir, "$(nprefix)_aux.nc")

    statenames = flattenednames(vars_state(bl, FT))
    #statenames = ("rho", "rhou", "rhov", "rhow", "rhoe")
    auxnames = flattenednames(vars_aux(bl, FT))

    # TODO: need to use `dgngrp.writer`
    write_interpolated_data(dgngrp.interpol, istate, statenames, statefilename)
    write_interpolated_data(dgngrp.interpol, iaux, auxnames, auxfilename)

    return nothing
end

function dump_state_and_aux_fini(dgngrp, currtime) end
