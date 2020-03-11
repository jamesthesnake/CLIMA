using FileIO
using JLD2

struct JLD2Writer <: AbstractWriter end

function write_data(
    jld::JLD2Writer,
    filename,
    dims,
    varvals,
)
    jldopen(filename, "a+") do ds
        dimnames = collect(keys(dims))
        for dn in dimnames
            ds["dim_$dn"] = dn
        end
        for (dn, dv) in dims
            ds[dn] = dv
        end
        for (vn, vv) in varvals
            defVar(ds, vn, vv, dimnames)
        end
    end
    return nothing
end
