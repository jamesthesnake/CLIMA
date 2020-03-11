module Writers

export AbstractWriter, NetCDFWriter, JLD2Writer, write_data

abstract type AbstractWriter end

"""
    write_data(
        writer,
        filename,
        dims,
        varvals,
    )

Writes the specified dimension names, dimensions, axes, variable names
and variable values to a file. Specialized by every `Writer` subtype.

# Arguments:
# - `writer`: instance of a subtype of `AbstractWriter`.
# - `filename`: into which to write data.
# - `dims`: Dict of dimension name to axis.
# - `varvals`: Dict of variable name to array of values.
"""
function write_data end

include("netcdf_writer.jl")
include("jld2_writer.jl")

end
