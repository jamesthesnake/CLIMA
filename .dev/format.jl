#!/bin/bash
#=
export JULIA_PROJECT="$(dirname ${BASH_SOURCE[0]})"
julia --color=yes --startup-file=no -e 'using Pkg; Pkg.instantiate()'
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#

using JuliaFormatter

format(
    ".";
    verbose = true,
    indent = 4,
    margin = 80,
    always_for_in = true,
    whitespace_typedefs = true,
    whitespace_ops_in_indices = true,
    remove_extra_newlines = false,
)
