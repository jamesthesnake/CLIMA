
@kernel function update!(dQ, Q, rka, rkb, dt, slow_δ, slow_dQ, slow_scaling)
    i = @index(Global, Linear)
    @inbounds begin
        if slow_δ !== nothing
            dQ[i] += slow_δ * slow_dQ[i]
        end
        Q[i] += rkb * dt * dQ[i]
        dQ[i] *= rka
        if slow_scaling !== nothing
            slow_dQ[i] *= slow_scaling
        end
    end
end

@kernel function lsrk_mri_update!(dQ, Q, rka, rkb, τ, dt, γs, Rs)
    i = @index(Global, Linear)
    @inbounds begin
        NΓ = length(γs)
        Ns = length(γs[1])
        dqi = dQ[i]

        for s in 1:Ns
            ri = Rs[s][i]
            sc = γs[NΓ][s]
            for k in (NΓ - 1):-1:1
                sc = sc * τ + γs[k][s]
            end
            dqi += sc * ri
        end

        Q[i] += rkb * dt * dqi
        dQ[i] = rka * dqi
    end
end
