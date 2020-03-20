using Test
using CLIMA
using CLIMA.ODESolvers
using CLIMA.LinearSolvers
using StaticArrays
using LinearAlgebra

CLIMA.init()
ArrayType = CLIMA.array_type()

mrigark_methods = [
    (MRIGARKERK33aSandu, 3)
    (MRIGARKERK45aSandu, 4)
]

fast_mrigark_methods = [
    (LSRK54CarpenterKennedy, 4)
    (LSRK144NiegemannDiehlBusch, 4)
]

@testset "3-rate ODE" begin
    ω1, ω2, ω3 =  100,  10,  1
    λ1, λ2, λ3 = -100, -10, -1
    β1, β2, β3 = 2, 2, 2

    ξ12 = λ2 / (λ1 + λ2)
    ξ13 = λ3 / (λ1 + λ3)
    ξ23 = λ3 / (λ2 + λ3)

    α12, α13, α23 = 1, 1, 1

    η12 = ((1 - ξ12) / α12) * (λ1 - λ2)
    η13 = 0#((1 - ξ13) / α13) * (λ1 - λ3)
    η23 = ((1 - ξ23) / α23) * (λ2 - λ3)

    η21 = ξ12 * α12 * (λ2 - λ1)
    η31 = 0#ξ13 * α13 * (λ3 - λ1)
    η32 = ξ23 * α23 * (λ3 - λ2)

    Ω = @SMatrix [
        λ1 η12 η13
        η21 λ2 η23
        η31 η32 λ3
    ]

    function rhs1!(dQ, Q, param, t; increment)
        @inbounds begin
            increment || (dQ .= 0)
            y1, y2, y3 = Q[1], Q[2], Q[3]
            g = @SVector [
                (-β1 + y1^2 - cos(ω1 * t)) / 2y1,
                (-β2 + y2^2 - cos(ω2 * t)) / 2y2,
                (-β3 + y3^2 - cos(ω3 * t)) / 2y3,
            ]
            dQ[1] += Ω[1, :]' * g - ω1 * sin(ω1 * t) / 2y1
        end
    end
    function rhs2!(dQ, Q, param, t; increment)
        @inbounds begin
            increment || (dQ .= 0)
            y1, y2, y3 = Q[1], Q[2], Q[3]
            g = @SVector [
                (-β1 + y1^2 - cos(ω1 * t)) / 2y1,
                (-β2 + y2^2 - cos(ω2 * t)) / 2y2,
                (-β3 + y3^2 - cos(ω3 * t)) / 2y3,
            ]
            dQ[2] += Ω[2, :]' * g - ω2 * sin(ω2 * t) / 2y2
        end
    end
    function rhs3!(dQ, Q, param, t; increment)
        @inbounds begin
            increment || (dQ .= 0)
            y1, y2, y3 = Q[1], Q[2], Q[3]
            g = @SVector [
                (-β1 + y1^2 - cos(ω1 * t)) / 2y1,
                (-β2 + y2^2 - cos(ω2 * t)) / 2y2,
                (-β3 + y3^2 - cos(ω3 * t)) / 2y3,
            ]
            dQ[3] += Ω[3, :]' * g - ω3 * sin(ω3 * t) / 2y3
        end
    end
    function rhs12!(dQ, Q, param, t; increment)
        rhs1!(dQ, Q, param, t; increment = increment)
        rhs2!(dQ, Q, param, t; increment = true)
    end
    function rhs_null!(dQ, Q, param, t; increment)
      increment || (dQ .= 0)
    end

    exactsolution(t) =
        [sqrt(β1 + cos(ω1 * t)), sqrt(β2 + cos(ω2 * t)), sqrt(β3 + cos(ω3 * t))]

    @testset "MRI-GARK method" begin
        finaltime = 1 / 2
        dts = [2.0^(-k) for k in 1:6]
        error = similar(dts)
        for (rate3_method, rate3_order) in mrigark_methods
            for (rate2_method, rate2_order) in mrigark_methods
                for (rate1_method, rate1_order) in fast_mrigark_methods
                  rate1_method = LSRKEulerMethod
                    for (n, dt) in enumerate(dts)
                        Q = exactsolution(0)
                        fastsolver = rate1_method(rhs_null!, Q; dt = dt)
                        midsolver =
                            rate2_method(rhs12!, fastsolver, Q, dt = dt / ω2)
                        slowsolver = rate3_method(rhs3!, midsolver, Q, dt = dt)
                        solve!(Q, slowsolver; timeend = finaltime)
                        error[n] = norm(Q - exactsolution(finaltime))
                    end

                    rate = log2.(error[1:(end - 1)] ./ error[2:end])
                    @show error
                    @show rate
                    @show rate3_order, rate2_order, rate1_order
                    # min_order = min(rate3_order, rate2_order, rate1_order)
                    # max_order = max(rate3_order, rate2_order, rate1_order)
                    # @test 2 <= rate[end]
                end
            end
        end
    end
end
