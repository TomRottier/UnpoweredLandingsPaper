using Test, PointMass, OrdinaryDiffEq

@testset "all" verbose = true begin
    # test parameters
    p = Parameters(
        TrigAerodynamics(1.5, 0.03, 2.0),
        ConstantAlpha(0.2)
    )
    tspan = (0.0, 5.0)

    # energy conservation
    @testset "energy conservation" verbose = true begin
        p = Parameters(
            TrigAerodynamics(1.5, 0.0, 0.0),
            ConstantAlpha(0.2)
        )
        for _ in 1:10
            u₀ = rand(4)
            sol = simulate(full_dynamics!, u₀, p, tspan)
            TE₀ = total_energy(sol(0.0), p, 0.0)
            TE₁ = total_energy(sol(sol.t[end]), p, sol.t[end])

            @test TE₀ ≈ TE₁
        end
    end

    # power conservation
    @testset "power conservation" verbose = true begin
        for _ in 1:10
            u₀ = rand(4)
            sol = simulate(dynamics!, u₀, p, tspan)
            net_power = map(sol.t) do t
                u = sol(t)
                dke_dt = u[1] * (-CD(alpha(u, p.control, t), p.aerodynamics) * u[1]^2 - sin(u[2]))
                aero_power = aerodynamic_power(u, p, t)
                grav_power = gravitational_power(u, p, t)

                return aero_power + grav_power - dke_dt
            end
            @test net_power ≈ fill(0.0, length(sol)) atol = 1e-5
        end
    end

    # test negative velocity
    # @testset "negative velocity" verbose = true begin
    #     @test dynamics([1.0, 0.0], p, 0.0) ≈ dynamics([-1.0, π], p, 0.0)
    # end

    # test OrdinaryDiffEq extension
    @testset "DiffEq extension" verbose = true begin
        for _ in 1:10
            sol = simulate(full_dynamics!, rand(4))
            @test sol.retcode === ReturnCode.Success
            sol2 = simulate(dynamics!, rand(4), p, tspan)
            @test sol2.retcode === ReturnCode.Success
        end
    end
end
