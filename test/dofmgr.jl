@testitem "DofManager" begin
    using AtomsBase
    using AtomsBuilder
    using AtomsCalculators
    using EmpiricalPotentials
    using GeometryOptimization
    using LinearAlgebra
    using StaticArrays
    using Unitful
    using UnitfulAtomic
    AC = AtomsCalculators
    GO = GeometryOptimization

    silicon_rattle = rattle!(bulk(:Si; cubic=true) * (2, 2, 1), 0.2u"Å")
    F = I + 0.01randn(3, 3)
    new_cell_vectors = tuple([F * v for v in cell_vectors(silicon_rattle)]...)
    silicon = AbstractSystem(silicon_rattle; cell_vectors=new_cell_vectors)

    @testset "Fixed cell getter / setter (no clamped)" begin
        dofmgr = GO.DofManager(silicon; variablecell=false)
        x = GO.get_dofs(silicon, dofmgr)
        # x = @inferred GO.get_dofs(silicon, dofmgr)
        @test length(x) == 3 * length(silicon)
        @test typeof(x) == Vector{Float64}
        @test all(iszero, x)

        u = 0.01 * randn(length(x))
        newsys = GO.set_dofs(silicon, dofmgr, u)

        X = dofmgr.X0 + reinterpret(SVector{3, Float64}, u) * dofmgr.r0
        @test maximum(x -> austrip(maximum(abs, x)), position(newsys, :) - X) < 1e-14
    end

    @testset "Variable cell getter / setter (no clamped)" begin
        dofmgr = GO.DofManager(silicon; variablecell=true)
        x = GO.get_dofs(silicon, dofmgr)
        # x = @inferred GO.get_dofs(silicon, dofmgr)
        @test length(x) == 3 * length(silicon) + 9
        @test typeof(x) == Vector{Float64}
        @test maximum(abs, x[1:end-9]) < 1e-14
        @test maximum(abs, reshape(x[end-8:end], 3, 3) - I) < 1e-14

        u = 0.01 * randn(3length(silicon))
        F = I + 0.01randn(3, 3)
        newsys = GO.set_dofs!(silicon, dofmgr, [u; vec(F)])

        X = Ref(F) .* ( dofmgr.X0 + reinterpret(SVector{3, Float64}, u) * dofmgr.r0 )
        X_diff  = position(newsys, :)   -  X
        cv_diff = cell_vectors(newsys) .- [F * b for b in cell_vectors(silicon)]

        @test maximum(x -> austrip(maximum(abs, x)),  X_diff) < 1e-14
        @test maximum(x -> austrip(maximum(abs, x)), cv_diff) < 1e-14
    end

    @testset "energy_dofs / gradient_dofs agrees with raw energy" begin
        sw = StillingerWeber()
        dofmgrs = (GO.DofManager(silicon; variablecell=false),
                   GO.DofManager(silicon; variablecell=true ))

        for dofmgr in dofmgrs
            x = GO.get_dofs(silicon, dofmgr)

            ps    = AC.get_parameters(sw)
            state = AC.get_state(sw)
            Edof  = GO.energy_dofs(silicon, sw, dofmgr, x, ps, state).energy_unitless
            @test Edof * u"hartree" ≈ AC.potential_energy(silicon, sw)

            gdof = GO.gradient_dofs(silicon, sw, dofmgr, x, ps, state).grad
            @test length(gdof) == length(x)
            @test typeof(gdof) == Vector{Float64}

            ε = 1e-5
            d = randn(size(x))
            get_ene(X) = GO.energy_dofs(silicon, sw, dofmgr, X, ps, state).energy_unitless
            gd_ref = (get_ene(x + ε * d) - get_ene(x - ε * d)) / 2ε
            @test abs(gd_ref - dot(d, gdof)) < 1e-8
        end
    end

    # TODO Test clamped version of getters / setters
    # TODO Test clamped version of energy_dofs / gradient_dofs
end
