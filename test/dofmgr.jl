@testitem "DofManager" begin
    using AtomsBase
    using AtomsBuilder
    using EmpiricalPotentials
    using GeometryOptimization
    using LinearAlgebra
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
        x = @inferred GO.get_dofs(silicon, dofmgr)
        @test length(x) == 3 * length(silicon)
        @test typeof(x) == Vector{Float64}
        @test all(iszero, x)

        u = 0.01 * randn(length(x))
        newsys = GO.set_dofs(silicon, dofmgr, u)

        X = dofmgr.X0 + reinterpret(SVector{3, Float64}, u) * dofmgr.r0
        @test position(newsys, :) == X
    end

    @testset "Variable cell getter / setter (no clamped)" begin
        dofmgr = GO.DofManager(silicon; variablecell=true)
        x = @inferred GO.get_dofs(silicon, dofmgr)
        @test length(x) == 3 * length(silicon) + 9
        @test typeof(x) == Vector{Float64}
        @test all(iszero, x[1:end-9])
        @test reshape(x[end-8:end], 3, 3) == Matrix(I, 3, 3)

        u = 0.01 * randn(3length(silicon))
        F = I + 0.01randn(3, 3)
        newsys = GO.set_dofs!(silicon, dofmgr, [u; vec(F)])

        X = Ref(F) .* ( dofmgr.X0 + reinterpret(SVector{3, Float64}, u) * dofmgr.r0 )
        @test position(newsys, :)  == X
        @test bounding_box(newsys) == tuple([F * b for b in cell_vectors(newsys)]...)
    end

    @testset "energy_dofs / gradient_dofs agrees with raw energy" begin
        sw = StillingerWeber()
        dofmgrs = (GO.DofManager(silicon; variablecell=false),
                   GO.DofManager(silicon; variablecell=true ))

        for dofmgr in dofmgrs
            x = GO.get_dofs(silicon, dofmgr)

            ps    = AC.get_parameters(sw)
            state = AC.get_state(sw)
            Edof  = GO.energy_dofs(silicon,   sw, dofmgr, x, ps, state)
            Eref  = potential_energy(silicon, sw)
            @test Edof.energy_unitless * u"hartree" ≈ Eref

            gdof  = GO.gradient_dofs(silicon, sw, dofmgr, x, ps, state)
            @test length(gdof) == length(x)
            @test typeof(gdof) == Vector{Float64}

            ε = 1e-4
            d = randn(size(x))
            gd_ref = (  GO.energy_dofs(silicon, sw, dofmgr, x + ε * d, ps, state)
                      - GO.energy_dofs(silicon, sw, dofmgr, x - ε * d, ps, state)) / 2ε
            @test abs(gd_ref - dot(d, gdof)) < 1e-8
        end
    end

    # TODO Test clamped version of getters / setters
    # TODO Test clamped version of energy_dofs / gradient_dofs
end
