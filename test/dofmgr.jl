@testitem "DofManager" begin




#=

using AtomsBase, DecoratedParticles, AtomsBuilder,
      GeomOpt, Test, StaticArrays, Unitful, LinearAlgebra

using AtomsCalculators: virial, forces, potential_energy

GO = GeomOpt
DP = DecoratedParticles

##

sys = AosSystem( rattle!(bulk(:Si, cubic=true) * 2, 0.1) )
dofmgr = GO.DofManager(sys)
x = GO.get_dofs(sys, dofmgr)
@test length(x) == 3 * length(sys)
@test typeof(x) == Vector{Float64}
@test all(iszero, x)
u = 0.01 * randn(length(x))
X = dofmgr.X0 + reinterpret(SVector{3, Float64}, u) * dofmgr.r0
GO.set_dofs!(sys, dofmgr, u)
@test position(sys, :) == X

##

sys = AosSystem( rattle!(bulk(:Si, cubic=true) * 2, 0.1) )
dofmgr = GO.DofManager(sys; variablecell = true)
x = GO.get_dofs(sys, dofmgr)
@test length(x) == 3 * length(sys) + 9
@test typeof(x) == Vector{Float64}
@test all(iszero, x[1:end-9])
@test x[end-8:end] == [1, 0, 0, 0, 1, 0, 0, 0, 1]
u = 0.01 * randn(length(x)-9)
F = SMatrix{3, 3}([1 0 0; 0 1 0; 0 0 1] + 0.01 * randn(3, 3))
x = [u; F[:]]
X = Ref(F) .* ( dofmgr.X0 + reinterpret(SVector{3, Float64}, u) * dofmgr.r0 )
bb_new = tuple([F * b for b in bounding_box(sys)]...)
GO.set_dofs!(sys, dofmgr, x)
@test position(sys, :) == X
@test bounding_box(sys) == bb_new

##

using EmpiricalPotentials, AtomsCalculators
using AtomsCalculators: potential_energy
using ACEbase

sys = AosSystem( rattle!(bulk(:Si, cubic=true) * (2,2,1), 0.1) )
dofmgr = GO.DofManager(sys; variablecell=false)

sw = StillingerWeber()
E1 = potential_energy(sys, sw)
x = GO.get_dofs(sys, dofmgr)
E2 = GO.energy_dofs(sys, sw, dofmgr, x)
@test E2 * u"eV" ≈ E1

##

g = GO.gradient_dofs(sys, sw, dofmgr, x)
@test length(g) == length(x)
@test length(g) == length(sys) * 3
@test typeof(g) == Vector{Float64}

_fd = ACEbase.Testing.fdtest( x -> GO.energy_dofs(sys, sw, dofmgr, x),
                        x -> GO.gradient_dofs(sys, sw, dofmgr, x),
                        x )
@test _fd

##

sys = AosSystem( rattle!(bulk(:Si, cubic=true) * (2,2,1), 0.1) )
dofmgr = GO.DofManager(sys; variablecell=true)

sw = StillingerWeber()
E1 = potential_energy(sys, sw)
x = GO.get_dofs(sys, dofmgr)
E2 = GO.energy_dofs(sys, sw, dofmgr, x)
@test E2 * u"eV" ≈ E1

g = GO.gradient_dofs(sys, sw, dofmgr, x)
@test length(g) == length(x)
@test length(g) == length(sys) * 3 + 9
@test typeof(g) == Vector{Float64}

_fd = ACEbase.Testing.fdtest( x -> GO.energy_dofs(sys, sw, dofmgr, x),
                        x -> GO.gradient_dofs(sys, sw, dofmgr, x),
                        x )
@test _fd

=#
end
