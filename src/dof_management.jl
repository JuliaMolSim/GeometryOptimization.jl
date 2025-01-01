"""
`DofManager`:

Constructor:
```julia
DofManager(sys; variablecell = false, r0 =..., free=..., clamp=..., mask=...)
```
* `variablecell` determines whether the cell is fixed or allowed to change
  during optimization
* `r0` is a reference length-scale, default set to one (in the unit of sys),
   this is used to non-dimensionalize the degrees of freedom. 

In addition set at most one of the kwargs:
* no kwarg: all atoms are free
* `free` : list of free atom indices (not dof indices)
* `clamp` : list of clamped atom indices (not dof indices)
* `mask` : 3 x N Bool array to specify individual coordinates to be clamped

### Meaning of dofs

On call to the constructor, `DofManager` stores positions and cell
`X0, C0`, dofs are understood *relative* to this initial configuration.
`get_dofs(sys, dm::DofManager)` returns a vector that represents the
non-dimensional displacement and a deformation matrix `(U, F)`. The new configuration extracted from a dof vector 
is understood as 
* The new cell: `C = F * C0` 
* The new positions: `𝐫[i] = F * (X0[i] + U[i] * r0)`
One aspect of this definition is that clamped atom positions still change via
the deformation `F`. This is natural in the context of optimizing the 
cell shape. 
"""
mutable struct DofManager{D,T}
   variablecell::Bool
   ifree::Vector{Int}   # extract the free position dofs
   r0::T
   X0::Vector{SVector{D,T}}       # reference positions
   C0::NTuple{D, SVector{D, T}}   # reference cell
end

# NOTES: 
#  - length units are implicitly given by the units in X0, C0, r0.
#    there should be no explicit length-unit stripping but this should be
#    implicitly through the reference length-scale r0
#  - at the moment energy-nondimensionalization is achieved simply by
#    stripping. A better approach would be to enforce this to happen in
#    the preconditioner, which could simply be a rescaling operation.

# ========================================================================
#  Constructors

function DofManager(sys::AbstractSystem{D};
                    variablecell=false,
                    r0=_auto_r0(sys),
                    free=nothing,
                    clamp=nothing,
                    mask=nothing) where {D}
    if D != 3
        error("this package assumes D = 3; please file an issue if you neeed a different use case.")
    end
    X0 = copy(position(sys, :))
    C0 = cell_vectors(sys)
    ifree = analyze_mask(sys, free, clamp, mask)
    DofManager(variablecell, ifree, r0, X0, C0)
end
function _auto_r0(sys)
    r = position(sys, 1)[1]
    one(ustrip(r)) * unit(r)
end


"""
`analyze_mask` : helper function to generate list of dof indices from
lists of atom indices indicating free and clamped atoms
"""
function analyze_mask(sys, free, clamp, mask)
    if sum(!isnothing, (free, clamp, mask)) > 1
        error("DofManager: At most one of `free`, `clamp`, `mask` may be provided")
    end
    if all(isnothing, (free, clamp, mask))
        # in this case (default) all atoms are free
        return collect(1:3length(sys))
    end

    # determine free dof indices
    n_atom = length(sys)
    if !isnothing(clamp)  # revert to setting free
        free = setdiff(1:n_atom, clamp)
    end
    if !isnothing(free)   # revert to setting mask
        mask = fill(false, 3, n_atom)
        if !isempty(free)
            mask[:, free] .= true
        end
    end
    return findall(mask[:])
end

# ========================================================================
#   DOF Conversions

length_unit(dm::DofManager)      = unit(dm.r0)
length_unit(sys::AbstractSystem) = unit(position(sys, 1)[1])
function check_length_units(sys, dm::DofManager)
    if length_unit(dm) != length_unit(sys) 
        error("System `sys` and DofManager have inconsistent units.")
    end
    if length(sys) != length(dm.X0)
        error("System `sys` and DofManager have inconsistent size.")
    end
end

variablecell(dofmgr::DofManager) = dofmgr.variablecell
fixedcell(dofmgr::DofManager)    = !variablecell(dofmgr)

# there is a type-instability here!! 
_posdofs(x, dofmgr::DofManager) = dofmgr.variablecell ? (@view x[1:end-9]) : x

function _pos2dofs(U::AbstractVector{SVector{3, T}}, dofmgr) where {T}
    @view(reinterpret(T, U)[dofmgr.ifree])
end

function _dofs2pos(x::AbstractVector{T}, dofmgr)  where {T} 
   u = zeros(T, 3 * length(dofmgr.X0))
   u[dofmgr.ifree] .= _posdofs(x, dofmgr)
   return reinterpret(SVector{3, T}, u)
end
_defm2dofs(F, dofmgr) = Matrix(F)[:]

function _dofs2defm(x::AbstractVector{T}, dofmgr) where {T}
    if dofmgr.variablecell
        SMatrix{3, 3, T}(x[end-8:end])
    else
        SMatrix{3, 3, T}([1 0 0; 0 1 0; 0 0 1])
    end
end

function get_dofs(sys::AbstractSystem, dofmgr::DofManager)
    check_length_units(sys, dofmgr)

    # obtain the positions and their underlying floating point type
    X = position(sys, :)
    if fixedcell(dofmgr)
        # there are allocations here that could maybe be avoided
        return collect(_pos2dofs((X - dofmgr.X0)/dofmgr.r0, dofmgr))
    else
        # variable cell case: note we already checked units and can strip
        #   (otherwise we would have problems inverting)
        bb = cell_vectors(sys)
        F = ustrip.(hcat(bb...)) / ustrip.(hcat(dofmgr.C0...))
        # Xi = F * (X0i + Ui * r0)  =>  Ui = (F \ Xi - X0i) / r0 
        U = [ (F \ X[i] - dofmgr.X0[i]) / dofmgr.r0 for i = 1:length(X) ]
        return [ _pos2dofs(U, dofmgr);
                 _defm2dofs(F, dofmgr) ]
    end
end


function set_dofs(system::AbstractSystem, dofmgr::DofManager,
                  x::AbstractVector{T} ) where {T <: AbstractFloat}
    check_length_units(system, dofmgr)

    # get the displacement from the dof vector
    U = _dofs2pos(x, dofmgr)
    F = _dofs2defm(x, dofmgr)

    # convert the displacements to positions
    positions = [ F * (dofmgr.X0[i] + U[i] * dofmgr.r0) for i = 1:length(U) ]
    bb_old = dofmgr.C0
    bb_new = ntuple(i -> F * bb_old[i], 3)

    # and update the system
    particles = [Atom(atom; position) for (atom, position) in zip(system, positions)]
    AbstractSystem(system; particles, cell_vectors=bb_new)
end



# ========================================================================
#   Compute the gradient with respect to dofs 
#   from forces and virials 

function energy_dofs(system, calculator, dofmgr, x::AbstractVector, ps, state)
    res = calculate(Energy(), set_dofs(system, dofmgr, x), calculator, ps, state)
    (; energy_unitless=ustrip(res.energy), res...)
end

function gradient_dofs(system, calculator, dofmgr, x::AbstractVector{T}, ps, state) where {T}
    # Compute and transform forces and virial into a gradient w.r.t. x
    if fixedcell(dofmgr)
        # fixed cell version
        # fi = - ∇_𝐫i E  [eV/A]
        # 𝐫i = X0[i] + r0 * U[i]
        # g_iα = - fiα * r0  [eV] => same unit as E so can strip

        res = calculate(Forces(), set_dofs(system, dofmgr, x), calculator, ps, state)
        g_pos = [ ustrip( - dofmgr.r0 * f ) for f in res.forces ]
        grad = collect(_pos2dofs(g_pos, dofmgr))::Vector{T}
    else
        # variable cell version
        # fi = - ∇_𝐫i E  [eV/A]     𝐫i = F * (X0[i] + r0 * U[i])
        # ∇_𝐮i' = - fi' * ∂𝐫i/∂𝐮i = - fi' * (r0 * F)   =>   ∇_𝐮i = - F' * r0 * fi
        # ∂F E |_{F = I} = - virial  => ∂F E = - virial / F'

        res = calculate((Forces(), Virial()), set_dofs(system, dofmgr, x),
                        calculator, ps, state)

        F = _dofs2defm(x, dofmgr)
        g_pos = [ - ustrip(dofmgr.r0 * F' * f) for f in res.forces ]

        grad = [ _pos2dofs(g_pos, dofmgr);
                 ( - ustrip.(res.virial) / F' )[:] ]::Vector{T}
    end

    (; grad, res...)
end
