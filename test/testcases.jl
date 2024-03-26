@testsetup module TestCases
    using AtomsBase
    using ASEconvert
    using Unitful
    using UnitfulAtomic

    a = 5.431u"angstrom"
    lattice = a / 2 * [[0, 1, 1.], [1, 0, 1.], [1, 1 ,0.]]
    atoms     = [:Si => ones(3)/8, :Si => -ones(3)/8]
    silicon = periodic_system(atoms, lattice; fractional=true)
    rep = 3
    supercell_ase = convert_ase(silicon) * pytuple((rep, 1, 1))
    silicon_supercell = pyconvert(AbstractSystem, supercell_ase)

    a = 335.9u"pm"
    lattice = a  * [[1., 0, 0], [0 , 1, 0.], [0, 0 ,1.]]
    atoms     = [:Po => zeros(3)]
    polonium = periodic_system(atoms, lattice; fractional=true)
    supercell_ase = convert_ase(polonium) * pytuple((rep, 1, 1))
    polonium_supercell = pyconvert(AbstractSystem, supercell_ase)
end
