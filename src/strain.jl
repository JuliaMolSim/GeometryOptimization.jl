export apply_voigt_strain

function voigt_to_full(v)
	[v[1] .5*v[6] .5*v[5];
	 .5*v[6] v[2] .5*v[4];
	 .5*v[5] .5*v[4] v[3]]
end

function full_to_voigt(A)
	[A[1, 1], A[2, 2], A[3, 3], 2. * A[2, 3], 2. * A[1, 3], 2. * A[1, 2]]
end

function apply_voigt_strain(system, strain)
	deformation = I + voigt_to_full(strain)
	new_lattice = eachrow(
			      reduce(hcat, bounding_box(system))' * deformation)
	# Also deform coordinates, since internally do not work in fractional coordinates.
	new_positions = [transpose(pos) * deformation for pos in position(system)]
	new_generalized_positions = ComponentVector(atoms=new_positions , bounding_box=new_lattice)
	update_positions(system, new_generalized_positions)
end
