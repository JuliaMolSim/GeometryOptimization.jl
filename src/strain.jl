export apply_voigt_strain, compute_voigt_strain

function voigt_to_full(v)
	[v[1] .5*v[6] .5*v[5];
	 .5*v[6] v[2] .5*v[4];
	 .5*v[5] .5*v[4] v[3]]
end

function full_to_voigt(A)
	[A[1, 1], A[2, 2], A[3, 3], 2. * A[2, 3], 2. * A[1, 3], 2. * A[1, 2]]
end

""" 
Transform a bounding box in list of vectors form in matrix form 
(box vectors in the rows).
"""
function bbox_to_matrix(bbox)
	reduce(hcat, bbox)
end
	
@doc raw"""
Deform system by applying given strain. Unit cell vectors as well as atomic 
positions are scaled. For a matrix of strains E, the matrix L of unit cell vectors 
is scaled accodring to L' = (I + E) L. 
Input strains should be given as a vector (Voigt notation).
"""
function apply_voigt_strain(system, strain)
	deformation_tensor = I + voigt_to_full(strain)
	new_lattice = eachcol(
			      deformation_tensor * bbox_to_matrix(bounding_box(system)))
	# Also deform coordinates, since internally do not work in fractional coordinates.
	new_positions = [deformation_tensor * pos for pos in position(system)]
	new_generalized_positions = ComponentVector(atoms=new_positions , bounding_box=new_lattice)
	update_positions(system, new_generalized_positions)
end

@doc raw"""
Compute the strain (in Voigt notation) needed to deform the unit cell `bouding_box` into 
`deformed_bounding_box`.
"""
function compute_voigt_strain(bbox, deformed_bbox)
	strain_full = (bbox_to_matrix(deformed_bbox .- bbox)) / bbox_to_matrix(bbox)
	full_to_voigt(strain_full)
end
