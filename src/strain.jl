
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
(box vectors in the columns).
"""
function bbox_to_matrix(bbox)
    reduce(hcat, bbox)
end

function matrix_to_bbox(matrix)
    collect.(eachcol(matrix))
end

@doc raw"""
Compute the strain (in Voigt notation) needed to deform the unit cell `bouding_box` into 
`deformed_bounding_box`.
"""
function compute_voigt_strain(bbox, deformed_bbox)
    strain_full = (bbox_to_matrix(deformed_bbox .- bbox)) / bbox_to_matrix(bbox)
    full_to_voigt(strain_full)
end
