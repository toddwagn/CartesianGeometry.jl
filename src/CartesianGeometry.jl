# Cartestian geometry functions

module CartesianGeometry

using GeometricalPredicates
using MatrixHelper

export  # Export functions
	find_bounding_polygon,
	generate_triangles,
	simple_convex_polygon_area,
	triangle_area

function find_bounding_polygon(pointarray)
	startpt = selectpermrows(pointarray,1)
	
	# closure
	close = pointarray[startpt, :]

	curpt = starpt+1 < size(pointarray, 1)+1 ? startpt+1 : 1

	while !isequal(curpt, startpt) && mapreduce(x->isequal(x,zero(x)), &, close)
	end

	return false
end

function generate_triangles(pointarray)
	startpt = selectpermrows(pointarray,1)

	vertextups = [ map(x->(x<size(pointarray,1)+1?x:mod(x,size(pointarray,1))), (startpt, startpt+i, startpt + i + 1) ) for i in 1:(size(pointarray,1)-2) ]

	return [ pointarray[union(x),:] for x in vertextups ]
end

function triangle_area{T<:Number}(verticies::Array{T,2})
	v = broadcast(-, verticies, verticies[1,:])

	return abs(v[2,1]*v[3,2] - v[2,2]*v[3,1])/2.
end

function simple_convex_polygon_area(pointlist)
	ptarray = convert_pointlist(pointlist)
	trianglelist = generate_trianges(ptarray)

	return mapreduce(triangle_area, +, trianglelist)
end

end # Module
