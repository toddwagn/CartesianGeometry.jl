# Cartestian geometry functions

module CartesianGeometry

using GeometricalPredicates
using MatrixHelper

import 
	Base.isless,
	Base.==,
	GeometricalPredicates.area,
	GeometricalPredicates.getx,
	GeometricalPredicates.gety

export  # Export functions
	find_bounding_polygon,
	generate_triangles,
	getx,
	gety,
	isless,
	line_isequal,
	simple_convex_polygon_area,
	triangle_area,

	# types
	Point2D_Unbounded,

	# constants
	unbounded_max_coord,
	unbounded_min_coord



# unscaled 2D point for use with GeometricalPredicates
const unbounded_max_coord = Float64(2^24)
const unbounded_min_coord = -unbounded_max_coord

immutable Point2D_Unbounded <: AbstractPoint2D
	_x::Float64
	_y::Float64
	_min_x::Float64
	_max_x::Float64
	_min_y::Float64
	_max_y::Float64
	Point2D_Unbounded(x::Float64,y::Float64,min_x::Float64,max_x::Float64,min_y::Float64,max_y::Float64) = new(x,y,min_x,max_x,min_y,max_y)
end
Point2D_Unbounded(x::Float64, y::Float64) = Point2D_Unbounded(x,y,unbounded_min_coord,unbounded_max_coord,unbounded_min_coord,unbounded_max_coord)

function rescale_coordinate(c::Float64, cmin::Float64, cmax::Float64)
	if c < cmin # Error?
		c = cmin
	end
	if c > cmax # Error?
		c = cmax
	end
	if c == cmax
		return max_coord
	end
	return min_coord + (c-cmin)/(cmax - cmin)
end

GeometricalPredicates.getx(p::Point2D_Unbounded) = rescale_coordinate(p._x, p._min_x, p._max_x)
GeometricalPredicates.gety(p::Point2D_Unbounded) = rescale_coordinate(p._y, p._min_y, p._max_y)

function isless{T<:AbstractPoint2D}(p1::T, p2::T)
	(getx(p1) < getx(p2)) || (getx(p1) == getx(p2) && gety(p1) < gety(p2))
end

# Line functions
function line_isequal(ax, ay, bx, by, cx, cy, dx, dy)
	return ((cy-ay)*(bx-ax) == (by-ay)*(cx-ax)) && ((dy-ay)*(bx-ax) == (by-ay)*(dx-ax))
end

function =={T<:AbstractPoint2D}(lineA::Line2D{T}, lineB::Line2D{T})
	return (orientation(lineA, geta(lineB)) == 0 && orientation(lineA, getb(lineB)) == 0)
end

function line2D_ysolve{T<:AbstractPoint2D}(line::Line2D{T}, x)
	return (gety(getb(line))-gety(geta(line)))*(x-getx(geta(line)))/(getx(getb(line))-getx(geta(line)))+gety(geta(line))
end

function line2D_xsolve{T<:AbstractPoint2D}(line::Line2D{T}, y)
	return (y-gety(geta(line)))*(getx(getb(line))-getx(geta(line)))/(gety(getb(line))-gety(geta(line)))+getx(geta(line))
end

function rand_line2D()
	coords = sortrows(rand(2,2))
	return Line2D(Point(coords[1],coords[3]), Point(coords[2],coords[4]))
end

function intersection{T<:AbstractPoint2D}(line::Line2D{T}, point::T)
	return orientation(line, point) == 0
end
intersection{T<:AbstractPoint2D}(point::T, line::Line2D{T}) = intersection(line, point)

function intersection{T<:AbstractPoint2D}(lineA::Line2D{T}, lineB::Line2D{T})
	# lines are the same
	if lineA == lineB
		return lineA
	end

	lineA_dx = getx(getb(lineA)) - getx(geta(lineA))
	lineB_dx = getx(getb(lineB)) - getx(geta(lineB))

	# Vertical lines
	if (lineA_dx == 0. && lineB_dx == 0.) # Both lines vertical & not identical (test above)
		return AbstractPoint2D
	elseif lineA_dx == 0. # lineA is vertical
		return Point(getx(geta(lineA)), line2D_ysolve(lineB, getx(geta(lineA))))
	elseif lineB_dx == 0. # lineB is vertical
		return Point(getx(geta(lineB)), line2D_ysolve(lineA, getx(geta(lineB))))
	end

	lineA_slope = (gety(getb(lineA)) - gety(geta(lineA))) / lineA_dx
	lineA_int = gety(geta(lineA)) - lineA_slope*getx(geta(lineA))

	lineB_slope = (gety(getb(lineB)) - gety(geta(lineB))) / lineB_dx
	lineB_int = gety(geta(lineB)) - lineB_slope*getx(geta(lineB))

	# Lines are parallel & not identical (test above)
	if lineA_slope == lineB_slope
		return AbstractPoint2D
	end

	return Point((lineB_int-lineA_int)/(lineA_slope-lineB_slope),
			lineA_slope*(lineB_int-lineA_int)/(lineA_slope-lineB_slope)+lineA_int)
end

function isbetween(x, a, b)
	if a >= b
		return b <= x <= a
	else
		return a <= x <= b
	end
end

function segment_intersection{T<:AbstractPoint2D}(line::Line2D{T}, point::T)
	pt = intersection(line, point)
	
	if (pt == point && 
		isbetween(getx(pt), getx(geta(line)), getx(getb(line))) && 
		isbetween(gety(pt), gety(geta(line)), gety(getb(line))) )
		return pt
	else
		return AbstractPoint2D
	end
end
segment_intersection{T<:AbstractPoint2D}(point::T, line::Line2D{T}) = segment_intersection(line,point)

function segment_intersection{T<:AbstractPoint2D}(lineA::Line2D{T}, lineB::Line2D{T})
	pt = intersection(lineA, lineB)

	if pt == segment_intersection(lineA, pt) && pt == segment_intersection(lineB, pt)
		return pt
	else
		return AbstractPoint2D
	end
end

# Polygon functions

function intersection{T<:AbstractPoint2D}(polygon::Polygon2D{T}, line::Line2D{T})
	intersection_points = AbstractPoint2D[]

	for l in getlines(polygon)
		if orientation(line, geta(l)) != orientation(line, getb(l))
			line_int = intersection(line, l)
			if line_int != AbstractPoint2D
				push!(intersection_points, line_int)
			end
		end
	end

	if length(intersection_points) > 0
		return intersection_points
	else
		return AbstractPoint2D
	end
end
intersection{T<:AbstractPoint2D}(line::Line2D{T}, polygon::Polygon2D{T}) = intersection(polygon, line)

function generate_triangles{T}(pointarray::Vector{T})
	startpt = selectperm(pointarray,1)

	vertextups = [ map(x->(x<size(pointarray,1)+1?x:mod(x,size(pointarray,1))), (startpt, startpt+i, startpt + i + 1) ) for i in 1:(size(pointarray,1)-2) ]

	return [ pointarray[union(x),:] for x in vertextups ]
end

function generate_triangles{T<:AbstractMatrix}(pointarray::T)
	startpt = selectpermrows(pointarray,1)

	vertextups = [ map(x->(x<size(pointarray,1)+1?x:mod(x,size(pointarray,1))), (startpt, startpt+i, startpt + i + 1) ) for i in 1:(size(pointarray,1)-2) ]

	return [ pointarray[union(x),:] for x in vertextups ]
end

function triangle_area{T<:Number}(verticies::Array{T,2})
	v = broadcast(-, verticies, verticies[1,:])

	return abs(v[2,1]*v[3,2] - v[2,2]*v[3,1])/2.
end

function triangle_area{T<:AbstractPoint2D}(m::Matrix{T})
	return triangle_area(vec(m))
end

function triangle_area{T<:AbstractPoint2D}(v::Vector{T})
	if length(v) != 3
		error("Not a triangle")
	end

	orig_x = getx(v[1])
	orig_y = gety(v[1])

	return abs((getx(v[2])-orig_x)*(gety(v[3])-orig_y) - (gety(v[2])-orig_y)*(getx(v[3])-orig_x))*0.5
end

function simple_convex_polygon_area(ptarray)
	trianglelist = generate_triangles(ptarray)

	return mapreduce(triangle_area, +, trianglelist)
end

end # Module
