using CartesianGeometry
using GeometricalPredicates
using Base.Test

# Test helper function
@test CartesianGeometry.isbetween(rand(), 0., 1.)
@test !CartesianGeometry.isbetween(2., 0., 1.)

# Test Rescaled 2D Point
@test CartesianGeometry.rescale_coordinate(0., -1., 1.) == 1.5
@test CartesianGeometry.rescale_coordinate(0., unbounded_min_coord, unbounded_max_coord) == 1.5
@test 1.0 <= CartesianGeometry.rescale_coordinate(exp(100*rand()^2), unbounded_min_coord, unbounded_max_coord) < 2.0

p = Point2D_Unbounded((rand()-0.5)*unbounded_max_coord, (rand()-0.5)*unbounded_max_coord)
@test isa(p, AbstractPoint2D)
@test 1.0 <= getx(p) < 2.0
@test 1.0 <= gety(p) < 2.0

# Test Point comparison
@test Point2D(1., 1.5) < Point2D(1.1, 1.5)
@test Point2D(1.2, 1.5) < Point2D(1.2, 1.6)
@test Point2D_Unbounded(-10., 3.) < Point2D_Unbounded(-9.9, 0.1)

# Line functions
@test line_isequal(0., 0., 1.,1., -1., -1., 2., 2.)
@test Line2D(Point(1., 1.), Point(1.5,1.5)) == Line2D(Point(1.7, 1.7), Point(prevfloat(2.), prevfloat(2.)))

@test line2D_ysolve(Line(Point(1., 1.), Point(2., 1.5)), 0.) == 0.5
@test line2D_xsolve(Line(Point(1., 1.), Point(1.5, 2.)), 0.) == 0.5

@test isa(rand_line2D(), Line2D{Point2D})

# Line - Point intersections
# This next test fails due to orientation bug in GeometricalPredicates
#@test intersection(Line(Point(1., 1.), Point(1.8, 1.5)), Point(1.4, 1.25)) == Point(1.4, 1.25)
# Picking coordinates outside of the restricted range functions correctly
@test intersection(Line(Point(0., 0.), Point(5., 5.)), Point(1., 1.)) == Point(1., 1.)
@test intersection(Line(Point(1., 1.), Point(2., 1.5)), Point(1.5, 1.)) == AbstractPoint2D
@test intersection(Line(Point(1., 1.), Point(2., 1.5)), Point(1.5, 1.5)) == AbstractPoint2D
@test intersection(Point(1., 1.), Line(Point(0.,0.), Point(5., 5.))) == Point(1., 1.)

# Line - Line intersections
@test intersection(Line(Point(1., 1.), Point(1.8, 1.8)), Line(Point(1., 1.8), Point(1.8, 1.))) == Point(1.4, 1.4)
@test intersection(Line(Point(1., 1.), Point(1.8, 1.8)), Line(Point(1.1, 1.1), Point(1.4, 1.4))) == Line(Point(1., 1.), Point(1.8, 1.8))
@test intersection(Line(Point(1., 1.), Point(1., 1.8)), Line(Point(1., 1.2), Point(1., 1.4))) == Line(Point(1., 1.), Point(1., 1.8))
@test intersection(Line(Point(1., 1.), Point(1., 1.8)), Line(Point(1.2, 1.), Point(1.2, 1.8))) == AbstractPoint2D
@test intersection(Line(Point(1.4, 1.), Point(1.4, 1.8)), Line(Point(1.0, 1.2), Point(1.8, 1.6))) == Point(1.4, 1.4)
@test intersection(Line(Point(1., 1.), Point(1.4, 1.2)), Line(Point(1., 1.2), Point(1.4, 1.4))) == AbstractPoint2D

# Line segment - Point intersections
@test segment_intersection(Line(Point(1.,1.), Point(1.4, 1.4)), Point(1.3, 1.3)) == Point(1.3, 1.3)
@test segment_intersection(Line(Point(1.,1.), Point(1.4, 1.4)), Point(1.4, 1.4)) == Point(1.4, 1.4)
@test segment_intersection(Line(Point(1.,1.), Point(1.4, 1.4)), Point(1.5, 1.5)) == AbstractPoint2D

# Line segment - line segment intersections
# Failure due to orientation issue in GeometricalPredicates
#@test segment_intersection(Line(Point(1., 1.), Point(1.4, 1.4)), Line(Point(1.2, 1.), Point(1.4, 1.6))) == Point(1.3, 1.3)

# Polygon functions
# intersection is found slightly misplaced from Point(1.2, 1.2)
@test intersection(Polygon(Point(1., 1.), Point(1.3, 1.), Point(1.3, 1.3)), Line(Point(1.1, 1.), Point(1.3, 1.4))) == AbstractPoint2D[Point(1.1, 1.), Point(nextfloat(1.2), nextfloat(1.2))]

#@test generate_triangles currently uncovered

# area functions
@test triangle_area([0. 0.; 3. 0.; 3. 4.]) == 6.

@test triangle_area([Point(0., 0.) Point(3., 0.) Point(3., 4.)]) == 6.
@test triangle_area([Point(0., 0.); Point(3., 0.); Point(3., 4.)]) == 6.

@test simple_convex_polygon_area([0. 0.; 0. 4.; 4. 4.; 0. 4.]) == 16.
@test simple_convex_polygon_area([-2. -2.; 2. -2.; 3. 0.; 2. 2.; -2. 2.; -3. 0.]) == 20.
