using CartesianGeometry
using GeometricalPredicates
using Base.Test

# Test Rescaled 2D Point
@test CartesianGeometry.rescale_coordinate(0., -1., 1.) == 1.5
@test CartesianGeometry.rescale_coordinate(0., unbounded_min_coord, unbounded_max_coord) == 1.5
@test 1.0 <= CartesianGeometry.rescale_coordinate(exp(100*rand()^2), unbounded_min_coord, unbounded_max_coord) < 2.0

p = Point2D_Unbounded((rand()-0.5)*unbounded_max_coord, (rand()-0.5)*unbounded_max_coord)
@test typeof(p) <: AbstractPoint2D
@test 1.0 <= getx(p) < 2.0
@test 1.0 <= gety(p) < 2.0

# Test Point comparison
@test Point2D(1., 1.5) < Point2D(1.1, 1.5)
@test Point2D(1.2, 1.5) < Point2D(1.2, 1.6)
@test Point2D_Unbounded(-10., 3.) < Point2D_Unbounded(-9.9, 0.1)

# Line functions
@test line_isequal(0., 0., 1.,1., -1., -1., 2., 2.)
@test Line2D(Point(1., 1.), Point(1.5,1.5)) == Line2D(Point(1.7, 1.7), Point(prevfloat(2.), prevfloat(2.)))



