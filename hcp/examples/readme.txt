# Example Programs

These programs are small examples to demonstrate what is required to perform
different calculations using this implementation of the homogeneous closest
point algorithm.

To run, run make in this directory.  A bin folder will be created with the
executables ex1.exe through ex8.exe.  Each demonstrates a particular
calculation.  The output of a program on a different platform, or
compiled using a different compiler, may vary slightly because of floating
point representation and display precision.  A reasonable output is shown
with each example:

## ex1.exe

Closest point example.  Finds the closest point on the tetrahedron to a query
point.

The tetrahedron has sides bounded by the planes:
  yz plane
  zx plane
  xy plane
  the plane with normal and point (1/sqrt(3), 1/sqrt(3), 1/sqrt(3)) 

The query point is (1, 2, 3).

The correct query return value is 1, signifying a feasible set.  The exact
closest point is (0, (sqrt(3)-1)/2, (sqrt(3)+1)/2).

The output should be similar to:
```
result = 1, closest point = (0 0.366025 1.36603 1)
```

## ex2.exe

Support point example.  Finds a point p in the given tetrahedron which
maximizes p dot q, where q is the query direction.

The tetrahedron is the same as that used in ex1.exe.

The query direction q is (-1, 1,-1).

The correct query return value is 1, signifying a feasible set.  The exact
maximal point is (0, sqrt(3), 0).

The output should be similar to:
```
result = 1, support point = (5.96046e-08 1.73205 0 1)
```

## ex3.exe

Unbounded support point example.  Finds a vector p in the given open pyramid
which maximizes p dot q.  In this case there is no bounded maximum, so the
point is at infinity.  This is represented by a direction vector with last
(projective) component equal to 0.

The pyramid is bounded by the planes:
  zx plane
  xy plane
  the plane with normal and point (1/sqrt(3), 1/sqrt(3), 1/sqrt(3))

The query direction q is (-1, 1,-1).

The correct query return value is 1, signifying a feasible set.  The exact
maximal direction is any positive multiple of (-1, 1, 0).

The output should be similar to:
```
result = 1, support direction = (-1 1 0 0)
```

## ex4.exe

Intersection example.  The closest point to the origin in the intersection
of a tetrhedron and a cube is found.

The tetrahedron is the same as that used in ex1.exe.

The cube is axis-aligned with unit width, and centered at (3/4, 3/4, 0).

The correct query return value is 1, signifying a feasible set.  The exact
closest point to the origin in the intersection is (1/4, 1/4, 0).

The output should be similar to:
```
result = 1, closest point = (0.25 0.25 0 1)
```

## ex5.exe

Separation example.  Attempt to find the closest point to the origin in the
intersection of a tetrhedron and a cube, however the shapes have no
intersection.

The tetrahedron is the same as that used in ex1.exe.

The cube is axis-aligned with unit width, and centered at (7/4, 7/4, 0).

The correct query return value is 0, signifying an infeasible set.

The output should be:
```
result = 0
```

## ex6.exe

Time of Intersection example.  Calculates the time when a moving cube
intersects a static tetrahedron, and returns the point of intersection.

The tetrahedron is the same as that used in ex1.exe.

The cube is axis-aligned with unit width, and centered at (5, 5, 5).  It has
velocity (-1,-1,-1).

The correct ToI return value is 1, signifying that there is a time when the
sets intesect.  The exact time of intersection is 4 - 1/sqrt(3).  The exact
point of intersection at that time is (1/sqrt(3), 1/sqrt(3), 1/sqrt(3)).

The output should be similar to:
```
result = 1, ToI = 3.42265, point of intersection = (0.57735 0.57735 0.57735 1)
```

## ex7.exe

Closest point example with exactly the same setup as ex1.exe.  Here, HCP
specialized for D = 3 is used.  The result should be the same as that of
ex1.exe.

## ex8.exe

Time of Intersection example with exactly the same setup as ex6.exe.  Here,
HCP specialized for D = 4 is used.  (A 4-dimensional calculation is used  for
ToI in 3 dimensions.)  The result should be the same as that of ex6.exe.
