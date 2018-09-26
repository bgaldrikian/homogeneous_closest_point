# Example Programs

These programs are small examples to demonstrate what is required to perform
different calculations using this implementation of the homogeneous closest
point algorithm.

To run, run make in this directory.  A bin folder will be created with the
executables ex1.exe through ex8.exe.  Each demonstrates a particular
calculation:

## ex1.exe

Closest-point example.  Finds the closest point on a tetrahedron to a query
point.

The output should be:
```
result = 1, closest point = (0 0.366025 1.36603 1)
```

## ex2.exe

Support point example.  Finds a point p in the given tetrahedron which
maximizes p dot q.

The output should be:
```
result = 1, support point = (5.96046e-08 1.73205 0 1)
```

## ex3.exe

Unbounded support point example.  Finds a point p in the given open pyramid
which maximizes p dot q.  In this case there is no bounded maximum, so the
point is at infinity.  This is represented by a direction vector with last
(projective) component equal to 0.

The output should be:
```
result = 1, support direction = (-1 1 0 0)
```

## ex4.exe

Intersection example.  The closest point to a query point in the intersection
of two shapes is found.

The output should be:
```
result = 1, closest point = (0.25 0.25 0 1)
```

## ex5.exe

Separation example.  Attempt to find the closest point to the intersection of
two shapes that have no intersection.  HCP returns 0 to indicate there is no
intersection.

The output should be:
```
result = 0
```

## ex6.exe

Time of Intersection example.  The time when two moving shapes intersect.
The point of intersection is also returned.

The output should be:
```
result = 1, ToI = 3.42265, point of intersection = (0.57735 0.57735 0.57735 1)
```

## ex7.exe

Closest-point example, specialized for D = 3.  Finds the closest point on a
tetrahedron to a query point.

The output should be:
```
result = 1, closest point = (0 0.366025 1.36603 1)
```

## ex8.exe

Time of Intersection example, specialized for D = 3.  The time when two moving
shapes intersect.  The point of intersection is also returned.

The output should be:
```
result = 1, ToI = 3.42265, point of intersection = (0.57735 0.57735 0.57735 1)
```
