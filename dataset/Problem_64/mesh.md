# Mesh

The physical domain is the disk of radius $R=0.3$ centered at the origin.
The solver uses mshr.generate_mesh(Circle(Point(0,0), R), 40), producing a
conforming triangular approximation of the curved boundary. The full circle
boundary is clamped in the scalar membrane model.
