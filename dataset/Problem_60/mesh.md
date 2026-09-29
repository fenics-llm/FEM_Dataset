# Mesh

The domain is the square $[-1,1]^2$ with an open circular hole of radius
0.25 centered at $(0.5,0.5)$. The solver constructs a conforming triangular
CSG mesh with mshr.generate_mesh(domain, 40). Outer and hole boundaries are
identified geometrically. The regenerated mesh contains 4,618 triangular
cells and 2,414 vertices.
