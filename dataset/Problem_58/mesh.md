# Mesh

The computational domain is the three-dimensional unit cube
$\Omega=(0,1)^3$. The solver creates a conforming tetrahedral mesh with
UnitCubeMesh(16, 16, 16). Boundary locations are identified geometrically:
$x=1$ is the velocity inflow, $x=0$ carries the pressure condition, and
$y=0,1$ are no-slip walls.
