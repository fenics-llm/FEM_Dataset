# Mesh

The canonical mesh is a uniform $128\times64$ crossed triangular subdivision
of $[0,2]\times[-0.5,0.5]$. The solver identifies $x=0$ and $x=2$
periodically and applies no slip on $y=\pm0.5$.

```python
mesh = RectangleMesh(Point(0.0, -0.5), Point(2.0, 0.5), 128, 64, "crossed")
```
