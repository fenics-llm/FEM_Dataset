# Mesh

The canonical mesh is a uniform $64\times160$ crossed triangular subdivision
of $[0,2]\times[-0.5,0.5]$.

```python
mesh = RectangleMesh(Point(0.0, -0.5), Point(2.0, 0.5), 64, 160, "crossed")
```
