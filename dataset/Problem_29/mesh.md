# Mesh

Use $440\times44$ uniform rectangular subdivisions of the
$1.0\ \mathrm{m}\times0.10\ \mathrm{m}$ channel. Each rectangle is split with
a crossed triangulation so the discrete mesh respects reflection about the
channel midline.

The spacings are equal:

$$h_x=h_y=\frac{1}{440}=\frac{0.10}{44}
=2.272727\times10^{-3}\ \mathrm{m}.$$

At $U_{\max}=0.01\ \mathrm{m/s}$ and $D=10^{-5}\ \mathrm{m^2/s}$, the
streamwise cell Peclet number is $Pe_h=1.13636$.

## FEniCS mesh code

```python
mesh = RectangleMesh(
    Point(0.0, 0.0), Point(L, H), nx, ny, diagonal="crossed"
)
```
