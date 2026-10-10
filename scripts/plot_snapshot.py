import sys
import yt

# usage: python plot_snapshot.py kh2d/kh_2d0015.dat
ds = yt.load(sys.argv[1])
print("time =", ds.current_time)

p = yt.SlicePlot(ds, "z", ("gas", "density"))
p.set_cmap(("gas", "density"), "inferno")
p.save("density.png")
