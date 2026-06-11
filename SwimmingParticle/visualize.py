dump_image = False

import os
import numpy as np
if dump_image:
    import matplotlib
    matplotlib.use("Agg")
from matplotlib import pyplot
from tqdm import tqdm

rc = np.load("output/xc.npy")
tc = np.load("output/yc.npy")

tc = np.concatenate([tc, tc[:1]])

tc, rc = np.meshgrid(tc, rc)

xc = rc * np.cos(tc)
yc = rc * np.sin(tc)

root = "output"
c_files = [f"{root}/{f}" for f in os.listdir(root) if f.startswith("concentration_") and f.endswith(".npy")]
c_files = sorted(c_files)

levels = np.linspace(0., 1., 11)

figure = pyplot.figure(figsize=(4, 4))
axes = [
    figure.add_axes([0, 0, 1, 1]),
]

for cnt, c_file in enumerate(tqdm(c_files)):
    c = np.load(c_file)
    axes[0].clear()
    axes[0].contourf(
        xc,
        yc,
        c[:, 1:],
        cmap="viridis",
        levels=levels,
        vmin=levels[0],
        vmax=levels[-1],
        extend="both"
    )
    axes[0].set_xlim([-10., +10.])
    axes[0].set_ylim([-10., +10.])
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].set_aspect("equal")
    if dump_image:
        figure.savefig(f"image/image{cnt:03d}.png")
    else:
        pyplot.show(block=False)
        pyplot.pause(1e-1)

pyplot.close()

