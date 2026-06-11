import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

DIRECTORY = "output"
INTERVAL_MS = 1000
CONTOUR_LEVELS = 51

g_files = sorted([f"{DIRECTORY}/{f}" for f in os.listdir(DIRECTORY) if f.startswith("g") and f.endswith(".npy")])

n_frames = len(g_files)

figure = plt.figure()
axes = [
    figure.add_subplot(),
]

def update(step):
    figure.suptitle(f"{step} / {n_frames}")
    axes[0].clear()
    t = np.sum(np.load(g_files[step]), axis=-1)
    axes[0].contourf(t, levels=CONTOUR_LEVELS, vmin=-0.5, vmax=0.5, cmap="seismic")
    config = {
        "aspect": "equal",
    }
    axes[0].set(**config)

update(0)

animation = FuncAnimation(
    figure,
    update,
    frames=n_frames,
    interval=INTERVAL_MS,
    blit=False,
    repeat=False,
)

plt.show()
plt.close()

