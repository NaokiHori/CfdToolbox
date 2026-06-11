import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

xlim = [0., 1.]
ylim = [0., 2.]

nx = 64
ny = 128

np.random.seed(0)

def compute_boundary_condition(ys: np.ndarray) -> [np.ndarray, np.ndarray]:
    uy_m = np.zeros(ys.shape)
    uy_p = np.zeros(ys.shape)
    for k in range(1, 8):
        amp = np.power(float(k), -1)
        phase = 2 * np.pi * np.random.random_sample()
        uy_m += amp * np.sin(2 * np.pi * (k * ys / (ylim[1] - ylim[0]) + phase))
    return uy_m, uy_p

# stream function in frequency domain
def get_f(k: float, x: float) -> float:
    return [
        1 * np.exp(+ k * x),
        x * np.exp(+ k * x),
        1 * np.exp(- k * x),
        x * np.exp(- k * x),
    ]

# differentiated stream function in frequency domain
def get_dfdx(k: float, x: float) -> float:
    return [
        (0 + k * 1) * np.exp(+ k * x),
        (1 + k * x) * np.exp(+ k * x),
        (0 - k * 1) * np.exp(- k * x),
        (1 - k * x) * np.exp(- k * x),
    ]

def compute_stream_function(xs: np.ndarray, ys: np.ndarray, bc: [np.ndarray, np.ndarray]) -> np.ndarray:
    # convert wall velocities to frequency domain
    uy_m = - np.fft.rfft(bc[0])
    uy_p = - np.fft.rfft(bc[1])
    # compute stream function in frequency domain
    psi_freq = np.zeros((ny // 2 + 1, nx), dtype="complex128")
    mat = np.zeros((4, 4), dtype="float64")
    vec = np.zeros((4), dtype="complex128")
    # zero-th mode cannot be determined: just assign zero
    psi_freq[0, :] = 0
    for k in range(1, ny // 2 + 1):
        # for each wavenumber, compute coefficients
        mat[0] = get_f(k, xlim[0])
        mat[1] = get_f(k, xlim[1])
        mat[2] = get_dfdx(k, xlim[0])
        mat[3] = get_dfdx(k, xlim[1])
        vec[0] = 0
        vec[1] = 0
        vec[2] = uy_m[k]
        vec[3] = uy_p[k]
        coefs = np.linalg.solve(mat, vec)
        # for each wall-normal position, compute stream function
        for i, x in enumerate(xs):
            psi_freq[k, i] = np.dot(coefs, get_f(k, x))
    # convert frequency domain to physical domain
    psi = np.zeros((ny, nx), dtype="float64")
    for i, x in enumerate(xs):
        psi[:, i] = np.fft.irfft(psi_freq[:, i])
    return psi

def show(xs: np.ndarray, ys: np.ndarray, bc: np.ndarray, psi: np.ndarray):
    def show_stream_function():
        figure = pyplot.figure()
        axis = figure.add_subplot()
        axis.contourf(xs, ys, psi, cmap="gray")
        keywords = {
            "title": "Stream function",
            "aspect": "equal",
            "xlim": xlim,
            "ylim": ylim,
        }
        axis.set(**keywords)
        axis.spines["left"].set_color("red")
        axis.spines["right"].set_color("blue")
        axis.spines["left"].set_linewidth(2)
        axis.spines["right"].set_linewidth(2)
        figure.savefig("cartesian_stream_function.jpg")
        pyplot.close()
    def show_wall_velocities():
        figure = pyplot.figure()
        axis = figure.add_subplot()
        axis.plot(ys / (ylim[1] - ylim[0]), bc[0], label="left", color="red")
        axis.plot(ys / (ylim[1] - ylim[0]), bc[1], label="right", color="blue")
        keywords = {
            "title": "Wall velocities",
            "xlabel": "$y / L_y$",
            "ylabel": "$u_y$",
        }
        axis.set(**keywords)
        axis.legend()
        figure.savefig("cartesian_wall_velocities.jpg")
        pyplot.close()
    show_stream_function()
    show_wall_velocities()

def main():
    xs = (xlim[1] - xlim[0]) * (np.arange(nx) + 0.5) / nx + xlim[0]
    ys = np.linspace(ylim[0], ylim[1], ny)
    bc = compute_boundary_condition(ys)
    psi = compute_stream_function(xs, ys, bc)
    show(xs, ys, bc, psi)

main()

