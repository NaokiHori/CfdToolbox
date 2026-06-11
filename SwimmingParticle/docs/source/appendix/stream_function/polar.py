import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

xlim = [1., 4.]
ylim = [0., 2. * np.pi]
nx = 32
ny = 256

np.random.seed(0)

def compute_boundary_condition(ys: np.ndarray) -> [np.ndarray, np.ndarray]:
    uy_m = np.zeros(ys.shape)
    uy_p = np.zeros(ys.shape)
    for k in range(1, 8):
        amp = np.power(float(k), -1)
        phase = 2 * np.pi * np.random.random_sample()
        uy_m += amp * np.sin(k * ys + phase)
    return uy_m, uy_p

# stream function in frequency domain
def get_f(k: int, x: float) -> float:
    if 0 == k:
        raise RuntimeError("This code should not be reachable")
    elif 1 == k:
        return [
            np.power(x, -1),
            x,
            x * np.log(x),
            np.power(x, 3),
        ]
    else:
        return [
            np.power(x, 0 - k),
            np.power(x, 0 + k),
            np.power(x, 2 - k),
            np.power(x, 2 + k),
        ]

# differentiated stream function in frequency domain
def get_dfdx(k: int, x: float) -> float:
    if 0 == k:
        raise RuntimeError("This code should not be reachable")
    elif 1 == k:
        return [
            - np.power(x, -2),
            1,
            np.log(x) + 1,
            3 * np.power(x, 2),
        ]
    else:
        return [
            (0 - k) * np.power(x, 0 - k - 1),
            (0 + k) * np.power(x, 0 + k - 1),
            (2 - k) * np.power(x, 2 - k - 1),
            (2 + k) * np.power(x, 2 + k - 1),
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

def show(rs: np.ndarray, ts: np.ndarray, bc: np.ndarray, psi: np.ndarray):
    def show_stream_function(rs: np.ndarray, ts: np.ndarray):
        figure = pyplot.figure()
        axis = figure.add_subplot()
        rs, ts = np.meshgrid(rs, ts)
        xs = rs * np.cos(ts)
        ys = rs * np.sin(ts)
        axis.contourf(xs, ys, psi, cmap="gray")
        keywords = {
            "title": "Stream function",
            "aspect": "equal",
            "xlim": [- xlim[1], + xlim[1]],
            "ylim": [- xlim[1], + xlim[1]],
            "xticks": [],
            "yticks": [],
        }
        axis.set(**keywords)
        figure.savefig("polar_stream_function.jpg")
        pyplot.close()
    def show_wall_velocities(rs: np.ndarray, ts: np.ndarray):
        figure = pyplot.figure()
        axis = figure.add_subplot()
        rs, ts = np.meshgrid(rs, ts)
        xs = rs * np.cos(ts)
        ys = rs * np.sin(ts)
        axis.plot(np.linspace(0, 2 * np.pi, ny), bc[0], label="inner", color="red")
        axis.plot(np.linspace(0, 2 * np.pi, ny), bc[1], label="outer", color="blue")
        keywords = {
            "title": "Wall velocities",
            "xlabel": "$\\theta$",
            "ylabel": "$u_{\\theta}$",
        }
        axis.set(**keywords)
        axis.legend()
        figure.savefig("polar_wall_velocities.jpg")
        pyplot.close()
    show_stream_function(rs, ts)
    show_wall_velocities(rs, ts)

def main():
    xs = (xlim[1] - xlim[0]) * (np.arange(nx) + 0.5) / nx + xlim[0]
    ys = np.linspace(0, 2 * np.pi, ny)
    bc = compute_boundary_condition(ys)
    psi = compute_stream_function(xs, ys, bc)
    show(xs, ys, bc, psi)

main()

