import numpy as np

DESTINATION = "grid/"


def save(file, dtype, data):
    np.save(f"{DESTINATION}/{file}", np.array(data, dtype=dtype))


class Coordinate:
    def __init__(self, length: float, nitems: int, grad: float):
        self.length = length
        self.nitems = nitems
        # Pre-allocate arrays
        self.faces = np.full(nitems + 2, np.nan)
        self.centers = np.empty(nitems + 2)
        self.face_scale_factors = np.full(nitems + 2, np.nan)
        self.center_scale_factors = np.full(nitems + 2, np.nan)
        # Calculate faces
        i_faces = np.arange(1, nitems + 2)
        x = (i_faces - 1) * length / nitems
        self.faces[1 : nitems + 1 + 1] = (
            (2.0 * grad - 2.0) * (x**3) + (-3.0 * grad + 3.0) * (x**2) + grad * x
        )
        # Calculate centers
        self.centers[0] = 0.0
        self.centers[1 : nitems + 1] = (
            0.5 * self.faces[1 : nitems + 1] + 0.5 * self.faces[2 : nitems + 2]
        )
        self.centers[nitems + 1] = length
        # Calculate scale factors
        self.face_scale_factors[1 : nitems + 2] = (
            self.centers[1 : nitems + 2] - self.centers[0 : nitems + 1]
        )
        self.center_scale_factors[1 : nitems + 1] = (
            self.faces[2 : nitems + 2] - self.faces[1 : nitems + 1]
        )


class PoissonOperator:
    def __init__(self, coordinate: Coordinate):
        nitems = coordinate.nitems
        center_scale_factors = coordinate.center_scale_factors[1:-1]
        face_scale_factors = coordinate.face_scale_factors[1:]
        lower_diagonal = 1.0 / center_scale_factors / face_scale_factors[:-1]
        upper_diagonal = 1.0 / center_scale_factors / face_scale_factors[1:]
        laplace_operator = np.zeros((nitems, nitems))
        for j in range(nitems):
            lower_component = lower_diagonal[j]
            upper_component = upper_diagonal[j]
            main_component = 0.0
            if 0 < j:
                main_component -= lower_component
                laplace_operator[j, j - 1] = lower_component
            if j < nitems - 1:
                main_component -= upper_component
                laplace_operator[j, j + 1] = upper_component
            laplace_operator[j, j] = main_component
        symmetrization_diagonal_components = np.sqrt(
            coordinate.center_scale_factors[1:-1]
        )
        symmetrization_operator = np.diag(symmetrization_diagonal_components)
        inverse_symmetrization_operator = np.diag(
            1.0 / symmetrization_diagonal_components
        )
        symmetric_laplace_operator = (
            symmetrization_operator @ laplace_operator @ inverse_symmetrization_operator
        )
        eigenvalues, eigenvectors = np.linalg.eigh(symmetric_laplace_operator)
        self.eigenvalues = eigenvalues
        self.forward = eigenvectors.T @ symmetrization_operator
        self.backward = inverse_symmetrization_operator @ eigenvectors


def main(lx, ly, nx, ny):
    save("lengths.npy", "<f8", [lx, ly])
    save("nitems.npy", "<u8", [nx, ny])
    x_coordinate = Coordinate(length=lx, nitems=nx, grad=0.50)
    y_coordinate = Coordinate(length=ly, nitems=ny, grad=0.25)
    save("x_centers.npy", "<f8", x_coordinate.centers)
    save("x_faces.npy", "<f8", x_coordinate.faces)
    save("x_center_scale_factors.npy", "<f8", x_coordinate.center_scale_factors)
    save("x_face_scale_factors.npy", "<f8", x_coordinate.face_scale_factors)
    save("y_centers.npy", "<f8", y_coordinate.centers)
    save("y_faces.npy", "<f8", y_coordinate.faces)
    save("y_center_scale_factors.npy", "<f8", y_coordinate.center_scale_factors)
    save("y_face_scale_factors.npy", "<f8", y_coordinate.face_scale_factors)
    x_poisson_operator = PoissonOperator(x_coordinate)
    save("x_eigenvalues.npy", "<f8", x_poisson_operator.eigenvalues)
    save("x_forward_operator.npy", "<f8", x_poisson_operator.forward)
    save("x_backward_operator.npy", "<f8", x_poisson_operator.backward)


if __name__ == "__main__":
    lx, ly = 1.0, 1.0
    nx, ny = 65, 65
    main(lx, ly, nx, ny)

