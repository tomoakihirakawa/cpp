import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from time import perf_counter_ns

# ---------------------------------------------------------------------------- #
# テスト
a = np.linalg.lstsq(
    [[2.11651903e+00, 8.31724217e-01, 2.19315279e-01, 4.27776241e-02, 3.05439267e-01],
     [9.08353866e-01, -4.87077837e-02, -
         1.71868677e-01, -5.69860472e-02, 2.35151844e-01],
     [0.00000000e+00, 6.53767046e-01, 1.12044981e-01, -
         1.84849052e-01, -1.00718344e-01],
     [0.00000000e+00, 0.00000000e+00, 9.54523812e-01,
         7.61932353e-02, -6.21088335e-01],
     [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.14245406e-01, 3.04759635e-01],
     [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.48336609e-16]],
    [7.41619849, 0., 0., 0., 0., 0.]
)

print("a", a)
# ---------------------------------------------------------------------------- #


class ArnoldiProcess:
    def __init__(self, linear_map, V0, n_iter):
        self.H = np.zeros((n_iter + 1, n_iter))
        self.V = np.zeros((n_iter + 1, V0.shape[0]))
        self.V[0] = V0
        for j in range(n_iter):
            # Compute next Krylov vector
            w = linear_map(self.V[j])
            # print(j,", w=",w,"Q[j]=",self.V[j],)
            # Gram-Schmidt orthogonalization
            for i in range(j + 1):
                self.H[i, j] = np.dot(w, self.V[i])
                w -= self.H[i, j] * self.V[i]
                # print((i,j),self.H[i, j] * self.V[i],",w=",w)
            self.H[j + 1, j] = np.linalg.norm(w)
            # Add new vector to basis
            self.V[j + 1] = w / self.H[j + 1, j]
            # print(j, ", ", self.V[j + 1], ", ", w, ", ", self.H[j + 1,j])


def gmres(linear_map, b, x0, n_iter):
    # Initialization
    r0 = b - linear_map(x0)
    beta = np.linalg.norm(r0)
    ap = ArnoldiProcess(linear_map, r0 / beta, n_iter)
    # Find best approximation in the basis V
    e1 = np.zeros(n_iter + 1)
    e1[0] = beta
    print("H", ap.H)
    print("e1", e1)
    print("V", ap.V)
    y = np.linalg.lstsq(ap.H, e1, rcond=None)[0]
    print("y", y)
    # Convert result back to full basis and return
    x_new = x0 + ap.V[:-1].T @ y
    return x_new


# -------------------------------------------------------- #
np.random.seed(179)

n = 5
N = n
shape = (n, n)

# Create random sparse (n, n) matrix with N non-zero entries
coords = np.random.choice(n * n, size=N, replace=False)
coords = np.unravel_index(coords, shape)
values = np.random.normal(size=N)
A_sparse = scipy.sparse.coo_matrix((values, coords), shape=shape)
A_sparse = A_sparse.toCRS()
A_sparse += scipy.sparse.eye(n)
A_dense = A_sparse.toarray()

# solution sould be
# {15.3565, -20.5184, -4.7769, 27.2011, -7.22571}

A_sparse = np.array([[0.0247911, 0.161413, 0.625419, 0.465341, 0.794249],
                     [0.895294,   0.215363, 0.280354, 0.0206005, 0.906597],
                     [0.457972, 0.76661,   0.590316, 0.535627, 0.00733951],
                     [0.315392, 0.925959, 0.412796,   0.825637, 0.322538],
                     [0.572894, 0.0998945, 0.738812, 0.30581,   0.904702]]
                    )

b = np.array([1, 2, 3, 4, 5])
# b = np.random.normal(size=N)
# n = 5
# -------------------------------------------------------- #

# # Solve using np.linalg.lstsq
# time_before = perf_counter_ns()
# x = np.linalg.lstsq(A_dense, b, rcond=None)[0]
# time_taken = (perf_counter_ns() - time_before) * 1e-6
# error = np.linalg.norm(A_dense @ x - b) ** 2
# print(f"Using dense solver: error: {error:.4e} in time {time_taken:.1f}ms")
#
# # Solve using inverse matrix
# time_before = perf_counter_ns()
# x = np.linalg.inv(A_dense) @ x
# time_taken = (perf_counter_ns() - time_before) * 1e-6
# error = np.linalg.norm(A_dense @ x - b) ** 2
# print(f"Using matrix inversion: error: {error:.4e} in time {time_taken:.1f}ms")
# Solve using GMRES
#
time_before = perf_counter_ns()
x = scipy.sparse.linalg.gmres(A_sparse, b, tol=1e-8)[0]
time_taken = (perf_counter_ns() - time_before) * 1e-6
error = np.linalg.norm(A_sparse @ x - b) ** 2
print(f"Using sparse solver: error: {error:.4e} in time {time_taken:.1f}ms")
#
# Try out the GMRES routine
#
time_before = perf_counter_ns()
x0 = np.zeros(n)
def linear_map(x): return A_sparse @ x


#
#
#
n_ite = 4
x = gmres(linear_map, b, x0, n_ite)
time_taken = (perf_counter_ns() - time_before) * 1e-6
error = np.linalg.norm(A_sparse @ x - b) ** 2
print("solution=", x)
print(f"Using GMRES: error: {error:.4e} in time {time_taken:.1f}ms")
