import Linear as L
import numpy as np

def TestInversion(epsilon = 1e-15, N = 100, LOW_n = 1, HIGH_n = 10, LOW_r = 1, HIGH_r = 1000):
    rng = np.random.default_rng()
    failed1 = 0
    failed2 = 0
    total = 0
    for i in range(N):
        n = rng.integers(low = LOW_n, high = HIGH_n)
        A = rng.integers(low = LOW_r, high = HIGH_r, size = (n, n))
        A = A.astype(float)
        try:
            Ainv1 = L.MatrixInverse(A)
            Ainv2 = L.MatrixInverseGEPP(A)
            Ainv3 = np.linalg.inv(A)
        except (np.linalg.LinAlgError, RuntimeWarning, ZeroDivisionError):
            continue
        total += 1
        for x, y in np.nditer([Ainv1, Ainv3]):
            if abs(x-y) < epsilon:
                continue
            else:
                failed1 += 1
                break
        for z, w in np.nditer([Ainv2, Ainv3]):
            if abs(z-w) < epsilon:
                continue
            else:
                failed2 += 1
                break
    print("Proportion of tests passed for inversion by Standard Gaussian Elimination :", 100*((total - failed1)/total), "%.")
    print("Proportion of tests passed for inversion by Gaussian Elimination with partial pivoting :", 100*((total - failed2)/total), "%.")
    return


TestInversion(1e-15, 10000)
