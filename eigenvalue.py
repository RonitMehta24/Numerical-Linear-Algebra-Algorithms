import numpy as np

def EuclideanNorm(x):
    sum = 0
    for elem in x:
        sum += elem*elem
    return np.sqrt(sum)


def DotProduct(x, y):
    sum = 0
    for z, w in zip(x, y):
        sum += z*w
    return sum

def EigenVector(A, x0, epsilon):
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    lambda_true = eigenvalues[-1]
    v_true = eigenvectors[:, -1]
    error = epsilon + 1
    x = x0/EuclideanNorm(x0)
    rprev = 0
    E = []
    D = []
    
    while error > epsilon:
        y = A@x
        r = DotProduct(x, y)
        #print(y)
        xnew = y/EuclideanNorm(y)
        #print(xnew)
        #print(x)
        error = EuclideanNorm(xnew - x)
        actual_error = EuclideanNorm(xnew - v_true) 
        
        E.append(actual_error)
        D.append(abs(r-rprev))
        #print(error)
        r = rprev
        x = xnew
    print(np.log(D[-1]/D[-2])/np.log(E[-2]/E[-3]))

    return x

A = np.array([[4, 1, 0],
              [1, 3, 1],
              [0, 1, 2]])
x0 = np.array([1, 1, 1])
epsilon = 1e-9
x = EigenVector(A, x0, epsilon)
print(x)
