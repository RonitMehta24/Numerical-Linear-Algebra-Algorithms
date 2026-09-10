import numpy as np
    

def isValidMatrix(A):
    return True

class LinearSystem:
    def __init__(self, A, b):
        if not isValidMatrix(A):
            raise ValueError("Please Enter a Valid Matrix")
        if len(A) != len(b):
            raise ValueError("Please Enter a Valid Linear System")
        self.rows = len(A)
        self.cols = len(A[0])
        self.A = A.astype(float)
        self.b = b.astype(float)

    def __str__(self):
        return(str(self.A) + "=" + '\n' + str(self.b))


    # Does LU Decomposition in case we want to use it later for different rhs, so the lower part of the matrix isnt wasted
    def GaussianElimination(self):
        for k in range(self.rows-1):
            den = self.A[k][k]
            for i in range(k+1, self.rows):
                m = self.A[i][k]/den
                self.A[i][k] = m
                for j in range(k+1, self.cols):
                    self.A[i][j] -= m*self.A[k][j]
                self.b[i] -= m*self.b[k]
        return

#I do not understand why the indexing for the back and forward substitution works or if it even works in the general case. I have no idea, and i need to figure out what is going on
# Doing indexing to the solution in the forwardsubstitution but not in the backwardsubstitution works for inverting





    
    def BackSubstitution(self, r = []):
        if len(r) == 0:
            r = [i for i in range(self.rows)]
        # Assumes A is upper Trianglular with respect to the given indexing r of rows
        x = np.zeros(self.rows)
        for i in range(self.rows - 1, -1, -1):
            sol = self.b[r[i]]
            for j in range(i + 1, self.cols):
                sol -= x[j]*self.A[r[i]][j]
            x[i] = sol/self.A[r[i]][i]
        return x

    def ForwardSubstitution(self, r = []):
        if len(r) == 0:
            r = [i for i in range(self.rows)]
        # Assumes A is a lower triangluar(with respect to the given indexing r of rows) with 1s on the main diagonal, but only accesses the lower triangular part of the matrix
        x = np.zeros(self.rows)
        for i in range(self.rows):
            sol = self.b[r[i]]
            for j in range(i):
                sol -= x[r[j]]*self.A[r[i]][j]
            x[r[i]] = sol
        return x




    
            
    # Since we are doing partial pivoting, we need to return the row indexing vector r so the user knows which order to do backsubstitution in
    def GEPP(self):
        r = [i for i in range(self.rows)] # Initialize the index vector
        for k in range(self.rows - 1):
            maxv = self.A[r[k]][k]
            maxi = k
            for l in range(k + 1, self.rows):
                if self.A[r[l]][k] > maxv:
                    maxi = l
                    maxv = self.A[r[l]][k]
            r[k], r[maxi] = r[maxi], r[k]
            den = self.A[r[k]][k]
            for i in range(k + 1, self.rows):
                m = self.A[r[i]][k]/den
                self.A[r[i]][k] = m
                for j in range(k + 1, self.cols):
                    self.A[r[i]][j] -= m*self.A[r[k]][j]
                self.b[r[i]] -= m*self.b[r[k]]
        return r

    def Solve_LU(self, rhs):
        # rhs is a list of vectors, so rhs[i] is the ith rhs
        # Returns a list of vectors
        rhs_copy = np.transpose(rhs)
        self.GaussianElimination()
        sol = []
        for r in rhs_copy:
            self.b = r
            self.b = self.ForwardSubstitution()
            sol.append(self.BackSubstitution())
        return np.transpose(np.array(sol))

    def Solve_LU_GEPP(self, rhs):
        rhs_copy = np.transpose(rhs)
        indexing = self.GEPP()
        sol = []
        for r in rhs_copy:
            self.b = r
            self.b = self.ForwardSubstitution(indexing)
            sol.append(self.BackSubstitution(indexing))
        return np.transpose(np.array(sol))

def MatrixInverse(A):
    if len(A) == 0:
        raise ValueError("The matrix to be inverted must not be empty.")
    if len(A) != len(A[0]):
        raise ValueError("The matrix to be inverted must not be a square matrix and cannot be inverted")
    Sys = LinearSystem(A, np.zeros(len(A)))
    return Sys.Solve_LU(np.eye(len(A)))

def MatrixInverseGEPP(A):
    if len(A) == 0:
        raise ValueError("The matrix to be inverted must not be empty.")
    if len(A) != len(A[0]):
        raise ValueError("The matrix to be inverted must not be a square matrix and cannot be inverted")
    Sys = LinearSystem(A, np.zeros(len(A)))
    return Sys.Solve_LU_GEPP(np.eye(len(A)))
    
                       
    
        
               
                    
        
            

                
                    
    
