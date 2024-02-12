def make_diagonally_dominant(Aaug):
    """
    Attempts to make the given matrix diagonally dominant by swapping rows.
    :param Aaug: The augmented matrix of the system of linear equations.
    :return: The possibly modified augmented matrix with better diagonal dominance.
    """
    n = len(Aaug)  # Number of equations (rows in Aaug)
    for i in range(n):
        row_with_max_diagonal = max(range(i, n), key=lambda k: abs(Aaug[k][i]))
        if i != row_with_max_diagonal:
            Aaug[i], Aaug[row_with_max_diagonal] = Aaug[row_with_max_diagonal], Aaug[i]
    return Aaug

def GaussSeidel(Aaug, x, Niter=15):
    """
    Solves the system of linear equations using the Gauss-Seidel iterative method.
    :param Aaug: The augmented matrix [A|b] of the system.
    :param x: Initial guess for the solution vector.
    :param Niter: Maximum number of iterations to perform.
    :return: The estimated solution vector after Niter iterations.
    """
    n = len(Aaug)  # Number of equations
    for _ in range(Niter):
        for i in range(n):
            sum_ax = sum(Aaug[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (Aaug[i][-1] - sum_ax) / Aaug[i][i]
    return x

def separateAugmented(Aaug):
    """
    Separates the augmented matrix [A|b] into its components A and b.
    :param Aaug: The augmented matrix [A|b].
    :return: A tuple (A, b) where A is the coefficient matrix and b is the vector of constants.
    """
    A = [row[:-1] for row in Aaug]
    b = [row[-1] for row in Aaug]
    return A, b

def checkMatrixSoln(Aaug, x, verbose=True):
    """
    Checks the solution of a linear system given the augmented matrix and the solution vector.
    :param Aaug: The augmented matrix [A|b].
    :param x: The solution vector.
    :param verbose: If True, prints the comparison of the original b and the calculated b.
    :return: True if the solution is correct within a tolerance, False otherwise.
    """
    A, b = separateAugmented(Aaug)
    b_calc = [sum(a * xi for a, xi in zip(row, x)) for row in A]
    if verbose:
        print("Original b:", b)
        print("Calculated b:", b_calc)
    return all(abs(bi - bci) < 1e-5 for bi, bci in zip(b, b_calc))

def matrixMult(A, B):
    """
    Multiplies two matrices A and B.
    :param A: Matrix A.
    :param B: Matrix B.
    :return: The product of matrices A and B.
    """
    zip_b = list(zip(*B))
    return [[sum(ai * bj for ai, bj in zip(row_a, col_b)) for col_b in zip_b] for row_a in A]
