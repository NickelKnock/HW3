# DoolittleMethod.py
from copy import deepcopy as dcpy
from math import sqrt
"""
I gutted and reworked a lot of things but started with Smay's code provided from his Github.
I used chatGPT to reformat several instances of my code for efficiency on this part
"""

def is_symmetric(A):
    """
    Check if the matrix A is symmetric.
    :param A: nxn matrix
    :return: Boolean indicating if A is symmetric
    """
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            if A[i][j] != A[j][i]:
                return False
    return True


def is_positive_definite(A):
    """
    Check if the matrix A is positive definite using Cholesky decomposition attempt.
    :param A: nxn matrix
    :return: Boolean indicating if A is positive definite
    """
    try:
        cholesky_decomposition(A)
        return True
    except:
        return False


def cholesky_decomposition(A):
    """
    Perform Cholesky decomposition on a positive definite matrix A.
    :param A: nxn positive definite matrix
    :return: Lower triangular matrix L such that A = L * L^T
    """
    n = len(A)
    L = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            sum_k = sum(L[i][k] * L[j][k] for k in range(j))

            if i == j:  # Diagonal elements
                L[i][j] = sqrt(max(A[i][i] - sum_k, 0))
            else:
                L[i][j] = (1.0 / L[j][j]) * (A[i][j] - sum_k) if L[j][j] != 0 else 0
    return L


def LUFactorization(A):
    """
    Lower-Upper factorization part of Doolittle's method.
    :param A: a nxn matrix
    :return: a tuple with (L, U)
    """
    n = len(A)
    U = [([0] * n if r != 0 else A[0][:]) for r in range(n)]
    L = [[(1 if r == c else 0) for c in range(n)] for r in range(n)]
    for j in range(1, n):
        for i in range(j, n):
            sum = 0
            for k in range(j):
                sum += L[i][k] * U[k][j]
            L[i][j] = (A[i][j] - sum) / U[j][j]
        for i in range(j, n):
            sum = 0
            for k in range(j):
                sum += L[j][k] * U[k][i]
            U[j][i] = A[j][i] - sum
    return L, U


def BackSolve(A, b, UT=True):
    """
    Backsolving algorithm for a matrix and b vector.
    :param A: A triangular matrix (Upper or Lower)
    :param b: the RHS of a matrix equation Ax=b
    :param UT: boolean for upper triangular (True) or lower triangular (False)
    :return: the solution vector x, from Ax=b
    """
    n = len(b)
    x = [0] * n
    if UT:
        for i in reversed(range(n)):
            s = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x[i] = (b[i] - s) / A[i][i]
    else:
        for i in range(n):
            s = sum(A[i][j] * x[j] for j in range(i))
            x[i] = (b[i] - s) / A[i][i]
    return x

def separateAugmented(Aaug):
    """
    Separates augmented matrix into A and b.
    It goes beyond the scope of our implementation.
    """
    A = [row[:-1] for row in Aaug]
    b = [row[-1] for row in Aaug]
    return A, b


def Doolittle(Aaug):
    """
    Solves the matrix equation Ax=b using Doolittle's LU factorization.
    :param Aaug: the augmented matrix [A|b]
    :return: the solution vector x
    """
    # A, b = separateAugmented(Aaug)  # Assume implementation of separateAugmented exists
    L, U = LUFactorization(A)
    y = BackSolve(L, b, UT=False)
    x = BackSolve(U, y, UT=True)
    return x
