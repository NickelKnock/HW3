#recurssive calls
#list comprehension
#turnary functions?

"""
These are notes and programs I used for reference during class but do not get called in the the homework
"""

def make_diagonally_dominant(Aaug):
    """
    Attempts to make the given matrix diagonally dominant by swapping rows.

    This function swaps rows in the matrix to move the largest (by absolute value) element
    in each column to the diagonal position, if possible.

    :param: Aaug (list of list): The augmented matrix of the system of linear equations.
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

    :param: Aaug (list of list of floats): The augmented matrix [A|b] of the system.
    :param: x (list of floats): Initial guess for the solution vector.
    :param: Niter (int): Maximum number of iterations to perform.

    :returns: list of floats: The estimated solution vector after Niter iterations.
    """
    n = len(Aaug)  # Number of equations
    for _ in range(Niter):
        for i in range(n):
            sum_ax = sum(Aaug[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (Aaug[i][-1] - sum_ax) / Aaug[i][i]
    return x


def main():
    """
    Main function to demonstrate the solving of systems of linear equations
    using the Gauss-Seidel method with preprocessing for diagonal dominance.

    It sets up two systems of linear equations, preprocesses their augmented matrices
    to enhance diagonal dominance, applies the Gauss-Seidel method, and prints the solutions.
    """
    # System definitions and initial guesses
    Aaug1 = [
        [3, 1, -1, 2],
        [1, 4, 1, 12],
        [2, 1, 2, 10]
    ]
    x_guess1 = [0, 0, 0]

    Aaug1 = make_diagonally_dominant(Aaug1)
    solution1 = GaussSeidel(Aaug1, x_guess1, Niter=15)
    print("Solution for preprocessed system 1:", solution1)

    Aaug2 = [
        [1, -10, 2, 4, 2],
        [3, 1, 4, 12, 12],
        [9, 2, 3, 4, 21],
        [-1, 2, 7, 3, 37]
    ]
    x_guess2 = [0, 0, 0, 0]

    Aaug2 = make_diagonally_dominant(Aaug2)
    solution2 = GaussSeidel(Aaug2, x_guess2, Niter=15)
    print("Solution for preprocessed system 2:", solution2)


if __name__ == "__main__":
    main()
