import DoolittleMethod as DM
import hw2c # I didn't end up using this exact program but used it for reference

def main():
    """
    It seems that both matrices are symmetric positive definite, so I'm not falling back to Doolittles method
    but it is there to check to verify the state of the matrix used.
    The Cholesky method is implemented in the Doolittle method code because I had trouble calling is from here
    to run in the Doolittle code, not sure why
    """

    # Define matrices and vectors for testing
    A1 = [[1, -1, 3, 2],
          [-1, 5, -5, -2],
          [3, -5, 19, 3],
          [2, -2, 3, 21]]
    b1 = [15, -35, 94, 1]

    A2 = [[4, 2, 4, 0],
          [2, 2, 3, 2],
          [4, 3, 6, 3],
          [0, 2, 3, 9]]
    b2 = [20, 36, 60, 122]

    # Process and solve using the appropriate method for each matrix
    for A, b in [(A1, b1), (A2, b2)]:
        if DM.is_symmetric(A) and DM.is_positive_definite(A):
            print("Matrix is symmetric and positive definite. Using Cholesky method.")
            L = DM.cholesky_decomposition(A)
            # BackSolve for y in Ly = b
            y = DM.BackSolve(L, b, UT=False)
            # Transpose L to get LT, then BackSolve for x in LT x = y
            LT = [list(i) for i in zip(*L)]
            x = DM.BackSolve(LT, y, UT=True)
        else:
            print("Matrix does not meet Cholesky criteria. Using Doolittle method.")
            aug = [A[i] + [b[i]] for i in range(len(A))]  # Create augmented matrix [A|b]
            x = DM.Doolittle(aug)

        x_rounded = [round(xi, 3) for xi in x]
        print("Solution:", x_rounded)

if __name__ == "__main__":
    main()
