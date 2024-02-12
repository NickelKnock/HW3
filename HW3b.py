import math

def t_distribution_pdf(t, v):
    """
    Probability density function for the t-distribution.
    """
    gamma_part = math.gamma((v + 1) / 2) / (math.sqrt(v * math.pi) * math.gamma(v / 2))
    return gamma_part * (1 + t ** 2 / v) ** (-(v + 1) / 2)

def Simpsons(f, a, b, n):
    """
    Numerical integration using Simpson's 1/3rd rule.
    """
    h = (b - a) / n
    s = f(a) + f(b)

    for i in range(1, n, 2):
        s += 4 * f(a + i * h)
    for i in range(2, n - 1, 2):
        s += 2 * f(a + i * h)

    return s * h / 3

def compute_cdf(v, t_value):
    """
    Computes the cumulative distribution function (CDF) for the t-distribution
    up to a given t-value.
    """
    # Integrate from a large negative number to t_value to simulate -infinity (actually -20) to t_value
    probability = Simpsons(lambda t: t_distribution_pdf(t, v), -20, t_value, 1000)
    return probability

def compute_t_value_from_Fz(v, desired_Fz):
    """
    Approximates the t-value for a given F(z) using numerical integration to compute the CDF.
    """
    # Approximating by looping through a range of t-values and find the one
    # whose computed CDF is closest to the desired F(z).
    min_diff = float('inf')
    closest_t_value = None
    # used chat to help me sor this section out
    for z in range(-500, 500):  # A range of t-values
        z /= 100  # Convert integer range to decimal t-values.
        cdf_value = compute_cdf(v, z)
        diff = abs(cdf_value - desired_Fz)
        if diff < min_diff:
            min_diff = diff
            closest_t_value = z

    return closest_t_value

def main():
    """
    Looping the input-output to ask the user as many times as needed for the degrees of freedom and the
    F(z) number needed. This was more elegant than hard coding the 7, 1, and 15 IMHO
    """
    while True:  # Start an infinite loop
        print("Enter the degrees of freedom (v):")
        v = int(input().strip())
        print("Enter the desired cumulative probability (F(z)):")
        desired_Fz = float(input().strip())

        closest_t_value = compute_t_value_from_Fz(v, desired_Fz)
        print(f"The t-value that corresponds to F(z)={desired_Fz} with {v} degrees of freedom is approximately {closest_t_value:.2f}")

        # Ask the user if they want to continue, this makes it so I user can use a capital Y because I kept
        # doing that and breaking my code...lol
        answer = input("Do you want to continue? (y or yes to continue): ").strip().lower()
        if answer not in ['y', 'yes']:
            break  # Exit the loop if the answer is not affirmative

if __name__ == "__main__":
    main()
