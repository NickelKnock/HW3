import math
"""
I used this program to test ideas and have a frame of reference for the things I changed.
This can be ignored
"""
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
    # Integrate from a large negative number to t_value to simulate -infinity to t_value
    probability = Simpsons(lambda t: t_distribution_pdf(t, v), -20, t_value, 1000)
    return probability

def compute_t_value_from_Fz(v, desired_Fz):
    """
    Approximates the t-value for a given F(z) using numerical integration to compute the CDF.
    """
    # We will use a simple approximation: loop through a range of t-values and find the one
    # whose computed CDF is closest to the desired F(z).
    min_diff = float('inf')
    closest_t_value = None
    for z in range(-500, 500):  # A range of t-values; in practice, you would use a finer step.
        z /= 100  # Convert integer range to decimal t-values.
        cdf_value = compute_cdf(v, z)
        diff = abs(cdf_value - desired_Fz)
        if diff < min_diff:
            min_diff = diff
            closest_t_value = z

    return closest_t_value

def main():
    print("Enter the degrees of freedom :")
    v = int(input().strip())
    print("Enter the desired cumulative probability (F(z)):")
    desired_Fz = float(input().strip())

    closest_t_value = compute_t_value_from_Fz(v, desired_Fz)
    print(f"The t-value that corresponds to F(z)={desired_Fz} with {v} degrees of freedom is approximately {closest_t_value:.2f}")

if __name__ == "__main__":
    main()
