
n = 0

# ------------------- Parameters and Setup --------------------- #

def get_B(n):
    """
    Obtain parameter B, the upper bound for small prime numbers to use
    in the quadratic sieve.

    From the paper, this should involve the following steps:
    - Choose parameter epsilon
    - Calculate bound X (and u?? Need to check paper regarding u)
    - Calculate bound B

    I imagine this part will be done first by trial and error, since
    we don't know what's a good way to pick epsilon. But as we actually
    finish up and run the code, we can do trial and error to pick the
    best epsilon that gives us the shortest runtime.

    :param n: Input integer
    :return: Bound B
    """
    return 0  # TODO

def get_primes(B):
    """
    Find all the primes up to B. There should be pi(B) of them.

    (This is probably one of the easier helper methods.)

    :param B: Bound for primes
    :return: List of all primes in [2,B]
    """
    return []  # TODO

# ------------------- Sieving --------------------- #

def solve_quad_congruence(n, p):
    """
    Helper method to solve the quadratic congruence x^2 - n = 0 (mod p).

    Since p<=B is small, a simple brute force should be enough for this
    function, but better algorithms may exist.
    Useful in sieving below.
    :param n: n
    :param p: p
    :return: List of solutions to x^2 - n = 0 (mod p).
    """
    return []  # TODO

def sieve(n, B):
    """
    Perform the Quadratic Sieve.

    The procedure should be as follows:
    - Find all primes up to B and solve all congruences x^2-n=0 (mod p)
      for each p.
        - At this point, if x=k1 or x=k2 (mod p) are solutions, we know
          that we can "mark" or "cross" out all such x's with regards to p.
    - Consider all integers x, starting from x=ceil(sqrt(n)).
        - Before we start iterating through all such x's, for each p,
          we have to find the first value to cross out, i.e. the
          smallest x such that x=k1 (mod p) and x>=ceil(sqrt(n)).
          Same for x=k2 (mod p).
        - When we cross out a number x^2-n, we keep dividing its value
          by p until the result is no longer divisible, to obtain the
          exponent of p in the factorization of x^2-n.
          We keep both the exponent and the new value (the quotient).
        - Although in theory, for each p, we should cross out all such
          x^2-n from the starting point to infinity, that has terrible
          runtime. Instead, we can do a "delayed crossing": For each p,
          cross out the first x first, and then when we reach that x,
          cross out the next one which would simply be x+p.
    - Now iterate through all x and consider each x^2-n.
    - For each x, check if x has been crossed out and reduced to 1.
      If yes, it's B-smooth. Add it to the list of B-smooth numbers.
      Then do the "delayed crossing".
    - When we have pi(B) B-smooth numbers, try the next step: Find a
      subset of them whose exponents sum up to an even number for each
      p, and attempt to find a solution from there. [Details in next
      section]

    This is probably the hardest part to implement. It might be the next
    thing on my TODO list.

    :param n: n
    :param B: Bound B for small primes
    :return: A factor of n, or -1 if there are none (Need to increase B)
    """
    return -1  # TODO

# ----- Gathering numbers and Producing Solution ----- #

def find_subset(nums, exp_vecs):
    """
    Given a list of numbers (x^2-n) and their corresponding exponent
    vectors, find a subset of them whose exponent vectors sum up to
    0 mod 2.

    Need to convert the exponent vectors to a matrix mod 2, and then
    use Gaussian Reduction.

    :param nums: List of numbers (of the form x^2-n), given as either
        x's or x^2-n's (depending on implementation)
    :param exp_vecs: Exponent vectors of each number (can be either
        list or dict, depending on implementation)
    :return: List of chosen numbers, or empty list if solution does
        not exist
    """
    return []  # TODO

def find_factor(nums, exp_vecs):
    """
    Gien a list of chosen numbers (x^2-n), and their corresponding
    exponent vectors that SUM UP TO 0 MOD 2, attempt to find a factor
    of n.

    This is done by:
    - Compute a, the square root of product of each x^2-n (mod n), by
      summing up all exponent vectors and dividing by 2.
    - Compute b, the product of all x's mod n.
    - We now have a^2=b^2 (mod n). Check if a=b (mod n), and if not,
      we have a factor gcd(a-b, n).

    This is the final attempt before successfully factorizing n.

    :param nums: List of numbers (of the form x^2-n), given as either
        x's or x^2-n's (depending on implementation)
    :param exp_vecs: Exponent vectors of each number (can be either
        list or dict, depending on implementation)
    :return: A factor of n, or -1 if unsuccessful
    """
    return -1  # TODO

# ------------------- Main --------------------- #

if __name__ == '__main__':
    n = int(input('Enter n:'))

    """
    Procedure:
    - Find bound B using get_B().
    - Call sieve() to hopefully get a factor of n.
    - If successful, repeat on both d and n/d (Repeat or continue??? 
      Wait to try it out). Otherwise, increase B and try again (restart or 
      build on existing structures??? Wait to try it out).
    """