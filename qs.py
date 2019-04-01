from functions import *
import math

N = 63787
B = 19

# B_list = get_prime(B)  # List of primes

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
	return 19  # TODO


# ------------------- Sieving --------------------- #


def sieve(N, B):
	"""
	Perform the Quadratic Sieve.

	:param N: N
	:param B: Bound B for small primes
	:return: A factor of N, or -1 if there are none (Need to increase B)
	"""
	x_start = math.ceil(math.sqrt(N))
	x_upper_bound = 2 * x_start  # Stop the search when we reach this x (TODO)

	factor_base = get_prime(B)  # All primes up to B
	min_nums = len(factor_base)  # Minimum number of B-smooth numbers required

	nums = []  # All x's such that x^2-n is B-smooth
	exp_vecs = []  # Exponent vectors (raw) of each x^2-n

	# Accumulate x's and B-smooth remainders (x^2-n)
	for x in range(x_start, x_upper_bound + 1):
		x_sqr = (x ** 2) % N  # x^2 - n
		vec = is_cand(x_sqr, factor_base)
		if not vec:  # Check B-smoothness
			continue
		nums.append(x)
		exp_vecs.append(vec)
		if len(nums) >= min_nums:
			# Check if there exists subset of exponent vectors that sum to 0
			chosen_nums_sets, chosen_vecs_sets = find_subset(nums, exp_vecs)
			for i in range(len(chosen_nums_sets)):
				chosen_nums = chosen_nums_sets[i]  # subset of x's
				chosen_vecs = chosen_vecs_sets[i]
				factor = find_factor(N, chosen_nums, chosen_vecs, factor_base)
				if factor != -1:
					return factor


print(sieve(N, B))