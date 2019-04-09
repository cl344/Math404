import time

from functions import *
import math
from math import *

N = 3744843080529615909019181510330554205500926021947  # 63787


def sieve(N, B):
	"""
	Perform the Quadratic Sieve, with the bound B already given.

	:param N: N
	:param B: Bound B for small primes
	:return: A factor of N, or -1 if there are none (Need to increase B)
	"""
	print('B = %d' % B)
	x_start = math.ceil(math.sqrt(N))
	x_upper_bound = int(x_start * 1.01)  # Stop the search when we reach this x

	factor_base = get_prime(B)  # All primes up to B
	print('%d primes' % len(factor_base))

	nums = []  # All x's such that x^2-n is B-smooth
	exp_vec_dicts = []  # Exponent vectors (as dict of p:exp) of each x^2-n

	def try_solve():
		"""
		Try to solve for a factor given the current exponent vectors.
		:return: A factor of N, or -1 if it doesn't exist
		"""
		# Convert dicts to exponent vectors
		exp_vecs = []
		for d in exp_vec_dicts:
			v = [0] * len(factor_base)
			for p, exp in d.items():
				v[p_to_index[p]] = exp
			exp_vecs.append(v)

		# Check if there exists subset of exponent vectors that sum to 0
		chosen_nums_sets, chosen_vecs_sets = find_subset(nums, exp_vecs)
		for i in range(len(chosen_nums_sets)):
			chosen_nums = chosen_nums_sets[i]  # subset of x's
			chosen_vecs = chosen_vecs_sets[i]
			factor = find_factor(N, chosen_nums, chosen_vecs, factor_base)
			if factor != -1:
				return factor
		return -1

	flags = {}  # crossed out number: [list of primes that crosses it out]
	# Pick initial numbers to cross out for each p
	effective_factor_base = []
	for p in factor_base:
		x1, x2 = solve_quad_congruence(N, p)
		if x1 == 0:
			return p  # n == 0 (mod p)
		if x1 == -1:
			continue
		effective_factor_base.append(p)
		# Pick the smallest x >= x_start s.t. x == x1 (mod p)
		for res in {x1, x2}:  # CAREFUL! x1=x2 when p=2
			x = x_start // p * p  # Largest multiple of p that's <= x_start
			x += res
			if x < x_start:
				x += p
			if x not in flags:
				flags[x] = []
			flags[x].append(p)
	factor_base = effective_factor_base
	min_nums = len(factor_base)  # Minimum number of B-smooth numbers required

	p_to_index = {}  # Lookup index of p in factor_base
	for i in range(len(factor_base)):
		p_to_index[factor_base[i]] = i
	print('%d useful primes' % len(factor_base))

	# Accumulate x's and B-smooth remainders (x^2-n)
	x_sqr = (x_start - 1) ** 2 % N
	for x in range(x_start, x_upper_bound + 1):
		x_sqr = (x_sqr + 2 * x - 1) % N  # x^2 = (x-1)^2 + 2x - 1
		if x_sqr == 0:
			return x  # x^2 == 0 (mod n)
		if x % 10000000 == 0:  # Print progress
			print('  x = %d' % x)
			print('  # of B-smooth numbers found = %d' % len(nums))

		# Check B-smoothness by taking the crossed out p's and divide x
		if x not in flags and x_sqr != 1:
			continue
		x_sqr_remain = x_sqr
		vec_dict = {}
		for p in flags[x]:
			exp = 0
			while x_sqr_remain % p == 0:
				exp += 1
				x_sqr_remain = x_sqr_remain // p
			vec_dict[p] = exp
			if x + p not in flags:  # Push back flags
				flags[x + p] = []
			flags[x + p].append(p)
		del flags[x]
		if x_sqr_remain != 1:  # Not B-smooth
			continue

		nums.append(x)
		exp_vec_dicts.append(vec_dict)
		if len(nums) >= min_nums:  # DEBUG
			print('Found %d B-smooth numbers' % len(nums))
			factor = try_solve()
			if factor != -1:
				return factor

	return try_solve()


def opt_bound(N):
	"""
	Given N, find the range of possible B values according to the
	formula given below.
	:param N: N
	:return: - Lower bound for B
			- Upper bound for B
	"""
	low = pow(exp(sqrt(log(N)*log(log(N)))), 1/2)
	up = low**3
	return int(low),int(up)


def sieve_auto(N):
	"""
	Perform the Quadratic Sieve without the bound B.

	This function first chooses B from the range given by opt_bound,
	then starting from the lower bound for B, pick a B value and try
	factorizing using sieve(N, B). If it fails to find a factor,
	increase B by a factor of 10 and try again, until the upper bound
	is reached.

	:param N: N
	:return: A factor of N, or -1 if there are none (Need to increase B)
	"""
	start_time = time.time()

	B, upper = opt_bound(N)
	result = sieve(N, B)
	while result == -1 and B <= upper:
		B = B * 10
		result = sieve(N, B)

	if B > upper:
		return -1

	end_time = time.time()
	print('Time taken (in seconds): %.3f' % (end_time - start_time))
	return result
	
		
factor_q = sieve_auto(N)
factor_p = N/factor_q
print('factors of %d are %d and %d' % (N, factor_q, factor_p))