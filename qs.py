from functions import *

N = 63787
B = 19

B_list = get_prime(B)  # List of primes

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
	# Placeholder variables (TODO)
	cand_num_target = 7  # Number of x^2-n candidates required
	cand_num = 0  # Number of x^2-n candidates found
	cand = B+1  # Current x

	num = []
	primes = []

	# Accumulate x^2-n's and B-smooth remainders
	while cand_num < cand_num_target:
		cand_sqr = (cand**2) % N  # x^2 - n
		cand_list = is_cand(cand_sqr, B_list)  # Check B-smoothness
		if cand_list:
			cand_num = cand_num + 1
			num.append(cand)
			primes.append(cand_list)
		cand = cand + 1

	# Choose a subset of x^2-n's whose product is also a square mod n
	# (i.e. exponent vectors sum up to 0 mod 2)
	select_list = select_cand(num, primes)

	# Compute x and y s.t. x^2 = y^2 (mod N) and try getting a factor
	# TODO: Back out if x = y (mod N)
	x = 1
	y = 1

	for i in range(len(num)):
		x = x*num[i]
		y = y*(B_list[i]**primes[i])

	y = sqrt(y)

	return gcd(x,y)


print(sieve(N, B))