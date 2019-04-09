import random
from itertools import chain


def get_prime(B):
	"""
	Find all the primes up to B. There should be pi(B) of them.

	:param B: Bound for primes
	:return: List of all primes in [2,B]
	"""
	# Use the same "delayed crossing" idea:
	#	- Initialize a list with all numbers
	#	- Examine each number x one by one:
	#		- If it's not crossed out, it's prime, leave a flag "x" on 2x
	#		- If it's crossed out, take all flags f and push the flax x to x+f
	nums = list(range(2, B+1))
	primes = []
	flags = {}  # crossed out number: [list of primes that crosses it out]
	for x in nums:
		if x not in flags.keys():
			primes.append(x)
			if x + x not in flags:
				flags[x + x] = []
			flags[x + x].append(x)
		else:
			for f in flags[x]:
				if x + f not in flags:
					flags[x + f] = []
				flags[x + f].append(f)
	return primes


def gcd(a, b):
	"""
	Calculate the Greatest Common Divisor.
	:param a: a
	:param b: b
	:return: gcd(a,b)
	"""
	if a < b:
		return gcd(b, a)
	while b != 0:
		a, b = b, a % b
	return a


def transpose(matrix):
	"""
	Transpose a matrix so columns become rows, makes list comp easier to
	work with.
	:param matrix: Original matrix
	:return: Transpose of matrix
	"""
	return [[row[j] for row in matrix] for j in range(len(matrix[0]))]


def pow(x, n):
	"""
	Calculate x^n.
	https://stackoverflow.com/questions/5246856/how-did-python-implement-the-built-in-function-pow
	:param x: x
	:param n: n
	:return: x^n
	"""
	ans = 1
	while n:
		if n % 2:
			ans *= x
		x *= x
		n //= 2
	return ans


def pow_mod(x, n, p):
	"""
	Calculate x^n mod p.
	https://stackoverflow.com/questions/5246856/how-did-python-implement-the-built-in-function-pow
	:param x: x
	:param n: n
	:param p: p
	:return: x^n mod p
	"""
	ans = 1
	while n:
		if n % 2:
			ans = ans * x % p
		x = x * x % p
		n //= 2
	return ans


def solve_quad_congruence(n, p):
	"""
	Helper method to solve the quadratic congruence x^2 - n = 0 (mod p),
	using the Tonelli–Shanks algorithm.

	https://en.wikipedia.org/wiki/Tonelli–Shanks_algorithm
	https://rosettacode.org/wiki/Tonelli-Shanks_algorithm

	Useful in sieving below.
	:param n: n
	:param p: p
	:return: Tuple of solutions to x^2 - n = 0 (mod p), or (-1,-1) if none
		exists
	"""
	n %= p
	if n == 0:
		return 0, 0  # x = 0 (We just factorized n!!)
	if p == 2:
		return 1, 1

	# First, check n is indeed a square
	# i.e. n^((p-1)/2) = 1 (mod p)
	if pow_mod(n, (p-1) // 2, p) != 1:
		return -1, -1

	# find Q and S such that p − 1 = Q*2^S with Q odd
	Q = p - 1
	S = 0
	while not Q % 2:
		Q //= 2
		S += 1
	if S == 1:
		# p = 3 (mod 4), use formula for solutions
		r = pow_mod(n, (p+1) // 4, p)
		return r, p-r

	# Search for z which is a quadratic non-residue
	z = 0
	for z_cand in range(2, p):
		if pow_mod(z_cand, (p-1) // 2, p) == p - 1:
			z = z_cand
			break
	if z == 0:  # Unknown error
		return -1, -1

	# Set variables
	M = S
	c = pow_mod(z, Q, p)
	t = pow_mod(n, Q, p)
	R = pow_mod(n, (Q+1) // 2, p)

	# Loop
	while True:
		if t == 0:  # n == 0 (mod p)
			return 0, 0
		if t == 1:
			return R, p-R
		i = 0
		t_pow = t
		while t_pow != 1 and i < M:
			t_pow = t_pow * t_pow % p
			i += 1
		if i == M:
			return -1, -1
		b = c  # b = c^(2^(M-i-1))
		for j in range(M-i-1):
			b = b * b % p
		M = i
		c = b * b % p
		t = t * b * b % p
		R = R * b % p

	# Unreachable, but just in case
	for x in range(1, p // 2 + 1):  # Check [1, p/2] since -x is a solution
		if x * x % p == n:
			return x, p-x

	return -1, -1


# ----- Gathering numbers and Producing Solution ----- #


def build_matrix_from_vecs(exp_vecs):
	"""
	Construct a matrix of exponent vectors (from B-smooth x^2-n) mod 2,
	with the raw exponent vectors provided.
	:param exp_vecs: Exponent vectors of each number as list
	:return: Matrix with rows as exponent vectors mod 2
	"""
	return [[i % 2 for i in v] for v in exp_vecs]


def gauss_elim(M):
	"""
	Perform Gaussian Elimination to solve Mx=0 (mod 2).
	(This method only works in Z/2Z.)

	:param M: Matrix
	:return: List of all non-zero solutions
	"""
	RANDOM_TRIALS = 5  # Random trials to assign free variables a value
	marks = [False]*len(M[0])  # If a pivot has been constructed in this column
	pivot_row = [-1] * len(M[0])  # For each pivot column, the corresponding
									# row for which value is 1

	# Find a pivot in all columns, and adjust the other rows accordingly
	for i in range(len(M)): #do for all rows
		row = M[i]
		for j in range(len(row)): #search for pivot
			if row[j] == 1:
				# Make this a pivot column
				marks[j] = True
				pivot_row[j] = i

				# Add this row to all other rows with this column being 1
				for k in chain(range(0, i), range(i+1, len(M))):
					if M[k][j] == 1:
						for l in range(len(M[k])):
							M[k][l] = M[k][l] ^ row[l]  # (M[k][l] + row[l]) % 2
				break

	# Find all free columns: Each column corresponds to a free variable t,
	# and values of the pivot column variables can then be deduced.
	# This currently takes O(2^(free columns)) time.
	free_cols = [i for i in range(len(M[0])) if not marks[i]]
	if not free_cols:
		return []  # Only has trivial solution: fail

	def find_solution(assignment):
		"""
		Compute a solution given an assignment of free variables.
		:param assignment: Value of each free variables (corresponding to
			free_cols)
		:return: Solution to Mx=0 (mod 2)
		"""
		sol = []
		free_index = 0
		for i in range(len(M[0])):
			if marks[i]:  # Deduce values of pivot variables
				sol.append(sum(M[pivot_row[i]][free_cols[j]] * assignments[j]
								for j in range(len(free_cols))) % 2)
			else:
				sol.append(assignment[free_index])
				free_index += 1
		return sol

	# Search over 5 random assignments of free variables
	solutions = []
	# assignments = [0] * (len(free_cols) - 1) + [1]  # Value of each free variable
	# while assignments[0] <= 1:
	# 	solutions.append(find_solution(assignments))
	# 	assignments[len(assignments) - 1] += 1
	# 	i = len(assignments) - 1
	# 	while i > 0 and assignments[i] > 1:
	# 		assignments[i - 1] += 1
	# 		assignments[i] -= 2
	# 		i -= 1
	for i in range(RANDOM_TRIALS):
		assignments = [0] * len(free_cols)
		while assignments == [0] * len(free_cols):
			for j in range(len(free_cols)):
				assignments[j] = random.randint(0, 1)
		solutions.append(find_solution(assignments))

	return solutions


def find_subset(nums, exp_vecs):
	"""
	Given a list of numbers (x^2-n) and their corresponding exponent
	vectors, find all possible subsets of them whose exponent vectors
	sum up to 0 mod 2.

	Need to convert the exponent vectors to a matrix mod 2, and then
	use Gaussian Reduction.

	:param nums: List of numbers (of the form x^2-n), given as x's
	:param exp_vecs: Exponent vectors of each number (can be either
		list or dict, depending on implementation)
	:return: - List of lists (subsets) of chosen numbers, or empty list
		if solution does not exist
				e.g. [set1, set2, ...] where set1 = [x1, x2, ...]
			 - List of lists (subsets) of exponent vectors corresponding
		to the chosen vectors, or empty list if solution does not exist
	"""
	M = build_matrix_from_vecs(exp_vecs)
	M = transpose(M)
	sols = gauss_elim(M)
	chosen_nums = []
	chosen_vecs = []
	for sol in sols:
		chosen_nums.append([nums[i] for i in range(len(sol)) if sol[i] == 1])
		chosen_vecs.append(
			[exp_vecs[i] for i in range(len(sol)) if sol[i] == 1])
	return chosen_nums, chosen_vecs


def find_factor(N, nums, exp_vecs, factor_base):
	"""
	Given a list of chosen numbers (x^2-n), and their corresponding
	exponent vectors that SUM UP TO 0 MOD 2, attempt to find a factor
	of n.

	This is done by:
	- Compute a, the product of all x's mod n.
	- Compute b, the square root of product of each x^2-n (mod n), by
	  summing up all exponent vectors and dividing by 2.
	- We now have a^2=b^2 (mod n). Check if a=b (mod n), and if not,
	  we have a factor gcd(a-b, n).

	This is the final attempt before successfully factorizing n.

	:param nums: List of numbers (of the form x^2-n), given as x's
	:param exp_vecs: Exponent vectors of each number (can be either
		list or dict, depending on implementation)
	:param factor_base: List of available prime factors
	:param N: n
	:return: A factor of n, or -1 if unsuccessful
	"""
	a = 1
	for k in nums:
		a = a * k % N

	sum_exps = [0] * len(exp_vecs[0])  # Sum of all exponent vectors (for product)
	for v in exp_vecs:
		sum_exps = [x+y for x, y in zip(sum_exps, v)]

	b = 1
	for i in range(len(sum_exps)):
		b = b * pow(factor_base[i], (sum_exps[i] // 2)) % N

	return -1 if (a - b) % N == 0 or (a + b) % N == 0 else gcd(a - b + N, N)

