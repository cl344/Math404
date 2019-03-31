from math import fabs, ceil, sqrt, exp, log
import random
from itertools import chain


def get_prime(B):
	"""
	Find all the primes up to B. There should be pi(B) of them.

	:param B: Bound for primes
	:return: List of all primes in [2,B]
	"""
	# TODO (Optional): Improve efficiency using sieving idea (??)
	prime_list = []
	for num in range(B+1):
		# prime numbers are greater than 1
		if num > 1:
			check_prime = True
			for i in range(2, num):
				if (num % i) == 0:
					check_prime = False
					break
			if check_prime:
				prime_list.append(num)
	return prime_list


def is_cand(n, prime_list):
	"""
	Checks if n can be factored as a product of all primes in
	prime_list, and returns the exponent vector
	:param n: Number to be factored
	:param prime_list: List of available primes
	:return: List of exponents in the factorization for each prime in
		prime_list (0 if n is not a multiple)
	TODO: Shouldn't there be some sort of indication if n can't be
		completely factorized?
	"""
	exp_list = []
	for p in prime_list:
		occur = 0
		while n % p == 0:
			occur = occur + 1
			n = n/p
		exp_list.append(occur)
	return exp_list


#def select_cand(mat):


def factorize(n, factor_base):
	"""
	Factorize n (mod N) into a product of vectors in factor_base.

	This method currently uses trial division.

	:param n: Number to be factorized
	:param factor_base: List of available prime factors
	:return: Prime factors of n, with repetition if necessary
	TODO: Shouldn't there be some sort of indication if n can't be
		completely factorized?
	"""
	N = 63787  # Placeholder (TODO)
	n = n % N
	print("function entered")
	factors = []
	for p in factor_base:
		while(n!=0 and n%p == 0):
			factors.append(p)
			n =n/p
	return factors


def gcd(a, b):
	"""
	Calculate the Greatest Common Divisor.
	:param a: a
	:param b: b
	:return: gcd(a,b)
	"""
	if(b==0):
		return a
	else:
		return gcd(b,a%b)


def transpose(matrix):
	"""
	Transpose a matrix so columns become rows, makes list comp easier to
	work with.
	:param matrix: Original matrix
	:return: Transpost of matrix
	"""
	return [[row[j] for row in matrix] for j in range(len(matrix[0]))]


def build_matrix(smooth_nums,factor_base):
	"""
	Construct a matrix of exponent vectors (from B-smooth x^2-n) mod 2.
	:param smooth_nums: List of B-smooth x^2-n's
	:param factor_base: List of available prime factors
	:return: - Whether one of the B-smooth numbers is a perfect square
			   mod n itself
			 - If first argument is True, the perfect square;
			   Otherwise, the entire matrix
	TODO: Pass in exponent vetor as a parameter instead of
		re-factorizing everything?
	"""

	# generates exponent vectors mod 2 from previously obtained smooth numbers, then builds matrix
	M = []
	#factor_base.insert(0,-1)


	for n in smooth_nums:
		exp_vector = [0]*(len(factor_base))
		n_factors = factorize(n,factor_base)
		print(n,n_factors)
		for i in range(len(factor_base)):
			if factor_base[i] in n_factors:
				exp_vector[i] = (exp_vector[i] + n_factors.count(factor_base[i])) % 2

		#print(n_factors, exp_vector)
		if 1 not in exp_vector: #search for squares
			return True, n
		else:
			pass

		M.append(exp_vector)

	#print("Matrix built:")
	#mprint(M)
	return(False, transpose(M))


def gauss_elim(M):
	"""
	Perform Gaussian Elimination to solve Mx=0 (mod 2).
	(This method only works in Z/2Z.)

	:param M: Matrix
	:return:
	"""
	marks = [False]*len(M[0])
	for i in range(len(M)): #do for all rows
		row = M[i]
		#print(row)

		for num in row: #search for pivot
			if num == 1:
				#print("found pivot at column " + str(row.index(num)+1))
				j = row.index(num) # column index
				marks[j] = True

				for k in chain(range(0,i),range(i+1,len(M))): #search for other 1s in the same column
					if M[k][j] == 1:
						for i in range(len(M[k])):
							M[k][i] = (M[k][i] + row[i])%2
				break

	print(marks)
	M = transpose(M)
	#mprint(M)

	sol_rows = []
	for i in range(len(marks)): #find free columns (which have now become rows)
		if marks[i]== False:
			free_row = [M[i],i]
			sol_rows.append(free_row)

	if not sol_rows:
		return("No solution found. Need more smooth numbers.")

	print("Found {} potential solutions".format(len(sol_rows)))
	return sol_rows,marks,M



N = 63787
# print(get_prime(19))
# print(is_cand((439**2)%N,get_prime(19)))
# print(gcd(25,85))


smooth_nums = [439**2, 441**2, 444**2, 445**2, 447**2, 449**2]
factor_base = [2,3,5,7,11,13,17]

is_sqr, M = build_matrix(smooth_nums, factor_base)
print(M)
# print([False]*len(M[0]))
# sol_rows,marks,matrix = gauss_elim(M)
# #print(sol_rows)
# print(M)





