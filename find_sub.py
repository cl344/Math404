import numpy as np


def find_subset(nums, exp_vecs):

	#nums - list of x^2-n that are B-smooth
	#exp_vecs - list (or dict) of exponent vector lists (assume dict here)
	# ie exp_vecs[x] is the list of exps for x^2-n

	#reduce the dict to mod2
	for n in nums :
		mod2_vec = [x%2 for x in exp_vecs[n]]
		exp_vecs.pop(n)
		exp_vecs.update({n : mod2_vec})

	print(exp_vecs);

	# find the solution to Ax=b where
	# we assume we have pi(B)+1 vectors
	# b is a random chosen exp_vec as a column
	# A is the square matrix of the remaining exp_vecs as columns
	# x is then the coef of each exp_vec (mod2)

	#TODO : find out how matricies work
	A = np.matrix()
	print(A)
	#create the matricies
	#for n in nums :
	
	return -1



n = 1649;
B = 5;

nums = [10, 11, 12, 13];
exp_vecs = {nums[0] : [2, 3, 2], nums[1] : [1, 0, 3], nums[2] : [0, 4, 1], nums[3] : [11, 2, 1]}

print(nums)
print(exp_vecs)

find_subset(nums, exp_vecs)



