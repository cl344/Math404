N = 63787
B = 19

cand_num_target = 7
cand_num = 0
cand = B+1

B_list = get_prime(B);


num = []
primes = []

while(cand_num < cand_num_target):
	cand_sqr = (cand**2) % N
	cand_list = is_cand(cand_sqr, B_list)
	if(cand_list):
		cand_num = cand_num + 1
		num.append(cand)
		primes.append(cand_list)
	

	cand = cand + 1


select_list = select_cand(num, primes)

x = 1
y = 1

for i = 1 in range(len(num)):
	x = x*num[i]
	y = y*(B_list[i]**primes[i])

y = sqrt(y)

return gcd(x,y)





def get_prime(B):
	prime_list = []
	for num in range(B+1):
    # prime numbers are greater than 1
   		if num > 1:
   			check_prime = True
   			for i in range(2, num):
   				if (num % i) == 0:
   					check_prime = False
   					break


   			if(check_prime):
   				prime_list.append(num)

	return prime_list 


def is_cand(n, prime_list):
	exp_list = []
	for p in prime_list:
		occur = 0
		while(n%p == 0):
			occur = occur+1
			n = n/p
		exp_list.append(occur)

	return exp_list

def gcd(a,b): 
    if(b==0): 
        return a 
    else: 
        return gcd(b,a%b) 







