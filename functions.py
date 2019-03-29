
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

#def select_cand(mat):


def gcd(a,b): 
    if(b==0): 
        return a 
    else: 
        return gcd(b,a%b) 


N = 63787
print(get_prime(19))
print(is_cand((439**2)%N,get_prime(19)))
print(gcd(25,85))

