from collections import deque

def test_replicate(n):
	'''
	Test the replication circuit

	Params:
		n (int): the array length

	Returns:
		(bool): True if the circuit is correct, False otherwise
	'''

	arr = deque([0 for i in range(2*n)]);
	arr[0] = 17;
	print("input array: ", arr); 

	orig = deepcopy(arr);

	K = ZZ(floor(log(n,2)));
	print("floored log_2(n) is: ", K);

	e = 1;

	for ind in range(K-1,-1,-1):
		tmp = deepcopy(arr);
		tmp.rotate(e);
		arr = deque([arr[i] + tmp[i] for i in range(2*n)]);
		e *= 2
		if ((ZZ(n).digits(2))[ind] == 1):
			arr.rotate(1);
			arr = deque([arr[i] + orig[i] for i in range(2*n)]);
			e += 1
	print(arr);

	for i in range(n):
		if arr[i] != arr[0]:
			return False;
		if arr[i+n] != 0:
			return False;

	return True;

def big_test_replicate(n):
	'''
	Test the rotation circuit for several vector length

	Params:
		n (int): the maximal vector length
	'''
	for i in range(1,n+1):
		if not test_replicate(i):
			print("Test has failed for i =", i);
			return;

	print("Passed!");

def ps_complexity(p):
	'''
	Outputs the number of multiplications needed for the Paterson-Stockmeyer algorithm to compute the comparison polynomial modulo p without optimization

	Params:
		p (int): prime number

	Returns:
		mul_num(int): number of multiplications
	'''
	if(not is_prime(p)):
		print("p must be prime");
		return 0;

	d = p-1;
	min_mul = p;
	min_k = 0;
	for k in range(1,d+1):
		mul = k-1;
		m = ceil(log(d/k,2));
		if m != 0:
			mul += (2^(m-1)-1 + m - 1);
			if((2^m-1) != (d/k)):
				mul += m;

		#multiplications to compute equality
		top_term = p-1;
		baby_index = top_term%k;
		giant_index = 0;
		if m != 0:
			giant_index = floor(top_term/k);
		if baby_index == 0:
			baby_index = k;
			giant_index -= 1;
		mul += sum(ZZ(giant_index).digits(2));
		if len(ZZ(giant_index).digits(2)) > m:
			mul += 1;

		if mul < min_mul:
			min_mul = mul;
			min_k = k;

	return min_k, min_mul;

def ps_complexity_optimized(p):
	'''
	Outputs the number of multiplications needed for the Paterson-Stockmeyer algorithm to compute the comparison polynomial modulo p with optimization

	Params:
		p (int): prime number

	Returns:
		mul_num(int): number of multiplications
	'''
	if(not is_prime(p)):
		print("p must be prime");
		return 0;

	if(p <= 3):
		print("p must be bigger than 3");
		return 0;

	d = (p-3)/2;
	min_mul = p;
	min_k = 0;
	for k in range(1,d+1):
		#optimization: add multiplication by x and squaring
		mul = 2
		
		mul += k-1;
		m = ceil(log(d/k,2));
		if m != 0:
			mul += (2^(m-1)-1 + m - 1);
			if((2^m-1) != (d/k)):
				mul += m;
		
		#multiplications to compute the top term
		top_term = (p-1)/2;
		baby_index = top_term%k;
		giant_index = floor(top_term/k);
		if baby_index == 0:
			baby_index = k;
			giant_index -= 1;
		mul += sum(ZZ(giant_index).digits(2));
		if m != 0:
			if len(ZZ(giant_index).digits(2)) > m:
				mul += 1;
		else:
			if len(ZZ(giant_index).digits(2)) > 1:
				mul += 1;
				print("Danger");

		#print(k, mul);
		if mul < min_mul:
			min_mul = mul;
			min_k = k;

	return min_k, min_mul;

def test_ps_complexity(n):
	'''Outputs the multiplicative complexity of the comparison circuit with and without optimizations for all primes up to n

	Params:
		n(int): range
	'''
	p = next_prime(3);
	while(p < n):
		min_k, min_mul = ps_complexity(p);
		min_k_optimized, min_mul_optimized = ps_complexity_optimized(p);
		print("p:", p, (min_k, min_mul), (min_k_optimized, min_mul_optimized), min_mul >= min_mul_optimized);
		p = next_prime(p);

def univariate_circuit(p):
	'''Returns the coefficients of the univariate comparison polynomial modulo p 

	Params:
		p(int): prime modulus

	Returns:
		poly(polynomial over F_p): univariate circuit as polynomial 
	'''
	if(not is_prime(p)):
		print("p must be prime");
		return 0;

	if(p < 3):
		print("p must be bigger or equal to 3");
		return 0;

	F = GF(p);
	R.<z> = PolynomialRing(F);
	poly = R(0);
	half_p = (p-1) / 2;
	for a in range(-half_p, 0):
		poly += (1 - (z-F(a))^(p-1));
	return poly;
