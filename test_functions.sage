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
	print "input array: ", arr; 

	orig = deepcopy(arr);

	K = ZZ(floor(log(n,2)));
	print "floored log_2(n) is: ", K;

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
	print arr;

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
			print "Test has failed for i =", i;
			return;

	print "Passed!"