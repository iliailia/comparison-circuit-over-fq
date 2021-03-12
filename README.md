# Faster homomorphic comparison operations for BGV and BFV

This repository contains the LaTeX and C++ code of the homomorphic pattern matching algorithm from the paper ["Faster homomorphic comparison operations for BGV and BFV"](https://eprint.iacr.org/2021/315) by Ilia Iliashenko and Vincent Zucca.

## Installation guide
To build the code, install [HElib](https://github.com/homenc/HElib) (>= 2.0.0), go to the `code` folder and run 

    cmake .

and finally

    make

## How to use
### Integer comparison
To test the basic comparison of integers, use the following command
  
    ./comparison_circuit circuit_type p d m q l runs print_debug_info
    
where
+ `circuit_type` takes one of three values `U`, `B` or `T` corresponding to our univariate, bivariate circuits and the circuit of [Tan et al.](https://eprint.iacr.org/2019/332).
+ `p`: the plaintext modulus, must be a prime number.
+ `d`: the dimension of a vector space over the slot finite field.
+ `m`: the order of the cyclotomic ring.
+ `q`: the minimal bitsize of the ciphertext modulus in ciphertexts. The actual size of the modulus is automatically chosen by HElib.
+ `l`: the length of finite field vectors to be compared.
+ `runs`: the number of experiments.
+ `print_debug_info`: type `y` or `n` to show/hide more details on computation.
More details on these parameters can be found in Section 5 of the paper.

The following line performs 10 tests of our bivariate comparison circuit comparing vectors of length 3 over the finite field of order 49.
  
    ./comparison_circuit B 7 2 300 90 3 10 y

### Sorting
The sorting algorithm can be executed by the following line

    ./sorting_circuit p d m q l N runs print_debug_info
    
where the most arguments are analogous to the comparison circuit above. The additional argument is
+ `N`: the number of values to be sorted (see Section 4.1 of the paper).

By default, the sorting algorithm uses the univariate circuit. 

### Minimum of an array
To test the minimum function, run this command

    ./min_max_circuit p d m q l N T runs print_debug_info
    
where
+ `N`: the number of values in the input array.
+ `T`: the number of tournament stages (see Section 4.2 of the paper).

Other arguments have the same meaning as above. The univariate circuit is used here by default.
