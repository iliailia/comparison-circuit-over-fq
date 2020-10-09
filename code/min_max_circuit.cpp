/* Copyright (C) 2019 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */
#include <iostream>
#include <time.h>
#include <random>

#include <helib/helib.h>
#include <helib/debugging.h>
#include <helib/Context.h>
#include <helib/polyEval.h>
#include "tools.h"
#include "comparator.h"

using namespace std;
using namespace NTL;
using namespace helib;

// the main function that takes 8 arguments (type in Terminal: ./sorting_circuit argv[1] argv[2] argv[3] argv[4] argv[5] argv[6] argv[7] argv[8])
// argv[1] - the plaintext modulus
// argv[2] - the dimension of a vector space over a finite field
// argv[3] - the order of the cyclotomic ring
// argv[4] - the bitsize of the ciphertext modulus in ciphertexts (HElib increases it to fit the moduli chain). The modulus used for public-key generation
// argv[5] - the length of vectors to be compared
// argv[6] - the number of values to be sorted
// argv[7] - the number of consecutive comparisons
// argv[8] - the number of experiment repetitions
// argv[9] - print debug info (y/n)

// some parameters for quick testing
// 7 1 75 90 1 4 1 10 y
// 7 1 300 90 1 6 2 10 y
// 17 1 145 120 1 7 2 10 y
int main(int argc, char *argv[]) {
  if(argc < 10)
  {
   throw invalid_argument("There should be exactly 9 arguments\n");
  }

  bool verbose = false;
  if (!strcmp(argv[9], "y"))
    verbose = true;

  //////////PARAMETER SET UP////////////////
  // Plaintext prime modulus
  unsigned long p = atol(argv[1]);
  // Field extension degree
  unsigned long d = atol(argv[2]);
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = atol(argv[3]);
  // Number of ciphertext prime bits in the modulus chain
  unsigned long nb_primes = atol(argv[4]);
  // Number of columns of Key-Switching matrix (default = 2 or 3)
  unsigned long c = 3;
  cout << "Initialising context object..." << endl;
  // Intialise context
  Context context(m, p, 1);
  context.scale = 6;
  // Modify the context, adding primes to the modulus chain
  cout  << "Building modulus chain..." << endl;
  buildModChain(context, nb_primes, c);

  // Print the security level
  cout << "Q size: " << context.logOfProduct(context.ctxtPrimes)/log(2.0) << endl;
  cout << "Q*P size: " << context.logOfProduct(context.fullPrimes())/log(2.0) << endl;
  cout << "Security: " << context.securityLevel() << endl;

  // Print the context
  context.zMStar.printout();
  cout << endl;

  //maximal number of digits in a number
  unsigned long expansion_len = atol(argv[5]);

  // Secret key management
  cout << "Creating secret key..." << endl;
  // Create a secret key associated with the context
  SecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();
  cout << "Generating key-switching matrices..." << endl;
  // Compute key-switching matrices that we need
  if (expansion_len > 1)
  {
    if (context.zMStar.numOfGens() == 1)
    {
      set<long> automVals;
      long e = 1;
      long ord = context.zMStar.OrderOf(0);
      bool native = context.zMStar.SameOrd(0);
      if(!native)
        automVals.insert(context.zMStar.genToPow(0, -ord));
      while (e < expansion_len){
        long atm = context.zMStar.genToPow(0, ord-e);
        //cout << "Automorphism " << -e << " is " << atm << endl;
        automVals.insert(atm);
        e <<=1;
      }
      addTheseMatrices(secret_key, automVals);
    }
    else
    {
      addSome1DMatrices(secret_key);
    }
  }

  if (d > 1)
    addFrbMatrices(secret_key); //might be useful only when d > 1

  // create Comparator (initialize after buildModChain)
  Comparator comparator(context, d, expansion_len, secret_key, verbose);

  // number of values to be sorted
  int input_len = atoi(argv[6]);

  // levels of consecutive comparisons
  long depth = atol(argv[7]); 

  //repeat experiments 'runs' times
  int runs = atoi(argv[8]);

  //test sorting
  comparator.test_array_min(input_len, depth, runs);

  printAllTimers(cout);

  return 0;
}