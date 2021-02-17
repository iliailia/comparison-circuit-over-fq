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
using namespace he_cmp;

// the main function that takes 7 arguments (type in Terminal: ./comparison_circuit argv[1] argv[2] argv[3] argv[4] argv[5] argv[6] argv[7] argv[8])
// argv[1] - circuit type (U, B or T)
// argv[2] - the plaintext modulus
// argv[3] - the dimension of a vector space over a finite field
// argv[4] - the order of the cyclotomic ring
// argv[5] - the bitsize of the ciphertext modulus in ciphertexts (HElib increases it to fit the moduli chain). The modulus used for public-key generation
// argv[6] - the length of vectors to be compared
// argv[7] - the number of experiment repetitions
// argv[8] - print debug info (y/n)

// some parameters for quick testing
// B 7 1 75 90 1 10 y
// B 7 1 300 90 1 10 y
// U 17 1 145 120 1 10 y
int main(int argc, char *argv[]) {
  if(argc < 9)
  {
   throw invalid_argument("There should be exactly 8 arguments\n");
  }


  CircuitType type = UNI;
  if (!strcmp(argv[1], "B"))
  {
    type = BI;
  }
  else if (!strcmp(argv[1], "T"))
  {
    type = TAN;
  }
  else if (!strcmp(argv[1], "U"))
  {
    type = UNI;
  }
  else
  {
    throw invalid_argument("Choose a valid circuit type (U for univariate, B for bivariate and T for Tan et al.\n");
  }

  bool verbose = false;
  if (!strcmp(argv[8], "y"))
    verbose = true;

  //////////PARAMETER SET UP////////////////
  // Plaintext prime modulus
  unsigned long p = atol(argv[2]);
  // Field extension degree
  unsigned long d = atol(argv[3]);
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = atol(argv[4]);
  // Number of ciphertext prime bits in the modulus chain
  unsigned long nb_primes = atol(argv[5]);
  // Number of columns of Key-Switching matix (default = 2 or 3)
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
  unsigned long expansion_len = atol(argv[6]);

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
      std::set<long> automVals;
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
  Comparator comparator(context, type, d, expansion_len, secret_key, verbose);

  //repeat experiments several times
  int runs = atoi(argv[7]);
  
  //test comparison circuit
  comparator.test_compare(runs);

  //test min/max circuit
  //comparator.test_min_max(runs);

  //printAllTimers(cout);

  return 0;
}
