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

// the main function that takes 5 arguments (type in Terminal: ./comparison_circuit argv[1] argv[2] argv[3] argv[4] argv[5] argv[6] argv[7])
// argv[1] - the plaintext modulus
// argv[2] - the extension degree of a finite field
// argv[3] - the order of the cyclotomic ring
// argv[4] - the bitsize of the ciphertext modulus in ciphertexts (HElib increases it to fit the moduli chain). The modulus used for public-key generation
// argv[5] - the length of vectors to be compared
// argv[6] - the number of experiment repetitions
// argv[7] - print debug info (y/n)
// argv[8] - exact or randomized circuit (e/r)

// some parameters for quick testing
// 7 75 90 1 10 y
// 7 300 90 1 10 y
// 17 145 120 1 10 y
int main(int argc, char *argv[]) {
  if(argc < 8)
  {
   throw invalid_argument("There should be exactly 7 arguments\n");
  }

  bool verbose = false;
  if (!strcmp(argv[7], "y"))
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
  unsigned long expansion_len = atol(argv[5]);

  // Secret key management
  cout << "Creating secret key..." << endl;
  // Create a secret key associated with the context
  SecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();
  cout << "Generating key-switching matrices..." << endl;
  // Compute key-switching matrices that we need
  addSome1DMatrices(secret_key);
  if (d > 1)
    addFrbMatrices(secret_key); //might be useful only when d > 1

  // create Comparator (initialize after buildModChain)
  Comparator comparator(context, d, expansion_len, secret_key, verbose);

  //repeat experiments several times
  int runs = atoi(argv[6]);
  
  comparator.test(runs);

  printAllTimers(cout);

  return 0;
}
