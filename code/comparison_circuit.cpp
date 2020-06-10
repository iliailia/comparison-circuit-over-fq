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

using namespace std;
using namespace NTL;
using namespace helib;

double intlog(unsigned int base, unsigned int input)
{
  return floor(log2(input)/log2(base));
}

void convertToBaset(vector<long>& decomp, unsigned int input, unsigned int base, int nslots)
{
  decomp.clear();
  decomp.resize(nslots,0);
  int power = static_cast<int>(intlog(base, input)) + 1;
  if (power > nslots)
  {
    cout << "Input character is too big to be converted" << endl;
    exit(1);
  }
  unsigned int rest = input;
  unsigned int coeff;

  int i = 0;
  while(i < power)
    {
      coeff = rest % base;
      decomp[i] = static_cast<long>(coeff);
      rest = (rest - coeff) / base;
      i++;
    }
}

void shift_and_add(Ctxt& x, long range_len, const EncryptedArray& ea)
{
  long e = 1;
  while (e < range_len){
    Ctxt tmp = x;
    ea.shift(tmp, -e);
    x += tmp;
    e *=2;
  }
}

void rotate_and_add(Ctxt& x, long range_len, const EncryptedArray& ea)
{
  long nb_slots = ea.size();
  if (range_len == 0)
  {
    vector<long> ptxt1(nb_slots,0);
    ea.encrypt(x,x.getPubKey(),ptxt1);
    return;
  }
  Ctxt y(x.getPubKey());
  vector<long> ptxt0(nb_slots,0);
  ea.encrypt(y,x.getPubKey(),ptxt0);
  int rot_x = 1;
  int rot_y = range_len;
  Ctxt temp(x.getPubKey());
  int tmp_len = range_len;
  while (tmp_len > 1)
    {
      if (tmp_len % 2 == 0)
        {
          temp = x;
          ea.rotate(temp, -rot_x);
          x += temp;
          rot_x *=2;
          tmp_len = tmp_len/2;
        }
      else //if (tmp_len % 2 == 1)
        {
          rot_y -= rot_x;
          temp = x;
          ea.rotate(temp, -rot_y);
          y += temp;
          temp = x;
          ea.rotate(temp, -rot_x);
          x += temp;
          rot_x *= 2;
          tmp_len = (tmp_len - 1)/2;
        }
    }
  x +=y;
}

//replicate function
//the input ciphertext is assumed to have a non-zero value in the first slot and 0 elsewhere
void replicate_n_times(Ctxt& ctxt, long slots, const EncryptedArray& ea)
{
  Ctxt ctxt_orig = ctxt; 

  long k = NTL::NumBits(slots);
  long e = 1;

  // now process bits k-2 down to 0
  for (long j = k-2; j >= 0; j--) {
    // e -> 2*e
    Ctxt tmp = ctxt;
    ea.rotate(tmp, e);
    ctxt += tmp;
    e = 2*e;
    
    long b = NTL::bit(slots, j); // bit j of sz
    // e -> e+b
    if (b) {
      ea.rotate(ctxt, 1); // "don't care"
      ctxt += ctxt_orig;
      e++;
    }
  }
}

//send non-zero elements of a field F_{p^d} to 1 and zero to 0
//if d = 1, this map operates on elements of the basic field F_p
void mapTo01_subfield(const EncryptedArray& ea, Ctxt& ctxt, long d)
{
  long p = ctxt.getPtxtSpace();
  if (p != ea.getPAlgebra().getP()) // ptxt space is p^r for r>1
    throw helib::LogicError("mapTo01 not implemented for r>1");

  if (p>2)
    ctxt.power(p-1); // set y = x^{p-1}

  if (d>1) { // compute the product of the d automorphisms
    std::vector<Ctxt> v(d, ctxt);
    for (long i=1; i<d; i++)
      v[i].frobeniusAutomorph(i);
    totalProduct(ctxt, v);
  }
}

void print_decrypted(Ctxt& ctxt, long ord_p, const EncryptedArray& ea, const SecKey& sk)
{
    long nSlots = ea.size();
    vector<ZZX> decrypted(nSlots);
    ea.decrypt(ctxt, sk, decrypted);

    for(int i = 0; i < nSlots; i++)
    {
      printZZX(cout, decrypted[i], ord_p);
      cout << endl;
    }
}

void less_simple(Ctxt& ctxt_x, Ctxt& ctxt_y, long p, long ord_p, const EncryptedArray& ea, const SecKey& sk, bool verbose)
{
  FHE_NTIMER_START(ComparisonCircuit);
  long nSlots = ea.size(); 
  //check that p <= number of slots
  assertInRange(p, 0l, nSlots, "replication failed (the plaintext modulus must be in [0, nSlots))");

  //replicate a value in the position pos, the rest must be zero
  cout << "Replication" << endl;
  replicate_n_times(ctxt_x, p-1, ea);
  replicate_n_times(ctxt_y, p-1, ea);

  if(verbose)
  {
    print_decrypted(ctxt_x, ord_p, ea, sk);
    cout << endl;
    print_decrypted(ctxt_y, ord_p, ea, sk);
  }
  //generate a plaintext with slots {0,...,p-2}
  cout << "Subtraction of field elements" << endl;
  vector<long> vec_field_elements(nSlots,0);
  for (long i = 0; i < p-1; i++)
  {
    vec_field_elements[i] = -i;
  }
  NTL::ZZX ptxt_field_elements;
  ea.encode(ptxt_field_elements, vec_field_elements);
  //subtract this plaintext from x
  ctxt_x.addConstant(ptxt_field_elements);
  //generate a plaintext with slots {1,...,p-1}
  for (long i = 0; i < p-1; i++)
  {
    vec_field_elements[i] = -i-1;
  }
  ea.encode(ptxt_field_elements, vec_field_elements);
  //subtract this plaintext from y
  ctxt_y.addConstant(ptxt_field_elements);

  if(verbose)
  {
    print_decrypted(ctxt_x, ord_p, ea, sk);
    cout << endl;
    print_decrypted(ctxt_y, ord_p, ea, sk);
  }
  //compute mapTo01
  cout << "Mapping to 0 and 1" << endl;
  mapTo01_subfield(ea, ctxt_x, 1);
  mapTo01_subfield(ea, ctxt_y, 1);

  if(verbose)
  {
    print_decrypted(ctxt_x, ord_p, ea, sk);
    cout << endl;
    print_decrypted(ctxt_y, ord_p, ea, sk);
  }
  //generate a plaintext with 1 in all the slots
  cout << "Computing NOT" << endl;
  vector<long> vec_ones(nSlots,0);
  for (int i = 0; i< p-1; i++)
    vec_ones[i] = 1;
  NTL::ZZX ptxt_ones;
  ea.encode(ptxt_ones, vec_ones);
  //compute 1 - mapTo01(x-i) and 1 - mapTo01(y-i)
  ctxt_x.negate();
  ctxt_x.addConstant(ptxt_ones,1);
  ctxt_y.negate();
  ctxt_y.addConstant(ptxt_ones,1);

  if(verbose)
  {
    print_decrypted(ctxt_x, ord_p, ea, sk);
    cout << endl;
    print_decrypted(ctxt_y, ord_p, ea, sk);
  }
  //compute shift_and_add that results in a ciphertext containing the partial sums of (1 - mapTo01(y-i)) in different slots
  //DEBUG HERE
  cout << "Rotating and adding y slots" << endl;
  shift_and_add(ctxt_y, p-1, ea);

  if(verbose)
    print_decrypted(ctxt_y, ord_p, ea, sk);

  //multiply the result by 1 - mapTo01(x-i)
  cout << "Multiplying" << endl;
  ctxt_x.multiplyBy(ctxt_y);

  if(verbose)
    print_decrypted(ctxt_x, ord_p, ea, sk);
  //rotate_and_add
  cout << "Final summation" << endl;
  shift_and_add(ctxt_x, p-1, ea);

  if(verbose)
    print_decrypted(ctxt_x, ord_p, ea, sk);
  FHE_NTIMER_STOP(ComparisonCircuit);
}

// the main function that takes 5 arguments (type in Terminal: ./comparison_circuit argv[1] argv[2] argv[3] argv[4] argv[5])
// argv[1] - the plaintext modulus
// argv[2] - the order of the cyclotomic ring
// argv[3] - the bitsize of the ciphertext modulus in ciphertexts (HElib increases it to fit the moduli chain). The modulus used for public-key generation
// argv[4] - the length of vectors to be compared
// argv[5] - the number of experiment repetitions
// argv[6] - print debug info or not (y/n)

// some parameters for quick testing
// 7 75 90 1 10 y
// 17 145 120 1 10 y
int main(int argc, char *argv[]) {
  if(argc < 7)
  {
   throw invalid_argument("There should be exactly 6 arguments\n");
  }

  bool verbose = false;
  if (!strcmp(argv[6], "y"))
    verbose = true;

  // initialize the random generator
  random_device rd;
  mt19937 eng(rd());
  uniform_int_distribution<unsigned int> distr_u;
  uniform_int_distribution<int> distr_i;

  //////////PARAMETER SET UP////////////////
  // Plaintext prime modulus
  unsigned int p = atol(argv[1]);
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = atol(argv[2]);
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of ciphertext prime bits in the modulus chain = depth of the computation
  unsigned long nb_primes = atol(argv[3]);
  // Number of columns of Key-Switching matix (default = 2 or 3)
  unsigned long c = 3;
  std::cout << "Initialising context object..." << std::endl;
  // Intialise context
  Context context(m, p, r);
  context.scale = 6;
  // Modify the context, adding primes to the modulus chain
  cout  << "Building modulus chain..." << endl;
  buildModChain(context, nb_primes, c);

  // Print the context
  context.zMStar.printout();
  cout << endl;

  //determine the order of p in (Z/mZ)*
  long ord_p = context.zMStar.getOrdP();

  // Print the security level
  cout << "Q size: " << context.logOfProduct(context.ctxtPrimes)/log(2.0) << endl;
  cout << "Q*P size: " << context.logOfProduct(context.fullPrimes())/log(2.0) << endl;
  cout << "Security: " << context.securityLevel() << endl;

  // Secret key management
  cout << "Creating secret key..." << endl;
  // Create a secret key associated with the context
  SecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();
  cout << "Generating key-switching matrices..." << endl;
  // Compute key-switching matrices that we need
  add1DMatrices(secret_key);
  addFrbMatrices(secret_key);

  // Public key management
  // Set the secret key (upcast: SecKey is a subclass of PubKey)
  const PubKey& public_key = secret_key;

  // Get the EncryptedArray of the context
  const EncryptedArray& ea = *(context.ea);

  // Get the number of slot (phi(m))
  long nslots = ea.size();
  cout << "Number of slots: " << nslots << endl;
  cout << "Extension degree of a slot:  " << ord_p << endl;

  //vector length
  long vec_len = atol(argv[4]);
  //pattern copies in one ciphertext
  //long pat_copies = nslots/pat_len;

  //timers
  setTimersOn();

  //repeat experiments several times
  int runs = atoi(argv[5]);
  long min_capacity = 1000;
  long capacity;
  for (int run = 0; run < runs; run++)
  {
    printf("Run %d started\n", run);

    vector<ZZX> expected_result(nslots);
    vector<ZZX> decrypted(nslots);

    // Create the plaintext polynomials for the text and for the pattern
    vector<ZZX> pol_x(nslots);
    vector<ZZX> pol_y(nslots);
    
    unsigned int input_coef_x;
    unsigned int input_coef_y;
    ZZX pol_slot;

    for (int i = 0; i < vec_len; i++)
    {
      input_coef_x = distr_u(eng) % p;
      input_coef_y = distr_u(eng) % p;

      if(verbose)
      {
        cout << "Input" << endl;
        cout << input_coef_x << endl;
        cout << input_coef_y << endl;
      }

      if (input_coef_x < input_coef_y)
      {
        expected_result[i] = ZZX(INIT_MONO, 0, 1);
      }
      else
      {
        expected_result[i] = ZZX(INIT_MONO, 0, 0);
      }

      vector<long> decomp_char;
      convertToBaset(decomp_char, input_coef_x, p, ord_p);
      if(verbose)
      {
        cout << "Input slots" << endl;
        cout << decomp_char << endl;
      }
      for (int j = 0; j < ord_p; j++)
      {
        SetCoeff(pol_slot, j, decomp_char[j]);
      }
      pol_x[i] = pol_slot;
      convertToBaset(decomp_char, input_coef_y, p, ord_p);
      if(verbose)
        cout << decomp_char << endl;
      for (int j = 0; j < ord_p; j++)
      {
        SetCoeff(pol_slot, j, decomp_char[j]);
      }
      pol_y[i] = pol_slot;
    }

    //printZZX(cout, pol_x[0], 1);
    //printZZX(cout, pol_y[0], 1);
    cout << endl;

    Ctxt ctxt_x(public_key);
    Ctxt ctxt_y(public_key);
    ea.encrypt(ctxt_x, public_key, pol_x);
    ea.encrypt(ctxt_y, public_key, pol_y);

    FHE_NTIMER_START(Comparison);
    less_simple(ctxt_x, ctxt_y, p, ord_p, ea, secret_key, verbose);
    printNamedTimer(cout, "ComparisonCircuit");
    FHE_NTIMER_STOP(Comparison);
    printNamedTimer(cout, "Comparison");

    // remove the line below if it gives bizarre results 
    ctxt_x.cleanUp();
    capacity = ctxt_x.bitCapacity();
    cout << "Final capacity: " << capacity << endl;
    if (capacity < min_capacity)
      min_capacity = capacity;
    cout << "Min. capacity: " << min_capacity << endl;
    cout << "Final size: " << ctxt_x.logOfPrimeSet()/log(2.0) << endl;
    ea.decrypt(ctxt_x, secret_key, decrypted);

    /*
    for(int i = 0; i < p-1; i++)
    {
      printZZX(cout, decrypted[i], ord_p);
      cout << endl;
    }
    */

    for(int i = 0; i < vec_len; i++)
    { 
      if (decrypted[i] != expected_result[i])
      {
        printf("Slot %d: ", i);
        printZZX(cout, decrypted[i], ord_p);
        cout << endl;
        cout << "Failure" << endl;
        return 1;
      }
    }
    cout << endl;
  }

  return 0;
}
