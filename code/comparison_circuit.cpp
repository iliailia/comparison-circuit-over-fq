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
//#include <NTL/ZZ_pX.h>

using namespace std;
using namespace NTL;
using namespace helib;

double intlog(unsigned int base, unsigned int input)
{
  return floor(log2(input)/log2(base));
}

void digit_decomp(vector<long>& decomp, unsigned int input, unsigned int base, int nslots)
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

ZZX getG(const EncryptedArray& ea)
{
  NTL::ZZX G;
  switch (ea.getTag()) {
  case PA_GF2_tag:
    G = NTL::conv<NTL::ZZX>(ea.getDerived(PA_GF2()).getG());
    break;
  case PA_zz_p_tag:
    convert(G, ea.getDerived(PA_zz_p()).getG());
    break;
  case PA_cx_tag:
    throw helib::LogicError("Cannot get polynomial modulus G when scheme is CKKS");
    break;
  default:
    throw helib::LogicError("No valid tag found in EncryptedArray");
  }
  return G;
}

void int_to_slot(ZZX& poly, unsigned long input, unsigned long p, unsigned long d, unsigned long ord_p, const EncryptedArray& ea)
{
    if (ord_p%d != 0)
    {
      cout << "Extension degree must divide the order" << endl;
      exit(1);
    }
    unsigned long p_d = power_long(p, d);
    unsigned long field_index = ord_p / d;
    unsigned long gen_exp = 1;
    for (unsigned long i = 1; i < field_index; i++)
    {
      gen_exp += power_long(p_d, i);
    }
    //cout << "Exponent of X: " << gen_exp << endl;
    //ZZX G = getG(ea);
    //printZZX(cout, G, ord_p+1);
    //cout << endl;

    vector<long> decomp;

    //decomposition of a digit
    digit_decomp(decomp, input, p, d);
    PolyMod poly_mod(ea.getContext().slotRing);
    poly_mod = ZZX(INIT_MONO, 0, 0);
    for (int k = 0; k < d; k++)
    {
        //TODO: wrong assumption that X is the generator of the finite field
        poly_mod+=ZZX(INIT_MONO, k * gen_exp, decomp[k]);
        //SetCoeff(poly, k * field_index, decomp[k]);
    }
    poly = poly_mod.getData();
    //ZZX def_poly = getG(ea);
    //NTL:rem(poly, poly, def_poly);
}

void batch_shift(Ctxt& ctxt, long start, long shift, long range_len, const EncryptedArray& ea)
{
  if(shift == 0)
    return;

  long nSlots = ea.size();
  long batch_size = nSlots / range_len;

  // left cyclic rotation
  ea.rotate(ctxt, shift);

  // create a mask
  vector<long> mask_vec(nSlots,1);
  // set zeros in the unused slots
  long nEndZeros = nSlots - batch_size * range_len;
  for (int i = 1; i <= nEndZeros; i++)
  {
    long indx = (start + nSlots - i)%nSlots;
    mask_vec[indx] = 0;
  }

  // masking values rotated outside their batches
  for (long i = 0; i < batch_size; i++)
  {
    if (shift < 0)
    {
      for (long j = 0;  j < -shift; j++)
      {
        long indx = (start + (i + 1) * range_len - j - 1)%nSlots;
        mask_vec[indx] = 0;
      }
    }
    else if (shift > 0)
    {
      for (long j = 0;  j < shift; j++)
      {
        long indx = (start + i * range_len + j)%nSlots;
        mask_vec[indx] = 0;
      }
    }
  }
  ZZX mask_ptxt;
  ea.encode(mask_ptxt, mask_vec);
  ctxt.multByConstant(mask_ptxt);
}

void batch_shift(Ctxt& ctxt, const vector<long>& batch_borders, long shift, long range_len, const EncryptedArray& ea)
{
  if(shift == 0)
    return;

  long nSlots = ea.size();
  size_t batch_size = batch_borders.size();

  // left cyclic rotation
  ea.rotate(ctxt, shift);

  // create a mask
  vector<long> mask_vec(nSlots,1);
  // masking values rotated outside their batches
  for (size_t i = 0; i < batch_size; i++)
  {
    if (shift < 0)
    {
      for (long j = 0;  j < -shift; j++)
      {
        long indx = (nSlots + batch_borders[i] - j - 1)%nSlots;
        mask_vec[indx] = 0;
      }
    }
    else if (shift > 0)
    {
      for (long j = 0;  j < shift; j++)
      {
        long indx = (batch_borders[i] + j)%nSlots;
        mask_vec[indx] = 0;
      }
    }
  }
  ZZX mask_ptxt;
  ea.encode(mask_ptxt, mask_vec);
  ctxt.multByConstant(mask_ptxt);
}

// shift_direction is false for the left shift and true for the right shift
void shift_and_add(Ctxt& x, const vector<long>& batch_borders, long range_len, const EncryptedArray& ea, const long shift_direction = false)
{
  long shift_sign = -1;
  if(shift_direction)
    shift_sign = 1;

  long e = 1;

  // shift and add
  while (e < range_len){
    Ctxt tmp = x;
    batch_shift(tmp, batch_borders, e * shift_sign, range_len, ea);
    x += tmp;
    e <<=1;
  }
}

void shift_and_add(Ctxt& x, long start, long range_len, const EncryptedArray& ea, const long shift_direction = false)
{
  long shift_sign = -1;
  if(shift_direction)
    shift_sign = 1;

  long e = 1;

  // shift and add
  while (e < range_len){
    Ctxt tmp = x;
    batch_shift(tmp, start, e * shift_sign, range_len, ea);
    x += tmp;
    e <<=1;
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
//the input ciphertext is assumed to have a non-zero value in the "start" slot and 0 elsewhere
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
void mapTo01_subfield(const EncryptedArray& ea, Ctxt& ctxt, unsigned long d)
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

void basic_less(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y, unsigned long p, unsigned long d, unsigned long ord_p, const EncryptedArray& ea, const SecKey& sk, bool verbose)
{
  FHE_NTIMER_START(ComparisonCircuit);

  //total number of slots
  unsigned long nSlots = ea.size();
  //size of a finite field keeping digits
  unsigned long field_size = power_long(p,d);
  //check that p <= number of slots
  assertInRange(field_size, 0UL, nSlots, "replication failed (the plaintext modulus must be in [0, nSlots))");

  //maximal number of slots to be compared per ciphertext
  long batch_cpcty = nSlots / (field_size - 1);
  //number of batches (ciphertexts) to compare all slots
  long batch_num = static_cast<long>(ceil(static_cast<float>(nSlots) / static_cast<float>(batch_cpcty)));
  //number of elements in the last batch
  long tail = nSlots%batch_cpcty;
  if(tail == 0)
    tail += batch_cpcty; 

  for (long batch_indx = 0; batch_indx < batch_num; batch_indx++)
  {
    //define batch split
    vector<long> batch_split;
    //fill the last subarray
    if (batch_indx == (batch_num - 1))
    {
      for(long indx = 0; indx < tail; indx++)
      {
        batch_split.push_back(batch_indx + indx * batch_num);
      }
    }
    else //fill other subarrays
    {
      for (long indx = 0; indx < min(tail+1, batch_cpcty); indx++)
      {
        batch_split.push_back(batch_indx + indx * batch_num);
      }
      for (long indx  = tail+1; indx < batch_cpcty; indx++)
      {
        batch_split.push_back(batch_split[indx-1] + batch_num - 1);
      }
    }

    if(verbose)
    {
      cout << "Batch split:" << endl;
      for (size_t i = 0; i < batch_split.size(); i++)
      {
        cout << batch_split[i] << endl;
      }
    }

    //extract a batch of digits
    //create a mask
    vector<long> mask_vec(nSlots,0);

    //encode 1's in the positions defined by batch_split
    for (size_t split_indx = 0; split_indx < batch_split.size(); split_indx++)
    {
        mask_vec[batch_split[split_indx]] = 1;
    }
    ZZX mask_ptxt;
    ea.encode(mask_ptxt, mask_vec);

    //multiply both inputs by a mask
    Ctxt ctxt_tmp_x = ctxt_x;
    ctxt_tmp_x.multByConstant(mask_ptxt);
    Ctxt ctxt_tmp_y = ctxt_y;
    ctxt_tmp_y.multByConstant(mask_ptxt);

    //replicate a value in the position pos, the rest must be zero
    cout << "Replication" << endl;
    replicate_n_times(ctxt_tmp_x, field_size-1, ea);
    replicate_n_times(ctxt_tmp_y, field_size-1, ea);

    if(verbose)
    {
      print_decrypted(ctxt_tmp_x, ord_p, ea, sk);
      cout << endl;
      print_decrypted(ctxt_tmp_y, ord_p, ea, sk);
    }

    //TODO: subtract elements of the field extension 
    cout << "Subtraction of field elements" << endl;
    vector<ZZX> vec_field_elements(nSlots, ZZX(INIT_MONO, 0, 0));
    ZZX slot;
    //generate a plaintext with slots {0,...,p-2}
    for (size_t split_indx = 0; split_indx < batch_split.size(); split_indx++)
    {
      for (long field_elem = 1; field_elem < field_size-1; field_elem++)
      {
        long indx = static_cast<long>((batch_split[split_indx] + field_elem)%nSlots);
        int_to_slot(slot, field_elem, p, d, ord_p, ea);
        NTL::negate(vec_field_elements[indx], slot);
      }
    }
    ZZX ptxt_field_elements;
    ea.encode(ptxt_field_elements, vec_field_elements);
    //subtract this plaintext from x
    ctxt_tmp_x.addConstant(ptxt_field_elements);
    //generate a plaintext with slots {1,...,p-1}
    for (size_t split_indx = 0; split_indx < batch_split.size(); split_indx++)
    {
      for (long field_elem = 0; field_elem < field_size-1; field_elem++)
      {
        long indx = static_cast<long>((batch_split[split_indx] + field_elem)%nSlots);
        int_to_slot(slot, field_elem + 1, p, d, ord_p, ea);
        NTL::negate(vec_field_elements[indx], slot);
      }
    }
    ea.encode(ptxt_field_elements, vec_field_elements);
    //subtract this plaintext from y
    ctxt_tmp_y.addConstant(ptxt_field_elements);

    if(verbose)
    {
      print_decrypted(ctxt_tmp_x, ord_p, ea, sk);
      cout << endl;
      print_decrypted(ctxt_tmp_y, ord_p, ea, sk);
    }

    //compute mapTo01
    cout << "Mapping to 0 and 1" << endl;
    mapTo01_subfield(ea, ctxt_tmp_x, d);
    mapTo01_subfield(ea, ctxt_tmp_y, d);

    if(verbose)
    {
      print_decrypted(ctxt_tmp_x, ord_p, ea, sk);
      cout << endl;
      print_decrypted(ctxt_tmp_y, ord_p, ea, sk);
    }

    //generate a plaintext with 1 in all the slots
    cout << "Computing NOT" << endl;
    vector<long> vec_ones(nSlots,0);
    for (size_t split_indx = 0; split_indx < batch_split.size(); split_indx++)
    {
      for (long slot_indx = 0; slot_indx < field_size-1; slot_indx++)
      {
        long indx = (batch_split[split_indx] + slot_indx)%nSlots;
        vec_ones[indx] = 1;
      }
    }
    NTL::ZZX ptxt_ones;
    ea.encode(ptxt_ones, vec_ones);
    //compute 1 - mapTo01(x-i) and 1 - mapTo01(y-i)
    ctxt_tmp_x.negate();
    ctxt_tmp_x.addConstant(ptxt_ones,1);
    ctxt_tmp_y.negate();
    ctxt_tmp_y.addConstant(ptxt_ones,1);

    if(verbose)
    {
      print_decrypted(ctxt_tmp_x, ord_p, ea, sk);
      cout << endl;
      print_decrypted(ctxt_tmp_y, ord_p, ea, sk);
    }
    
    //compute shift_and_add that results in a ciphertext containing the partial sums of (1 - mapTo01(y-i)) in different slots
    cout << "Rotating and adding y slots" << endl;
    shift_and_add(ctxt_tmp_y, batch_split, field_size-1, ea);

    if(verbose)
      print_decrypted(ctxt_tmp_y, ord_p, ea, sk);

    //multiply the result by 1 - mapTo01(x-i)
    cout << "Multiplying" << endl;
    ctxt_tmp_x.multiplyBy(ctxt_tmp_y);

    if(verbose)
      print_decrypted(ctxt_tmp_x, ord_p, ea, sk);
    //rotate_and_add
    cout << "Final summation" << endl;
    shift_and_add(ctxt_tmp_x, batch_split, field_size-1, ea);
    ctxt_tmp_x.multByConstant(mask_ptxt);
    //extract the results

    if(verbose)
      print_decrypted(ctxt_tmp_x, ord_p, ea, sk);

    ctxt_res += ctxt_tmp_x;
  }

  cout << "Final results of all digit comparisons" << endl;
  if(verbose)
      print_decrypted(ctxt_res, ord_p, ea, sk);
  FHE_NTIMER_STOP(ComparisonCircuit);
}

//randomized equality circuit
void random_equality(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y, long expansion_len, unsigned long ord_p, const EncryptedArray& ea, const SecKey& sk, bool verbose)
{
  FHE_NTIMER_START(RandomEqualityCircuit);
  long nSlots = ea.size();
  long numbers_size = nSlots / expansion_len; 

  Ctxt ctxt_tmp_x = ctxt_x;
  Ctxt ctxt_tmp_y = ctxt_y;

  //Remove the least significant digit and shift to the left
  //TODO: test the lines below
  cout << "Remove the least significant digit" << endl;
  batch_shift(ctxt_tmp_x, 0, -1, expansion_len, ea);
  batch_shift(ctxt_tmp_y, 0, -1, expansion_len, ea);

  if(verbose)
  {
    print_decrypted(ctxt_tmp_x, ord_p, ea, sk);
    cout << endl;
    print_decrypted(ctxt_tmp_y, ord_p, ea, sk);
    cout << endl;
  }

  // Subtraction (x_i - y_i)
  cout << "Subtraction" << endl;
  ctxt_res = ctxt_tmp_x;
  ctxt_res -= ctxt_tmp_y;

  if(verbose)
  {
    print_decrypted(ctxt_res, ord_p, ea, sk);
    cout << endl;
  }

  //compute the multiplication with the random polynomial: r_i*(x_i - y_i)
  cout << "Multiplication by a random element" << endl;
  Ptxt<BGV> poly_r(ea.getContext());
  poly_r.random();
  ctxt_res.multByConstant(poly_r);

  if(verbose)
  {
    print_decrypted(ctxt_res, ord_p, ea, sk);
    cout << endl;
  }

  //compute running sums: sum_i r_i*(x_i - y_i)
  cout << "Rotating and adding slots" << endl;
  shift_and_add(ctxt_res, 0, expansion_len, ea);

  if(verbose)
  {
    print_decrypted(ctxt_res, ord_p, ea, sk);
    cout << endl;
  }

  //compute mapTo01
  cout << "Mapping to 0 and 1" << endl;
  mapTo01_subfield(ea, ctxt_res, ord_p);

  if(verbose)
  {
    print_decrypted(ctxt_res, ord_p, ea, sk);
    cout << endl;
  }
  //generate a plaintext with 1 in all the slots
  cout << "Computing NOT" << endl;
  vector<long> vec_ones(nSlots,0);
  for (long i = 0; i < numbers_size; i++)
  {
    for (long j = 0; j < expansion_len; j++)
    {
      vec_ones[i * (expansion_len) + j] = 1;
    }
  }
  ZZX ptxt_ones;
  ea.encode(ptxt_ones, vec_ones);

  //compute 1 - mapTo01(r_i*(x_i - y_i))
  ctxt_res.negate();
  ctxt_res.addConstant(ptxt_ones,1);

  if(verbose)
  {
    print_decrypted(ctxt_res, ord_p, ea, sk);
    cout << endl;
  }

  FHE_NTIMER_STOP(RandomEqualityCircuit);
}

bool test_equality(const PubKey& public_key, const EncryptedArray& ea, const SecKey& sk, long p, long ord_p, long field_size, long expansion_len, long nslots, long numbers_size, bool verbose)
{
    long min_capacity = 1000;
    long capacity;
    // initialize the random generator
    random_device rd;
    mt19937 eng(rd());
    uniform_int_distribution<unsigned int> distr_u;
    uniform_int_distribution<int> distr_i;

    //generate input
    long input_range = static_cast<long>(pow(field_size, expansion_len));
    vector<ZZX> expected_result(nslots);
    vector<ZZX> decrypted(nslots);

    // Create the plaintext polynomials for the text and for the pattern
    vector<ZZX> pol_x(nslots);
    vector<ZZX> pol_y(nslots);
    
    unsigned int input_x;
    unsigned int input_y;
    ZZX pol_slot;

    for (int i = 0; i < numbers_size; i++)
    {
      input_x = distr_u(eng) % input_range;
      input_y = distr_u(eng) % input_range;

      if(verbose)
      {
        cout << "Input integers:" << endl;
        cout << input_x << endl;
        cout << input_y << endl;
      }

      vector<long> decomp_int_x;
      vector<long> decomp_int_y;
      vector<long> decomp_char;

      //decomposition of input integers
      digit_decomp(decomp_int_x, input_x, field_size, expansion_len);
      digit_decomp(decomp_int_y, input_y, field_size, expansion_len);

      //expected output
      expected_result[(i + 1) * expansion_len - 1] = ZZX(INIT_MONO, 0, 1);

      for (int j = expansion_len - 2; j >= 0; j--)
      {
        if (decomp_int_x[j+1] == decomp_int_y[j+1])
        {
          expected_result[i * expansion_len + j] = ZZX(INIT_MONO, 0, 1) * expected_result[i * expansion_len + j + 1];
        }
        else
        {
          expected_result[i * expansion_len + j] = ZZX(INIT_MONO, 0, 0);
        }
      }

      //encoding of slots
      for (int j = 0; j < expansion_len; j++)
      {
          //decomposition of a digit
          digit_decomp(decomp_char, decomp_int_x[j], p, ord_p);
          for (int k = 0; k < ord_p; k++)
          {
            SetCoeff(pol_slot, k, decomp_char[k]);
          }
          pol_x[i * expansion_len + j] = pol_slot;
      }

      for (int j = 0; j < expansion_len; j++)
      {
          //decomposition of a digit
          digit_decomp(decomp_char, decomp_int_y[j], p, ord_p);
          for (int k = 0; k < ord_p; k++)
          {
            SetCoeff(pol_slot, k, decomp_char[k]);
          }
          pol_y[i * expansion_len + j] = pol_slot;
      }
    }

    if(verbose)
    {
      cout << "Input" << endl;
      for(int i = 0; i < nslots; i++)
      {
          printZZX(cout, pol_x[i], 1);
          printZZX(cout, pol_y[i], 1);
          cout << endl;
      }
    }

    Ctxt ctxt_x(public_key);
    Ctxt ctxt_y(public_key);
    ea.encrypt(ctxt_x, public_key, pol_x);
    ea.encrypt(ctxt_y, public_key, pol_y);
    Ctxt ctxt_res(public_key);

    random_equality(ctxt_res, ctxt_x, ctxt_y, expansion_len, ord_p, ea, sk, verbose);
    printNamedTimer(cout, "RandomEqualityCircuit");

    // remove the line below if it gives bizarre results 
    ctxt_res.cleanUp();
    capacity = ctxt_res.bitCapacity();
    cout << "Final capacity: " << capacity << endl;
    if (capacity < min_capacity)
      min_capacity = capacity;
    cout << "Min. capacity: " << min_capacity << endl;
    cout << "Final size: " << ctxt_res.logOfPrimeSet()/log(2.0) << endl;
    ea.decrypt(ctxt_res, sk, decrypted);

    for(int i = 0; i < numbers_size; i++)
    { 
      for (int j = 0; j < expansion_len; j++)
      {
        long indx = i * expansion_len + j;
        if (decrypted[indx] != expected_result[indx])
        {
          printf("Slot %ld: ", indx);
          printZZX(cout, decrypted[indx], ord_p);
          cout << endl;
          cout << "Failure" << endl;
          return false;
        }
      }
    }
    cout << endl;
    cout << "Random equality test passed" << endl;
    return true;
}
// the main function that takes 5 arguments (type in Terminal: ./comparison_circuit argv[1] argv[2] argv[3] argv[4] argv[5] argv[6] argv[7])
// argv[1] - the plaintext modulus
// argv[2] - the extension degree of a finite field
// argv[3] - the order of the cyclotomic ring
// argv[4] - the bitsize of the ciphertext modulus in ciphertexts (HElib increases it to fit the moduli chain). The modulus used for public-key generation
// argv[5] - the length of vectors to be compared
// argv[6] - the number of experiment repetitions
// argv[7] - print debug info (y/n)

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

  // initialize the random generator
  random_device rd;
  mt19937 eng(rd());
  uniform_int_distribution<unsigned long> distr_u;
  uniform_int_distribution<long> distr_i;

  //////////PARAMETER SET UP////////////////
  // Plaintext prime modulus
  unsigned long p = atol(argv[1]);
  // Field extension degree
  unsigned long d = atol(argv[2]);
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = atol(argv[3]);
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of ciphertext prime bits in the modulus chain = depth of the computation
  unsigned long nb_primes = atol(argv[4]);
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
  unsigned long ord_p = context.zMStar.getOrdP();
  if (ord_p%d != 0)
  {
    cout << "Field extension must divide the order of the plaintext modulus" << endl;
    return 0;
  }
  //index of the finite field of digits in the finite field of a slot
  unsigned long field_index = ord_p / d;
  //size of a subfield
  unsigned long subfield_size = power_long(p,d);

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

  //find the generator of a slot
  /*ZZX gen;

  for (int i = 0; i < p; i++)
  {
    ZZ_pX genP = ZZ_pX(INIT_MONO, 0, i);
    for (int j = 0; j < p; j++)
    {
      genP += ZZ_pX(INIT_MONO, 1, j);
      ZZ_pX G = conv<ZZ_pX>(getG(ea));
      PowerMod(genP, genP, p, G);
    }
  }*/

  // Get the number of slot (phi(m))
  unsigned long nslots = ea.size();
  cout << "Number of slots: " << nslots << endl;
  cout << "Extension degree of a slot:  " << ord_p << endl;

  //size of the finite field used to keep digits
  unsigned long field_size = power_long(p,d);
  //maximal number of digits in a number
  unsigned long expansion_len = atol(argv[5]);
  //amount of numbers in one ciphertext
  unsigned long numbers_size = nslots / expansion_len;
  unsigned long occupied_slots = numbers_size * expansion_len;

  //check that field_size^expansion_len fits into 64-bits
  int space_bit_size = static_cast<int>(ceil(expansion_len * log2(field_size)));
  unsigned long input_range = LONG_MAX;
  if(space_bit_size < 64)
  {
    input_range = power_long(field_size, expansion_len);
  }
  cout << "Maximal input: " << input_range << endl; 

  //timers
  setTimersOn();

  //repeat experiments several times
  int runs = atoi(argv[6]);
  long min_capacity = 1000;
  long capacity;
  for (int run = 0; run < runs; run++)
  {
    printf("Run %d started\n", run);

    /*
    if (!test_equality(public_key, ea, secret_key, p, ord_p, field_size, expansion_len, nslots, numbers_size, verbose))
      return 1;

    continue;
    */

    vector<ZZX> expected_result(occupied_slots);
    vector<ZZX> decrypted(occupied_slots);

    // Create the plaintext polynomials for the text and for the pattern
    vector<ZZX> pol_x(nslots);
    vector<ZZX> pol_y(nslots);
    
    unsigned long input_x;
    unsigned long input_y;
    ZZX pol_slot;

    for (int i = 0; i < numbers_size; i++)
    {
      input_x = distr_u(eng) % input_range;
      input_y = distr_u(eng) % input_range;

      if(verbose)
      {
        cout << "Input" << endl;
        cout << input_x << endl;
        cout << input_y << endl;
      }

      if (input_x < input_y)
      {
        expected_result[i * expansion_len] = ZZX(INIT_MONO, 0, 1);
      }
      else
      {
        expected_result[i * expansion_len] = ZZX(INIT_MONO, 0, 0);
      }

      vector<long> decomp_int_x;
      vector<long> decomp_int_y;
      vector<long> decomp_char;

      //decomposition of input integers
      digit_decomp(decomp_int_x, input_x, field_size, expansion_len);
      digit_decomp(decomp_int_y, input_y, field_size, expansion_len);

      //encoding of slots
      for (int j = 0; j < expansion_len; j++)
      {
          //decomposition of a digit
          int_to_slot(pol_slot, decomp_int_x[j], p, d, ord_p, ea);
          pol_x[i * expansion_len + j] = pol_slot;
      }

      for (int j = 0; j < expansion_len; j++)
      {
          //decomposition of a digit
          int_to_slot(pol_slot, decomp_int_y[j], p, d, ord_p, ea);
          pol_y[i * expansion_len + j] = pol_slot;
      }
    }

    if(verbose)
    {
      cout << "Input" << endl;
      for(int i = 0; i < nslots; i++)
      {
          printZZX(cout, pol_x[i], ord_p);
          printZZX(cout, pol_y[i], ord_p);
          cout << endl;
      }
    }

    Ctxt ctxt_x(public_key);
    Ctxt ctxt_y(public_key);
    ea.encrypt(ctxt_x, public_key, pol_x);
    ea.encrypt(ctxt_y, public_key, pol_y);
    Ctxt ctxt_less(public_key);
    Ctxt ctxt_res(public_key);

    FHE_NTIMER_START(Comparison);
    basic_less(ctxt_less, ctxt_x, ctxt_y, p, d, ord_p, ea, secret_key, verbose);
    random_equality(ctxt_res, ctxt_x, ctxt_y, expansion_len, ord_p, ea, secret_key, verbose);

    ctxt_res.multiplyBy(ctxt_less);
    shift_and_add(ctxt_res, 0, expansion_len, ea);

    if(verbose)
    {
      cout << "Input" << endl;
      for(int i = 0; i < nslots; i++)
      {
          printZZX(cout, pol_x[i], ord_p);
          printZZX(cout, pol_y[i], ord_p);
          cout << endl;
      }

      cout << "Output" << endl;
      print_decrypted(ctxt_res, ord_p, ea, secret_key);
      cout << endl;
    }
    FHE_NTIMER_STOP(Comparison);
    printNamedTimer(cout, "ComparisonCircuit");
    printNamedTimer(cout, "Comparison");

    // remove the line below if it gives bizarre results 
    ctxt_res.cleanUp();
    capacity = ctxt_res.bitCapacity();
    cout << "Final capacity: " << capacity << endl;
    if (capacity < min_capacity)
      min_capacity = capacity;
    cout << "Min. capacity: " << min_capacity << endl;
    cout << "Final size: " << ctxt_res.logOfPrimeSet()/log(2.0) << endl;
    ea.decrypt(ctxt_res, secret_key, decrypted);

    for(int i = 0; i < numbers_size; i++)
    { 
      if (decrypted[i * expansion_len] != expected_result[i * expansion_len])
      {
        printf("Slot %ld: ", i * expansion_len);
        printZZX(cout, decrypted[i * expansion_len], ord_p);
        cout << endl;
        cout << "Failure" << endl;
        return 1;
      }
    }
    cout << endl;
  }

  return 0;
}
