//TODO: finish this function
ZZX get_subfield_gen(const EncryptedArray& ea)
{
  long p = ea.getPAlgebra().getP();
  long ord_p = ea.getPAlgebra().getOrdP();

  //find the generator of a slot
  ZZX gen;
  ZZ_p modP;
  modP.init(ZZ(p));
  modP = 0;
  ZZ_pX genP(modP);

  //polynomial modulus of a slot
  ZZ_pXModulus G(conv<ZZ_pX>(getG(ea)));
  cout << "G: " << G << endl;

  //size of the slot finite field
  ZZ slot_group_size = power_ZZ(p, ord_p) - ZZ(1);
  cout << "Slot field size: " << slot_group_size << endl;

  //factorize the order - 1
  vector<ZZ> factors;
  factorize(factors, slot_group_size);

  //factor orders
  vector<ZZ> factor_orders(factors.size(), ZZ(0));
  for (int i = 0; i < factors.size(); i++)
  {
    ZZ tmp = slot_group_size;
    ZZ factor_order(0);
    while(true)
    {
      if(tmp%factors[i]!=0)
        break;
      tmp/=factors[i];
      factor_order++;
    }
    factor_orders[i] = factor_order;
  }
  for (int i = 0; i < factor_orders.size(); i++)
    cout << factors[i] << ' ' << factor_orders[i] << endl;
  
  //compute the number of divisors
  ZZ orders_prod(1);
  for (int i = 0; i < factor_orders.size(); i++)
    orders_prod*=(factor_orders[i]+1);
  cout << "Number of divisors: " << orders_prod << endl;

  //TODO: exponentiate random polynomials to divisors of slot_field_size
  //generate a random polynomial of a field
  genP = random_ZZ_pX(ord_p);
  cout << genP << endl;

  genP = PowerMod(genP, 2, G);
  cout << genP << endl;

  return gen;
}
*/

/*
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

void batch_shift_for_mul(Ctxt& ctxt, const vector<long>& batch_borders, long shift, long range_len, const EncryptedArray& ea)
{
  FHE_NTIMER_START(BatchShiftForMul);
  if(shift == 0)
    return;

  long nSlots = ea.size();
  size_t batch_size = batch_borders.size();

  // left cyclic rotation
  ea.rotate(ctxt, shift);

  // create a mask
  vector<long> mask_vec(nSlots,1);
  ZZX mask_ptxt;
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
  ea.encode(mask_ptxt, mask_vec);
  ctxt.multByConstant(mask_ptxt);

  mask_ptxt = 1 - mask_ptxt;

  ctxt.addConstant(mask_ptxt);
  FHE_NTIMER_STOP(BatchShiftForMul);
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

// shift_direction is false for the left shift and true for the right shift
void shift_and_mul(Ctxt& x, const vector<long>& batch_borders, long range_len, const EncryptedArray& ea, const long shift_direction = false)
{
  long shift_sign = -1;
  if(shift_direction)
    shift_sign = 1;

  long e = 1;

  // shift and add
  while (e < range_len){
    Ctxt tmp = x;
    batch_shift_for_mul(tmp, batch_borders, e * shift_sign, range_len, ea);
    x.multiplyBy(tmp);
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
*/

/*
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
    ZZX ptxt_ones;
    ea.encode(ptxt_ones, vec_ones);
    //compute 1 - mapTo01(x-i) and 1 - mapTo01(y-i)
    ctxt_tmp_x.negate();
    ctxt_tmp_x.addConstant(ptxt_ones);
    ctxt_tmp_y.negate();
    ctxt_tmp_y.addConstant(ptxt_ones);

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
*/

/*
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