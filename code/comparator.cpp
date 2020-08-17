#include "comparator.h"
#include "tools.h"
#include <helib/debugging.h>
#include <helib/polyEval.h>
#include <random>
#include <map> 

void Comparator::create_shift_mask(ZZX& mask_ptxt, long shift)
{
	cout << "Mask for shift " << shift << " is being created" << endl;
	// get EncryptedArray
  	const EncryptedArray& ea = *(m_context.ea);

  	//extract slots
	long nSlots = ea.size();

	//number of batches in one slot
	long batch_size = nSlots / m_expansionLen;

	// create a mask vector
	vector<long> mask_vec(nSlots,1);

	//starting position of all batches
	long start = 0;

	// set zeros in the unused slots
	long nEndZeros = nSlots - batch_size * m_expansionLen;
	for (int i = 1; i <= nEndZeros; i++)
	{
	long indx = (start + nSlots - i) % nSlots;
	mask_vec[indx] = 0;
	}

	// masking values rotated outside their batches
	for (long i = 0; i < batch_size; i++)
	{
	if (shift < 0)
	{
	  for (long j = 0;  j < -shift; j++)
	  {
	    long indx = (start + (i + 1) * m_expansionLen - j - 1) % nSlots;
	    mask_vec[indx] = 0;
	  }
	}
	else if (shift > 0)
	{
	  for (long j = 0;  j < shift; j++)
	  {
	    long indx = (start + i * m_expansionLen + j) % nSlots;
	    mask_vec[indx] = 0;
	  }
	}
	}
	ea.encode(mask_ptxt, mask_vec);
}

void Comparator::create_all_shift_masks()
{
	long shift = 1;
	while (shift < m_expansionLen)
	{
	    ZZX mask_ptxt;
	    create_shift_mask(mask_ptxt, -shift);
	    m_mulMasks.push_back(mask_ptxt); 
	    shift <<=1;
	}
}

void Comparator::create_poly()
{
	// get p
	unsigned long p = m_context.zMStar.getP();;

	// polynomial coefficient
	ZZ_p coef;
	coef.init(ZZ(p));

	// field element
	ZZ_p field_elem;
	field_elem.init(ZZ(p));

	// initialization of the comparison polynomial with x^{p-1} * (p+1)/2
	m_poly = ZZX(INIT_MONO, p-1, ZZ((p+1) >> 1));

	// loop over all odd coefficient indices
	for (long indx = 1; indx < p - 1; indx+=2)
	{ 
		// coefficient f_i = sum_a a^{p-1-indx} where a runs over [1,...,(p-1)/2]
		coef = 1;
		for(long a = 2; a <= ((p-1) >> 1); a++)
		{
		  field_elem = a;
		  coef += power(field_elem, p - 1 - indx);
		}

		m_poly += ZZX(INIT_MONO, indx, rep(coef));
	}

	if (m_verbose)
	{
		cout << "Comparison polynomial: " << endl;
		printZZX(cout, m_poly, p-1);
		cout << endl;
	}
}

Comparator::Comparator(const Context& context, unsigned long d, unsigned long expansion_len, const SecKey& sk, bool verbose): m_context(context), m_slotDeg(d), m_expansionLen(expansion_len), m_sk(sk), m_pk(sk), m_verbose(verbose)
{
	//determine the order of p in (Z/mZ)*
	unsigned long ord_p = context.zMStar.getOrdP();
	//check that the extension degree divides the order of p
	if (ord_p % d != 0)
	{
		throw invalid_argument("Field extension must divide the order of the plaintext modulus\n");
	}

	create_all_shift_masks();
	create_poly();

}

const ZZX& Comparator::get_mask(long index) const
{
	return m_mulMasks[index];
}

const ZZX& Comparator::get_poly() const
{
	return m_poly;
}

void Comparator::print_decrypted(const Ctxt& ctxt) const
{
	// get EncryptedArray
	const EncryptedArray& ea = *(m_context.ea);

	// get order of p
	unsigned long ord_p = m_context.zMStar.getOrdP();

    long nSlots = ea.size();
    vector<ZZX> decrypted(nSlots);
    ea.decrypt(ctxt, m_sk, decrypted);

    for(int i = 0; i < nSlots; i++)
    {
      printZZX(cout, decrypted[i], ord_p);
      cout << endl;
    }
}

void Comparator::batch_shift(Ctxt& ctxt, long start, long shift) const
{
	FHE_NTIMER_START(BatchShift);
	// get EncryptedArray
	const EncryptedArray& ea = *(m_context.ea);
	
	// if shift is zero, do nothing
	if(shift == 0)
		return;

	// left cyclic rotation
	ea.rotate(ctxt, shift);

	// masking elements shifted out of batch
	long index = static_cast<long>(intlog(2, -shift));
	//cout << "Mask index: " << index << endl;
	ZZX mask = get_mask(index);
	ctxt.multByConstant(mask);
	FHE_NTIMER_STOP(BatchShift);
}

void Comparator::batch_shift_for_mul(Ctxt& ctxt, long start, long shift) const
{
	FHE_NTIMER_START(BatchShiftForMul);
	// get EncryptedArray
	const EncryptedArray& ea = *(m_context.ea);
	
	// if shift is zero, do nothing
	if(shift == 0)
		return;

	// left cyclic rotation
	ea.rotate(ctxt, shift);

	long index = static_cast<long>(intlog(2, -shift));
	//cout << "Mask index: " << index << endl;
	ZZX mask = get_mask(index);
	ctxt.multByConstant(mask);

	// add 1 to masked slots
	ctxt.addConstant(1-mask);

	FHE_NTIMER_STOP(BatchShiftForMul);
}

void Comparator::shift_and_add(Ctxt& x, long start, long shift_direction) const
{
  FHE_NTIMER_START(ShiftAdd);
  long shift_sign = -1;
  if(shift_direction)
    shift_sign = 1;

  long e = 1;

  // shift and add
  while (e < m_expansionLen){
    Ctxt tmp = x;
    batch_shift(tmp, start, e * shift_sign);
    x += tmp;
    e <<=1;
  }
  FHE_NTIMER_STOP(ShiftAdd);
}

void Comparator::shift_and_mul(Ctxt& x, long start, long shift_direction) const
{
  FHE_NTIMER_START(ShiftMul);
  long shift_sign = -1;
  if(shift_direction)
    shift_sign = 1;

  long e = 1;

  // shift and add
  while (e < m_expansionLen){
    Ctxt tmp = x;
    batch_shift_for_mul(tmp, start, e * shift_sign);
    x.multiplyBy(tmp);
    e <<=1;
  }
  FHE_NTIMER_STOP(ShiftMul);
}

void Comparator::mapTo01_subfield(Ctxt& ctxt, long d) const
{
  FHE_NTIMER_START(MapTo01);	
  // get EncryptedArray
  const EncryptedArray& ea = *(m_context.ea);

  // get p
  long p = ctxt.getPtxtSpace();
  if (p != ea.getPAlgebra().getP()) // ptxt space is p^r for r>1
    throw helib::LogicError("mapTo01 not implemented for r>1");

  if (p > 2)
    ctxt.power(p - 1); // set y = x^{p-1}

  if (d > 1) { // compute the product of the d automorphisms
    vector<Ctxt> v(d, ctxt);
    for (long i = 1; i < d; i++)
      v[i].frobeniusAutomorph(i);
    totalProduct(ctxt, v);
  }
  FHE_NTIMER_STOP(MapTo01);
}

void Comparator::less_than(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
  FHE_NTIMER_START(ComparisonCircuit);

  unsigned long p = m_context.zMStar.getP();

  // Subtraction z = x - y
  cout << "Subtraction" << endl;
  Ctxt ctxt_z = ctxt_x;
  ctxt_z -= ctxt_y;

  if(m_verbose)
  {
    print_decrypted(ctxt_z);
    cout << endl;
  }

  map<unsigned long, unsigned long> pat_stock_ks{
  	{13, 4}, //4 (15), 2..5
  	{47, 7}, //7 (26), 3..8
  	{61, 6}, //6 (26), 3..6
  	{67, 6}, //6 (26), 4..8
  	{131, 14}, //10 (33), 10..13
  	{167, 10}, //10 (38), 8..12
  	{173, 12} //11 (39), 7..12
  };

  long k = -1;
  if(pat_stock_ks.count(p) > 0)
  	k = pat_stock_ks[p];

  // compute polynomial function for 'z < 0'
  cout << "Compute comparison polynomial" << endl;
  polyEval(ctxt_res, m_poly, ctxt_z, k);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  FHE_NTIMER_STOP(ComparisonCircuit);
}

void Comparator::random_equality(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
  FHE_NTIMER_START(RandomEqualityCircuit);

  // get order of p
  unsigned long ord_p = m_context.zMStar.getOrdP(); 

  // Subtraction (x_i - y_i)
  cout << "Subtraction" << endl;
  ctxt_res = ctxt_x;
  ctxt_res -= ctxt_y;

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  //compute the multiplication with the random polynomial: r_i*(x_i - y_i)
  cout << "Multiplication by a random element" << endl;
  Ptxt<BGV> poly_r(ctxt_res.getContext());
  poly_r.random();
  ctxt_res.multByConstant(poly_r);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  //compute running sums: sum_i r_i*(x_i - y_i)
  cout << "Rotating and adding slots" << endl;
  shift_and_add(ctxt_res, 0);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  //compute mapTo01
  cout << "Mapping to 0 and 1" << endl;
  mapTo01_subfield(ctxt_res, ord_p);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }
  
  cout << "Computing NOT" << endl;
  //generate a plaintext with 1 in all the data slots
  ZZX ptxt_ones(INIT_MONO, 0, 1);

  //compute 1 - mapTo01(r_i*(x_i - y_i))
  ctxt_res.negate();
  ctxt_res.addConstant(ptxt_ones);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  //Remove the least significant digit and shift to the left
  cout << "Remove the least significant digit" << endl;
  batch_shift_for_mul(ctxt_res, 0, -1);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  FHE_NTIMER_STOP(RandomEqualityCircuit);
}

void Comparator::exact_equality(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
  FHE_NTIMER_START(EqualityCircuit);

  // Subtraction (x_i - y_i)
  cout << "Subtraction" << endl;
  ctxt_res = ctxt_x;
  ctxt_res -= ctxt_y;

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  //compute mapTo01: (x_i - y_i)^{p^d-1}
  cout << "Mapping to 0 and 1" << endl;
  mapTo01_subfield(ctxt_res, m_slotDeg);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  cout << "Computing NOT" << endl;
  //generate a plaintext with 1 in all the data slots
  ZZX ptxt_ones(INIT_MONO, 0, 1);

  //compute 1 - mapTo01(r_i*(x_i - y_i))
  ctxt_res.negate();
  ctxt_res.addConstant(ptxt_ones);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  //compute running products: prod_i 1 - (x_i - y_i)^{p^d-1}
  cout << "Rotating and multiplying slots" << endl;
  shift_and_mul(ctxt_res, 0);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  //Remove the least significant digit and shift to the left
  cout << "Remove the least significant digit" << endl;
  batch_shift_for_mul(ctxt_res, 0, -1);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  FHE_NTIMER_STOP(EqualityCircuit);
}

void Comparator::compare(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y, bool is_randomized) const
{
	FHE_NTIMER_START(Comparison);

    if(m_expansionLen == 1) // if you compare slots one by one, use only less_than function
    {
      less_than(ctxt_res, ctxt_x, ctxt_y);
    }
    else // if you compare vectors of slots lexicographically, use equality function
    {
      Ctxt ctxt_less(ctxt_x.getPubKey());

      // compare all digits by less_than
      less_than(ctxt_less, ctxt_x, ctxt_y);

      // compute equality function on consecutive leading digits
      if(is_randomized)
        random_equality(ctxt_res, ctxt_x, ctxt_y);
      else
        exact_equality(ctxt_res, ctxt_x, ctxt_y);

      // combine above results
      ctxt_res.multiplyBy(ctxt_less);
      shift_and_add(ctxt_res, 0);
    }

    FHE_NTIMER_STOP(Comparison);
}

//TODO: finish this function
void Comparator::int_to_slot(ZZX& poly, unsigned long input) const
{ 
	unsigned long p = m_context.zMStar.getP();
	unsigned long ord_p = m_context.zMStar.getOrdP();

    unsigned long p_d = power_long(p, m_slotDeg);
    unsigned long field_index = ord_p / m_slotDeg;
    unsigned long gen_exp = 1;
    for (unsigned long i = 1; i < field_index; i++)
    {
      gen_exp += power_long(p_d, i);
    }

    vector<long> decomp;

    //decomposition of a digit
    digit_decomp(decomp, input, p, m_slotDeg);
    PolyMod poly_mod(m_context.slotRing);
    poly_mod = ZZX(INIT_MONO, 0, 0);
    //TODO: this loop works correctly only when d = 1
    for (int k = 0; k < m_slotDeg; k++)
    {
        //TODO: wrong assumption that X is the generator of the finite field
        poly_mod+=ZZX(INIT_MONO, k * gen_exp, decomp[k]);
    }
    poly = poly_mod.getData();
}

void Comparator::test(long runs, bool is_randomized) const
{
  //reset timers
  setTimersOn();
  
  // initialize the random generator
  random_device rd;
  mt19937 eng(rd());
  uniform_int_distribution<unsigned long> distr_u;
  uniform_int_distribution<long> distr_i;

  // get EncryptedArray
  const EncryptedArray& ea = *(m_context.ea);

  //extract number of slots
  long nslots = ea.size();

  //get p
  unsigned long p = m_context.zMStar.getP();

  //order of p
  unsigned long ord_p = m_context.zMStar.getOrdP();

  //amount of numbers in one ciphertext
  unsigned long numbers_size = nslots / m_expansionLen;

  // number of slots occupied by encoded numbers
  unsigned long occupied_slots = numbers_size * m_expansionLen;

  //encoding base, (p+1)/2
  //if 2-variable comparison polynomial is used, it must be p
  unsigned long enc_base = (p + 1) >> 1;

  //check that field_size^expansion_len fits into 64-bits
  int space_bit_size = static_cast<int>(ceil(m_expansionLen * log2(enc_base)));
  unsigned long input_range = LONG_MAX;
  if(space_bit_size < 64)
  {
    //input_range = power_long(field_size, expansion_len);
    input_range = power_long(enc_base, m_expansionLen);
  }
  cout << "Maximal input: " << input_range << endl;

  long min_capacity = 1000;
  long capacity;
  for (int run = 0; run < runs; run++)
  {
    printf("Run %d started\n", run);

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

      if(m_verbose)
      {
        cout << "Input" << endl;
        cout << input_x << endl;
        cout << input_y << endl;
      }

      if (input_x < input_y)
      {
        expected_result[i * m_expansionLen] = ZZX(INIT_MONO, 0, 1);
      }
      else
      {
        expected_result[i * m_expansionLen] = ZZX(INIT_MONO, 0, 0);
      }

      vector<long> decomp_int_x;
      vector<long> decomp_int_y;
      vector<long> decomp_char;

      //decomposition of input integers
      digit_decomp(decomp_int_x, input_x, enc_base, m_expansionLen);
      digit_decomp(decomp_int_y, input_y, enc_base, m_expansionLen);

      //encoding of slots
      for (int j = 0; j < m_expansionLen; j++)
      {
          //decomposition of a digit
          int_to_slot(pol_slot, decomp_int_x[j]);
          pol_x[i * m_expansionLen + j] = pol_slot;
      }

      for (int j = 0; j < m_expansionLen; j++)
      {
          //decomposition of a digit
          int_to_slot(pol_slot, decomp_int_y[j]);
          pol_y[i * m_expansionLen + j] = pol_slot;
      }
    }

    if(m_verbose)
    {
      cout << "Input" << endl;
      for(int i = 0; i < nslots; i++)
      {
          printZZX(cout, pol_x[i], ord_p);
          printZZX(cout, pol_y[i], ord_p);
          cout << endl;
      }
    }

    Ctxt ctxt_x(m_pk);
    Ctxt ctxt_y(m_pk);
    ea.encrypt(ctxt_x, m_pk, pol_x);
    ea.encrypt(ctxt_y, m_pk, pol_y);
    
    Ctxt ctxt_res(m_pk);

    // comparison function
    compare(ctxt_res, ctxt_x, ctxt_y, is_randomized);

    if(m_verbose)
    {
      cout << "Input" << endl;
      for(int i = 0; i < nslots; i++)
      {
          printZZX(cout, pol_x[i], ord_p);
          printZZX(cout, pol_y[i], ord_p);
          cout << endl;
      }

      cout << "Output" << endl;
      print_decrypted(ctxt_res);
      cout << endl;
    }
    printNamedTimer(cout, "ComparisonCircuit");
    printNamedTimer(cout, "RandomEqualityCircuit");
    printNamedTimer(cout, "EqualityCircuit");
    printNamedTimer(cout, "ShiftMul");
    printNamedTimer(cout, "ShiftAdd");
    printNamedTimer(cout, "Comparison");

    // remove the line below if it gives bizarre results 
    ctxt_res.cleanUp();
    capacity = ctxt_res.bitCapacity();
    cout << "Final capacity: " << capacity << endl;
    if (capacity < min_capacity)
      min_capacity = capacity;
    cout << "Min. capacity: " << min_capacity << endl;
    cout << "Final size: " << ctxt_res.logOfPrimeSet()/log(2.0) << endl;
    ea.decrypt(ctxt_res, m_sk, decrypted);

    for(int i = 0; i < numbers_size; i++)
    { 
      if (decrypted[i * m_expansionLen] != expected_result[i * m_expansionLen])
      {
        printf("Slot %ld: ", i * m_expansionLen);
        printZZX(cout, decrypted[i * m_expansionLen], ord_p);
        cout << endl;
        cout << "Failure" << endl;
        return;
      }
    }
    cout << endl;
  }
}