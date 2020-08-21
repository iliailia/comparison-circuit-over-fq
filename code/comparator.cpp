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

void Comparator::compute_poly_params()
{
	// get p
	ZZ p = ZZ(m_context.zMStar.getP());
	long p_long = conv<long>(p);

	// hardcoded babystep numbers
	map<unsigned long, unsigned long> bs_nums
	{
  		{13, 3}, // 3 (16), 1..5
  		{47, 5}, // 5 (26), 2..11 
  		{61, 6}, // 6 (27), 4..8 
  		{67, 5}, // 5 (27), 4..8
  		{131, 8}, // 8 (32), 4..11
  		{167, 10}, // 10 (36), 8..12
  		{173, 10}  // 10 (36), 8..12
  	};

  	m_bs_num = -1;
  	if(bs_nums.count(p_long) > 0)
  		m_bs_num = bs_nums[p_long];

  	long d = deg(m_poly);

  	// How many baby steps: set sqrt(n/2), rounded up/down to a power of two

	// FIXME: There may be some room for optimization here: it may be possible to choose this number as something other than a power of two and still maintain optimal depth, in principle we can try all possible values of m_babystep_num between two consecutive powers of two and choose the one that gives the least number of multiplies, conditioned on minimum depth.

  	if (m_bs_num <= 0) 
	{
		long kk = static_cast<long>(sqrt(d/2.0));
		m_bs_num = 1L << NextPowerOfTwo(kk);

    	// heuristic: if k>>kk then use a smaler power of two
    	if ((m_bs_num==16 && d>167) || (m_bs_num>16 && m_bs_num>(1.44*kk)))
      		m_bs_num /= 2;
  	}

	if(m_verbose)
	{
		cout << "Number of baby steps: " << m_bs_num << endl;
	}

	// n = ceil(deg(p)/k), deg(p) >= k*n
	m_gs_num = divc(d,m_bs_num);      

	// If n is not a power of two, ensure that poly is monic and that
	// its degree is divisible by k, then call the recursive procedure

	// top coefficient is equal to (p^2 - 1)/8 mod p
	// its inverse is equal to -8 mod p
	m_top_coef = LeadCoeff(m_poly);
	ZZ topInv = ZZ(-8) % p; // the inverse mod p of the top coefficient of poly (if any)
	bool divisible = (m_gs_num * m_bs_num == d); // is the degree divisible by k?
	//long nonInvertibe = InvModStatus(topInv, top, p);
	   // 0 if invertible, 1 if not

	// FIXME: There may be some room for optimization below: instead of
	// adding a term X^{n*k} we can add X^{n'*k} for some n'>n, so long
	// as n' is smaller than the next power of two. We could save a few
	// multiplications since giantStep[n'] may be easier to compute than
	// giantStep[n] when n' has fewer 1's than n in its binary expansion.

	m_extra_coef = ZZ::zero();    // extra!=0 denotes an added term extra*X^{n*k}

	if (m_gs_num != (1L << NextPowerOfTwo(m_gs_num)))
	{
		if (!divisible) 
		{  // need to add a term
	    	m_top_coef = NTL::to_ZZ(1);  // new top coefficient is one
	    	topInv = m_top_coef;    // also the new inverse is one
	    	// set extra = 1 - current-coeff-of-X^{n*k}
	    	m_extra_coef = SubMod(m_top_coef, coeff(m_poly, m_gs_num * m_bs_num), p);
	    	SetCoeff(m_poly, m_gs_num * m_bs_num); // set the top coefficient of X^{n*k} to one
		}

		if (!IsOne(m_top_coef)) 
		{
	    	m_poly *= topInv; // Multiply by topInv to make into a monic polynomial
	    	for (long i = 0; i <= m_gs_num * m_bs_num; i++) rem(m_poly[i], m_poly[i], p);
	    	m_poly.normalize();
		}
	}

	long top_deg = conv<long>(p-1) >> 1;
	m_baby_index = top_deg % m_bs_num;
	m_giant_index = top_deg / m_bs_num;
	if(m_baby_index == 0)
	{
		m_baby_index = m_bs_num;
		m_giant_index -= 1;
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
	m_poly = ZZX(INIT_MONO, 0, 0);

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

		m_poly += ZZX(INIT_MONO, (indx-1) >> 1, rep(coef));
	}

	if (m_verbose)
	{
		cout << "Comparison polynomial: " << endl;
		printZZX(cout, m_poly, (p-1)>>1);
		cout << endl;
	}

	compute_poly_params();
}


void Comparator::slot_init()
{
	/*
	// get p
	long p = m_context.zMStar.getP();

	// coefficient vector of a candidate generator
	vector<ZZ_p> gen_coefs(m_slotDeg);

	// get the generator of the entire slot
	ZZX full_gen = m_context.slotRing->G;

	//ZZX cand_gen;
	while(true)
	{
		for (int iCoef = 0; iCoef < m_slotDeg; iCoef++)
		{
			gen_coefs[iCoef].init(ZZ(p));
			random(gen_coefs[iCoef]);
		}
		//cand_gen = 
	}
	*/
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
	slot_init();
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

void Comparator::less_than_mod_2(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	//Comp(x,y) = y(x+1)
	cout << "Compute comparison polynomial" << endl;

	// x + 1
	Ctxt x_plus_1 = ctxt_x;
	x_plus_1.addConstant(ZZ(1));

	// y(x+1)
	ctxt_res = ctxt_y;
	ctxt_res.multiplyBy(x_plus_1);

	if(m_verbose)
	  {
	    print_decrypted(ctxt_res);
	    cout << endl;
	  }
}

void Comparator::less_than_mod_3(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	//Comp(x,y) = -y(x-y)(x+1)
	cout << "Compute comparison polynomial" << endl;

	// x + 1
	Ctxt x_plus_1 = ctxt_x;
	x_plus_1.addConstant(ZZ(1));

	// y(x - y)
	ctxt_res = ctxt_x;
	ctxt_res -= ctxt_y;
	ctxt_res.multiplyBy(ctxt_y);

	// -y(x-y)(x+1)
	ctxt_res.multiplyBy(x_plus_1);
	ctxt_res.negate();

	if(m_verbose)
	  {
	    print_decrypted(ctxt_res);
	    cout << endl;
	  }
}

void Comparator::less_than_mod_5(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	//Comp(x,y)=−(x+1) y(x−y) (x (x+1) − y(x−y)).
	cout << "Compute comparison polynomial" << endl;

	// y(x - y)
	Ctxt y_x_min_y = ctxt_x;
	y_x_min_y -= ctxt_y;
	y_x_min_y.multiplyBy(ctxt_y);

	// x + 1
	Ctxt x_plus_1 = ctxt_x;
	x_plus_1.addConstant(ZZ(1));

	// x * (x+1)
	ctxt_res = ctxt_x;
	ctxt_res.multiplyBy(x_plus_1);

	// x * (x+1) - y * (x-y)
	ctxt_res -= y_x_min_y;

	// y * (x-y) * (x * (x+1) - y * (x-y))
	ctxt_res.multiplyBy(y_x_min_y);

	// -(x+1) * y * (x-y) * (x * (x+1) - y * (x-y))
	ctxt_res.multiplyBy(x_plus_1);
	ctxt_res.negate();

	if(m_verbose)
	  {
	    print_decrypted(ctxt_res);
	    cout << endl;
	  }
}

void Comparator::less_than_mod_7(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	// Comp(x,y) = -y(x-y)(x+1)(x(x+1)(x(x+1)+3) + 5y(x-y)(x(x+1)+2x+3y(x-y)))
	cout << "Compute comparison polynomial" << endl;

	// x
	Ctxt y_x_min_y = ctxt_x;
	// x - y
	y_x_min_y -= ctxt_y;
	// y(x-y)
	y_x_min_y.multiplyBy(ctxt_y);

	// x + 1
	Ctxt x_plus_1 = ctxt_x;
	x_plus_1.addConstant(ZZ(1));

	// x(x+1)
	Ctxt x_x_plus_1 = x_plus_1;
	x_x_plus_1.multiplyBy(ctxt_x);

	// x(x+1)(x(x+1)+3)
	Ctxt tmp = x_x_plus_1;
	tmp.addConstant(ZZ(3));
	tmp.multiplyBy(x_x_plus_1);
	ctxt_res = tmp;

	// x(x+1) + 2x + 3y(x-y)
	tmp = y_x_min_y;
	tmp.multByConstant(ZZ(3));
	tmp += x_x_plus_1;
	tmp += ctxt_x;
	tmp += ctxt_x;

	// 5y(x-y)(x(x+1) + 2x + 3y(x-y))
	tmp.multiplyBy(y_x_min_y);
	tmp.multByConstant(ZZ(5));

	// (x^2+x)(x^2+x+3) + 5y(x-y)(x^2+x+2x+3y(x-y))
	ctxt_res += tmp;

	// -y(x-y)(x+1)((x^2+x)(x^2+x+3) + 5y(x-y)(x^2+x+2x+3y(x-y)))
	ctxt_res.multiplyBy(y_x_min_y);
	ctxt_res.multiplyBy(x_plus_1);
	ctxt_res.negate();

	if(m_verbose)
	{
		print_decrypted(ctxt_res);
		cout << endl;
	}
}

void Comparator::less_than_mod_11(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
        // Y = y(x-y)
        // X = x(x+1)
	// Comp(x,y) = Y(x+1)[(7X + 8X^2 + 8X^3 - X^4) + (3 + (7x+4)(X+1)^2 + 4(X+1)^3)Y + (x^2+6x + 4(X +5x)^2)Y^2 + 5(x^2+7x)Y^3 - Y^4]
	cout << "Compute comparison polynomial" << endl;

	Ctxt Y = ctxt_x, x_plus_1 = ctxt_x;

	// x - y
	Y -= ctxt_y;
	// Y = y(x-y)
	Y.multiplyBy(ctxt_y);
	// Y(x+1)
	x_plus_1.addConstant(ZZ(1));
	x_plus_1.multiplyBy(Y);

	// X = x^2 + x
	Ctxt X = ctxt_x;
	X.multiplyBy(ctxt_x);
	X += ctxt_x;
	
	
	// f0 = 7X + 8X^2 + 8X^3 - X^4
	Ctxt f0 = X, tmp = X, X_2 = X, X_3 = X;
	// 7X
	f0.multByConstant(ZZ(7));
	
	// 8(X^2 + X^3)
	X_2.multiplyBy(X);
	X_3.multiplyBy(X_2);
	
	tmp = X_2;
	tmp += X_3;
	tmp += tmp;
	tmp +=tmp;
	tmp += tmp;
	f0 += tmp;

	//X^4
	Ctxt X_4 = X_2;
	X_4.multiplyBy(X_2);
	f0 -= X_4;

	// X' = X + 1
	//f1 = Y(3 + (7x+4)X'^2 + 4X'^3)
	Ctxt f1 = X_2;
	f1 += X;
	f1 += X;
	f1.addConstant(ZZ(1));

	Ctxt X1_2 = f1;
        
	// 7x + 4 = 4(1-x) mod 11 
	tmp = ctxt_x;
	tmp.negate();
	tmp.addConstant(ZZ(1));
	tmp += tmp;
	tmp+=tmp;

	tmp.multiplyBy(Y);
	f1.multiplyBy(tmp);

	// 4X'^3*Y
	tmp = X;
	tmp.addConstant(ZZ(1));
	tmp += tmp;
	tmp+= tmp;
	tmp.multiplyBy(Y);
	tmp.multiplyBy(X1_2);

	f1 += tmp;
	f1 += Y;
	f1 += Y;
	f1 += Y;

	//X' = x^2 + 6x = X + 5x
	//f2 = 3X' + 4X'^2
	X1_2 = ctxt_x;
	X1_2 += X1_2;
	X1_2 += X1_2;
	X1_2 += ctxt_x;
	X1_2 += X;

	tmp = X1_2;
	tmp += tmp;
	tmp += tmp;
	tmp.addConstant(ZZ(3));
	tmp.multiplyBy(X1_2);

	Ctxt Y2 = Y;
	Y2.multiplyBy(Y);
	Ctxt f2 = Y2;
	
	f2.multiplyBy(tmp);

	// X' = x^2 + 7x = X1_2 + x
	// f3 = 5X
	X1_2 += ctxt_x;
	Ctxt f3 = X1_2;
	f3 += f3;
	f3 += f3;
	f3 += X1_2;
	f3.multiplyBy(Y);
	f3.multiplyBy(Y2);

	f0 += f1;
	f2 += f3;
	f0 += f2;

	Y2.multiplyBy(Y2);

	f0 -= Y2;

	ctxt_res = f0;
	ctxt_res.multiplyBy(x_plus_1);
	
	if(m_verbose)
	{
		print_decrypted(ctxt_res);
		cout << endl;
	}
}



void Comparator::evaluate_poly(Ctxt& ret, const Ctxt& x, const Ctxt& x2) const
{
	// get p
	ZZ p = ZZ(m_context.zMStar.getP());

	DynamicCtxtPowers babyStep(x2, m_bs_num);
	const Ctxt& x2k = babyStep.getPower(m_bs_num);

	DynamicCtxtPowers giantStep(x2k, m_gs_num);

	// Special case when deg(p)>k*(2^e -1)
	if (m_gs_num == (1L << NextPowerOfTwo(m_gs_num))) 
	{ // n is a power of two
		cout << "I'm computing degPowerOfTwo" << endl;
    	degPowerOfTwo(ret, m_poly, m_bs_num, babyStep, giantStep);
    }
    else
    {
	  	recursivePolyEval(ret, m_poly, m_bs_num, babyStep, giantStep);

	  	if (!IsOne(m_top_coef)) 
	  	{
	    	ret.multByConstant(m_top_coef);
		}

		if (!IsZero(m_extra_coef)) 
		{ // if we added a term, now is the time to subtract back
	    	Ctxt topTerm = giantStep.getPower(m_gs_num);
	    	topTerm.multByConstant(m_extra_coef);
	    	ret -= topTerm;
		}
	}
	ret.multiplyBy(x);

	Ctxt top_term = babyStep.getPower(m_baby_index);
	top_term.multiplyBy(giantStep.getPower(m_giant_index));
	top_term.multByConstant(ZZ((p+1)>> 1));

	ret += top_term;
}


void Comparator::less_than(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
  FHE_NTIMER_START(ComparisonCircuit);

  unsigned long p = m_context.zMStar.getP();

  if(p == 2)
    {
      less_than_mod_2(ctxt_res, ctxt_x, ctxt_y);
      return;
    }

  
  if(p == 3)
  {
  	less_than_mod_3(ctxt_res, ctxt_x, ctxt_y);
  	return;
  }

  if(p == 5)
  {
  	less_than_mod_5(ctxt_res, ctxt_x, ctxt_y);
  	return;
  }

  if(p == 7)
  {
  	less_than_mod_7(ctxt_res, ctxt_x, ctxt_y);
  	return;
  }

  if(p == 11)
    {
      less_than_mod_11(ctxt_res, ctxt_x, ctxt_y);
      return;
    }
  
  // Subtraction z = x - y
  cout << "Subtraction" << endl;
  Ctxt ctxt_z = ctxt_x;
  ctxt_z -= ctxt_y;

  if(m_verbose)
  {
    print_decrypted(ctxt_z);
    cout << endl;
  }

  // compute polynomial function for 'z < 0'
  cout << "Compute comparison polynomial" << endl;

  // z^2
  Ctxt ctxt_z2 = ctxt_z;
  ctxt_z2.square();

  evaluate_poly(ctxt_res, ctxt_z, ctxt_z2);

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
  if (p == 3 || p == 5 || p == 7)
  {
  	enc_base = p;
  }

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

    const FHEtimer* comp_timer = getTimerByName("Comparison");

    cout << "Avg. time per 64-bit integer: " << 1000.0 * comp_timer->getTime()/static_cast<double>(run+1)/static_cast<double>(numbers_size) << " ms" << endl;

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
