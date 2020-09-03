#include "comparator.h"
#include "tools.h"
#include <helib/debugging.h>
#include <helib/polyEval.h>
#include <random>
#include <map> 
#include <NTL/ZZ_pE.h>
#include <NTL/mat_ZZ_pE.h>
#include <helib/Ptxt.h>

DoubleCRT Comparator::create_shift_mask(double& size, long shift)
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
	ZZX mask_zzx;
	ea.encode(mask_zzx, mask_vec);

	size = conv<double>(embeddingLargestCoeff(mask_zzx, m_context.zMStar));

	DoubleCRT mask_crt = DoubleCRT(mask_zzx, m_context, m_context.allPrimes());
	return mask_crt;
}

void Comparator::create_all_shift_masks()
{
	long shift = 1;
	while (shift < m_expansionLen)
	{
		double size;
	    DoubleCRT mask_ptxt = create_shift_mask(size, -shift);
	    m_mulMasks.push_back(mask_ptxt);
	    m_mulMasksSize.push_back(size);

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
		{11, 3}, // 3 (19), 1..4
  		{13, 3}, // 3 (16), 1..5
  		{17, 4}, // 4 (13), 1..5
  		{19, 3}, // 3 (14), 1..2
  		{29, 5}, // 5 (16), 1..6
  		{31, 5}, // 5 (15), 4..6
  		{47, 5}, // 5 (18), 2..11 
  		{61, 6}, // 6 (19), 4..8 
  		{67, 5}, // 5 (20), 4..8
  		{71, 4}, // 4 (20), 3..7
  		{101, 7}, // 7 (21), 4..8
  		{131, 8}, // 8 (24), 4..11
  		{167, 10}, // 10 (26), 8..12
  		{173, 10},  // 10 (26), 8..12
  		{271, 9},  // 9 (31), 9..10
  		{401, 12}  // 12 (33), 9..14
  	};

  	m_bs_num = -1;
  	if(bs_nums.count(p_long) > 0)
  		m_bs_num = bs_nums[p_long];

  	long d = deg(m_univar_poly);

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
	m_top_coef = LeadCoeff(m_univar_poly);
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
	    	m_extra_coef = SubMod(m_top_coef, coeff(m_univar_poly, m_gs_num * m_bs_num), p);
	    	SetCoeff(m_univar_poly, m_gs_num * m_bs_num); // set the top coefficient of X^{n*k} to one
		}

		if (!IsOne(m_top_coef)) 
		{
	    	m_univar_poly *= topInv; // Multiply by topInv to make into a monic polynomial
	    	for (long i = 0; i <= m_gs_num * m_bs_num; i++) rem(m_univar_poly[i], m_univar_poly[i], p);
	    	m_univar_poly.normalize();
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

	if(m_isUnivar)
	{
		// polynomial coefficient
		ZZ_p coef;
		coef.init(ZZ(p));

		// field element
		ZZ_p field_elem;
		field_elem.init(ZZ(p));

		// initialization of the univariate comparison polynomial
		m_univar_poly = ZZX(INIT_MONO, 0, 0);

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

			m_univar_poly += ZZX(INIT_MONO, (indx-1) >> 1, rep(coef));
		}

		compute_poly_params();
	}
	else
	{
		// computing the coefficients of the bivariate polynomial
		m_bivar_coefs.SetDims(p,p);

		// y^{p-1}
		m_bivar_coefs[0][p-1] = ZZ(1);

		// (p+1)/2 * x^{(p-1)/2} * y^{(p-1)/2}
		m_bivar_coefs[(p-1) >> 1][(p-1) >> 1] = ZZ((p+1) >> 1);

		// iterator
		ZZ_p field_elem;
		field_elem.init(ZZ(p));

		// inner sum
		ZZ_p inner_sum;
		inner_sum.init(ZZ(p));
		
		// outer sum
		ZZ_p outer_sum;
		outer_sum.init(ZZ(p));

		for (long i = 1; i < p; i++)
		{
			for (long j = 1; j < p; j++)
			{
				// x^i * y^i have the zero coefficient except for i = (p-1)/2
				if (i == j)
					continue;

				outer_sum = 0;
				// sum_{a=1}^{p-1} a^{p-1-i} sum_{b = a+1}^{p-1} b^{p-1-j} 
				for (long a = 1; a < p; a++)
				{
					inner_sum = 0;
					// sum_{b = a+1}^{p-1} b^{p-1-j} 
					for (long b = a+1; b < p; b++)
					{
						// b^{p-1-j}
						field_elem = b;
						field_elem = power(field_elem, p - 1 - j);

						inner_sum += field_elem;
					}
					// a^{p-1-i}
					field_elem = a;
					field_elem = power(field_elem, p - 1 - i);

					inner_sum *= field_elem;
					outer_sum += inner_sum;
				}
				m_bivar_coefs[i][j] = rep(outer_sum);
			}
		}

		cout << "Bivariate coefficients" << endl << m_bivar_coefs << endl;

		if (m_verbose)
		{
			cout << "Comparison polynomial: " << endl;
			printZZX(cout, m_univar_poly, (p-1)>>1);
			cout << endl;
		}
	}
}


void Comparator::extraction_init()
{
	// get the total number of slots
	const EncryptedArray& ea = *(m_context.ea);
	long nslots = ea.size();

	// get p
	long p = m_context.zMStar.getP();

	// get the order of p
	long d = m_context.zMStar.getOrdP();

	// get the defining polynomial of a slot mod p
	ZZX def_poly = m_context.slotRing->G;

	//cout << "Def. poly" << def_poly << endl;

	ZZ_pX def_poly_p;
	for (long iCoef = 0; iCoef <= deg(def_poly); iCoef++)
	{
		ZZ_p coef;
		coef.init(ZZ(p));
		coef = conv<long>(def_poly[iCoef]);
		SetCoeff(def_poly_p, iCoef, coef);
	}
	//cout << "Def. poly mod p " << def_poly_p << endl;
	
	// build the trace matrix
	mat_ZZ_pE trace_matrix;
	trace_matrix.SetDims(d, d);
	
	ZZ_pE prim_elem;
	prim_elem.init(def_poly_p);
	prim_elem = conv<ZZ_pE>(ZZ_pX(INIT_MONO, 1, ZZ_p(1)));

	//cout << "Primitive element " << prim_elem << endl;

	ZZ_pE coef;
	coef.init(def_poly_p);

	for (long iRow = 0; iRow < d; iRow++)
	{
		for (long iCol = 0; iCol < d; iCol++)
		{
			// x^(iRow * p^iCol)
			coef = power(prim_elem, iRow * power_long(p, iCol));
			trace_matrix[iRow][iCol] = coef;
		}
	}
	//cout << "Trace matrix" << trace_matrix << endl;

	mat_ZZ_pE inv_trace_matrix = inv(trace_matrix);
	//cout << "Inverse of trace matrix" << inv_trace_matrix << endl;

	//cout << "Extraction consts: " << endl;
	for (long iCoef = 0; iCoef < d; iCoef++)
	{
		vector<DoubleCRT> tmp_crt_vec;
		vector<double> size_vec;

		for (long iFrob = 0; iFrob < d; iFrob++)
		{
			ZZX tmp = conv<ZZX>(rep(inv_trace_matrix[iFrob][iCoef]));

			//cout << tmp << endl;

			vector<ZZX> vec_const(nslots, tmp);
			ea.encode(tmp, vec_const);

			DoubleCRT tmp_crt(tmp, m_context, m_context.allPrimes());

			double const_size = conv<double>(embeddingLargestCoeff(tmp, m_context.zMStar));
			size_vec.push_back(const_size);

			tmp_crt_vec.push_back(tmp_crt);
		}
		m_extraction_const.push_back(tmp_crt_vec);
		m_extraction_const_size.push_back(size_vec);
	}
}

void Comparator::extract_mod_p(vector<Ctxt>& mod_p_coefs, const Ctxt& ctxt_x) const
{
	FHE_NTIMER_START(Extraction);
	mod_p_coefs.clear();

	if (m_slotDeg == 1)
	{
		mod_p_coefs.push_back(ctxt_x);
		return;
	}

	const EncryptedArray& ea = *(m_context.ea);
	long nslots = ea.size();

	// get max slot degree
	long d = m_context.zMStar.getOrdP();

	// TODO: how to use key switching hoisting from CRYPTO'18?
	vector<Ctxt> ctxt_frob(d-1, ctxt_x);
	for(long iFrob = 1; iFrob < d; iFrob++)
	{
		ctxt_frob[iFrob-1].frobeniusAutomorph(iFrob);
	} 

	for(long iCoef = 0; iCoef < m_slotDeg; iCoef++)
	{
		cout << "Extract coefficient " << iCoef << endl;
		Ctxt mod_p_ctxt = ctxt_x;

		mod_p_ctxt.multByConstant(m_extraction_const[iCoef][0], m_extraction_const_size[iCoef][0]);

		for(long iFrob = 1; iFrob < d; iFrob++)
		{
			Ctxt tmp = ctxt_frob[iFrob-1];
			tmp.multByConstant(m_extraction_const[iCoef][iFrob], m_extraction_const_size[iCoef][iFrob]);
			mod_p_ctxt += tmp; 
		}
		mod_p_coefs.push_back(mod_p_ctxt);
	}
	FHE_NTIMER_STOP(Extraction);
}

Comparator::Comparator(const Context& context, unsigned long d, unsigned long expansion_len, const SecKey& sk, bool verbose): m_context(context), m_slotDeg(d), m_expansionLen(expansion_len), m_sk(sk), m_pk(sk), m_verbose(verbose)
{
	//determine the order of p in (Z/mZ)*
	m_isUnivar = true;
	unsigned long p = context.zMStar.getP();
	if (p == 2 || p == 3 || p == 5 || p == 7 || p == 11)
		m_isUnivar = false;
	unsigned long ord_p = context.zMStar.getOrdP();
	//check that the extension degree divides the order of p
	if (ord_p < d != 0)
	{
		throw invalid_argument("Field extension must be larger than the order of the plaintext modulus\n");
	}

	create_all_shift_masks();
	create_poly();
	extraction_init();
}

const DoubleCRT& Comparator::get_mask(double& size, long index) const
{
	size = m_mulMasksSize[index];
	return m_mulMasks[index];
}

const ZZX& Comparator::get_poly() const
{
	return m_univar_poly;
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
	double size;
	DoubleCRT mask = get_mask(size, index);
	ctxt.multByConstant(mask, size);
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
	double mask_size;
	DoubleCRT mask = get_mask(mask_size, index);
	ctxt.multByConstant(mask, mask_size);

	// add 1 to masked slots
	ctxt.addConstant(ZZ(1));
	mask.Negate();
	ctxt.addConstant(mask, mask_size);

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

void Comparator::mapTo01_subfield(Ctxt& ctxt, long pow) const
{
  FHE_NTIMER_START(MapTo01);	
  // get EncryptedArray
  const EncryptedArray& ea = *(m_context.ea);

  // get p
  long p = ctxt.getPtxtSpace();
  if (p != ea.getPAlgebra().getP()) // ptxt space is p^r for r>1
    throw helib::LogicError("mapTo01 not implemented for r>1");

  if (p % pow != 0)
  	throw helib::LogicError("Exponent must divide p");

  if (p > 2)
    ctxt.power((p - 1) / pow); // set y = x^{p-1}

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

void Comparator::evaluate_poly(Ctxt& ret, Ctxt& ctxt_p_1, const Ctxt& x) const
{
	FHE_NTIMER_START(ComparisonCircuitUnivar);
	// get p
	ZZ p = ZZ(m_context.zMStar.getP());

	if (p > ZZ(3)) //if p > 3, use the generic Paterson-Stockmeyer strategy
	{
	  // z^2
	  	Ctxt x2 = x;
	  	x2.square();

		DynamicCtxtPowers babyStep(x2, m_bs_num);
		const Ctxt& x2k = babyStep.getPower(m_bs_num);

		DynamicCtxtPowers giantStep(x2k, m_gs_num);

		// Special case when deg(p)>k*(2^e -1)
		if (m_gs_num == (1L << NextPowerOfTwo(m_gs_num))) 
		{ // n is a power of two
			cout << "I'm computing degPowerOfTwo" << endl;
	    	degPowerOfTwo(ret, m_univar_poly, m_bs_num, babyStep, giantStep);
	    }
	    else
	    {
		  	recursivePolyEval(ret, m_univar_poly, m_bs_num, babyStep, giantStep);

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

		// TODO: depth here is not optimal
		Ctxt top_term = babyStep.getPower(m_baby_index);
		top_term.multiplyBy(giantStep.getPower(m_giant_index));

		ctxt_p_1 = top_term; 

		top_term.multByConstant(ZZ((p+1)>> 1));

		ret += top_term;
	}
	else //circuit for p=3
	{
		ret = x;

		ctxt_p_1 = x;
		ctxt_p_1.square();

		Ctxt top_term = ctxt_p_1;
		top_term.multByConstant(ZZ(2));

		ret += top_term;
	}
	FHE_NTIMER_STOP(ComparisonCircuitUnivar);
}


void Comparator::less_than_bivar(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
  FHE_NTIMER_START(ComparisonCircuitBivar);

  // uncomment if you want to compare with Tan et al.
  less_than_bivar_tan(ctxt_res, ctxt_x, ctxt_y);
  return;

  unsigned long p = m_context.zMStar.getP();

  if(p == 2)
  {
    less_than_mod_2(ctxt_res, ctxt_x, ctxt_y);
  }
  
  if(p == 3)
  {
  	less_than_mod_3(ctxt_res, ctxt_x, ctxt_y);
  }

  if(p == 5)
  {
  	less_than_mod_5(ctxt_res, ctxt_x, ctxt_y);
  }

  if(p == 7)
  {
  	less_than_mod_7(ctxt_res, ctxt_x, ctxt_y);
  }

  if(p == 11)
  {
    less_than_mod_11(ctxt_res, ctxt_x, ctxt_y);
  }

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  FHE_NTIMER_STOP(ComparisonCircuitBivar);
}

void Comparator::less_than_bivar_tan(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	long p = m_context.zMStar.getP();

	DynamicCtxtPowers x_powers(ctxt_x, p-1);
	DynamicCtxtPowers y_powers(ctxt_y, p-1);

	ctxt_res = y_powers.getPower(p-1);

	for (long i = 1; i < p; i++)
	{
		// zero ciphertext
		Ctxt sum = Ctxt(ctxt_x.getPubKey());
		for (long j = 1; j < p; j++)
		{
			if (m_bivar_coefs[i][j] == ZZ(0))
				continue;
			Ctxt tmp = y_powers.getPower(j);
			tmp.multByConstant(m_bivar_coefs[i][j]);
			sum += tmp;
		}
		sum.multiplyBy(x_powers.getPower(i));
		ctxt_res += sum;
	}
}

void Comparator::is_zero(Ctxt& ctxt_res, const Ctxt& ctxt_z, long pow) const
{
  FHE_NTIMER_START(EqualityCircuit);

  ctxt_res = ctxt_z;

  //compute mapTo01: (z_i)^{p^d-1}
  cout << "Mapping to 0 and 1" << endl;
  mapTo01_subfield(ctxt_res, pow);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  cout << "Computing NOT" << endl;
  //compute 1 - mapTo01(z_i)
  ctxt_res.negate();
  ctxt_res.addConstant(ZZ(1));

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  FHE_NTIMER_STOP(EqualityCircuit);
}

void Comparator::compare(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	FHE_NTIMER_START(Comparison);

	vector<Ctxt> ctxt_less_p;
	vector<Ctxt> ctxt_eq_p;

	// bivariate circuit
	if (!m_isUnivar)
	{
		cout << "Extraction" << endl;
		// extract mod p coefficients
		vector<Ctxt> ctxt_x_p;
		extract_mod_p(ctxt_x_p, ctxt_x);

		if(m_verbose)
	    {
	    	for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
	    	{
	    		cout << "Ctxt x with coefficient " << iCoef << endl;
		    	print_decrypted(ctxt_x_p[iCoef]);
		    	cout << endl;
		    }
		}

		vector<Ctxt> ctxt_y_p;
		extract_mod_p(ctxt_y_p, ctxt_y);

		if(m_verbose)
	    {
	    	for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
	    	{
	    		cout << "Ctxt y with coefficient " << iCoef << endl;
		    	print_decrypted(ctxt_y_p[iCoef]);
		    	cout << endl;
		    }
		}

		cout << "Compute the less-than function modulo p" << endl;
		for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
		{
			Ctxt ctxt_tmp = Ctxt(ctxt_x.getPubKey());
			less_than_bivar(ctxt_tmp, ctxt_x_p[iCoef], ctxt_y_p[iCoef]);
			ctxt_less_p.push_back(ctxt_tmp);
		}

		cout << "Compute the equality function modulo p" << endl;
		for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
		{
			// Subtraction z = x - y
			cout << "Subtraction" << endl;
			Ctxt ctxt_z = ctxt_x_p[iCoef];
			ctxt_z -= ctxt_y_p[iCoef];
			Ctxt ctxt_tmp = Ctxt(ctxt_z.getPubKey());
			is_zero(ctxt_tmp, ctxt_z);
			ctxt_eq_p.push_back(ctxt_tmp);
		}
	}
	else // univariate circuit
	{
		// Subtraction z = x - y
		cout << "Subtraction" << endl;
		Ctxt ctxt_z = ctxt_x;
		ctxt_z -= ctxt_y;

		if(m_verbose)
		{
			print_decrypted(ctxt_z);
			cout << endl;
		}

		// extract mod p coefficients
		cout << "Extraction" << endl;
		vector<Ctxt> ctxt_z_p;
		extract_mod_p(ctxt_z_p, ctxt_z);

		if(m_verbose)
	    {
	    	for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
	    	{
	    		cout << "Ctxt x with coefficient " << iCoef << endl;
		    	print_decrypted(ctxt_z_p[iCoef]);
		    	cout << endl;
		    }
		}

		cout << "Compute the less-than and equality functions modulo p" << endl;
		for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
		{
			Ctxt ctxt_tmp = Ctxt(ctxt_z.getPubKey());
			Ctxt ctxt_tmp_eq = Ctxt(ctxt_z.getPubKey());

			// compute polynomial function for 'z < 0'
			cout << "Compute univariate comparison polynomial" << endl;
			evaluate_poly(ctxt_tmp, ctxt_tmp_eq, ctxt_z_p[iCoef]);

			if(m_verbose)
			{
			  cout << "Result of the less-than function" << endl;
			  print_decrypted(ctxt_tmp);
			  cout << endl;
			}

			ctxt_less_p.push_back(ctxt_tmp);

			cout << "Computing NOT" << endl;
			//compute 1 - mapTo01(r_i*(x_i - y_i))
			ctxt_tmp_eq.negate();
			ctxt_tmp_eq.addConstant(ZZ(1));

			if(m_verbose)
			{
			  cout << "Result of the equality function" << endl;
			  print_decrypted(ctxt_tmp_eq);
			  cout << endl;
			}

			ctxt_eq_p.push_back(ctxt_tmp_eq);
		}	
	}

	cout << "Compare digits" << endl;
	Ctxt ctxt_less = ctxt_less_p[m_slotDeg-1];
	Ctxt ctxt_eq = ctxt_eq_p[m_slotDeg-1];

	for (long iCoef = m_slotDeg-2; iCoef >= 0; iCoef--)
	{
		Ctxt tmp = ctxt_eq;
		tmp.multiplyBy(ctxt_less_p[iCoef]);
		ctxt_less += tmp;

		ctxt_eq.multiplyBy(ctxt_eq_p[iCoef]);
	}

	if(m_verbose)
	{
		cout << "Comparison results" << endl;
		print_decrypted(ctxt_less);
		cout << endl;

		cout << "Equality results" << endl;
		print_decrypted(ctxt_eq);
		cout << endl;
	}

	if(m_expansionLen == 1)
	{
		ctxt_res = ctxt_less;
		return;
	}


	//compute running products: prod_i 1 - (x_i - y_i)^{p^d-1}
	cout << "Rotating and multiplying slots with equalities" << endl;
	shift_and_mul(ctxt_eq, 0);

	if(m_verbose)
	{
		print_decrypted(ctxt_eq);
		cout << endl;
	}

	//Remove the least significant digit and shift to the left
	cout << "Remove the least significant digit" << endl;
	batch_shift_for_mul(ctxt_eq, 0, -1);

	if(m_verbose)
	{
		print_decrypted(ctxt_eq);
		cout << endl;
	}

	cout << "Final result" << endl;

	ctxt_res = ctxt_eq;
	ctxt_res.multiplyBy(ctxt_less);
	shift_and_add(ctxt_res, 0);

	if(m_verbose)
	{
		print_decrypted(ctxt_res);
		cout << endl;
	}

	if(m_verbose)
    {
      cout << "Input x: " << endl;
      print_decrypted(ctxt_x);
      cout << endl;
      cout << "Input y: " << endl;
      print_decrypted(ctxt_y);
      cout << endl;
    }	

    FHE_NTIMER_STOP(Comparison);
}

void Comparator::int_to_slot(ZZX& poly, unsigned long input, unsigned long enc_base) const
{ 
    vector<long> decomp;

    //decomposition of a digit
    digit_decomp(decomp, input, enc_base, m_slotDeg);
    poly = ZZX(INIT_MONO, 0, 0);
    for (int iCoef = 0; iCoef < m_slotDeg; iCoef++)
    {
        poly+=ZZX(INIT_MONO, iCoef, decomp[iCoef]);
    }
}

void Comparator::test(long runs) const
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

  //encoding base, ((p+1)/2)^d
  //if 2-variable comparison polynomial is used, it must be p^d
  unsigned long enc_base = (p + 1) >> 1;
  if (!m_isUnivar)
  {
  	enc_base = p;
  }

  unsigned long digit_base = power_long(enc_base, m_slotDeg);

  //check that field_size^expansion_len fits into 64-bits
  int space_bit_size = static_cast<int>(ceil(m_expansionLen * log2(digit_base)));
  unsigned long input_range = LONG_MAX;
  if(space_bit_size < 64)
  {
    //input_range = power_long(field_size, expansion_len);
    input_range = power_long(digit_base, m_expansionLen);
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
      digit_decomp(decomp_int_x, input_x, digit_base, m_expansionLen);
      digit_decomp(decomp_int_y, input_y, digit_base, m_expansionLen);

      //encoding of slots
      for (int j = 0; j < m_expansionLen; j++)
      {
          //decomposition of a digit
          int_to_slot(pol_slot, decomp_int_x[j], enc_base);
          pol_x[i * m_expansionLen + j] = pol_slot;
      }

      for (int j = 0; j < m_expansionLen; j++)
      {
          //decomposition of a digit
          int_to_slot(pol_slot, decomp_int_y[j], enc_base);
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
    cout << "Start of comparison" << endl;
    compare(ctxt_res, ctxt_x, ctxt_y);

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
    printNamedTimer(cout, "Extraction");
    printNamedTimer(cout, "ComparisonCircuitBivar");
    printNamedTimer(cout, "ComparisonCircuitUnivar");
    printNamedTimer(cout, "EqualityCircuit");
    printNamedTimer(cout, "ShiftMul");
    printNamedTimer(cout, "ShiftAdd");
    printNamedTimer(cout, "Comparison");

    const FHEtimer* comp_timer = getTimerByName("Comparison");

    cout << "Avg. time per 64-bit integer: " << 1000.0 * comp_timer->getTime()/static_cast<double>(run+1)/static_cast<double>(numbers_size) << " ms" << endl;
    cout << "Number of 64-bit integers in one ciphertext "<< numbers_size << endl;

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
