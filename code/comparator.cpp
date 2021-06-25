#include "comparator.h"
#include "tools.h"
#include <helib/debugging.h>
#include <helib/polyEval.h>
#include <random>
#include <map> 
#include <NTL/ZZ_pE.h>
#include <NTL/mat_ZZ_pE.h>
#include <helib/Ptxt.h>

using namespace he_cmp;

// polynomial coefficients of bivariate polynomial decomposition as in Theorem 2 for different plaintext moduli
map<unsigned long, vector<vector<long>>> fcoefs {
  {11,
  	{
  		{-4, 4, 2, -2, -2, 2, -4, -1},
  		{5, -5, 2, -2, -3, 4},
  		{-4, 4, 4, 4},
  		{2, 5},
  		{-1}
  	}
  },
  {13,
  	{
  		{-5, 5, 5, -5, -4, 4, 6, -6, -5, -1},
  		{3, -3, -5, 5, 4, -4, -4, 5},
  		{4, -4, -6, 6, -5, 1},
  		{3, -3, 6, 1},
  		{2, 6},
  		{1}
  	}
  },
  {17, 
    {
    	{-4, 4, -8, 8, 3, -3, 7, -7, -8, 8, -4, 4, -7, -1},
    	{4, -4, -1, 1, 4, -4, -4, 4, -7, 7, -6, 7},
    	{-6, 6, -4, 4, -7, 7, 7, -7, 3, 8},
    	{-4, 4, 3, -3, -2, 2, 3, 4},
    	{2, -2, -4, 4, 2, 2},
    	{-4, 4, 7, 8},
    	{-3, 5},
    	{1}
    }
  },
  {19,
  	{
  		{4, -4, 6, -6, 8, -8, -1, 1, 6, -6, 1, -1, 8, -8, -8, -1},
		{-2, 2, -6, 6, 2, -2, -7, 7, -1, 1, -5, 5, -7, 8},
    	{2, -2, -9, 9, 2, -2, -1, 1, 7, -7, 9, 3},
    	{9, -9, -7, 7, 6, -6, 5, -5, -8, -4},
    	{7, -7, -9, 9, 4, -4, -7, 9},
    	{5, -5, 3, -3, -9, -1},
    	{7, -7, 5, -9},
    	{-3, -4},
    	{-1}
  	}
  },
  {23,
  	{
  		{8, -8, -1, 1, -11, 11, 7, -7, -3, 3, -9, 9, 1, -1, 7, -7, -6, 6, -10, -1},
		{9, -9, 6, -6, -2, 2, -8, 8, -1, 1, 8, -8, -11, 11, 1, -1, -9, 10},
		{-10, 10, 4, -4, 3, -3, 0, 0, 7, -7, 2, -2, 7, -7, 2, -11},
		{5, -5, 9, -9, 1, -1, 0, 0, -8, 8, 2, -2, 10, -3},
		{-10, 10, -6, 6, -7, 7, -2, 2, 2, -2, -9, 7},
		{6, -6, -3, 3, -2, 2, -11, 11, 5, -8},
		{10, -10, -10, 10, 6, -6, -2, -2},
		{-6, 6, 5, -5, 10, -8},
		{10, -10, -2, -5},
		{4, -1},
		{-1}
  	}
  },
  {29,
  	{
  		{2,-2,-6,6,-2,2,14,-14,5,-5,-4,4,-10,10,-4,4,-1,1,0,0,-9,9,-8,8,-13,-1},
		{8,-8,0,0,10,-10,-12,12,4,-4,11,-11,-2,2,3,-3,-6,6,-3,3,-14,14,-12,13},
		{4,-4,9,-9,4,-4,12,-12,-11,11,0,0,-11,11,-14,14,-5,5,7,-7,1,-13},
		{13,-13,0,0,10,-10,-13,13,5,-5,-7,7,-2,2,7,-7,12,-12,-6,13},
		{-9,9,-14,14,3,-3,-7,7,-13,13,-1,1,7,-7,5,-5,-6,-2},
		{9,-9,7,-7,-13,13,-3,3,-8,8,-10,10,12,-12,-2,10},
		{-6,6,8,-8,3,-3,12,-12,-14,14,14,-14,6,-9},
		{-6,6,7,-7,-5,5,-1,1,-10,10,-14,4},
		{3,-3,2,-2,-9,9,10,-10,12,12},
		{-1,1,-6,6,-9,9,-2,-10},
		{-9,9,-12,12,-14,1},
		{-1,1,-4,-13},
		{-5,-6},
		{1}
  	}
  },
  {31,
  	{
  		{-14,14,9,-9,-5,5,4,-4,-9,9,-1,1,-15,15,-11,11,9,-9,0,0,-1,1,13,-13,12,-12,-14,-1},
		{8,-8,5,-5,9,-9,-1,1,-10,10,12,-12,15,-15,10,-10,-2,2,8,-8,9,-9,-10,10,-13,14},
		{-14,14,-15,15,-13,13,-10,10,2,-2,13,-13,0,0,-5,5,0,0,-12,12,7,-7,11,7},
		{-9,9,-6,6,6,-6,-8,8,-11,11,-2,2,-13,13,5,-5,14,-14,-4,4,8,-1},
		{-13,13,12,-12,-6,6,10,-10,-13,13,10,-10,-1,1,8,-8,-11,11,9,12},
		{-8,8,9,-9,-4,4,-9,9,-13,13,2,-2,5,-5,-15,15,-12,-15},
		{-14,14,3,-3,-10,10,-2,2,-4,4,10,-10,9,-9,-12,-6},
		{-4,4,6,-6,-2,2,-7,7,-1,1,9,-9,12,-10},
		{-11,11,10,-10,-9,9,-12,12,8,-8,-7,-11},
		{9,-9,0,0,12,-12,9,-9,3,-6},
		{-1,1,-9,9,-3,3,2,3},
		{-14,14,6,-6,15,-14},
		{-1,1,-12,-11},
		{-5,9},
		{-1}
  	}
  }
};

DoubleCRT Comparator::create_shift_mask(double& size, long shift)
{
	cout << "Mask for shift " << shift << " is being created" << endl;
	// get EncryptedArray
  	const EncryptedArray& ea = m_context.getEA();

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

	size = conv<double>(embeddingLargestCoeff(mask_zzx, m_context.getZMStar()));

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
	cout << "All masks are created" << endl;
}

void Comparator::compute_poly_params()
{
	// get p
	ZZ p = ZZ(m_context.getP());
	long p_long = conv<long>(p);

	// hardcoded babysteps sizes
	map<unsigned long, unsigned long> bs_nums
	{
		{5, 1},
		{7, 2}, // 4
		{11, 3}, // 3 (6), 1..4
  		{13, 3}, // 3 (6), 1..5
  		{17, 4}, // 4 (7), 1..5
  		{19, 3}, // 3 (8), 1..4
  		{23, 5}, // 5 (9), 3..6
  		{29, 5}, // 5 (10), 1..6
  		{31, 5}, // 5 (10), 4..6
  		{37, 5}, // 5 (12)
  		{47, 5}, // 5 (13), 2..11 
  		{61, 6}, // 6 (14), 4..8 
  		{67, 5}, // 5 (15), 4..8
  		{71, 4}, // 4 (15), 3..7
  		{101, 7}, // 7 (16), 4..8
  		{109, 7}, // 7 (19)
  		{131, 8}, // 8 (19), 4..11
  		{167, 10}, // 10 (21), 8..12
  		{173, 10},  // 10 (21), 8..12
  		{271, 9},  // 9 (26), 9..10
  		{401, 12},  // 12 (28), 9..14
  		{659, 11}	// 11 (41), 11..12
  	};

  	m_bs_num_comp = -1;
  	m_bs_num_min = -1;
  	if(bs_nums.count(p_long) > 0)
  	{
  		m_bs_num_comp = bs_nums[p_long];
  		m_bs_num_min = bs_nums[p_long];
  	}

  	// if p > 3, d = (p-3)/2
  	long d_comp = deg(m_univar_less_poly);
  	// if p > 3, d = (p-1)/2
  	long d_min = deg(m_univar_min_max_poly);

  	// How many baby steps: set sqrt(d/2), rounded up/down to a power of two

	// FIXME: There may be some room for optimization here: it may be possible to choose this number as something other than a power of two and still maintain optimal depth, in principle we can try all possible values of m_babystep_num between two consecutive powers of two and choose the one that gives the least number of multiplies, conditioned on minimum depth.

  	if (m_bs_num_comp <= 0) 
	{
		long kk = static_cast<long>(sqrt(d_comp/2.0)); //sqrt(d/2)
		m_bs_num_comp = 1L << NextPowerOfTwo(kk);

    	// heuristic: if #baby_steps >> kk then use a smaler power of two
    	if ((m_bs_num_comp==16 && d_comp>167) || (m_bs_num_comp>16 && m_bs_num_comp>(1.44*kk)))
      		m_bs_num_comp /= 2;
  	}
  	if (m_bs_num_min <= 0) 
	{
		long kk = static_cast<long>(sqrt(d_min/2.0)); //sqrt(d/2)
		m_bs_num_min = 1L << NextPowerOfTwo(kk);

    	// heuristic: if #baby_steps >> kk then use a smaler power of two
    	if ((m_bs_num_min==16 && d_min>167) || (m_bs_num_min>16 && m_bs_num_min>(1.44*kk)))
      		m_bs_num_min /= 2;
  	}

	if(m_verbose)
	{
		cout << "Number of baby steps for comparison: " << m_bs_num_comp << endl;
		cout << "Number of baby steps for min/max: " << m_bs_num_min << endl;
	}

	// #giant_steps = ceil(d/#baby_steps), d >= #giant_steps * #baby_steps
	m_gs_num_comp = divc(d_comp,m_bs_num_comp);
	m_gs_num_min = divc(d_min,m_bs_num_min);

	if(m_verbose)
	{
		cout << "Number of giant steps for comparison: " << m_bs_num_comp << endl;
		cout << "Number of giant steps for min/max: " << m_bs_num_min << endl;
	}      

	// If #giant_steps is not a power of two, ensure that poly is monic and that
	// its degree is divisible by #baby_steps, then call the recursive procedure

	// top coefficient is equal to (p^2 - 1)/8 mod p
	// its inverse is equal to -8 mod p
	m_top_coef_comp = LeadCoeff(m_univar_less_poly);
	m_top_coef_min = LeadCoeff(m_univar_min_max_poly);
	ZZ topInv_comp = ZZ(-8) % p; // the inverse mod p of the top coefficient of poly (if any)
	ZZ topInv_min = ZZ(-8) % p; // the inverse mod p of the top coefficient of poly (if any)
	bool divisible_comp = (m_gs_num_comp * m_bs_num_comp == d_comp); // is the degree divisible by #baby_steps?
	bool divisible_min = (m_gs_num_min * m_bs_num_min == d_min); // is the degree divisible by #baby_steps?

	// FIXME: There may be some room for optimization below: instead of
	// adding a term X^{n*k} we can add X^{n'*k} for some n'>n, so long
	// as n' is smaller than the next power of two. We could save a few
	// multiplications since giantStep[n'] may be easier to compute than
	// giantStep[n] when n' has fewer 1's than n in its binary expansion.

	m_extra_coef_comp = ZZ::zero();    // extra!=0 denotes an added term extra*X^{#giant_steps * #baby_steps}
	m_extra_coef_min = ZZ::zero();    // extra!=0 denotes an added term extra*X^{#giant_steps * #baby_steps}

	if (m_gs_num_comp != (1L << NextPowerOfTwo(m_gs_num_comp)))
	{
		if (!divisible_comp) 
		{  // need to add a term
	    	m_top_coef_comp = NTL::to_ZZ(1);  // new top coefficient is one
	    	topInv_comp = m_top_coef_comp;    // also the new inverse is one
	    	// set extra = 1 - current-coeff-of-X^{n*k}
	    	m_extra_coef_comp = SubMod(m_top_coef_comp, coeff(m_univar_less_poly, m_gs_num_comp * m_bs_num_comp), p);
	    	SetCoeff(m_univar_less_poly, m_gs_num_comp * m_bs_num_comp); // set the top coefficient of X^{n*k} to one
		}

		if (!IsOne(m_top_coef_comp)) 
		{
	    	m_univar_less_poly *= topInv_comp; // Multiply by topInv to make into a monic polynomial
	    	for (long i = 0; i <= m_gs_num_comp * m_bs_num_comp; i++) rem(m_univar_less_poly[i], m_univar_less_poly[i], p);
	    	m_univar_less_poly.normalize();
		}
	}

	/*
	cout << "Less-than poly: ";
	printZZX(cout, m_univar_less_poly, conv<long>(p));
	cout << endl;
	*/

	if (m_gs_num_min != (1L << NextPowerOfTwo(m_gs_num_min)))
	{
		if (!divisible_min) 
		{  // need to add a term
	    	m_top_coef_min = NTL::to_ZZ(1);  // new top coefficient is one
	    	topInv_min = m_top_coef_min;    // also the new inverse is one
	    	// set extra = 1 - current-coeff-of-X^{n*k}
	    	m_extra_coef_min = SubMod(m_top_coef_min, coeff(m_univar_min_max_poly, m_gs_num_min * m_bs_num_min), p);
	    	SetCoeff(m_univar_min_max_poly, m_gs_num_min * m_bs_num_min); // set the top coefficient of X^{n*k} to one
		}

		if (!IsOne(m_top_coef_min)) 
		{
	    	m_univar_min_max_poly *= topInv_min; // Multiply by topInv to make into a monic polynomial
	    	for (long i = 0; i <= m_gs_num_min * m_bs_num_min; i++) rem(m_univar_min_max_poly[i], m_univar_min_max_poly[i], p);
	    	m_univar_min_max_poly.normalize();
		}
	}

	/*
	cout << "Min-max poly: ";
	printZZX(cout, m_univar_min_max_poly, conv<long>(p));
	cout << endl;
	*/

	long top_deg = conv<long>(p-1) >> 1;
	m_baby_index = top_deg % m_bs_num_comp;
	m_giant_index = top_deg / m_bs_num_comp;
	if(m_baby_index == 0)
	{
		m_baby_index = m_bs_num_comp;
		m_giant_index -= 1;
	}
}

void Comparator::create_poly()
{
	cout << "Creating comparison polynomial" << endl;
	// get p
	unsigned long p = m_context.getP();;

	if(m_type == UNI)
	{
		// polynomial coefficient
		ZZ_p coef;
		coef.init(ZZ(p));

		// field element
		ZZ_p field_elem;
		field_elem.init(ZZ(p));

		// initialization of the univariate comparison polynomial
		m_univar_less_poly = ZZX(INIT_MONO, 0, 0);

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

			m_univar_less_poly += ZZX(INIT_MONO, (indx-1) >> 1, rep(coef));
		}

		/*
		cout << "Less-than poly: ";
		printZZX(cout, m_univar_less_poly, p);
		cout << endl;
		*/

		m_univar_min_max_poly = m_univar_less_poly * ZZX(INIT_MONO, 1, 1);

		/*
		cout << "Min-max poly: ";
		printZZX(cout, m_univar_min_max_poly, p);
		cout << endl;
		*/

		compute_poly_params();
	}
	else if (m_type == TAN)
	{
		// computing the coefficients of the bivariate polynomial of Tan et al.
		m_bivar_less_coefs.SetDims(p,p);

		// y^{p-1}
		m_bivar_less_coefs[0][p-1] = ZZ(1);

		// (p+1)/2 * x^{(p-1)/2} * y^{(p-1)/2}
		m_bivar_less_coefs[(p-1) >> 1][(p-1) >> 1] = ZZ((p+1) >> 1);

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
				m_bivar_less_coefs[i][j] = rep(outer_sum);
			}
		}

		cout << "Bivariate coefficients" << endl << m_bivar_less_coefs << endl;

		if (m_verbose)
		{
			cout << "Comparison polynomial: " << endl;
			printZZX(cout, m_univar_less_poly, (p-1)>>1);
			cout << endl;
		}
	}

	cout << "Comparison polynomial is created" << endl;
}

void Comparator::find_prim_root(ZZ_pE& root) const
{
	ZZ qm1 = root.cardinality() - 1;

	cout << "Slot order: " << qm1 << endl;
	cout << "Slot poly: " << root.modulus() << endl;

	vector<ZZ> facts;
	factorize(facts, qm1); // factorization of slot order

	NTL::set(root);

	for (unsigned long i = 0; i < facts.size(); i++) 
	{
		ZZ p = facts[i];
		ZZ pp = p;
		ZZ ee = qm1 / p;
		while (ee % p == 0) 
		{
	  		ee = ee / p;
	  		pp = pp * p;
		}
		// so now we have e = pp * ee, where pp is
		// the power of p that divides e.
		// Our goal is to find an element of order pp

		NTL::PrimeSeq s;
		ZZ_pE q = root;
		ZZ_pE qq = root;
		ZZ_pE qq1 = root;
		long iter = 0;
		do 
		{
	  		iter++;
	  		if (iter > 1000000)
	    		throw RuntimeError("FindPrimitiveRoot: possible infinite loop?");
	  		random(q);
	  		NTL::conv(qq, q);
	  		power(qq1, qq, qm1 / p);
		} 
		while (IsOne(qq1));
		power(qq1, qq, qm1 / pp); // qq1 has order pp

		mul(root, root, qq1);
	}

	// independent check that we have an e-th root of unity
	{
		ZZ_pE s;

		power(s, root, qm1);
		if (!IsOne(s))
	  		throw RuntimeError("FindPrimitiveRoot: internal error (1)");

		// check that s^{e/p} != 1 for any prime divisor p of e
		for (unsigned long i = 0; i < facts.size(); i++) 
		{
	  		ZZ e2 = qm1 / facts[i];
	  		power(s, root, e2); // s = root^{e/p}
	  		if (IsOne(s))
	    		throw RuntimeError("FindPrimitiveRoot: internal error (2)");
		}
	}
}

void Comparator::extraction_init()
{
	// get the total number of slots
	const EncryptedArray& ea = m_context.getEA();
	long nslots = ea.size();

	// get p
	long p = m_context.getP();
	//cout << "p: " << p << endl;

	// get the order of p
	long d = m_context.getOrdP();

	// get the defining polynomial of a slot mod p
	ZZX def_poly = m_context.getSlotRing()->G;

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
	mat_ZZ_pE inv_trace_matrix;
	trace_matrix.SetDims(d, d);
	
	ZZ_pE prim_elem;
	prim_elem.init(def_poly_p);
	prim_elem = conv<ZZ_pE>(ZZ_pX(INIT_MONO, 1, 1));

	//find_prim_root(prim_elem);
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
	//cout << "Trace matrix: " << trace_matrix << endl;
	//cout << "Modulus: " << trace_matrix[0][0].modulus() << endl; 

	inv_trace_matrix = NTL::inv(trace_matrix);
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

			double const_size = conv<double>(embeddingLargestCoeff(tmp, m_context.getZMStar()));
			size_vec.push_back(const_size);

			tmp_crt_vec.push_back(tmp_crt);
		}
		m_extraction_const.push_back(tmp_crt_vec);
		m_extraction_const_size.push_back(size_vec);
	}
}

void Comparator::extract_mod_p(vector<Ctxt>& mod_p_coefs, const Ctxt& ctxt_x) const
{
	HELIB_NTIMER_START(Extraction);
	mod_p_coefs.clear();

	if (m_slotDeg == 1)
	{
		mod_p_coefs.push_back(ctxt_x);
		return;
	}

	const EncryptedArray& ea = m_context.getEA();
	long nslots = ea.size();

	// get max slot degree
	long d = m_context.getOrdP();

	// TODO: how to use key switching hoisting from CRYPTO'18?
	vector<Ctxt> ctxt_frob(d-1, ctxt_x);
	for(long iFrob = 1; iFrob < d; iFrob++)
	{
		ctxt_frob[iFrob-1].frobeniusAutomorph(iFrob);
	} 

	for(long iCoef = 0; iCoef < m_slotDeg; iCoef++)
	{
		//cout << "Extract coefficient " << iCoef << endl;
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
	HELIB_NTIMER_STOP(Extraction);
}

Comparator::Comparator(const Context& context, CircuitType type, unsigned long d, unsigned long expansion_len, const SecKey& sk, bool verbose): m_context(context), m_type(type), m_slotDeg(d), m_expansionLen(expansion_len), m_sk(sk), m_pk(sk), m_verbose(verbose)
{
	//determine the order of p in (Z/mZ)*
	unsigned long ord_p = context.getOrdP();
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

const ZZX& Comparator::get_less_than_poly() const
{
	return m_univar_less_poly;
}

const ZZX& Comparator::get_min_max_poly() const
{
	return m_univar_min_max_poly;
}

void Comparator::print_decrypted(const Ctxt& ctxt) const
{
	// get EncryptedArray
	const EncryptedArray& ea = m_context.getEA();

	// get order of p
	unsigned long ord_p = m_context.getOrdP();

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
	HELIB_NTIMER_START(BatchShift);
	// get EncryptedArray
	const EncryptedArray& ea = m_context.getEA();
	
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
	HELIB_NTIMER_STOP(BatchShift);
}

void Comparator::batch_shift_for_mul(Ctxt& ctxt, long start, long shift) const
{
	HELIB_NTIMER_START(BatchShiftForMul);
	// get EncryptedArray
	const EncryptedArray& ea = m_context.getEA();
	
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

	HELIB_NTIMER_STOP(BatchShiftForMul);
}

void Comparator::shift_and_add(Ctxt& x, long start, long shift_direction) const
{
  HELIB_NTIMER_START(ShiftAdd);
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
  HELIB_NTIMER_STOP(ShiftAdd);
}

void Comparator::shift_and_mul(Ctxt& x, long start, long shift_direction) const
{
  HELIB_NTIMER_START(ShiftMul);
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
  HELIB_NTIMER_STOP(ShiftMul);
}

void Comparator::mapTo01_subfield(Ctxt& ctxt, long pow) const
{
  HELIB_NTIMER_START(MapTo01);	
  // get EncryptedArray
  const EncryptedArray& ea = m_context.getEA();

  // get p
  long p = ctxt.getPtxtSpace();
  if (p != ea.getPAlgebra().getP()) // ptxt space is p^r for r>1
    throw helib::LogicError("mapTo01 not implemented for r>1");

  if (p % pow != 0)
  	throw helib::LogicError("Exponent must divide p");

  if (p > 2)
    ctxt.power((p - 1) / pow); // set y = x^{p-1}

  HELIB_NTIMER_STOP(MapTo01);
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

void Comparator::less_than_mod_any(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	cout << "Compute comparison polynomial" << endl;
	
	Ctxt Y = ctxt_x;
	// x - y
	Y -= ctxt_y;
	// Y = y(x-y)
	Y.multiplyBy(ctxt_y);

	Ctxt x_plus_1 = ctxt_x;
	// x+1
	x_plus_1.addConstant(ZZ(1));

	unsigned long p = m_context.getP();

	unsigned long y_powers = ((p-3) >> 1);
	//powers of x
	DynamicCtxtPowers x_powers(ctxt_x, p-3);
	//powers of Y
	DynamicCtxtPowers Y_powers(Y, y_powers);
	Ctxt Ypow(m_pk);

	Ctxt fx(m_pk);

	vector<ZZX> fpolys(y_powers);
	for (size_t iPoly = 0; iPoly < y_powers; iPoly++)
	{
		for (size_t iCoef = 0; iCoef < fcoefs[p][iPoly].size(); iCoef++)
		{
			SetCoeff(fpolys[iPoly], iCoef+1, fcoefs[p][iPoly][iCoef]);	
		}
		if(iPoly == 0)
		{
			simplePolyEval(ctxt_res, fpolys[iPoly], x_powers);
		}
		else
		{
			simplePolyEval(fx, fpolys[iPoly], x_powers);
			Ypow = Y_powers.getPower(iPoly);
			fx.multiplyBy(Ypow);
			ctxt_res += fx;
		}
	}

	// c*Y^y_powers
	fx = Y_powers.getPower(y_powers);
	fx.multByConstant(ZZ(fcoefs[p][y_powers][0]));
	ctxt_res += fx;
	
	// (x+1)*f(x)
	ctxt_res.multiplyBy(x_plus_1);
	// Y*(x+1)*f(x)
	ctxt_res.multiplyBy(Y);
	
	if(m_verbose)
	{
		print_decrypted(ctxt_res);
		cout << endl;
	}
}

void Comparator::evaluate_univar_less_poly(Ctxt& ret, Ctxt& ctxt_p_1, const Ctxt& x) const
{
	HELIB_NTIMER_START(ComparisonCircuitUnivar);
	// get p
	ZZ p = ZZ(m_context.getP());

	if (p > ZZ(3)) //if p > 3, use the generic Paterson-Stockmeyer strategy
	{
	  // z^2
	  	Ctxt x2 = x;
	  	x2.square();

		DynamicCtxtPowers babyStep(x2, m_bs_num_comp);
		const Ctxt& x2k = babyStep.getPower(m_bs_num_comp);

		DynamicCtxtPowers giantStep(x2k, m_gs_num_comp);

		// Special case when #giant_steps is a power of two
		if (m_gs_num_comp == (1L << NextPowerOfTwo(m_gs_num_comp))) 
		{
			//cout << "I'm computing degPowerOfTwo" << endl;
	    	degPowerOfTwo(ret, m_univar_less_poly, m_bs_num_comp, babyStep, giantStep);
	    }
	    else
	    {
		  	recursivePolyEval(ret, m_univar_less_poly, m_bs_num_comp, babyStep, giantStep);

		  	if (!IsOne(m_top_coef_comp)) 
		  	{
		    	ret.multByConstant(m_top_coef_comp);
			}

			if (!IsZero(m_extra_coef_comp)) 
			{ // if we added a term, now is the time to subtract back
		    	Ctxt topTerm = giantStep.getPower(m_gs_num_comp);
		    	topTerm.multByConstant(m_extra_coef_comp);
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

		/*
		cout << "Computed baby steps" << endl;
		for(int i = 0; i < babyStep.size(); i++)
		{
			cout << i + 1 << ' ' << babyStep.isPowerComputed(i+1) << endl; 
		}

		cout << "Computed giant steps" << endl;
		for(int i = 0; i < giantStep.size(); i++)
		{
			cout << i + 1 << ' ' << giantStep.isPowerComputed(i+1) << endl;
		}
		*/
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
	HELIB_NTIMER_STOP(ComparisonCircuitUnivar);
}

void Comparator::evaluate_min_max_poly(Ctxt& ctxt_min, Ctxt& ctxt_max, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	HELIB_NTIMER_START(MinMaxCircuitUnivar);
	// get p
	ZZ p = ZZ(m_context.getP());

	// Subtraction z = x - y
	cout << "Subtraction" << endl;
	Ctxt ctxt_z = ctxt_x;
	ctxt_z -= ctxt_y;

	if(m_verbose)
	{
		print_decrypted(ctxt_z);
		cout << endl;
	}
  	
  	// z^2
  	Ctxt ctxt_z2 = ctxt_z;
  	ctxt_z2.square(); 

	if (p > ZZ(3)) //if p > 3, use the generic Paterson-Stockmeyer strategy
	{
		DynamicCtxtPowers babyStep(ctxt_z2, m_bs_num_min);
		const Ctxt& ctxt_z2k = babyStep.getPower(m_bs_num_min);

		DynamicCtxtPowers giantStep(ctxt_z2k, m_gs_num_min);

		// compute g(z^2)
		Ctxt g_z2 = Ctxt(ctxt_z2.getPubKey());;
		// Special case when #giant_steps is a power of two
		if (m_gs_num_min == (1L << NextPowerOfTwo(m_gs_num_min))) 
		{
			//cout << "I'm computing degPowerOfTwo" << endl;
	    	degPowerOfTwo(g_z2, m_univar_min_max_poly, m_bs_num_min, babyStep, giantStep);
	    }
	    else
	    {
		  	recursivePolyEval(g_z2, m_univar_min_max_poly, m_bs_num_min, babyStep, giantStep);

		  	if (!IsOne(m_top_coef_min)) 
		  	{
		    	g_z2.multByConstant(m_top_coef_min);
			}

			if (!IsZero(m_extra_coef_min)) 
			{ // if we added a term, now is the time to subtract back
		    	Ctxt topTerm = giantStep.getPower(m_gs_num_min);
		    	topTerm.multByConstant(m_extra_coef_min);
		    	g_z2 -= topTerm;
			}
		}

		// last term: ((p+1)/2) * (x + y) 
		Ctxt last_term = ctxt_x;
		last_term += ctxt_y; 
		last_term.multByConstant(ZZ((p+1)>> 1));

		ctxt_min = last_term;
		ctxt_min += g_z2;
		ctxt_max = last_term;
		ctxt_max -= g_z2;

		/*
		cout << "Computed baby steps" << endl;
		for(int i = 0; i < babyStep.size(); i++)
		{
			cout << i + 1 << ' ' << babyStep.isPowerComputed(i+1) << endl; 
		}

		cout << "Computed giant steps" << endl;
		for(int i = 0; i < giantStep.size(); i++)
		{
			cout << i + 1 << ' ' << giantStep.isPowerComputed(i+1) << endl;
		}
		*/
	}
	else //circuit for p=3
	{
		// last term: ((p+1)/2) * (x + y) 
		Ctxt last_term = ctxt_x;
		last_term += ctxt_y; 
		last_term.multByConstant(ZZ((p+1)>> 1));

		ctxt_min = last_term;
		ctxt_min += ctxt_z2;
		ctxt_max = last_term;
		ctxt_max -= ctxt_z2;
	}
	HELIB_NTIMER_STOP(MinMaxCircuitUnivar);
}

void Comparator::less_than_bivar(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
  HELIB_NTIMER_START(ComparisonCircuitBivar);

  //compare with the circuit of Tan et al.
  if (m_type == TAN)
  {
  	less_than_bivar_tan(ctxt_res, ctxt_x, ctxt_y);
  	return;
  }

  unsigned long p = m_context.getP();

  if(p > 31)
  {
  	throw helib::LogicError("Bivariate circuit is not implemented for p > 31");
  }

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

  if(p > 7)
  {
    less_than_mod_any(ctxt_res, ctxt_x, ctxt_y);
  }

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  HELIB_NTIMER_STOP(ComparisonCircuitBivar);
}

void Comparator::less_than_bivar_tan(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	cout << "Compute Tan's comparison polynomial" << endl;

	long p = m_context.getP();

	DynamicCtxtPowers x_powers(ctxt_x, p-1);
	DynamicCtxtPowers y_powers(ctxt_y, p-1);

	ctxt_res = y_powers.getPower(p-1);

	for (long i = 1; i < p; i++)
	{
		// zero ciphertext
		Ctxt sum = Ctxt(ctxt_x.getPubKey());
		for (long j = 1; j < p; j++)
		{
			if (m_bivar_less_coefs[i][j] == ZZ(0))
				continue;
			Ctxt tmp = y_powers.getPower(j);
			tmp.multByConstant(m_bivar_less_coefs[i][j]);
			sum += tmp;
		}
		sum.multiplyBy(x_powers.getPower(i));
		ctxt_res += sum;
	}
}

void Comparator::is_zero(Ctxt& ctxt_res, const Ctxt& ctxt_z, long pow) const
{
  HELIB_NTIMER_START(EqualityCircuit);

  ctxt_res = ctxt_z;

  //compute mapTo01: (z_i)^{p^d-1}
  //cout << "Mapping to 0 and 1" << endl;
  mapTo01_subfield(ctxt_res, pow);

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  //cout << "Computing NOT" << endl;
  //compute 1 - mapTo01(z_i)
  ctxt_res.negate();
  ctxt_res.addConstant(ZZ(1));

  if(m_verbose)
  {
    print_decrypted(ctxt_res);
    cout << endl;
  }

  HELIB_NTIMER_STOP(EqualityCircuit);
}

void Comparator::compare(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	HELIB_NTIMER_START(Comparison);

	vector<Ctxt> ctxt_less_p;
	vector<Ctxt> ctxt_eq_p;

	// bivariate circuit
	if (m_type == BI || m_type == TAN)
	{
		//cout << "Extraction" << endl;
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

		//cout << "Compute the less-than function modulo p" << endl;
		for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
		{
			Ctxt ctxt_tmp = Ctxt(ctxt_x.getPubKey());
			less_than_bivar(ctxt_tmp, ctxt_x_p[iCoef], ctxt_y_p[iCoef]);
			ctxt_less_p.push_back(ctxt_tmp);
		}

		//cout << "Compute the equality function modulo p" << endl;
		for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
		{
			// Subtraction z = x - y
			//cout << "Subtraction" << endl;
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
		//cout << "Subtraction" << endl;
		Ctxt ctxt_z = ctxt_x;
		ctxt_z -= ctxt_y;

		if(m_verbose)
		{
			print_decrypted(ctxt_z);
			cout << endl;
		}

		// extract mod p coefficients
		//cout << "Extraction" << endl;
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

		//cout << "Compute the less-than and equality functions modulo p" << endl;
		for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
		{
			Ctxt ctxt_tmp = Ctxt(ctxt_z.getPubKey());
			Ctxt ctxt_tmp_eq = Ctxt(ctxt_z.getPubKey());

			// compute polynomial function for 'z < 0'
			//cout << "Compute univariate comparison polynomial" << endl;
			evaluate_univar_less_poly(ctxt_tmp, ctxt_tmp_eq, ctxt_z_p[iCoef]);

			if(m_verbose)
			{
			  cout << "Result of the less-than function" << endl;
			  print_decrypted(ctxt_tmp);
			  cout << endl;
			}

			ctxt_less_p.push_back(ctxt_tmp);

			//cout << "Computing NOT" << endl;
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

	//cout << "Compare digits" << endl;
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
	//cout << "Rotating and multiplying slots with equalities" << endl;
	shift_and_mul(ctxt_eq, 0);

	if(m_verbose)
	{
		print_decrypted(ctxt_eq);
		cout << endl;
	}

	//Remove the least significant digit and shift to the left
	//cout << "Remove the least significant digit" << endl;
	batch_shift_for_mul(ctxt_eq, 0, -1);

	if(m_verbose)
	{
		print_decrypted(ctxt_eq);
		cout << endl;
	}

	//cout << "Final result" << endl;

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

    HELIB_NTIMER_STOP(Comparison);
}

void Comparator::min_max_digit(Ctxt& ctxt_min, Ctxt& ctxt_max, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	HELIB_NTIMER_START(MinMaxDigit);
	if(m_type != UNI)
		throw helib::LogicError("Min/Max is not implemented with the bivariate circuit");

	if(m_expansionLen != 1 || m_slotDeg != 1)
		throw helib::LogicError("Min/Max is not implemented for vectors over F_p");

	// get EncryptedArray
  	const EncryptedArray& ea = m_context.getEA();
  	//extract slots
	long nSlots = ea.size();

	vector<Ctxt> ctxt_min_p;
	vector<Ctxt> ctxt_max_p;

	// extract mod p coefficients
	cout << "Extraction" << endl;
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

	cout << "Compute min/max functions modulo p" << endl;
	for (long iCoef = 0; iCoef < m_slotDeg; iCoef++)
	{
		Ctxt ctxt_tmp_min = Ctxt(ctxt_x.getPubKey());
		Ctxt ctxt_tmp_max = Ctxt(ctxt_x.getPubKey());

		// compute polynomial function for 'z < 0'
		cout << "Compute univariate min/max polynomial" << endl;
		evaluate_min_max_poly(ctxt_tmp_min, ctxt_tmp_max, ctxt_x_p[iCoef], ctxt_y_p[iCoef]);

		if(m_verbose)
		{
		  cout << "Result of the min function" << endl;
		  print_decrypted(ctxt_tmp_min);
		  cout << endl;
		}

		if(m_verbose)
		{
		  cout << "Result of the max function" << endl;
		  print_decrypted(ctxt_tmp_max);
		  cout << endl;
		}

		ctxt_min_p.push_back(ctxt_tmp_min);
		ctxt_max_p.push_back(ctxt_tmp_max);
	}

	ctxt_min = ctxt_min_p[0];
	ctxt_max = ctxt_max_p[0];

	for (long iCoef = 1; iCoef < m_slotDeg; iCoef++)
	{
		vector<ZZX> x_power(nSlots, ZZX(INIT_MONO, iCoef, 1));
		ZZX x_power_ptxt;
		ea.encode(x_power_ptxt, x_power);

		// agregate minimum values
		Ctxt tmp = ctxt_min_p[iCoef];
		tmp.multByConstant(x_power_ptxt);
		ctxt_min += tmp;

		// agregate maximum values
		tmp = ctxt_max_p[iCoef];
		tmp.multByConstant(x_power_ptxt);
		ctxt_max += tmp;
	}

	HELIB_NTIMER_STOP(MinMaxDigit);
}

void Comparator::min_max(Ctxt& ctxt_min, Ctxt& ctxt_max, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const
{
	HELIB_NTIMER_START(MinMax);
	if(m_type == UNI && m_expansionLen == 1 && m_slotDeg == 1)
	{
		min_max_digit(ctxt_min, ctxt_max, ctxt_x, ctxt_y);
		return;
	}

	Ctxt ctxt_z = ctxt_x;
	ctxt_z -= ctxt_y;

	Ctxt ctxt_tmp = Ctxt(ctxt_z.getPubKey());
	compare(ctxt_tmp, ctxt_x, ctxt_y);
	ctxt_tmp.multiplyBy(ctxt_z);

	ctxt_min = ctxt_y;
	ctxt_min += ctxt_tmp;

	ctxt_max = ctxt_x;
	ctxt_max -= ctxt_tmp;

	if(m_verbose)
	{
		cout << "Minimum" << endl;
		print_decrypted(ctxt_min);
		cout << endl;
	}

	if(m_verbose)
	{
		cout << "Maximum" << endl;
		print_decrypted(ctxt_max);
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
	HELIB_NTIMER_STOP(MinMax);
}

void Comparator::array_min(Ctxt& ctxt_res, const vector<Ctxt>& ctxt_in, long depth) const
{
	HELIB_NTIMER_START(ArrayMin);

	if (depth < 0)
		throw helib::LogicError("depth parameter must be non-negative");

	cout << "Computing the minimum of an array" << endl;

	size_t input_len = ctxt_in.size();

	vector<Ctxt> ctxt_res_vec;
	for (size_t i  = 0; i < input_len; i++)
	{
		ctxt_res_vec.push_back(ctxt_in[i]);
	}

	size_t cur_len = input_len;
	long level = depth;

	while(cur_len > 1 && level > 0)
	{
		cout << "Comparison level: " << depth-level << endl;
		// compare x[i] and x[n-1-i] where n is the length of ctxt_res_vec
		for (size_t i  = 0; i < (cur_len >> 1); i++)
		{
			if(i != cur_len -  1 - i)
			{
				cout << "Comparing ciphertexts " << i << " and " << cur_len -  1 - i << endl;
				min_max(ctxt_res_vec[i], ctxt_res_vec[cur_len -  1 - i], ctxt_res_vec[i], ctxt_res_vec[cur_len -  1 - i]);
			}
		}
		cur_len = (cur_len >> 1) + (cur_len % 2);
		ctxt_res_vec.resize(cur_len, Ctxt(m_pk));
		level--;
	}

	if(cur_len > 1)
	{
		// plaintext modulus
  		long p = m_context.getP();
		// multiplications in the equality circuit
		long eq_mul_num = static_cast<long>(floor(log2(p-1))) + weight(ZZ(p-1)) - 1;
		long eq_depth = static_cast<long>(ceil(log2(p-1)));
		long prod_depth = static_cast<long>(ceil(log2(cur_len-1)));

		if ((((cur_len - 2 > eq_mul_num) && (eq_depth == prod_depth)) || (eq_depth < prod_depth)) && cur_len <= p)
		{
			cout << "Computing minimum via equality" << endl;
			cout << "Mult. of equality: " << eq_mul_num << endl;
			cout << "Depth of equality: " << eq_depth << endl;
			cout << "Depth of product: " << prod_depth << endl;
			// create a table with all pairwise comparisons and compute the Hamming weight of every row
			vector<Ctxt> ham_weights;
			get_sorting_index(ham_weights, ctxt_res_vec);

			cout << "Computing the minimum" << endl;
			ctxt_res = Ctxt(m_pk);
			for(size_t i = 0; i < ctxt_res_vec.size(); i++)
			{
				//compare the Hamming weight of the jth row with i
				Ctxt tmp_prod = ham_weights[i];
				tmp_prod.addConstant(ZZX(-(cur_len-1)));
				mapTo01_subfield(tmp_prod, 1);
				tmp_prod.negate();
				tmp_prod.addConstant(ZZX(1));

				//multiply by the jth input ciphertext
				tmp_prod.multiplyBy(ctxt_res_vec[i]);
				if(i == 0)
					ctxt_res = tmp_prod;
				else
					ctxt_res += tmp_prod;
			}
		}
		else 
		{
			cout << "Computing minimum via punctured products" << endl;
			cout << "Mult. of equality: " << eq_mul_num << endl;
			cout << "Depth of equality: " << eq_depth << endl;
			cout << "Depth of product: " << prod_depth << endl;
			vector<vector<Ctxt>> ctxt_products;
			// compute the product of every row
			for(size_t i = 0; i < ctxt_res_vec.size(); i++)
			{
				vector<Ctxt> ctxt_vec;
				ctxt_products.push_back(ctxt_vec);
			}

			cout << "Computing the comparison table" << endl;
			for (size_t i = 0; i < ctxt_res_vec.size() - 1; i++)
			{
				cout << "Computing Row " << i << endl;
				for(size_t j = i + 1; j < ctxt_res_vec.size(); j++)
				{
					cout << "Computing Column " << j << endl;
					// compute upper diagonal entries of the comparison table and multiply them
					Ctxt comp_col = Ctxt(m_pk);
					compare(comp_col, ctxt_res_vec[i], ctxt_res_vec[j]);

					if (ctxt_products[i].empty())
					{
						ctxt_products[i].push_back(comp_col);
					}
					else
					{
						long wt = weight(ZZ(j));
						int len_i = ctxt_products[i].size();
						if (wt > len_i)
						{
							ctxt_products[i].push_back(comp_col);
						}
						else
						{
							ctxt_products[i][len_i - 1].multiplyBy(comp_col);
							for (int k = len_i - 2; k >= (wt-1); k--)
							{
								ctxt_products[i][k].multiplyBy(ctxt_products[i][k+1]);
								ctxt_products[i].pop_back();	
							}
						}
					}
					
					// compute lower diagonal entries of the comparison table by transposition and logical negation of upper diagonal entries
					//NOT the result to multiply to the jth row
					comp_col.negate();
					comp_col.addConstant(ZZ(1));

					if (ctxt_products[j].empty())
					{
						ctxt_products[j].push_back(comp_col);
					}
					else
					{
						long wt = weight(ZZ(i+1));
						int len_j = ctxt_products[j].size();
						if (wt > len_j)
						{
							ctxt_products[j].push_back(comp_col);
						}
						else
						{
							ctxt_products[j][len_j - 1].multiplyBy(comp_col);
							for (int k = len_j - 2; k >= (wt-1); k--)
							{
								ctxt_products[j][k].multiplyBy(ctxt_products[j][k+1]);
								ctxt_products[j].pop_back();	
							}
						}
					}
				}
			}

			cout << "Computing the minimum" << endl;
			ctxt_res = Ctxt(m_pk);
			for(size_t i = 0; i < ctxt_res_vec.size(); i++)
			{
				int len_i = ctxt_products[i].size();
				for (int k = len_i - 2; k >= 0; k--)
				{
					ctxt_products[i][k].multiplyBy(ctxt_products[i][k+1]);
					ctxt_products[i].pop_back();	
				}
				Ctxt tmp_prod = ctxt_products[i][0];

				//multiply by the ith input ciphertext
				tmp_prod.multiplyBy(ctxt_res_vec[i]);

				//add to the result
				ctxt_res += tmp_prod;
			}
		}
	}
	else
	{
		ctxt_res = ctxt_res_vec[0];
	}

	HELIB_NTIMER_STOP(ArrayMin);
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

void Comparator::get_sorting_index(vector<Ctxt>& ctxt_out, const vector<Ctxt>& ctxt_in) const
{
	ctxt_out.clear();

	// length of the input vector
	size_t input_len = ctxt_in.size();

	// plaintext modulus
  	long p = m_context.getP();

	if (input_len > p)
		throw helib::LogicError("The number of ciphertexts cannot be larger than the plaintext modulus");

	// compute the Hamming weight of every row
	for(size_t i = 0; i < input_len; i++)
	{
		//initialize Hamming weights to zero
		Ctxt ctxt_tmp = Ctxt(ctxt_in[0].getPubKey());
		ctxt_out.push_back(ctxt_tmp);
	}

	cout << "Computing the comparison table" << endl;
	for (size_t i = 0; i < input_len - 1; i++)
	{
		cout << "Computing Row " << i << endl;
		for(size_t j = i + 1; j < input_len; j++)
		{
			cout << "Computing Column " << j << endl;
			// compute upper diagonal entries of the comparison table and sum them
			Ctxt comp_col_j = Ctxt(ctxt_in[0].getPubKey());
			compare(comp_col_j, ctxt_in[i], ctxt_in[j]);
			ctxt_out[i] += comp_col_j;

			// compute lower diagonal entries of the comparison table by transposition and logical negation of upper diagonal entries
			//NOT the result to add to the jth row
			comp_col_j.negate();
			comp_col_j.addConstant(ZZ(1));

			// add lower diagonal entries to Hamming weight accumulators of related rows
			ctxt_out[j] += comp_col_j;
		}
	}
}

void Comparator::sort(vector<Ctxt>& ctxt_out, const vector<Ctxt>& ctxt_in) const
{
	HELIB_NTIMER_START(Sorting);

	ctxt_out.clear();

	// length of the input vector
	size_t input_len = ctxt_in.size();

	// plaintext modulus
  	long p = m_context.getP();

	if (input_len > p)
		throw helib::LogicError("The number of ciphertexts cannot be larger than the plaintext modulus");

	// multiplications in the equality circuit
	long eq_mul_num = static_cast<long>(floor(log2(p-1))) + weight(ZZ(p-1)) - 1;
	cout << "Multiplications in the equality circuit: " << eq_mul_num << endl;

	// create a table with all pairwise comparisons and compute the Hamming weight of every row
	vector<Ctxt> ham_weights;

	get_sorting_index(ham_weights, ctxt_in);
	/*
	for(size_t i = 0; i < input_len; i++)
	{
		//initialize Hamming weights to zero
		Ctxt ctxt_tmp = Ctxt(ctxt_in[0].getPubKey());
		ham_weights.push_back(ctxt_tmp);
	}

	cout << "Computing the comparison table" << endl;
	for (size_t i = 0; i < input_len - 1; i++)
	{
		cout << "Computing Row " << i << endl;
		for(size_t j = i + 1; j < input_len; j++)
		{
			cout << "Computing Column " << j << endl;
			// compute upper diagonal entries of the comparison table and sum them
			Ctxt comp_col_j = Ctxt(ctxt_in[0].getPubKey());
			compare(comp_col_j, ctxt_in[i], ctxt_in[j]);
			ham_weights[i] += comp_col_j;

			// compute lower diagonal entries of the comparison table by transposition and logical negation of upper diagonal entries
			//NOT the result to add to the jth row
			comp_col_j.negate();
			comp_col_j.addConstant(ZZ(1));

			// add lower diagonal entries to Hamming weight accumulators of related rows
			ham_weights[j] += comp_col_j;
		}
	}
	*/

	// print Hamming weights
	/*
	for(size_t i = 0; i < input_len; i++)
	{
		cout << i << " Row:" << endl;
		print_decrypted(ham_weights[i]);
      	cout << endl;
	}
	*/

	if(eq_mul_num * input_len <= p - 2)
	//if(true)
	{
		for(size_t i = 0; i < input_len; i++)
		{
			cout << "Computing Element " << i << endl;
			Ctxt tmp_sum = Ctxt(ctxt_in[i].getPubKey());
			for(size_t j = 0; j < input_len; j++)
			{
				//compare the Hamming weight of the jth row with i
				Ctxt tmp_prod = ham_weights[j];
				tmp_prod.addConstant(ZZX(-i));
				mapTo01_subfield(tmp_prod, 1);
				tmp_prod.negate();
				tmp_prod.addConstant(ZZX(1));

				//multiply by the jth input ciphertext
				tmp_prod.multiplyBy(ctxt_in[j]);
				tmp_sum += tmp_prod;
			}
			ctxt_out.push_back(tmp_sum);
		}
	}
	else
	{
		// equality sums
		vector<Ctxt> eq_sums;

		// fill ctxt_out with zeros
		for(size_t i = 0; i < input_len; i++)
		{
			ctxt_out.push_back(Ctxt(ctxt_in[0].getPubKey()));
			eq_sums.push_back(Ctxt(ctxt_in[0].getPubKey()));
		}

		for (size_t i = 0; i < input_len; i++)
		{
			cout << "Adding element " << i << endl;

			// hw_i^j, j in [1,p-1]
			DynamicCtxtPowers hw_powers(ham_weights[i], p-1);
			
			// eq_sums[0] = 1 - hw_i^(p-1)
			eq_sums[0].clear();
			eq_sums[0] = hw_powers.getPower(p-1);
			eq_sums[0].negate();
			eq_sums[0].addConstant(ZZX(1));

			// eq_sums[0] * ctxt_in[i]
			eq_sums[0].multiplyBy(ctxt_in[i]);

			// sum_i eq_sums[0] * ctxt_in[i]
			ctxt_out[0] += eq_sums[0];	

			for (size_t k = 1; k < input_len; k++)
			{
				// zeroize
				eq_sums[k].clear();

				// current sorting index
				ZZ_p k_zzp;
				k_zzp.init(ZZ(p));
				k_zzp = k;

				ZZ_p k_power;
				k_power.init(ZZ(p));

				for (int j = 1; j < p; j++)
				{
					// k^(p-1-j) mod p
					k_power = power(k_zzp, p - 1 - j);
					// hw_i^j
					Ctxt tmp = hw_powers.getPower(j);
					// hw_i^j * k^(p-1-j)
					tmp.multByConstant(rep(k_power));
					// sum hw_i^j * k^(p-1-j)
					eq_sums[k] += tmp;
				}
				// add k^(p-1) to eq_sums
				eq_sums[k].addConstant(rep(power(k_zzp, p-1)));

				// 1 - sum_(j=0)^(p-1) hw_i^j * k^(p-1-j)
				eq_sums[k].negate();
				eq_sums[k].addConstant(ZZX(1));

				// eq_sums[k] * ctxt_in[i]
				eq_sums[k].multiplyBy(ctxt_in[i]);

				// sum_i eq_sums[k] * ctxt_in[i]
				ctxt_out[k] += eq_sums[k];
			}
		}
	}
	

	// print output ciphertexts
	/*
	for(size_t i = 0; i < input_len; i++)
	{
		cout << i << " ctxt:" << endl;
		print_decrypted(ctxt_out[i]);
      	cout << endl;
	}
	*/
	HELIB_NTIMER_STOP(Sorting);
}

void Comparator::test_sorting(int num_to_sort, long runs) const
{
	//reset timers
  setTimersOn();
  
  // initialize the random generator
  random_device rd;
  mt19937 eng(rd());
  uniform_int_distribution<unsigned long> distr_u;
  uniform_int_distribution<long> distr_i;

  // get EncryptedArray
  const EncryptedArray& ea = m_context.getEA();

  //extract number of slots
  long nslots = ea.size();

  //get p
  unsigned long p = m_context.getP();

  //order of p
  unsigned long ord_p = m_context.getOrdP();

  //amount of numbers in one ciphertext
  unsigned long numbers_size = nslots / m_expansionLen;

  // number of slots occupied by encoded numbers
  unsigned long occupied_slots = numbers_size * m_expansionLen;

  //encoding base, ((p+1)/2)^d
  //if 2-variable comparison polynomial is used, it must be p^d
  unsigned long enc_base = (p + 1) >> 1;
  if (m_type == BI || m_type == TAN)
  {
  	enc_base = p;
  }

  unsigned long digit_base = power_long(enc_base, m_slotDeg);

  //check that field_size^expansion_len fits into 64-bits
  int space_bit_size = static_cast<int>(ceil(m_expansionLen * log2(digit_base)));
  unsigned long input_range = ULONG_MAX;
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

    // all slots contain the same value
    vector<vector<ZZX>> expected_result;
    for (int i = 0; i < num_to_sort; i++)
    {
    	vector<ZZX> tmp_vec(occupied_slots,ZZX(INIT_MONO,0,0));
    	expected_result.push_back(tmp_vec);
    }
    
    // vector of input longs
    vector<vector<unsigned long>> input_xs;
    for (int i = 0; i < numbers_size; i++)
    {
    	vector<unsigned long> tmp_vec(num_to_sort,0);
    	input_xs.push_back(tmp_vec);
    }

    ZZX pol_slot;

    //ciphertexts to sort
    vector<Ctxt> ctxt_in;

    //sorted ciphertexts
    vector<Ctxt> ctxt_out;

    for (int i = 0; i < num_to_sort; i++)
    {
		// the plaintext polynomials
		vector<ZZX> pol_x(nslots);

		//encoding of slots
		for (int k = 0; k < numbers_size; k++)
		{
			unsigned long input_x = distr_u(eng) % input_range;

			input_xs[k][i] = input_x;
		
			if(m_verbose)
			{
				cout << "Input" << endl;
				cout << input_x << endl;
			}

			vector<long> decomp_int_x;

			//decomposition of input integers
			digit_decomp(decomp_int_x, input_x, digit_base, m_expansionLen);
			for (int j = 0; j < m_expansionLen; j++)
			{
			    //decomposition of a digit
			    int_to_slot(pol_slot, decomp_int_x[j], enc_base);
			    pol_x[k * m_expansionLen + j] = pol_slot;
			}
		}

      	if(m_verbose)
	    {
	      cout << "Input" << endl;
	      for(int j = 0; j < nslots; j++)
	      {
	          printZZX(cout, pol_x[j], ord_p);
	          cout << endl;
	      }
	    }

	    Ctxt ctxt_x(m_pk);
    	ea.encrypt(ctxt_x, m_pk, pol_x);

    	ctxt_in.push_back(ctxt_x);
    }

    //cout << "Input" << endl;
    for (int i = 0; i < numbers_size; i++)
    {
    	//for (int j = 0; j < num_to_sort; j++)
    	//	cout << input_xs[i][j] << " ";
    	//cout << endl;
    	std::sort(input_xs[i].begin(), input_xs[i].end());
    }

    //cout << "Expected results" << endl;
    for (int i = num_to_sort-1; i >= 0; i--)
    {
		for (int k = 0; k < numbers_size; k++)
		{
			vector<long> decomp_int_x;

			//decomposition of input integers
			digit_decomp(decomp_int_x, input_xs[k][i], digit_base, m_expansionLen);
			for (int j = 0; j < m_expansionLen; j++)
			{
			    //decomposition of a digit
			    int_to_slot(pol_slot, decomp_int_x[j], enc_base);
			    expected_result[num_to_sort-1-i][k * m_expansionLen + j] = pol_slot;
			}
		}

		/*
		cout << num_to_sort-1-i << endl;
		for(int j = 0; j < nslots; j++)
		{
			printZZX(cout, expected_result[num_to_sort-1-i][j], ord_p);
        	cout << endl;
		}
		*/	
    }

    // comparison function
    cout << "Start of sorting" << endl;
    this->sort(ctxt_out, ctxt_in);

    printNamedTimer(cout, "Extraction");
    printNamedTimer(cout, "ComparisonCircuitBivar");
    printNamedTimer(cout, "ComparisonCircuitUnivar");
    printNamedTimer(cout, "EqualityCircuit");
    printNamedTimer(cout, "ShiftMul");
    printNamedTimer(cout, "ShiftAdd");
    printNamedTimer(cout, "Comparison");
    printNamedTimer(cout, "Sorting");

    const FHEtimer* sort_timer = getTimerByName("Sorting");

    cout << "Avg. time per batch: " << 1000.0 * sort_timer->getTime()/static_cast<double>(run+1)/static_cast<double>(numbers_size) << " ms" << endl;
    cout << "Number of integers in one ciphertext "<< numbers_size << endl;

    // remove the line below if it gives bizarre results 
    ctxt_out[0].cleanUp();
    capacity = ctxt_out[0].bitCapacity();
    cout << "Final capacity: " << capacity << endl;
    if (capacity < min_capacity)
      min_capacity = capacity;
    cout << "Min. capacity: " << min_capacity << endl;
    cout << "Final size: " << ctxt_out[0].logOfPrimeSet()/log(2.0) << endl;
    
    for (int i = 0; i < num_to_sort; i++)
    {
    	vector<ZZX> decrypted(nslots);
	    ea.decrypt(ctxt_out[i], m_sk, decrypted);

	    for(int j = 0; j < numbers_size; j++)
	    { 
	    	for(int k = 0; k < m_expansionLen; k++)
	    	{
	    		if (decrypted[j * m_expansionLen + k] != expected_result[i][j * m_expansionLen + k])
			    {
			    	printf("Slot %ld: ", j * m_expansionLen + k);
			    	printZZX(cout, decrypted[j * m_expansionLen + k], ord_p);
			        cout << endl;
			        cout << "Failure" << endl;
			        return;
			    }
	    	}
	    }
    }
  }
}

void Comparator::test_compare(long runs) const
{
  //reset timers
  setTimersOn();
  
  // initialize the random generator
  random_device rd;
  mt19937 eng(rd());
  uniform_int_distribution<unsigned long> distr_u;
  uniform_int_distribution<long> distr_i;

  // get EncryptedArray
  const EncryptedArray& ea = m_context.getEA();

  //extract number of slots
  long nslots = ea.size();

  //get p
  unsigned long p = m_context.getP();

  //order of p
  unsigned long ord_p = m_context.getOrdP();

  //amount of numbers in one ciphertext
  unsigned long numbers_size = nslots / m_expansionLen;

  // number of slots occupied by encoded numbers
  unsigned long occupied_slots = numbers_size * m_expansionLen;

  //encoding base, ((p+1)/2)^d
  //if 2-variable comparison polynomial is used, it must be p^d
  unsigned long enc_base = (p + 1) >> 1;
  if (m_type == BI || m_type == TAN)
  {
  	enc_base = p;
  }

  unsigned long digit_base = power_long(enc_base, m_slotDeg);

  //check that field_size^expansion_len fits into 64-bits
  int space_bit_size = static_cast<int>(ceil(m_expansionLen * log2(digit_base)));
  cout << "Space bit size " << space_bit_size << endl;
  unsigned long input_range = ULONG_MAX;
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
        cout << "Input " << i << endl;
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

      if(m_verbose)
      {
      	cout << "Input decomposition into digits" << endl;
      	for(int j = 0; j < m_expansionLen; j++)
      	{
      		cout << decomp_int_x[j] << " " << decomp_int_y[j] << endl;
      	}
      }

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
      for(int j = 0; j < nslots; j++)
      {
          printZZX(cout, pol_x[j], ord_p);
          printZZX(cout, pol_y[j], ord_p);
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

    cout << "Avg. time per integer: " << 1000.0 * comp_timer->getTime()/static_cast<double>(run+1)/static_cast<double>(numbers_size) << " ms" << endl;
    cout << "Number of integers in one ciphertext "<< numbers_size << endl;

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

void Comparator::test_min_max(long runs) const
{
	//reset timers
  setTimersOn();
  
  // initialize the random generator
  random_device rd;
  mt19937 eng(rd());
  uniform_int_distribution<unsigned long> distr_u;
  uniform_int_distribution<long> distr_i;

  // get EncryptedArray
  const EncryptedArray& ea = m_context.getEA();

  //extract number of slots
  long nslots = ea.size();

  //get p
  unsigned long p = m_context.getP();

  //order of p
  unsigned long ord_p = m_context.getOrdP();

  //amount of numbers in one ciphertext
  unsigned long numbers_size = nslots / m_expansionLen;

  // number of slots occupied by encoded numbers
  unsigned long occupied_slots = numbers_size * m_expansionLen;

  //encoding base, ((p+1)/2)^d
  //if 2-variable comparison polynomial is used, it must be p^d
  unsigned long enc_base = (p + 1) >> 1;
  if (m_type == BI || m_type == TAN)
  {
  	enc_base = p;
  }

  unsigned long digit_base = power_long(enc_base, m_slotDeg);

  //check that field_size^expansion_len fits into 64-bits
  int space_bit_size = static_cast<int>(ceil(m_expansionLen * log2(digit_base)));
  unsigned long input_range = ULONG_MAX;
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

    vector<ZZX> expected_result_min(occupied_slots);
    vector<ZZX> expected_result_max(occupied_slots);
    vector<ZZX> decrypted_min(occupied_slots);
    vector<ZZX> decrypted_max(occupied_slots);

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

      if (input_x < input_y)
      {
      	for (int j = 0; j < m_expansionLen; j++)
      	{
	        expected_result_min[i * m_expansionLen + j] = pol_x[i * m_expansionLen + j];
	        expected_result_max[i * m_expansionLen + j] = pol_y[i * m_expansionLen + j];
	    }
      }
      else
      {
        for (int j = 0; j < m_expansionLen; j++)
      	{
	        expected_result_min[i * m_expansionLen + j] = pol_y[i * m_expansionLen + j];
	        expected_result_max[i * m_expansionLen + j] = pol_x[i * m_expansionLen + j];
	    }
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
    
    Ctxt ctxt_min(m_pk);
    Ctxt ctxt_max(m_pk);

    // comparison function
    cout << "Start of Min/Max" << endl;
    min_max(ctxt_min, ctxt_max, ctxt_x, ctxt_y);

    if(m_verbose)
    {
      cout << "Input" << endl;
      for(int i = 0; i < nslots; i++)
      {
          printZZX(cout, pol_x[i], ord_p);
          printZZX(cout, pol_y[i], ord_p);
          cout << endl;
      }

      cout << "Output min" << endl;
      print_decrypted(ctxt_min);
      cout << endl;

      cout << "Output max" << endl;
      print_decrypted(ctxt_max);
      cout << endl;
    }
    printNamedTimer(cout, "Extraction");
    printNamedTimer(cout, "MinMax");

    const FHEtimer* min_max_timer = getTimerByName("MinMax");

    cout << "Avg. time per integer: " << 1000.0 * min_max_timer->getTime()/static_cast<double>(run+1)/static_cast<double>(numbers_size) << " ms" << endl;
    cout << "Number of integers in one ciphertext "<< numbers_size << endl;

    // remove the line below if it gives bizarre results 
    ctxt_min.cleanUp();
    capacity = ctxt_min.bitCapacity();
    ctxt_max.cleanUp();
    cout << "Final capacity: " << capacity << endl;
    if (capacity < min_capacity)
      min_capacity = capacity;
    cout << "Min. capacity: " << min_capacity << endl;
    cout << "Final size: " << ctxt_min.logOfPrimeSet()/log(2.0) << endl;
    ea.decrypt(ctxt_min, m_sk, decrypted_min);
    ea.decrypt(ctxt_max, m_sk, decrypted_max);

    for(int i = 0; i < numbers_size; i++)
    { 
      if (decrypted_min[i * m_expansionLen] != expected_result_min[i * m_expansionLen])
      {
        printf("Slot %ld: ", i * m_expansionLen);
        printZZX(cout, decrypted_min[i * m_expansionLen], ord_p);
        cout << endl;
        cout << "Failure" << endl;
        return;
      }
    }
    cout << endl;
    for(int i = 0; i < numbers_size; i++)
    { 
      if (decrypted_max[i * m_expansionLen] != expected_result_max[i * m_expansionLen])
      {
        printf("Slot %ld: ", i * m_expansionLen);
        printZZX(cout, decrypted_max[i * m_expansionLen], ord_p);
        cout << endl;
        cout << "Failure" << endl;
        return;
      }
    }
    cout << endl;
  }
}

void Comparator::test_array_min(int input_len, long depth, long runs) const
{
	//reset timers
  setTimersOn();
  
  // initialize the random generator
  random_device rd;
  mt19937 eng(rd());
  uniform_int_distribution<unsigned long> distr_u;
  uniform_int_distribution<long> distr_i;

  // get EncryptedArray
  const EncryptedArray& ea = m_context.getEA();

  //extract number of slots
  long nslots = ea.size();

  //get p
  unsigned long p = m_context.getP();

  //order of p
  unsigned long ord_p = m_context.getOrdP();

  //amount of numbers in one ciphertext
  unsigned long numbers_size = nslots / m_expansionLen;

  // number of slots occupied by encoded numbers
  unsigned long occupied_slots = numbers_size * m_expansionLen;

  //encoding base, ((p+1)/2)^d
  //if 2-variable comparison polynomial is used, it must be p^d
  unsigned long enc_base = (p + 1) >> 1;
  if (m_type == BI || m_type == TAN)
  {
  	enc_base = p;
  }

  unsigned long digit_base = power_long(enc_base, m_slotDeg);

  //check that field_size^expansion_len fits into 64-bits
  int space_bit_size = static_cast<int>(ceil(m_expansionLen * log2(digit_base)));
  unsigned long input_range = ULONG_MAX;
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

    vector<ZZX> expected_result(occupied_slots,ZZX(INIT_MONO,0,0));
    
    // vector of input longs
    vector<vector<unsigned long>> input_xs;
    for (int i = 0; i < numbers_size; i++)
    {
    	vector<unsigned long> tmp_vec(input_len,0);
    	input_xs.push_back(tmp_vec);
    }

    ZZX pol_slot;

    //ciphertexts to sort
    vector<Ctxt> ctxt_in;

    //sorted ciphertexts
    Ctxt ctxt_out(m_pk);

    for (int i = 0; i < input_len; i++)
    {
		// the plaintext polynomials
		vector<ZZX> pol_x(nslots);

		//encoding of slots
		for (int k = 0; k < numbers_size; k++)
		{
			unsigned long input_x = distr_u(eng) % input_range;

			input_xs[k][i] = input_x;
		
			if(m_verbose)
			{
				cout << "Input" << endl;
				cout << input_x << endl;
			}

			vector<long> decomp_int_x;

			//decomposition of input integers
			digit_decomp(decomp_int_x, input_x, digit_base, m_expansionLen);
			for (int j = 0; j < m_expansionLen; j++)
			{
			    //decomposition of a digit
			    int_to_slot(pol_slot, decomp_int_x[j], enc_base);
			    pol_x[k * m_expansionLen + j] = pol_slot;
			}
		}

      	if(m_verbose)
	    {
	      cout << "Input" << endl;
	      for(int j = 0; j < nslots; j++)
	      {
	          printZZX(cout, pol_x[j], ord_p);
	          cout << endl;
	      }
	    }

	    Ctxt ctxt_x(m_pk);
    	ea.encrypt(ctxt_x, m_pk, pol_x);

    	ctxt_in.push_back(ctxt_x);
    }

    //cout << "Input" << endl;
    vector<unsigned long> output_xs(numbers_size, 0);
    for (int i = 0; i < numbers_size; i++)
    {
    	/*
    	for (int j = 0; j < input_len; j++)
    		cout << input_xs[i][j] << " ";
    	cout << endl;
    	*/
    	output_xs[i] = *std::min_element(input_xs[i].begin(), input_xs[i].end());
    	//cout << "Output: " << output_xs[i] << endl;
    }

    //cout << "Expected results" << endl;
	for (int k = 0; k < numbers_size; k++)
	{
		vector<long> decomp_int_x;

		//decomposition of input integers
		digit_decomp(decomp_int_x, output_xs[k], digit_base, m_expansionLen);
		for (int j = 0; j < m_expansionLen; j++)
		{
		    //decomposition of a digit
		    int_to_slot(pol_slot, decomp_int_x[j], enc_base);
		    expected_result[k * m_expansionLen + j] = pol_slot;
		}
	}

	/*
	cout << input_len-1-i << endl;
	for(int j = 0; j < nslots; j++)
	{
		printZZX(cout, expected_result[input_len-1-i][j], ord_p);
    	cout << endl;
	}
	*/

    // comparison function
    cout << "Start of array minimum" << endl;
    this->array_min(ctxt_out, ctxt_in, depth);

    printNamedTimer(cout, "Extraction");
    printNamedTimer(cout, "ComparisonCircuitBivar");
    printNamedTimer(cout, "ComparisonCircuitUnivar");
    printNamedTimer(cout, "EqualityCircuit");
    printNamedTimer(cout, "ShiftMul");
    printNamedTimer(cout, "ShiftAdd");
    printNamedTimer(cout, "Comparison");
    printNamedTimer(cout, "ArrayMin");

    const FHEtimer* sort_timer = getTimerByName("ArrayMin");

    cout << "Avg. time per batch: " << 1000.0 * sort_timer->getTime()/static_cast<double>(run+1)/static_cast<double>(numbers_size) << " ms" << endl;
    cout << "Number of integers in one ciphertext "<< numbers_size << endl;

    // remove the line below if it gives bizarre results 
    ctxt_out.cleanUp();
    capacity = ctxt_out.bitCapacity();
    cout << "Final capacity: " << capacity << endl;
    if (capacity < min_capacity)
      min_capacity = capacity;
    cout << "Min. capacity: " << min_capacity << endl;
    cout << "Final size: " << ctxt_out.logOfPrimeSet()/log(2.0) << endl;
    
	vector<ZZX> decrypted(nslots);
    ea.decrypt(ctxt_out, m_sk, decrypted);

    for(int j = 0; j < numbers_size; j++)
    { 
    	for(int k = 0; k < m_expansionLen; k++)
    	{
    		if (decrypted[j * m_expansionLen + k] != expected_result[j * m_expansionLen + k])
		    {
		    	printf("Slot %ld: ", j * m_expansionLen + k);
		    	printZZX(cout, decrypted[j * m_expansionLen + k], ord_p);
		        cout << endl;
		        cout << "Failure" << endl;
		        return;
		    }
    	}
    }
  }
}
