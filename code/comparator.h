/*
Comparator class: creates all auxiliary information necessary to compare integers
*/

#ifndef COMPARATOR_H
#define COMPARATOR_H

#include <helib/helib.h>

using namespace std;
using namespace NTL;
using namespace helib;


class Comparator{
	const Context& m_context;

    // field extension degree
  	unsigned long m_slotDeg;

  	// expansion length
  	unsigned long m_expansionLen;

  	// vector of multiplicative masks
  	vector<ZZX> m_mulMasks;

    // comparison polynomial
    ZZX m_poly;

    // polynomial evaluation parameters of the Patterson-Stockmeyer algorithm
    // number of baby steps
    long m_bs_num;
    // number of giant steps
    long m_gs_num;
    // leading coefficient
    ZZ m_top_coef;
    // extra coefficient
    ZZ m_extra_coef;
    
    // indexes to compute x^{p-1}
    long m_baby_index;
    long m_giant_index;

    // slot generator
    ZZX m_slot_gen;

    // secret key
    SecKey m_sk;

    // public key
    PubKey m_pk;

  	// print/hide flag for debugging
  	bool m_verbose;

    // create multiplicative masks for shifts
  	void create_shift_mask(ZZX& mask_ptxt, long shift);
  	void create_all_shift_masks();

    // compute Patterson-Stockmeyer parameters to evaluate the comparison polynomial
    void compute_poly_params();

    // create the comparison polynomial
    void create_poly();

    // initialize slot generator
    void slot_init();

    // shifts ciphertext slots to the left by shift within batches of size m_expansionLen starting at start. Slots shifted outside their respective batches are zeroized.
    void batch_shift(Ctxt& ctxt, long start, long shift) const;
    
    // shifts ciphertext slots to the left by shift within batches of size m_expansionLen starting at start. Slots shifted outside their respective batches filled with 1.
    void batch_shift_for_mul(Ctxt& ctxt, long start, long shift) const;

    // running sums of slot batches
    void shift_and_add(Ctxt& x, long start, long shift_direction = false) const;

    // running products of slot batches
    void shift_and_mul(Ctxt& x, long start, long shift_direction = false) const;

    // send non-zero elements of a field F_{p^d} to 1 and zero to 0
    // if d = 1, this map operates on elements of the prime field F_p
    void mapTo01_subfield(Ctxt& ctxt, long d) const;

    // comparison polynomial evaluation
    void evaluate_poly(Ctxt& ret, const Ctxt& x, const Ctxt& x2) const;

    // less than function comparing slots one by one
    void less_than(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // less than function comparing slots one by one in F_2
    void less_than_mod_2(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;
  
    // less than function comparing slots one by one in F_3
    void less_than_mod_3(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // less than function comparing slots one by one in F_5
    void less_than_mod_5(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // less than function comparing slots one by one in F_7
    void less_than_mod_7(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;  

    // less than function comparing slots one by one in F_11
    void less_than_mod_11(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;
  
    // randomized equality
    void random_equality(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // exact equality 
    void exact_equality(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // conversion to slots
    void int_to_slot(ZZX& poly, unsigned long input) const; 

public:
  // constructor
	Comparator(const Context& context, unsigned long d, unsigned long expansion_len, const SecKey& sk, bool verbose);

	const ZZX& get_mask(long index) const;
  const ZZX& get_poly() const;

  // decrypt and print ciphertext
  void print_decrypted(const Ctxt& ctxt) const;

  // comparison function
  void compare(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y, bool is_randomized = false) const;

  // test compare function 'runs' times
  void test(long runs, bool is_randomized = false) const;
};

#endif // #ifndef COMPARATOR_H
