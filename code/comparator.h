/*
Comparator class: creates all auxiliary information necessary to compare integers
*/

#ifndef COMPARATOR_H
#define COMPARATOR_H

#include <helib/helib.h>
#include <helib/Ptxt.h>
#include <helib/norms.h>
#include <NTL/mat_ZZ.h>

using namespace std;
using namespace NTL;
using namespace helib;

namespace he_cmp{
enum CircuitType{UNI, BI, TAN};

class Comparator{
    const Context& m_context;

    // field extension degree (d)
  	unsigned long m_slotDeg;

  	// expansion length
  	unsigned long m_expansionLen;

  	// vector of multiplicative masks
  	vector<DoubleCRT> m_mulMasks;

    // vector of multiplicative mask size
    vector<double> m_mulMasksSize;

    // univariate or bivariate circuit
    CircuitType m_type;

    // univariate comparison polynomial of the less-than function
    ZZX m_univar_less_poly;

    // univariate comparison polynomial of the less-than function
    ZZX m_univar_min_max_poly;

    // bivariate comparison polynomial coefficients of the less-than function
    mat_ZZ m_bivar_less_coefs; 

    // polynomial evaluation parameters of the Patterson-Stockmeyer algorithm
    // number of baby steps
    long m_bs_num_comp;
    long m_bs_num_min;
    // number of giant steps
    long m_gs_num_comp;
    long m_gs_num_min;
    // leading coefficient
    ZZ m_top_coef_comp;
    ZZ m_top_coef_min;
    // extra coefficient
    ZZ m_extra_coef_comp;
    ZZ m_extra_coef_min;
    
    // indexes to compute x^{p-1}
    long m_baby_index;
    long m_giant_index;

    // slot generator
    ZZX m_slot_gen;

    // secret key
    SecKey m_sk;

    // public key
    PubKey m_pk;

    // elements of F_{p^d} for extraction of F_p elements
    vector<vector<DoubleCRT>> m_extraction_const;
    vector<vector<double>> m_extraction_const_size;

  	// print/hide flag for debugging
  	bool m_verbose;

    // create multiplicative masks for shifts
  	DoubleCRT create_shift_mask(double& size, long shift);
  	void create_all_shift_masks();

    // compute Patterson-Stockmeyer parameters to evaluate the comparison polynomial
    void compute_poly_params();

    // create the comparison polynomial
    void create_poly();

    // initialize extraction constants
    void extraction_init();

    // extract F_p elements from slots
    void extract_mod_p(vector<Ctxt>& mod_p_coefs, const Ctxt& ctxt_x) const; 

    // shifts ciphertext slots to the left by shift within batches of size m_expansionLen starting at start. Slots shifted outside their respective batches are zeroized.
    void batch_shift(Ctxt& ctxt, long start, long shift) const;
    
    // shifts ciphertext slots to the left by shift within batches of size m_expansionLen starting at start. Slots shifted outside their respective batches filled with 1.
    void batch_shift_for_mul(Ctxt& ctxt, long start, long shift) const;

    // running sums of slot batches
    void shift_and_add(Ctxt& x, long start, long shift_direction = false) const;

    // running products of slot batches
    void shift_and_mul(Ctxt& x, long start, long shift_direction = false) const;

    // send non-zero elements of a field F_{p^d} to 1 and zero to 0
    // if pow = 1, this map operates on elements of the prime field F_p
    void mapTo01_subfield(Ctxt& ctxt, long pow) const;

    // univariate comparison polynomial evaluation
    void evaluate_univar_less_poly(Ctxt& ret, Ctxt& ctxt_p_1, const Ctxt& x) const;

    // univariate min/max polynomial evaluation
    void evaluate_min_max_poly(Ctxt& ctxt_min, Ctxt& ctxt_max, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // bivariate less than function comparing slots one by one
    void less_than_bivar(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // bivariate less than function comparing slots one by one as in Tan et al.
    void less_than_bivar_tan(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // less than function comparing slots one by one in F_2
    void less_than_mod_2(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;
  
    // less than function comparing slots one by one in F_3
    void less_than_mod_3(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // less than function comparing slots one by one in F_5
    void less_than_mod_5(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

    // less than function comparing slots one by one in F_7
    void less_than_mod_7(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;  

    // less than function comparing slots one by one in F_17
    void less_than_mod_any(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;    
  
    // exact equality 
    void is_zero(Ctxt& ctxt_res, const Ctxt& ctxt_z, long pow = 1) const;

    // minimum/maximum function for digits (vectors of dimension 1 over F_p)
    void min_max_digit(Ctxt& ctxt_min, Ctxt& ctxt_max, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const; 

    // conversion to slots
    void int_to_slot(ZZX& poly, unsigned long input, unsigned long enc_base) const;

    // compute an array of positions of ciphertexts in ctxt_in when sorted
    void get_sorting_index(vector<Ctxt>& ctxt_out, const vector<Ctxt>& ctxt_in) const;

    // find the primitive root of a SIMD slot
    void find_prim_root(ZZ_pE& root) const; 

public:
  // constructor
	Comparator(const Context& context, CircuitType type, unsigned long d, unsigned long expansion_len, const SecKey& sk, bool verbose);

	const DoubleCRT& get_mask(double& size, long index) const;
  const ZZX& get_less_than_poly() const;
  const ZZX& get_min_max_poly() const;

  // decrypt and print ciphertext
  void print_decrypted(const Ctxt& ctxt) const;

  // comparison function
  void compare(Ctxt& ctxt_res, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

  // minimum/maximum function for general vectors
  void min_max(Ctxt& ctxt_min, Ctxt& ctxt_max, const Ctxt& ctxt_x, const Ctxt& ctxt_y) const;

  // minimum/maximum of an array
  void array_min(Ctxt& ctxt_res, const vector<Ctxt>& ctxt_in, long depth = 0) const;

  // sorting
  void sort(vector<Ctxt>& ctxt_out, const vector<Ctxt>& ctxt_in) const;

  // test compare function 'runs' times
  void test_compare(long runs) const;

  // test min/max function 'runs' times
  void test_min_max(long runs) const;

  // test compare function 'runs' times
  void test_sorting(int num_to_sort, long runs) const;

  // test array_minn function
  void test_array_min(int input_len, long depth, long runs) const;
};
}

#endif // #ifndef COMPARATOR_H
