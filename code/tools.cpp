#include "tools.h"

void digit_decomp(vector<long>& decomp, unsigned long input, unsigned long base, int nslots)
{
  decomp.clear();
  decomp.resize(nslots,0);
  int power = intlog(base, input) + 1;
  if (power > nslots)
  {
    cout << "Input character is too big to be converted" << endl;
    exit(1);
  }
  unsigned long rest = input;
  unsigned long coeff;

  int i = 0;
  while(i < power)
    {
      coeff = rest % base;
      decomp[i] = coeff;
      rest = (rest - coeff) / base;
      i++;
    }
}

// Simple evaluation sum f_i * X^i, assuming that babyStep has enough powers
void simplePolyEval(Ctxt& ret, const NTL::ZZX& poly, DynamicCtxtPowers& babyStep)
{
  ret.clear();
  if (deg(poly)<0) return;       // the zero polynomial always returns zero

  //OLD: assert(deg(poly)<=babyStep.size()); // ensure that we have enough powers
  helib::assertTrue(deg(poly)<=babyStep.size(), "BabyStep has not enough powers (required more than deg(poly))");

  NTL::ZZ coef;
  NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
  for (long i=1; i<=deg(poly); i++) {
    rem(coef, coeff(poly,i),p);
    if (coef > p/2) coef -= p;

    Ctxt tmp = babyStep.getPower(i); // X^i
    tmp.multByConstant(coef);        // f_i X^i
    ret += tmp;
  }
  // Add the free term
  rem(coef, ConstTerm(poly), p);
  if (coef > p/2) coef -= p;
  ret.addConstant(coef);
  //  if (verbose) checkPolyEval(ret, babyStep[0], poly);
}

// The recursive procedure in the Paterson-Stockmeyer
// polynomial-evaluation algorithm from SIAM J. on Computing, 1973.
// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains >= k+delta powers
void PatersonStockmeyer(Ctxt& ret, const NTL::ZZX& poly, long k, long t, long delta,
		   DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep)
{
  if (deg(poly)<=babyStep.size()) { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }
  NTL::ZZX r = trunc(poly, k*t);      // degree <= k*2^e-1
  NTL::ZZX q = RightShift(poly, k*t); // degree == k(2^e-1) +delta

  const NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
  const NTL::ZZ& coef = coeff(r,deg(q));
  SetCoeff(r, deg(q), coef-1);  // r' = r - X^{deg(q)}

  NTL::ZZX c,s;
  DivRem(c,s,r,q); // r' = c*q + s
  // deg(s)<deg(q), and if c!= 0 then deg(c)<k-delta

  helib::assertTrue(deg(s)<deg(q), "Degree of s is not less than degree of q");
  helib::assertTrue(IsZero(c) || deg(c)<k - delta, "Nonzero c has not degree smaller than k - delta");
  SetCoeff(s,deg(q)); // s' = s + X^{deg(q)}, deg(s)==deg(q)

  // reduce the coefficients modulo p
  for (long i=0; i<=deg(c); i++) rem(c[i],c[i], p);
  c.normalize();
  for (long i=0; i<=deg(s); i++) rem(s[i],s[i], p);
  s.normalize();

  // Evaluate recursively poly = (c+X^{kt})*q + s'
  PatersonStockmeyer(ret, q, k, t/2, delta, babyStep, giantStep);

  Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
  simplePolyEval(tmp, c, babyStep);
  tmp += giantStep.getPower(t);
  ret.multiplyBy(tmp);

  PatersonStockmeyer(tmp, s, k, t/2, delta, babyStep, giantStep);
  ret += tmp;
}

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains >= k + (deg(poly) mod k) powers
void degPowerOfTwo(Ctxt& ret, const NTL::ZZX& poly, long k,
	      DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep)
{
  if (deg(poly)<=babyStep.size()) { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }
  long n = divc(deg(poly),k);        // We assume n=2^e or n=2^e -1
  n = 1L << NTL::NextPowerOfTwo(n); // round up to n=2^e
  NTL::ZZX r = trunc(poly, (n-1)*k);      // degree <= k(2^e-1)-1
  NTL::ZZX q = RightShift(poly, (n-1)*k); // 0 < degree < 2k
  SetCoeff(r, (n-1)*k);              // monic, degree == k(2^e-1)
  q -= 1;

  PatersonStockmeyer(ret, r, k, n/2, 0,	babyStep, giantStep);

  Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
  simplePolyEval(tmp, q, babyStep); // evaluate q

  // multiply by X^{k(n-1)} with minimum depth
  for (long i=1; i<n; i*=2) {  
    tmp.multiplyBy(giantStep.getPower(i));
  }
  ret += tmp;
}

void recursivePolyEval(Ctxt& ret, const NTL::ZZX& poly, long k,
		  DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep)
{
  if (deg(poly)<=babyStep.size()) { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }

  long delta = deg(poly) % k; // deg(poly) mod k
  long n = divc(deg(poly),k); // ceil( deg(poly)/k )
  long t = 1L<<(NTL::NextPowerOfTwo(n)); // t >= n, so t*k >= deg(poly)

  // Special case for deg(poly) = k * 2^e +delta
  if (n==t) {
    degPowerOfTwo(ret, poly, k, babyStep, giantStep);
    return;
  }

  // When deg(poly) = k*(2^e -1) we use the Paterson-Stockmeyer recursion
  if (n == t-1 && delta==0) {
    PatersonStockmeyer(ret, poly, k, t/2, delta, babyStep, giantStep);
    return;
  }

  t = t/2;

  // In any other case we have kt < deg(poly) < k(2t-1). We then set 
  // u = deg(poly) - k*(t-1) and poly = q*X^u + r with deg(r)<u
  // and recurse on poly = (q-1)*X^u + (X^u+r)

  long u = deg(poly) - k*(t-1);
  NTL::ZZX r = trunc(poly, u);      // degree <= u-1
  NTL::ZZX q = RightShift(poly, u); // degree == k*(t-1)
  q -= 1;
  SetCoeff(r, u);              // degree == u

  PatersonStockmeyer(ret, q, k, t/2, 0, babyStep, giantStep);

  Ctxt tmp = giantStep.getPower(u/k);
  if (delta!=0) { // if u is not divisible by k then compute it
    tmp.multiplyBy(babyStep.getPower(delta));
  }
  ret.multiplyBy(tmp);

  recursivePolyEval(tmp, r, k, babyStep, giantStep);
  ret += tmp;
}