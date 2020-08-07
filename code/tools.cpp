#include "tools.h"

void digit_decomp(vector<long>& decomp, unsigned long input, unsigned long base, int nslots)
{
  decomp.clear();
  decomp.resize(nslots,0);
  int power = static_cast<int>(intlog(base, input)) + 1;
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