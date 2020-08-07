/*
Auxiliary functions for integer encoding
*/

#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

inline double intlog(unsigned long base, unsigned long input)
{
  return floor(log2(input)/log2(base));
}

void digit_decomp(vector<long>& decomp, unsigned long input, unsigned long base, int nslots);

#endif // #ifndef TOOLS_H