#ifndef _BITOPS_H
#define _BITOPS_H

#include <iostream>

using namespace std;

inline unsigned long ls(unsigned long n) {
  if(n > sizeof(unsigned long)*8) return 0;
  else return (unsigned long)(1) << n;
}

#endif
