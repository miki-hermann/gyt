// Compile with 'g++ -O4 -o young-all-gmp young-all-gmp.cpp -lgmpxx -lgmp'

#include <iostream>
#include <vector>
#include <stack>
#include <set>
#include <gmpxx.h>
#include "gyt-common-gmp.hpp"

using namespace std;

const string header    = "Sequential All-Solution Multiprecision Multidimensional Generalized Young Tableaux";
const string underline = "==================================================================================";

//////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  mpz_class B;
  polynomial p;
  stack<val_tuple> stck;
  set<val_tuple> memo;
  set<val_tuple> sols;
  mpz_class nres = 0;
  mpz_class maxstack = 1;
  bool flip = true;
  mpz_class nback = 0;
  mpz_class split = 0;

  read_input(p, B);
  val_tuple bound = get_bounds(p, B);

  val_tuple val(p.k, 0);
  val[p.k-1] = bound[p.k-1];
  stck.push(val);
  memo.insert(val);
  mpz_class dbl = 0;
  while (!stck.empty()) {
    val = stck.top();
    stck.pop();
    nback += !flip;
    flip = false;

    while (test_bound(val, bound, p.k) && val[p.k-1] >= 0) {
      mpz_class result = eval(p, val);
      if (result == B) {
	if (sols.find(val) == sols.cend()) {
	  nres++;
	  cout << endl << "*** solution for values:" << endl;
	  for (unsigned int i = 0; i < p.k; ++i)
	    cout << "    x_" << i+1 << " = " << val[i] << endl;
	  sols.insert(val);
	  val[p.k-1]--;
	} else
	  break;
      } else if (result > B)
	val[p.k-1]--;
      else if (result < B) {
	mpz_class put = 0;
	for (unsigned int i = 0; i < p.k-1; ++i) {
	  val_tuple valx = val;
	  valx[i]++;
	  if (valx[i] <= bound[i] && memo.find(valx) == memo.cend()) {
	    stck.push(valx);
	    if (p.k > 2)
	      memo.insert(valx);
	    put++;
	    flip = true;
	  } else if (valx[i] <= bound[i])
	    dbl++;
	}
	split += put > 1;
	maxstack = max(maxstack, mpz_class(stck.size()));
	break;
      }
    }
  }

  cout << endl;
  cout << "+++ number of solutions = " << nres << endl;
  unsigned long long msize = memo.size();

  statistics(memo.size(), "stack", maxstack, split, nback, dbl);
}
//////////////////////////////////////////////////////////////////////////////
