// Compile with 'g++ -O4 -o young2d-gmp young2d-gmp.cpp -lgmpxx -lgmp'

#include <iostream>
#include <vector>
#include <array>
#include <gmpxx.h>
#include <climits>

using namespace std;

typedef array<mpz_class, 2> val_tuple;
typedef struct {
  vector<mpz_class> coeffs;	// coefficients
  vector<val_tuple> monomials;	// exponents
} polynomial;

void read_input (polynomial &p, mpz_class &B) {
  cout << "Young Tableaux 2D Multiprecision" << endl;
  cout << "================================" << endl;
  cout << endl;

  cerr << "+++ Input B: ";
  cin >> B;
  cout << endl << "*** B = " << B << endl;
  int k;
  cin >> k;
  if (k != 2) {
    cerr << "*** Example has not 2 variables" << endl;
    exit(1);
  }
  cerr << "+++ Monomials in the form 'coefficient exponent exponent', ";
  cerr << "one per line" << endl;
  cerr << "+++ Terminate with CTRL-D" << endl;
  
  mpz_class temp;
  while (cin >> temp) {
    p.coeffs.push_back(temp);
    val_tuple exponents;
    cin >> exponents[0];
    cin >> exponents[1];
    p.monomials.push_back(exponents);
    cerr << "*** monomial = " << temp;
    if (exponents[0] > 0)
      cerr << " (x_1)^" << exponents[0];
    if (exponents[1] > 0)
    cerr << " (x_2)^" << exponents[1];
    cerr << endl;
  }

  cout << "*** equation: ";
  bool plus = false;
  for (unsigned int i = 0; i < p.coeffs.size(); ++i) {
    if (plus)
      cout << " + ";
    else
      plus = true;
    cout << p.coeffs[i];
    if (p.monomials[i][0] > 0)
      cout << " (x_1)^" << p.monomials[i][0];
    if (p.monomials[i][1] > 0)
      cout << " (x_2)^" << p.monomials[i][1];
  }
  cout <<  " = " << B << endl;
}

mpz_class power(mpz_class x, mpz_class n) {
  if (n == 0)
    return 1;
  mpz_class y = 1;
  while (n > 1)
    if (n % 2 == 0) {
      x *= x;
      n >>= 1;
    } else {
      y *= x;
      x *= x;
      n = (n-1)/2;
    }
  return x * y;
}

mpf_class power(mpf_class x, mpz_class n) {
  if (n == 0)
    return 1;
  mpf_class y = 1.0;
  while (n > 1)
    if (n % 2 == 0) {
      x *= x;
      n >>= 1;
    } else {
      y *= x;
      x *= x;
      n = (n-1)/2;
    }
  return x * y;
}

mpz_class eval(const polynomial &p, const val_tuple &val) {
  mpz_class add = 0;
  for (unsigned int i = 0; i < p.coeffs.size(); ++i) {
    mpz_class mult = p.coeffs[i];
    mult *= power(val[0], p.monomials[i][0]);
    mult *= power(val[1], p.monomials[i][1]);
    add += mult;
  }
  return add;
}

// mpz_class Horner1 (const polynomial &p,
// 			    const mpz_class value) {
//   // Horner scheme for univariate polynomials
//   if (p.k != 1) {
//     cerr << "+++ not an univariate polynomial";
//     exit(1);
//   }
//   mpz_class res = 0;
//   for (unsigned int i = p.coeffs.size()-1; i > 0; --i) {
//     res += p.coeffs[i];
//     res *= power(value, p.monomials[i][1] - p.monomials[i-1][1]);
//   }
//   res += p.coeffs[0];
//   res *= power(value, p.monomials[0][1]);
//   return res;
// }

mpf_class nthroot (const mpf_class A, const mpz_class n) {
  mpf_class oldx, x = A/n;
  do {
    oldx = x;
    x = ((n-1)*oldx + A/power(oldx, n-1))/n;
  } while (abs(oldx-x) > 0.01);
  return x;
}

vector<mpz_class> get_bounds (const polynomial &p,
			      const mpz_class &B) {
  vector<mpz_class> bound(2);
  vector<mpz_class> minexp(2, UINT_MAX);
  vector<unsigned long long> minpos(2);
  for (int j = 0; j < p.monomials.size(); ++j)
    for (int i = 0; i <= 1; ++i)
      if (p.monomials[j][i] > 0 && p.monomials[j][i] < minexp[i]) {
	minexp[i] = p.monomials[j][i];
	minpos[i] = j;
      }
  for (int i = 0; i <= 1; ++i) {
    mpf_class bpc = B/p.coeffs[minpos[i]];
    bound[i] = nthroot(bpc, minexp[i]) + 1;
  }
  return bound;
}

//////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  mpz_class B;
  polynomial p;

  read_input(p, B);
  vector<mpz_class> bound = get_bounds(p, B);
  mpz_class row = 0;
  mpz_class column = bound[1];
  while (row < bound[0] && column >= 0) {
    val_tuple val = {row, column};
    mpz_class result = eval(p, val);
    if (result == B) {
      cout << endl << "+++ YES +++" << endl;
      cout << "*** for values:" << endl;
      cout << "    x_1 = " << val[0] << endl;
      cout << "    x_2 = " << val[1] << endl;
      exit(0);
    } else if (result < B)
      row++;
    else if (result > B)
      column--;
  }

  cout << endl << "+++ NO +++" << endl;
}
//////////////////////////////////////////////////////////////////////////////
