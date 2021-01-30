#include <iostream>
#include <vector>
#include <array>
#include <climits>

using namespace std;

typedef array<unsigned long long, 2> val_tuple;
typedef struct {
  vector<unsigned long long> coeffs;	// coefficients
  vector<val_tuple> monomials;		// exponents
} polynomial;

void read_input (polynomial &p, unsigned long long &B) {
  cout << "Young Tableaux 2D" << endl;
  cout << "=================" << endl;
  cout << endl;

  cerr << "+++ Input B: ";
  cin >> B;
  cout << "*** B = " << B << endl;
  cerr << "+++ Number of variables: ";
  int k;
  cin >> k;
  if (k != 2) {
    cerr << "*** Example has not 2 variables" << endl;
    exit(1);
  }
  cerr << "+++ Monomials in the form 'coefficient exponent exponent', ";
  cerr << "one per line" << endl;
  cerr << "+++ Terminate with CTRL-D" << endl;
  
  unsigned long long temp;
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

unsigned long long power(unsigned long long x, unsigned long long n) {
  if (n == 0)
    return 1;
  unsigned long long y = 1;
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

long double power(long double x, unsigned long long n) {
  if (n == 0)
    return 1;
  long double y = 1.0;
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

unsigned long long eval(const polynomial &p, const val_tuple &val) {
  unsigned long long add = 0;
  for (unsigned int i = 0; i < p.coeffs.size(); ++i) {
    unsigned long long mult = p.coeffs[i];
    mult *= power(val[0], p.monomials[i][0]);
    mult *= power(val[1], p.monomials[i][1]);
    add += mult;
  }
  return add;
}

// unsigned long long Horner1 (const polynomial &p,
// 			    const unsigned long long value) {
//   // Horner scheme for univariate polynomials
//   if (p.k != 1) {
//     cerr << "+++ not an univariate polynomial";
//     exit(1);
//   }
//   unsigned long long res = 0;
//   for (unsigned int i = p.coeffs.size()-1; i > 0; --i) {
//     res += p.coeffs[i];
//     res *= power(value, p.monomials[i][1] - p.monomials[i-1][1]);
//   }
//   res += p.coeffs[0];
//   res *= power(value, p.monomials[0][1]);
//   return res;
// }

long double nthroot (const long double A, const unsigned long long n) {
  long double oldx, x = A/n;
  do {
    oldx = x;
    x = ((n-1)*oldx + A/power(oldx, n-1))/n;
  } while (abs(oldx-x) > 0.01);
  return x;
}

vector<unsigned long long> get_bounds (const polynomial &p,
				       const unsigned long long &B) {
  vector<unsigned long long> bound(2);
  vector<unsigned long long> minexp(2, UINT_MAX);
  vector<unsigned long long> minpos(2);
  for (int j = 0; j < p.monomials.size(); ++j)
    for (int i = 0; i <= 1; ++i)
      if (p.monomials[j][i] > 0 && p.monomials[j][i] < minexp[i]) {
	minexp[i] = p.monomials[j][i];
	minpos[i] = j;
      }
  for (int i = 0; i <= 1; ++i)
    bound[i] = nthroot(B/p.coeffs[minpos[i]], minexp[i]) + 1;
  return bound;
}

//////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  unsigned long long B;
  polynomial p;

  read_input(p, B);
  vector<unsigned long long> bound = get_bounds(p, B);
  unsigned long long row = 0;
  unsigned long long column = bound[1];
  while (row <= bound[0] && column >= 0) {
    val_tuple val = {row, column};
    unsigned long long result = eval(p, val);
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
