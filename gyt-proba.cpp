#include <iostream>
#include <vector>
#include <climits>
#include <random>

using namespace std;

typedef vector<unsigned long long> val_tuple;
typedef struct {
  unsigned int k;		// number of variables
  val_tuple coeffs;		// coefficients
  vector<val_tuple> monomials;	// exponents
} polynomial;

void read_input (polynomial &p, unsigned long long &B) {
  cout << "Probabilistic Multidimensional Generalized Young Tableaux"
       << endl;
  cout << "========================================================="
       << endl;
  cout << endl;

  cerr << "+++ Input B: ";
  cin >> B;
  cout << "*** B = " << B << endl;
  cerr << "+++ Number of variables: ";
  cin >> p.k;
  cerr << "+++ Monomials in the form 'coefficient exponent ... exponent', ";
  cerr << "one per line" << endl;
  cerr << "+++ Terminate with CTRL-D" << endl;
  
  unsigned long long temp;
  while (cin >> temp) {
    p.coeffs.push_back(temp);
    val_tuple exponents;
    unsigned long long exp;
    for (unsigned int i = 0; i < p.k; ++i) {
      cin >> exp;
      exponents.push_back(exp);
    }
    p.monomials.push_back(exponents);
    cerr << "*** monomial = " << temp;
    for (unsigned int i = 0; i < p.k; ++i)
      if (exponents[i] > 0)
	cerr << " (x_" << i+1 << ")^" << exponents[i];
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
    for (unsigned int j = 0; j < p.k; ++j)
      if (p.monomials[i][j] > 0)
	cout << " (x_" << j+1 << ")^" << p.monomials[i][j];
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
    return 1.0;
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
    for (unsigned j = 0; j < p.k; ++j)
      mult *= power(val[j], p.monomials[i][j]);
    add += mult;
  }
  return add;
}

long double nthroot (const long double A, const unsigned long long n) {
  long double oldx, x = A/n;
  do {
    oldx = x;
    x = ((n-1)*oldx + A/power(oldx, n-1))/n;
  } while (abs(oldx-x) > 0.01);
  return x;
}

val_tuple get_bounds (const polynomial &p,
		      const unsigned long long &B) {
  val_tuple bound;
  val_tuple minpos(p.k);
  val_tuple minexp(p.k, ULLONG_MAX);
  for (unsigned int j = 0; j < p.monomials.size(); ++j)
    for (unsigned int i = 0; i < p.k; ++i)
      if (p.monomials[j][i] > 0 && p.monomials[j][i] < minexp[i]) {
	minexp[i] = p.monomials[j][i];
	minpos[i] = j;
      }
  for (unsigned int i = 0; i < p.k; ++i)
    bound.push_back(nthroot(B/p.coeffs[minpos[i]], minexp[i]) + 1);
  
  return bound;
}

bool test_bound (const val_tuple &val, const val_tuple &bound, int k) {
  for (unsigned int i = 0; i < k; ++i)
    if (val[i] > bound[i])
      return false;
  return true;
}

//////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  unsigned long long B;
  polynomial p;

  read_input(p, B);
  val_tuple bound = get_bounds(p, B);

  random_device rd;
  static uniform_int_distribution<int> uni_dist(0,p.k-2);
  static default_random_engine dre(rd());

  val_tuple val(p.k, 0);
  val[p.k-1] = bound[p.k-1];
  bool solution = false;
  unsigned long long choice = 0;

  while (test_bound(val, bound, p.k) &&
	 val[p.k-1] >= 0) {
    unsigned long long result = eval(p, val);
    if (result == B) {
      solution = true;
      break;
    } else if (result > B)
      val[p.k-1]--;
    else if (result < B) {
      val[uni_dist(dre)]++;
      choice++;
    }
  }

  if (solution) {
    cout << endl << "+++ YES +++" << endl;
    cout << "*** solution for values:" << endl;
    for (unsigned int i = 0; i < p.k; ++i)
      cout << "    x_" << i+1 << " = " << val[i] << endl;
  } else
    cout << endl << "+++ NO solution +++" << endl;
  cout << "*** # of choices = " << choice << endl;
}
//////////////////////////////////////////////////////////////////////////////
