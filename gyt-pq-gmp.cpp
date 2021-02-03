// Compile with 'g++ -O4 -o young-pq-gmp young-pq-gmp.cpp -lgmpxx -lgmp'

#include <iostream>
#include <vector>
#include <climits>
#include <queue>
#include <set>
#include <utility>
#include <gmpxx.h>

using namespace std;

typedef vector<mpz_class> val_tuple;
typedef struct {
  unsigned int k;		// number of variables
  val_tuple coeffs;		// coefficients
  vector<val_tuple> monomials;	// exponents
} polynomial;
typedef pair<val_tuple, mpz_class> val_res;
struct cmp_vr {
  bool operator()(const val_res &vr1, const val_res &vr2)
  {
    return vr1.second > vr2.second;
  }
};

void read_input (polynomial &p, mpz_class &B) {
  cout << "Priority-Driven Multidimensional Multiprecision Generalized Young Tableaux"
       << endl;
  cout << "=========================================================================="
       << endl;
  cout << endl;

  cerr << "+++ Input B: ";
  cin >> B;
  cout << endl << "*** B = " << B << endl;
  cerr << "+++ Number of variables: ";
  cin >> p.k;
  cout << endl << "*** k = " << p.k << endl;
  cerr << "+++ Monomials in the form 'coefficient exponent ... exponent', ";
  cerr << "one per line" << endl;
  cerr << "+++ Terminate with CTRL-D" << endl;
  
  mpz_class temp;
  while (cin >> temp) {
    p.coeffs.push_back(temp);
    val_tuple exponents;
    mpz_class exp;
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
    return 1.0;
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
    for (unsigned j = 0; j < p.k; ++j)
      mult *= power(val[j], p.monomials[i][j]);
    add += mult;
  }
  return add;
}

mpf_class nthroot (const mpf_class A, const mpz_class n) {
  mpf_class oldx, x = A/n;
  do {
    oldx = x;
    x = ((n-1)*oldx + A/power(oldx, n-1))/n;
  } while (abs(oldx-x) > 0.01);
  return x;
}

val_tuple get_bounds (const polynomial &p,
		      const mpz_class &B) {
  val_tuple bound(p.k);
  vector<unsigned long long> minpos(p.k);
  val_tuple minexp(p.k, UINT_MAX);
  for (unsigned int j = 0; j < p.monomials.size(); ++j)
    for (unsigned int i = 0; i < p.k; ++i)
      if (p.monomials[j][i] > 0 && p.monomials[j][i] < minexp[i]) {
	minexp[i] = p.monomials[j][i];
	minpos[i] = j;
      }
  for (unsigned int i = 0; i < p.k; ++i) {
    mpf_class bpc = B/p.coeffs[minpos[i]];
    bound[i] = nthroot(bpc, minexp[i]) + 1;
  }

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
  mpz_class B;
  polynomial p;
  priority_queue<val_res, vector<val_res>, cmp_vr> pq;
  set<val_tuple> memo;
  unsigned int maxstack = 1;
  bool flip = true;
  unsigned int nback = 0;
  unsigned int split = 0;

  read_input(p, B);
  val_tuple bound = get_bounds(p, B);

  val_tuple val(p.k, 0);
  val[p.k-1] = bound[p.k-1];
  pq.push(make_pair(val, 0));
  memo.insert(val);
  bool solution = false;
  unsigned int dbl = 0;
  while (!solution && !pq.empty()) {
    val = pq.top().first;
    pq.pop();
    nback += !flip;
    flip = false;

    while (test_bound(val, bound, p.k) &&
	   val[p.k-1] >= 0) {
      mpz_class result = eval(p, val);
      if (result == B) {
	solution = true;
	break;
      } else if (result > B)
	val[p.k-1]--;
      else if (result < B) {
	unsigned int put = 0;
	for (unsigned int i = 0; i < p.k-1; ++i) {
	  val_tuple valx = val;
	  valx[i]++;
	  mpz_class resx = eval(p, valx);
	  if (valx[i] <= bound[i] && memo.find(valx) == memo.cend()) {
	    pq.push(make_pair(valx, abs(resx-B)));
	    if (p.k > 2)
	      memo.insert(valx);
	    put++;
	    flip = true;
	  } else if (valx[i] <= bound[i])
	    dbl++;
	}
	split += put > 1;
	maxstack = max(maxstack, unsigned(pq.size()));
	break;
      }
    }
  }

  if (solution) {
    cout << endl << "+++ YES +++" << endl;
    cout << "*** solution for values:" << endl;
    for (unsigned int i = 0; i < p.k; ++i)
      cout << "    x_" << i+1 << " = " << val[i] << endl;
  } else
    cout << endl << "+++ NO solution +++" << endl;
  cout << "*** memo size       = " << memo.size() << endl;
  cout << "    max queue size  = " << maxstack << endl;
  cout << "    queue residue   = " << pq.size() << endl;
  cout << "    # of splits     = " << split << endl;
  cout << "    # of backtracks = " << nback << endl;
  cout << "    doubles reached = " << dbl << endl;
}
//////////////////////////////////////////////////////////////////////////////
