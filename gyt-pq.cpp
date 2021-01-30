#include <iostream>
#include <vector>
#include <climits>
#include <queue>
#include <set>
#include <utility>

using namespace std;

typedef vector<unsigned long long> val_tuple;
typedef struct {
  unsigned int k;		// number of variables
  val_tuple coeffs;		// coefficients
  vector<val_tuple> monomials;	// exponents
} polynomial;
typedef pair<val_tuple, unsigned long long> val_res;
struct cmp_vr {
  bool operator()(const val_res &vr1, const val_res &vr2)
  {
    return vr1.second > vr2.second;
  }
};

void read_input (polynomial &p, unsigned long long &B) {
  cout << "Priority-Driven Multidimensional Generalized Young Tableaux"
       << endl;
  cout << "==========================================================="
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
  // stack<val_tuple> stck;
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
      unsigned long long result = eval(p, val);
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
	  unsigned long long resx = eval(p, valx);
	  if (valx[i] <= bound[i] && memo.find(valx) == memo.cend()) {
	    long long rxB = resx - B;
	    pq.push(make_pair(valx, abs(rxB)));
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
