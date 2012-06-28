/* The following piece of C++ source code calculates the natural logarithm of 
 * the n-th Bell Number, i.e., the number of ways to partition a set of n 
 * elements. The function is based on the triangle algorithm, and requires 
 * O(n^2) time, and O(n) memory. The computations are done in double precision,
 * so no additional libraries are necessary. The results have been verified to
 * have a relative error of at most O(10^-15) for n up to one thousand. To get
 * the actual Bell number (in double precision) simply exponentiate the result.
 * -- code by Michael Mampaey, 2011
 */

#include <algorithm>
#include <math.h>

using namespace std;

double log_bell_number(int n)
{
  if (n == 0)
    return 0.0;

  double * t0 = new double[n];
  double * t1 = new double[n];

  t0[0] = 0.0;

  for (int i = 1; i < n; ++i) {
    t1[0] = t0[i-1];
    for (int j = 0; j < i; ++j)
      t1[j+1] = log(exp(t1[j]-t0[j]) + 1.0) + t0[j];
    swap(t0, t1);
  }

  double result = t0[n-1];

  delete [] t0;
  delete [] t1;

  return result;
}

