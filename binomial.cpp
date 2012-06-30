/* Computing binomial coefficients or probabilities according to a binomial 
 * distribution may seem like a trivial thing to do. However, browsing for 
 * code snippets around the web, it is surprising to see how often this is 
 * done naively by explicitly computing factorials. Needless to say, such 
 * code is not only inefficient, but unless you use some numeric library, it 
 * is also unusable beyond small numbers, say, 13 choose 7.
 * The following pieces of C++ code calculate binomial coefficients and 
 * distributions in logarithmic scale. The code only needs to include 
 * <math.h> for the log and exp functions, and simply uses double variables. 
 * Due to the log-scale, you can work with extremely large coefficients (or 
 * extremely small probabilities), albeit within the fixed precision of a 
 * double, but for many practical applications this is often sufficient. If 
 * you need the actual number rather than its natural logarithm, just 
 * exponentiate it (provided that it will fit in the target variable).
 * The time complexity is O(min(k,n-k)) for the first and third function, 
 * O(1) for the second function, and O(k) for the last one. Space complexity
 * is O(1) for all functions.
 * -- code by Michael Mampaey, 2011
 */

#include <math.h>
using namespace std;

/* compute the natural logarithm of the binomial coefficient
 * (n choose k) = n! / (k! (n-k)!) */
double log_binomial_coefficient(int n, int k)
{
  if (n < 0 || k < 0 || k > n)
    return log(0.0);
  if (k > n/2)
    k = n - k;
  double log_bc = 0.0;
  for (int i = 1; i <= k; ++i) {
    log_bc += log(n-i+1) - log(i);
  }
  return log_bc;
}


/* approximate the natural logarithm of the binomial
 * coefficient (n choose k) = n!/ (k! (n-k!))
 * using Stirling's approximation */
double log_binomial_coefficient_stirling(int n, int k)
{
  if (n < 0 || k < 0 || k > n)
    return log(0.0);
  else if (n == 0 || k == 0 || n == k)
    return log(1.0);
  const double pi = 3.1415926535897932384626433832795;
  return - 0.5*log(2*pi) - (k+0.5)*log(k) - (n-k+0.5)*log(n-k) + (n+0.5)*log(n);
}


/* compute the natural logarithm of the binomial probability
 * B(n;p)(k) = n! / (k! (n-k)!) p^k (1-p)^(n-k) */
double log_binomial_distribution(int n, double p, int k)
{
  if (n < 0 || k < 0 || k > n || p < 0.0 || p > 1.0)
    return log(0.0);
  double log_bp = k * log(p) + (n-k) * log(1.0-p);
  if (k > n/2)
    k = n - k;
  for (int i = 1; i <= k; ++i) {
    log_bp += log(n-i+1) - log(i);
  }
  return log_bp;
}


/* compute the logarithm of the cumulative binomial probability
 * Pr(B(n;p)(X) <= k) = \sum_{i=0}^{k} B(n;p)(i)*/
double log_cumulative_binomial(int n, double p, int k)
{
  if (n < 0 || k < 0 || k > n || p < 0.0 || p > 1.0)
    return log(0.0);
  const double log_success = log(p), log_fail = log(1.0-p);
  double log_term = n * log_fail;
  double log_cumulative = log_term;
  for (int i = 1; i <= k; ++i) {
    log_term = log_term + log(n-i+1) - log(i);
    log_term = log_term - log_fail + log_success;
    log_cumulative = log(exp(log_cumulative-log_term) + 1.0) + log_term;
  }
  return log_cumulative;
}
