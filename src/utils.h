#pragma once

#include <chrono>
#include <cmath>
#include <functional>
#include <vector>

// Util to benchmark code
struct benchmark {
  using clock = std::chrono::high_resolution_clock;
  clock::time_point start_point;

  // Reset the chronograph
  void reset() {
    start_point = clock::now();
  }

  // Return elapsed time in seconds
  double elapsed() {
    return std::chrono::duration<double>(clock::now() - start_point).count();
  }
};

// Type agnostic polynomial evaluation with Horner's rule.
//It's assumed that the coefficients are in increasing order.
template <typename R, typename T>
inline R polyeval(const T& c, const size_t order, const R x) {
  R result = c[order];
  for(int n = order - 1; n >=0; --n)
    //result = std::fma(result, x, c[n]);
    result = result * x + c[n];
  return result;
}

// ditto for a vector struct implementing ::size().
template <typename R, typename T>
inline R polyeval(const T& c, const R x) {
  return polyeval(c, c.size() - 1, x);
}

// Square a polynomial.
// Note that c2 must be twice the size of c.
template <typename T, typename R>
void polysquare(const T& c, const size_t order, R& c2) {
  const size_t n = order + 1;
  // Initialize to zero
  for (size_t i = 0; i < n*2; ++i)
    c2[i] = 0;
  // Compute square
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j)
      c2[j+i] += c[i]*c[j];
  }
}

// Simultaneously evaluate a polynomial and its derivatives.
// It returns pd = {c(x), dc(x), ddc(x), ...} upto the size of pd.
template <typename T>
void ddpoly(const T &c, const double x, std::vector<double> &pd) {
  const int nc = c.size() - 1;
  const int nd = pd.size() - 1;
  pd[0] = c[nc];
  for (int i = 1; i < nd; i++)
    pd[i] = 0.0;
  for (int i = nc - 1; i >= 0; --i) {
    int nnd = (nd < (nc - i) ? nd : nc-i);
    for (int j = nnd; j > 0; j--)
      pd[j] = std::fma(pd[j], x, pd[j-1]);
    pd[0] = std::fma(pd[0], x, c[i]);
  }
  double k = 2.0;
  for (int i = 2; i < nd+1; i++, k *= i)
    pd[i] *= k;
}

template <class T>
inline int sign(const T& z) {
  return (T(0) < z) - (z < T(0));
}

// Compute the crosstrack error to a polynomial curve of given order.
// In car coordinates, the CTE is the shortest distance of the polynomial
// to the origin. We compute it minimizing the squared distance by finding
// the root of the derivative closest to the initial point `x0`.
template <typename T>
double polyCTE(const T& c,
               const size_t order,
               double x0=0.0,
               const double tol=1e-4,
               const double eps=1e-14,
               const size_t max_iterations=20) {

  const size_t n = order + 1;

  // squared dist = x^2 + poly(c)^2
  std::vector<double> dist(n*2);
  polysquare(c, order, dist);
  dist[2] += 1;

  // This is the Newton algorithm to find the root of derivative of
  // the polynomial. Note that we need to compute the first and second
  // derivatives.
  std::vector<double> pd(3);
  for(size_t i = 0; i < max_iterations; i++) {
    // pd[1] is first derivative
    // pd[2] is second derivative
    ddpoly(dist, x0, pd);

    // Check division by 0.
    if (std::fabs(pd[2]) < eps)
      break;

    auto dx = pd[1] / pd[2];
    x0 -= dx;

    if (fabs(dx) <= tol * fabs(x0))
      break;
  }

  // It could be that the algorithm doesn't converge, for simplicity, I'm
  // going to ignore that case and work with the current estimate.
  return sign(polyeval(c, x0)) * std::sqrt(polyeval(dist, x0));
}
