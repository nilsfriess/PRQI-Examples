#pragma once

#include <blaze/Blaze.h>

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

using namespace std::complex_literals;

struct RQI {
  using complex = std::complex<double>;
  using Matrix = blaze::DynamicMatrix<complex>;
  using Vector = blaze::DynamicVector<complex>;

  RQI(const Matrix &mat, double tolerance = 10e-8)
      : A(mat), identity(A.rows()), tolerance(tolerance) {}

  inline double sqrResidual() {
    return blaze::sqrNorm((A - lambda * identity) * x).real();
  }

  std::tuple<double, Vector, size_t> iterate(const Vector &initVec,
                                             bool debug = false) {
    x = blaze::normalize(initVec);
    lambda = blaze::ctrans(x) * A * x;

    size_t it = 0;
    while (std::sqrt(sqrResidual()) >= tolerance && it < 100) {
      try {
        x = blaze::solve(A - lambda * identity, x);
      } catch (const std::runtime_error &e) {
        // If system is singular, we're done
        return {lambda.real(), x, it};
      }

      x = blaze::normalize(x);
      lambda = blaze::ctrans(x) * A * x;

      if (debug) {
        std::cout << std::setprecision(4);
        std::cout << "Iteration: " << std::setw(3) << (it + 1);
        std::cout << ", res^2 = " << std::scientific << std::setw(10)
                  << sqrResidual();
        auto eigval = 1. / (blaze::ctrans(x) * x) * (blaze::ctrans(x) * A * x);
        std::cout << ", lambda = " << std::fixed << std::setprecision(8)
                  << std::setw(9) << eigval.real() << "\n";
      }

      ++it;
    }

    return {lambda.real(), x, it};
  }

  const Matrix A;
  const blaze::IdentityMatrix<double> identity;
  Vector x;
  complex lambda;
  double tolerance;
  bool classic;
};
