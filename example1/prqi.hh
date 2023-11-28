#pragma once

#include <blaze/Blaze.h>

#include <complex>
#include <iomanip>
#include <iostream>

using namespace std::complex_literals;

struct PRQI {
  using complex = std::complex<double>;
  using Matrix = blaze::DynamicMatrix<complex>;
  using Vector = blaze::DynamicVector<complex>;

  PRQI(const Matrix &mat, double tolerance = 10e-8)
      : A(mat), identity(A.rows()), tolerance(tolerance) {}

  inline double residual(const Vector &x, double lambda) {
    return blaze::norm((A - lambda * identity) * x).real();
  }

  std::tuple<double, Vector, size_t> iterate(const Vector &initVec,
                                             bool debug = false) {
    Vector x = blaze::normalize(initVec);
    auto lambda = (blaze::ctrans(x) * A * x).real();

    auto gamma = residual(x, lambda);

    size_t it = 0;
    while (gamma >= tolerance && it < 100) {
      try {
        x = blaze::solve(A - (lambda + gamma * 1i) * identity, x);
      } catch (const std::runtime_error &e) {
        // If system is singular, solve last step
        x = blaze::normalize(blaze::real(x));
        return {lambda, x, it};
      }

      x = blaze::normalize(x);
      lambda = (blaze::ctrans(x) * A * x).real();
      gamma = residual(x, lambda);

      if (debug) {
        std::cout << std::setprecision(4);
        std::cout << "Iteration: " << std::setw(3) << (it + 1);
        std::cout << ", res^2 = " << std::scientific << std::setw(10) << gamma;
        auto eigval = 1. / (blaze::ctrans(x) * x) * (blaze::ctrans(x) * A * x);
        std::cout << ", lambda = " << std::fixed << std::setprecision(8)
                  << std::setw(9) << eigval.real() << "\n";
      }

      ++it;
    }

    x = blaze::normalize(blaze::real(x));
    return {lambda, x, it};
  }

  const Matrix A;
  const blaze::IdentityMatrix<double> identity;
  double tolerance;
  bool classic;
};
