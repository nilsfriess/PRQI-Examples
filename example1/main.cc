#include "prqi.hh"
#include "rqi.hh"
#include "simplex.hh"

#include <blaze/math/TransposeFlag.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/dense/StaticVector.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>
using namespace std::complex_literals;

constexpr double pi = M_PI;

blaze::StaticVector<double, 3> spherical_to_cartesian(double r, double theta,
                                                      double phi) {
  double x = r * std::cos(phi) * std::sin(theta);
  double y = r * std::sin(phi) * std::sin(theta);
  double z = r * std::cos(theta);

  return {x, y, z};
}

template <class Method>
void unit_sphere_test(size_t N, double g, const std::string &filename) {
  typename Method::Matrix mat(3, 3);
  mat(0, 0) = -1;
  mat(1, 1) = g;
  mat(2, 2) = 1;

  Method rqi(mat, 10e-12);

  auto thetas = blaze::linspace(N, 0., pi);
  auto phis = blaze::linspace(N, 0., pi);

  std::ofstream output;
  output.open(filename, std::ofstream::out | std::ofstream::trunc);

  for (auto theta : thetas) {
    for (auto phi : phis) {
      auto initialVec = spherical_to_cartesian(1, theta, phi);

      auto res = rqi.iterate(initialVec);

      auto t = simplexIntersection(initialVec);

      if (!t)
        continue;

      auto intersVec = *t * initialVec;
      auto vec2d = simplexTo2D(intersVec);

      output << std::setprecision(12);
      output << vec2d[0] << " ";
      output << vec2d[1] << " ";
      output << std::setprecision(4);
      output << std::get<0>(res) << "\n";
    }
  }

  output.close();
}

void test_problem_classic(double eigenvalue, size_t runs) {
  std::string filename = "classic_" + std::to_string(eigenvalue) + ".txt";
  unit_sphere_test<RQI>(runs, eigenvalue, filename);
}

void test_problem_complex(double eigenvalue, size_t runs) {
  std::string filename = "complex_" + std::to_string(eigenvalue) + ".txt";
  unit_sphere_test<PRQI>(runs, eigenvalue, filename);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Provide number of points as command line argument\n";
    return -1;
  }
  
  int runs = std::stoi(argv[1]);

  test_problem_classic(0.1, runs);
  test_problem_classic(0.98, runs);

  test_problem_complex(0.1, runs);
  test_problem_complex(0.98, runs);
}
