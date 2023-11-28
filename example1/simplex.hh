#pragma once

#include <blaze/math/Vector.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/dense/StaticVector.h>
#include <blaze/math/expressions/DVecScalarMultExpr.h>
#include <iostream>
#include <optional>

#include <blaze/Blaze.h>

// Source:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution

std::optional<double>
simplexIntersection(const blaze::DynamicVector<double> &vec) {

  const blaze::DynamicVector<double> v0{1, 0, 0};
  const blaze::DynamicVector<double> v1{0, 1, 0};
  const blaze::DynamicVector<double> v2{0, 0, 1};

  blaze::StaticVector<double, 3> N{1, 1, 1}; // normal vector to simplex

  auto angleWithN = blaze::dot(N, vec);
  if (std::abs(angleWithN) < 10e-10)
    return {};

  auto d = -N[0];
  auto t = -d / angleWithN;

  if (t < 0)
    return {};

  // Compute intersection point
  auto P = t * vec;

  auto edge0 = v1 - v0;
  auto vp0 = P - v0;
  auto C0 = blaze::cross(edge0, vp0);

  if (blaze::dot(N, C0) < 0)
    return {};

  auto edge1 = v2 - v1;
  auto vp1 = P - v1;
  auto C1 = blaze::cross(edge1, vp1);

  if (blaze::dot(N, C1) < 0)
    return {};

  auto edge2 = v0 - v2;
  auto vp2 = P - v2;
  auto C2 = blaze::cross(edge2, vp2);

  if (blaze::dot(N, C2) < 0)
    return {};

  return t;
}

// Source:
// https://stackoverflow.com/questions/44553886/how-to-create-2d-plot-of-arbitrary-coplanar-3d-curve/44559920#44559920
blaze::StaticVector<double, 2>
simplexTo2D(const blaze::StaticVector<double, 3> &vec) {
  const blaze::DynamicVector<double> v0{1, 0, 0};
  const blaze::DynamicVector<double> v1{0, 1, 0};
  const blaze::DynamicVector<double> v2{0, 0, 1};

  blaze::DynamicVector<double> u = v1 - v0;
  blaze::DynamicVector<double> v = v2 - v0;

  u = blaze::normalize(u);
  v = blaze::normalize(v);

  auto w = blaze::cross(u, v);
  u = blaze::cross(v, w);

  auto xp = blaze::dot(u, vec);
  auto yp = blaze::dot(v, vec);

  return {xp, yp};
}
