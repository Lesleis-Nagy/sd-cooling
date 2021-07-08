//
// Created by L. Nagy on 22/10/2020.
//

#ifndef SD_COOLING_BASIC_TYPES_HPP
#define SD_COOLING_BASIC_TYPES_HPP

#include <complex>

#include <Eigen/Dense>

typedef size_t Index;

typedef size_t ISize;

typedef long double Real;

typedef double RSize;

typedef std::complex<Real> Complex;

typedef Eigen::Matrix<Complex, 4, 4> ComplexMatrix4x4;

typedef Eigen::Matrix<Real, 3, 3> RealMatrix3x3;

typedef Eigen::Matrix<Real, 2, 2> RealMatrix2x2;

typedef Eigen::Matrix<Real, 2, 1> Vector2D;

typedef Eigen::Matrix<Real, 3, 1> Vector3D;

#endif //SD_COOLING_BASIC_TYPES_HPP
