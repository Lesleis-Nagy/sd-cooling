//
// Created by L. Nagy on 21/10/2020.
//

#ifndef SD_COOLING_UTILITIES_HPP
#define SD_COOLING_UTILITIES_HPP

#include <array>
#include <locale>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

#include <basic_types.hpp>
#include <constants.hpp>

/**
 * Sign function will return the sign of the input value
 *
 * @tparam Real the basic type used for real values.
 * @param val the value for which the sign is required.
 * @return -1 if the input is negative, 1 if the value is positive and 0 if the value is zero.
 */
template<typename Real>
int
sgn(Real val) {
  return (Real(0) < val) - (val < Real(0));
}

/**
 * Sign change function, tests whether an input function will change sign when evaluated at two points of an interval.
 *
 * @tparam Functor the function type - should be a single parameter scalar value function.
 * @tparam Real the basic type used for real values.
 * @param fun the function for which root bins should be calculated.
 * @param x0 the first point of the interval.
 * @param x1 the second point of the interval.
 * @return true if the sign of fun(x0) and fun(x1) differ, otherwise false.
 */
template<typename Function, typename Real>
bool
sign_change(Function fun, Real x0, Real x1) {
  if (sgn(fun(x0)) == sgn(fun(x1))) {
	return false;
  }
  return true;
}

/**
 * Calculate the iteration index start and end for an array of length n, divided among np processes/threads.
 *
 * @param n length or array / number of data items.
 * @param np numper of processes / threads
 * @return an array that gives the start/end indices of an array distributed over np processes.
 */
std::vector<std::array<size_t, 2> >
compute_iter_start_end(size_t n, size_t np) {
  using namespace std;

  // For each thread p, calculate the number of iterations and offsets.
  vector<array<size_t, 2> > start_end(np);

  // First calculate the number of iterations per thread.
  size_t n_per_p = (size_t)n / np;

  // The remaining/outstanding iterations
  size_t rem_n = n % np;

  // Deal with the first value.
  if (rem_n != 0) {
	// If there are any remaining iterations make sure to add an additional iteration for this process/thread.
	start_end[0] = {0, (n_per_p + 1)};
	// Count down the number of iterations remaining by one.
	rem_n = rem_n - 1;
  } else {
	start_end[0] = {0, n_per_p};
  }

  // Deal with the other values.
  for (size_t p = 1; p < np; ++p) {
	if (rem_n != 0) {
	  // If there are any remaining iterations make sure to add an additional iteration for this process/thread.
	  start_end[p] = {start_end[p - 1][1], start_end[p - 1][1] + n_per_p + 1};
	  // Count down the number of iterations remaining by one.
	  rem_n = rem_n - 1;
	} else {
	  start_end[p] = {start_end[p - 1][1], start_end[p - 1][1] + n_per_p};
	}
  }

  return start_end;
}

/**
 * Convert the input string to a lower case version.
 * @param str the input string.
 * @return a lower case version of the input string.
 */
std::string
to_lower_case(const std::string &str) {
  using namespace std;

  locale loc;
  string str_lower;
  for (auto elem : str) {
	str_lower.push_back(tolower(elem, loc));
  }
  return str_lower;
}

/**
 * Return the three dimensional rotation matrix.
 * @param ax the axis of rotation.
 * @param theta the amount to rotate by (in radians).
 * @return 3D rotation matrix.
 */
RealMatrix3x3
rotation_matrix(const Vector3D &ax, Real theta) {
  RealMatrix3x3 R;

  Vector3D nax = ax.normalized();

  // Diagonal
  R(0, 0) = cos(theta) + nax(0) * nax(0) * (1 - cos(theta));
  R(1, 1) = cos(theta) + nax(1) * nax(1) * (1 - cos(theta));
  R(2, 2) = cos(theta) + nax(2) * nax(2) * (1 - cos(theta));

  // Off diagonal (orthogonal) elements.
  R(0, 1) = nax(0) * nax(1) * (1 - cos(theta)) - nax(2) * sin(theta);
  R(1, 0) = nax(1) * nax(0) * (1 - cos(theta)) + nax(2) * sin(theta);

  R(0, 2) = nax(0) * nax(2) * (1 - cos(theta)) + nax(1) * sin(theta);
  R(2, 0) = nax(2) * nax(0) * (1 - cos(theta)) - nax(1) * sin(theta);

  R(1, 2) = nax(1) * nax(2) * (1 - cos(theta)) - nax(0) * sin(theta);
  R(2, 1) = nax(2) * nax(1) * (1 - cos(theta)) + nax(0) * sin(theta);

  return R;
}

/**
 * Sort a list of values by nearness to a key value, NOTE: `lst` is sorted in place.
 * @param lst the list of values to sort in order of nearness to the `key_value`.
 * @param key the key value, values in `lst` are sorted relative to how close they are to this value.
 */
void nearness_sort(std::vector<Real> &lst, Real key) {
  using namespace std;

  sort(lst.begin(), lst.end(),
	   [key](const Real &theta1, const Real &theta2) -> bool {
		 if (abs(theta1 - key) < abs(theta2 - key)) {
		   return true;
		 } else {
		   return false;
		 }
	   }
  );
}

/**
 * Enumeration to hold the three ways of projecting a vector a on to a vector u. These are 1) not projecting
 * a on to u, 2) the parallel component of projecting a on to u and 3) the perpendicular component resulting
 * from the projection of a on to u.
 */
enum ProjectionType {
  NONE,
  PARALLEL,
  PERPENDICULAR
};

/**
 * Operation to find the component of the vector v projected on to the vector u (i.e. proj_u[v])
 * @param u the vector that the projections are taken with respect to.
 * @param v the vector being projected onto u.
 * @param ptype the projection type: PARALLEL is the component of v parallel to u, PERPENDICULAR is the
 *           component of v perpendicular to u, NONE: no projection, i.e. v is returned.
 * @return the component of the vector v in the direction of u.
 */
Vector3D proj(const Vector3D &u,
			  const Vector3D &v,
			  ProjectionType ptype = NONE) {
  switch (ptype) {
  case NONE:return v;
  case PARALLEL: {
	Vector3D u_normed = u.normalized();
	return u_normed.dot(v) * u_normed;
  }
  case PERPENDICULAR: {
	Vector3D u_normed = u.normalized();
	Vector3D para = u_normed.dot(v) * u_normed;
	return Vector3D(
		v[0] - para[0],
		v[1] - para[1],
		v[2] - para[2]
	);
  }
  default:return v;
  }
}

/**
 * Operation to find scalar component of the vector v projected on to the
 * vector u, if no projection is defined, return the magnitude of v.
 */
Real scalar_proj(const Vector3D &u,
				 const Vector3D &v,
				 ProjectionType ptype = PARALLEL) {
  switch (ptype) {
  case PARALLEL: {
	Vector3D u_normed = u.normalized();
	return u_normed.dot(v);
  }
  case PERPENDICULAR: {
	auto u_norm = u.norm();
	auto v_norm = v.norm();
	Real theta = acos(u.dot(v) / (u_norm * v_norm));
	return sin(theta) * v_norm;
  }
  default: return v.norm();
  }
}

/**
 * Split a line of comma separated string in to a vector of strings.
 */
std::vector<std::string>
str_split(const std::string &str, const std::string &delim) {
  using namespace std;

  vector<string> tokens;
  size_t prev = 0, pos = 0;
  do {
	pos = str.find(delim, prev);
	if (pos == string::npos)
	  pos = str.length();
	string token = str.substr(prev, pos - prev);
	if (!token.empty())
	  tokens.push_back(token);
	prev = pos + delim.length();
  } while (pos < str.length() && prev<str.length());
  return tokens;
}

/**
 * Split a line of comma separated numbers in to a vector of numbers.
 */
std::vector<Real>
str_split_to_number(const std::string &line) {
  using namespace std;
  vector<Real> numbers;
  auto tokens = str_split(line, ",");
  for (const auto &tok : tokens) {
	numbers.push_back(boost::lexical_cast<Real>(tok));
  }
  return numbers;
}

Real convert_angle_range(Real theta) {
  return (2 * pi + theta) * (theta < 0) + theta * (theta >= 0);
}

#endif //SD_COOLING_UTILITIES_HPP
