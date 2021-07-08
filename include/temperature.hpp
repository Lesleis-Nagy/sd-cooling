//
// Created by L. Nagy on 30/10/2020.
//

#ifndef SD_COOLING_TEMPERATURE_HPP
#define SD_COOLING_TEMPERATURE_HPP

#include <cmath>
#include <exception>
#include <functional>
#include <locale>
#include <string>
#include <iostream>

#include "basic_types.hpp"
#include "utilities.hpp"

/**
 * This exception class is thrown if a user specifies an unknown type of temperature regime.
 */
class UnknownTemperatureRegime: public std::exception
{
    [[nodiscard]] const char* what() const noexcept override {
        return "An unknown temperature regime was specified.";
    }
};

/**
 * Enumeration to hold types of temperature regime.
 */
enum TemperatureRegime {
    NEWTON,
    LINEAR
};


/**
 * Retrieves the TemperatureRegime enumeration value corresponding to an input string.
 * @param name the input string temperature regime.
 * @return an enumeration corresponding to the remperature regime name input.
 */
TemperatureRegime
name_to_temperature_regime(const std::string & name) {
    using namespace std;

    string lower_name = to_lower_case(name);

    if (lower_name == "newton") {
        return NEWTON;
    }

    if (lower_name == "linear") {
        return LINEAR;
    }

    throw UnknownTemperatureRegime();
}

std::function<Real(Real)>
linear_temperature_function(Real t0, Real T0, Real t1, Real T1)
{
    Real m = (T1 - T0) / (t1 - t0);
    return [t0, T0, m] (Real t) -> Real {
        return m*(t - t0) + T0;
    };
}

std::function<Real(Real)>
d_linear_temperature_function_dt(Real t0, Real T0, Real t1, Real T1)
{
  Real m = (T1 - T0) / (t1 - t0);
  return [m] (Real t) -> Real {
    return m;
  };
}

/**
 * Return a function that will compute a linear temperature curve.
 * @param T0 a known temperature at t0.
 * @param dT the rate at which the temperature changes.
 * @param t0 the start time of the temperature curve (default t0 = 0.0s).
 * @return a function that gives the temperature as a function of time for a linear
 *         temperature curve.
 */
std::function<Real(Real)>
linear_temperature_function(Real T0, Real dT, Real t0 = 0.0)
{
    return [T0, t0, dT] (Real t) -> Real {
        return T0 + dT*(t - t0);
    };
}

/**
 * Returns function that will compute a Newtonian temperature curve.
 * @param Tamb The ambient temperature.
 * @param T0 The initial temperature of a sample.
 * @param T1 The temperature of a sample after an elapsed time of t1.
 * @param t1 The amount of time elapsed from when the sample was at T0 to T1.
 * @return a function that gives the temperature as a function of time for a Newtonian
 *         temperature curve.
 */
std::function<Real(Real)>
newtonian_temperature_function(Real Tamb, Real T0, Real T1, Real t1)
{
    using namespace std;

    Real k = (1.0/t1) * log((T0 - Tamb)/(T1 - Tamb));

    return [Tamb, T0, k] (Real t) -> Real {
        return Tamb + (T0 - Tamb) * exp(-1.0 * k * t);
    };
}

/**
 * Returns a function that will compute the first derivative of a Newtonian
 * temperature curve.
 * @param Tamb The ambient temperature.
 * @param T0 The initial temperature of a sample.
 * @param T1 The temperature of a sample after an elapsed time of t1.
 * @param t1 The amount of time elapsed from when the sample was at T0 to T1.
 * @return a function that gives the temperature rate of change as a function
 * 		   of time for a Newtonian temperature curve.
 */
std::function<Real(Real)>
d_newtonian_temperature_function_dt(Real Tamb, Real T0, Real T1, Real t1)
{
  using namespace std;
  using namespace std;

  Real k = (1.0/t1) * log((T0 - Tamb)/(T1 - Tamb));

  return [Tamb, T0, k] (Real t) -> Real {
	return -1.0 * k * (T0 - Tamb) * exp(-1.0 * k * t);
  };
}

/**
 * Returns function that will compute a Newtonian temperature curve.
 * @param Tamb The ambient temperature.
 * @param T0 The initial temperature of a sample.
 * @param k The the cooling rate.
 * @return a function that gives the temperature as a function of time for a Newtonian
 *         temperature curve.
 */
std::function<Real(Real)>
newtonian_temperature_function(Real Tamb, Real T0, Real k)
{
    return [Tamb, T0, k] (Real t) -> Real {
        return Tamb + (T0 - Tamb) * exp(-1.0 * t / k);
    };
}

/**
 * Returns a function that will compute the first derivative of a Newtonian
 * temperature curve.
 * @param Tamb The ambient temperature.
 * @param T0 The initial temperature of a sample.
 * @param k The the cooling rate.
 * @return a function that gives the temperature rate of change as a function
 * 		   of time for a Newtonian temperature curve.
 */
std::function<Real(Real)>
d_newtonian_temperature_function_dt(Real Tamb, Real T0, Real k)
{
  return [Tamb, T0, k] (Real t) -> Real {
	return -1.0 * k * (T0 - Tamb) * exp(-1.0 * k * t);
  };
}

/**
 * Returns function that will compute a Newtonian temperature curve.
 */
std::function<Real(Real)>
newtonian_temperature_function(Real Tamb, Real t0, Real T0, Real t1, Real T1) {
    return [Tamb, t0, T0, t1, T1] (Real t) -> Real {
        return (T0 - Tamb) * std::exp(
                1 / (t0 - t1) * std::log((T1 - Tamb) / (T0 - Tamb)) * (t0 - t)
        ) + Tamb;
    };
}

/**
 * Returns a function that will compute the first derivative of a Newtonian temperature curve.
 */
std::function<Real(Real)>
d_newtonian_temperature_function_dt(Real Tamb, Real t0, Real T0, Real t1, Real T1) {
    return [Tamb, t0, T0, t1, T1] (Real t) -> Real {
        return -(((T0 - Tamb) * pow((T1 - Tamb) / (T0 - Tamb),
                                  (-t + t0) / (t0 - t1)) *
                  log((T1 - Tamb) / (T0 - Tamb))) / (t0 - t1));
    };
}

#endif //SD_COOLING_TEMPERATURE_HPP
