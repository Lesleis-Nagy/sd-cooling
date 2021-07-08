//
// Created by L. Nagy on 26/10/2020.
//

#ifndef SD_COOLING_CONSTANTS_HPP
#define SD_COOLING_CONSTANTS_HPP

#include <cmath>

#include "basic_types.hpp"

/*
 * Mathematical constants.
 */

// Pi
const Real pi = acos(-1.0);

// The golden ratio.
const Real PHI = (1.0 + sqrt(5.0)) / 2.0;

// Square root of 2*pi
const Real sqrt2pi = sqrt(2.0 * acos(-1.0));

/*
 * Physics constants.
 */

// Permeability of free space
const Real mu0 = 4.0 * pi * 1E-7; // ( m kg ) / ( s^2 A^2 )

// Boltzmann's constant
const Real kb = 1.38064852E-23;   // ( m^2 kg ) / ( s K )

#endif //SD_COOLING_CONSTANTS_HPP
