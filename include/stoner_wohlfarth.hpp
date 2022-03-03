//
// Created by L. Nagy on 22/10/2020.
//

#ifndef SD_COOLING_STONER_WOHLFARTH_HPP
#define SD_COOLING_STONER_WOHLFARTH_HPP

#include <cmath>

#include <algorithm>
#include <functional>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <vector>
#include <unordered_map>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/eigen.hpp>

#include "constants.hpp"
#include "materials.hpp"
#include "grain_population.hpp"
#include "state_transition.hpp"
#include "formatter.hpp"
#include "debug.hpp"

/**
 * Long axis demagnetizing factor using Dunlop & Ozdemir's convention, (Note: Dunlop & Ozdemir denote the short axis `b`
 * and the long axis `a` [pp. 90]; Cullity & Graham denote the short axis `a` and the long axis `c` [pp. 54]).
 * @param m the aspect ratio of the prolate spheroid, m = a/b.
 * @return the long axis demagnetizing factor.
 */
inline Real
Na(Real m) {
    return (1.0 / (m * m - 1.0))
           * ((m / sqrt(m * m - 1.0)) * log(m + sqrt(m * m - 1.0)) - 1.0);
}

/**
 * Short axis demagnetizing factor using Dunlop & Ozdemir's convention, (Note: Dunlop & Ozdemir denote the short axis
 * `b` and the long axis `a` [pp. 90]; Cullity & Graham denote the short axis `a` and the long axis `c` [pp. 54]).
 * @param m the aspect ratio of the prolate spheroid, m = a/b.
 * @return the short axis demagnetizing factor.
 */
inline Real
Nb(Real m) {
    return 0.5 * (1.0 - Na(m));
}

/**
 * Return the C1 constant for the Stoner-Wohlfarth energy functional, this is the constant associated with the
 * demagnetizing part.
 * @param esvd equivalent spherical volume diameter [m].
 * @param elong elongation as a percentage (100% elongation is twice length to height) [%].
 * @param material the name of the material for which C1 is required.
 * @return the first S-W energy functional constant, as a function of temperature [J].
 */
std::function<Real(Real)>
C1_function(Real esvd, Real elong, const std::string &material) {
    Real m = (100.0 + elong) / 100.0;
    Real v = (1.0 / 6.0) * pi * esvd * esvd * esvd;
    Real C = 0.5 * (Nb(m) - Na(m)) * v * mu0;
    auto Ms = saturation_magnetization_function(name_to_material(material));

    return [C, Ms](Real T) -> Real {
        return C * Ms(T) * Ms(T);
    };
}

/**
 * Return the volume normalized C1 constant for the Stoner-Wohlfarth energy functional, this is the constant associated
 * with the demagnetizing part.
 * @param elong elongation as a percentage (100% elongation is twice length to height) [%].
 * @param material the name of the material for which C1 is required.
 * @return the first S-W energy functional constant, as a function of temperature [J].
 */
std::function<Real(Real)>
C1VN_function(Real elong, const std::string &material) {
    Real m = (100.0 + elong) / 100.0;
    Real C = 0.5 * (Nb(m) - Na(m)) * mu0;
    auto Ms = saturation_magnetization_function(name_to_material(material));

    return [C, Ms](Real T) -> Real {
        return C * Ms(T) * Ms(T);
    };
}

/**
 * Return the C2 constant for the Stoner-Wohlfarth energy functional, this is the constant associated with the
 * parallel component of the applied field part.
 * @param esvd equivalent spherical volume diameter [m].
 * @param H the strength of the applied external field [A/m].
 * @param phi the direction of the applied external field (relative to the axis of elongation) [rad].
 * @param material the name of the material for which C1 is required.
 * @return the second S-W energy functional constant, as a function of temperature.
 */
std::function<Real(Real)>
C2_function(Real esvd, Real H, Real phi, const std::string &material) {
    Real v = (1.0 / 6.0) * pi * esvd * esvd * esvd;
    Real C = H * v * mu0 * cos(phi);
    auto Ms = saturation_magnetization_function(name_to_material(material));

    return [C, Ms](Real T) -> Real {
        return C * Ms(T);
    };
}

/**
 * Return the volume normalized C2 constant for the Stoner-Wohlfarth energy functional, this is the constant associated
 * with the parallel component of the applied field part.
 * @param H the strength of the applied external field [A/m].
 * @param phi the direction of the applied external field (relative to the axis of elongation) [rad].
 * @param material the name of the material for which C1 is required.
 * @return the second S-W energy functional constant, as a function of temperature.
 */
std::function<Real(Real)>
C2VN_function(Real H, Real phi, const std::string &material) {
    Real C = H * mu0 * cos(phi);
    auto Ms = saturation_magnetization_function(name_to_material(material));

    return [C, Ms](Real T) -> Real {
        return C * Ms(T);
    };
}

/**
 * Return the C3 constant for the Stoner-Wohlfarth energy functional, this is the constant associated with the
 * perpendicular component of the applied field part.
 * @param esvd equivalent spherical volume diameter [m].
 * @param H the strength of the applied external field [A/m].
 * @param phi the direction of the applied external field (relative to the axis of elongation) [rad].
 * @param material the name of the material for which C1 is required.
 * @return the third S-W energy functional constant, as a function of temperature.
 */
std::function<Real(Real)>
C3_function(Real esvd, Real H, Real phi, const std::string &material) {
    Real v = (1.0 / 6.0) * pi * esvd * esvd * esvd;
    Real C = H * v * mu0 * sin(phi);
    auto Ms = saturation_magnetization_function(name_to_material(material));

    return [C, Ms](Real T) -> Real {
        return C * Ms(T);
    };
}

/**
 * Return the volume normalized C3 constant for the Stoner-Wohlfarth energy functional, this is the constant associated
 * with the perpendicular component of the applied field part.
 * @param H the strength of the applied external field [A/m].
 * @param phi the direction of the applied external field (relative to the axis of elongation) [rad].
 * @param material the name of the material for which C1 is required.
 * @return the third S-W energy functional constant, as a function of temperature.
 */
std::function<Real(Real)>
C3VN_function(Real H, Real phi, const std::string &material) {
    Real C = H * mu0 * sin(phi);
    auto Ms = saturation_magnetization_function(name_to_material(material));

    return [C, Ms](Real T) -> Real {
        return C * Ms(T);
    };
}

std::function<Complex(Real)>
cc1_function(Real elong, Real H, Real phi, const std::string &material) {
    Real m = (100.0 + elong) / 100.0;
    Real C = (2.0 * H) / ((Nb(m) - Na(m)));

    auto Ms = saturation_magnetization_function(name_to_material(material));

    return [C, phi, Ms](Real T) -> Complex {
        return Complex(C * cos(phi) / Ms(T), C * sin(phi) / Ms(T));
    };
}

std::function<Complex(Real)>
cc2_function(Real elong, Real H, Real phi, const std::string &material) {
    Real m = (100.0 + elong) / 100.0;
    Real C = (2.0 * H) / ((Nb(m) - Na(m)));

    auto Ms = saturation_magnetization_function(name_to_material(material));

    return [C, phi, Ms](Real T) -> Complex {
        return Complex(C * cos(phi) / Ms(T), (-1.0) * C * sin(phi) / Ms(T));
    };
}

/**
 * Return the Stoner-Wohlfarth energy functional for the given set of input parameters describing a grain.
 * @param esvd equivalent spherical volume diameter [m].
 * @param elong the elongation of the ellipsoid as a percentage [%].
 * @param H the strength of the applied field [A/m].
 * @param phi the direction of the applied field [rad].
 * @param material the material that the grain is made of.
 * @return the S-W energy functional (as a function of magnetization direction, theta and temperature in centigrade) for
 *         the grain described by the input parameters
 */
std::function<Real(Real, Real)>
stoner_wohlfarth_energy_function(Real esvd,
                                 Real elong,
                                 Real H,
                                 Real phi,
                                 const std::string &material) {
    auto C1 = C1_function(esvd, elong, material);
    auto C2 = C2_function(esvd, H, phi, material);
    auto C3 = C3_function(esvd, H, phi, material);

    return [C1, C2, C3](Real theta, Real T) -> Real {
        return C1(T) * sin(theta) * sin(theta) - C2(T) * cos(theta)
               - C3(T) * sin(theta);
    };
}

/**
 * Return the Stoner-Wohlfarth energy density functional for the given set of input parameters describing a grain.
 * @param elong the elongation of the ellipsoid as a percentage [%].
 * @param H the strength of the applied field [A/m].
 * @param phi the direction of the applied field [rad].
 * @param material the material that the grain is made of.
 * @return the S-W energy density functional (as a function of magnetization direction, theta and temperature in
 *         centigrade) for the grain described by the input parameters
 */
std::function<Real(Real, Real)>
stoner_wohlfarth_energy_density_function(Real elong,
                                         Real H,
                                         Real phi,
                                         const std::string &material) {
    auto C1VN = C1VN_function(elong, material);
    auto C2VN = C2VN_function(H, phi, material);
    auto C3VN = C3VN_function(H, phi, material);

    return [C1VN, C2VN, C3VN](Real theta, Real T) -> Real {
        return C1VN(T) * sin(theta) * sin(theta) - C2VN(T) * cos(theta)
               - C3VN(T) * sin(theta);
    };
}

/**
 * Return the first derivative of the Stoner-Wohlfarth energy functional for the given set of input parameters
 * describing a grain.
 * @param esvd equivalent spherical volume diameter [m].
 * @param elong the elongation of the ellipsoid as a percentage [%].
 * @param H the strength of the applied field [A/m].
 * @param phi the direction of the applied field [rad].
 * @param material the material that the grain is made of.
 * @return the first derivative S-W energy functional (as a function of magnetization direction, theta and temperature
 *         in centigrade) for the grain described by the input parameters
 */
std::function<Real(Real, Real)>
d_stoner_wohlfarth_energy_function(Real esvd,
                                   Real elong,
                                   Real H,
                                   Real phi,
                                   const std::string &material) {
    auto C1 = C1_function(esvd, elong, material);
    auto C2 = C2_function(esvd, H, phi, material);
    auto C3 = C3_function(esvd, H, phi, material);

    return [C1, C2, C3](Real theta, Real T) -> Real {
        return 2.0 * C1(T) * sin(theta) * cos(theta) + C2(T) * sin(theta)
               - C3(T) * cos(theta);
    };
}

/**
 * Return the first derivative of the Stoner-Wohlfarth energy-density functional for the given set of input parameters
 * describing a grain.
 * @param elong the elongation of the ellipsoid as a percentage [%].
 * @param H the strength of the applied field [A/m].
 * @param phi the direction of the applied field [rad].
 * @param material the material that the grain is made of.
 * @return the first derivative S-W energy-density functional (as a function of magnetization direction, theta and
 *         temperature in centigrade) for the grain described by the input parameters
 */
std::function<Real(Real, Real)>
d_stoner_wohlfarth_energy_density_function(Real elong,
                                           Real H,
                                           Real phi,
                                           const std::string &material) {
    auto C1VN = C1VN_function(elong, material);
    auto C2VN = C2VN_function(H, phi, material);
    auto C3VN = C3VN_function(H, phi, material);

    return [C1VN, C2VN, C3VN](Real theta, Real T) -> Real {
        return 2.0 * C1VN(T) * sin(theta) * cos(theta) + C2VN(T) * sin(theta)
               - C3VN(T) * cos(theta);
    };
}

/**
 * Return the second derivative of the Stoner-Wohlfarth energy functional for the given set of input parameters
 * describing a grain.
 * @param esvd equivalent spherical volume diameter [m].
 * @param elong the elongation of the ellipsoid as a percentage [%].
 * @param H the strength of the applied field [A/m].
 * @param phi the direction of the applied field [rad].
 * @param material the material that the grain is made of.
 * @return the second derivative S-W energy functional (as a function of magnetization direction, theta and temperature
 *         in centigrade) for the grain described by the input parameters
 */
std::function<Real(Real, Real)>
dd_stoner_wohlfarth_energy_function(Real esvd,
                                    Real elong,
                                    Real H,
                                    Real phi,
                                    const std::string &material) {
    auto C1 = C1_function(esvd, elong, material);
    auto C2 = C2_function(esvd, H, phi, material);
    auto C3 = C3_function(esvd, H, phi, material);

    return [C1, C2, C3](Real theta, Real T) -> Real {
        return 2.0 * C1(T) * cos(theta) * cos(theta)
               - 2.0 * C1(T) * sin(theta) * sin(theta) + C2(T) * cos(theta)
               + C3(T) * sin(theta);
    };
}

/**
 * Return the second derivative of the Stoner-Wohlfarth energy density functional for the given set of input parameters
 * describing a grain.
 * @param elong the elongation of the ellipsoid as a percentage [%].
 * @param H the strength of the applied field [A/m].
 * @param phi the direction of the applied field [rad].
 * @param material the material that the grain is made of.
 * @return the second derivative S-W energy density functional (as a function of magnetization direction, theta and
 *         temperature in centigrade) for the grain described by the input parameters
 */
std::function<Real(Real, Real)>
dd_stoner_wohlfarth_energy_density_function(Real elong,
                                            Real H,
                                            Real phi,
                                            const std::string &material) {
    auto C1VN = C1VN_function(elong, material);
    auto C2VN = C2VN_function(H, phi, material);
    auto C3VN = C3VN_function(H, phi, material);

    return [C1VN, C2VN, C3VN](Real theta, Real T) -> Real {
        return 2.0 * C1VN(T) * cos(theta) * cos(theta)
               - 2.0 * C1VN(T) * sin(theta) * sin(theta) + C2VN(T) * cos(theta)
               + C3VN(T) * sin(theta);
    };
}

/**
 * This exception class is thrown if a grain population is invalid.
 */
class InvalidGrainPopulation : public std::exception {
public:
    InvalidGrainPopulation() = default;

    explicit InvalidGrainPopulation(const std::string &msg) : _msg(msg) {}

    [[nodiscard]] const char *what() const noexcept override {
        return _msg.c_str();
    }

private:
    std::string _msg;
};

/**
 * Class to encapsulate a population of Stoner-Wohlfarth particles.
 */
class StonerWohlfarthGrainPopulation : public GrainPopulation {
public:
    /**
     * Default constructor.
     * @param esvd the equivalent spherical volume diameter of the grain population in meters.
     * @param elong the elongation of the particle as a %-age
     * @param H the magnitude of the applied field in A/m
     * @param hx the x component of the direction of the applied field (will be normalized)
     * @param hy the y component of the direction fo the applied field (will be normalized)
     * @param hz the z component of the direction of the applied field (will be normalized)
     * @param ux the x component of the grain orientation (will be normalized)
     * @param uy the y component of the grain orientation (will be normalized)
     * @param uz the z component of the grain orientation (will be normalized)
     * @param str_material the name of the material, either "magnetite" or "iron"
     * @param initial_temperature the initial temperature of the grain population in degrees centigrade
     * @param rho the initial distribution of states (will be normalized to sum to 1)
     * @param tau0 the switching frequency in Hz (i.e. per second), see Dunlop & Ozdemir (2001), pp. 202 (eqn. 8.3)
     * @param eps the epsilon value used to represent the models 'zero' quantity
     * @param n_polish the number of Newton-Raphson polishing steps
     */
    StonerWohlfarthGrainPopulation(Real esvd,
                                   Real elong,
                                   Real H,
                                   Real hx, Real hy, Real hz,
                                   Real ux, Real uy, Real uz,
                                   const std::string &str_material,
                                   Real initial_temperature,
                                   const std::vector<Real> &rho,
                                   Real tau0,
                                   Real eps = 1E-12,
                                   ISize n_polish = 0) :
            _tau0(tau0),
            _eps(eps),
            _n_polish(n_polish) {

        // Set the volume.
        _set_volume(esvd);

        // Set the elongation.
        _set_elongation(elong);

        // Set the orientation.
        _set_orientation(ux, uy, uz);

        // Set the material.
        _set_material(str_material);

        // Set the field.
        _set_field(H, hx, hy, hz);

        // Set temperature dependent constant functions.
        _set_cc1_cc2();

        // Set temperature depenent Stoner-Wohlfarth functions and their 1st/2nd derivatives.
        _set_sw();

        // Set the initial temperature.
        _set_initial_temperature(initial_temperature);

        // Set the population proportions.
        _set_rho(rho);

        // Set the current temperature to the inital temperature
        _initialize_model();

        // Validate the object.
        _check_object();
    }

    void update_temperature(Real T) override {
        // Update the current temperature.
        _temperature = T;
        _compute_lems_and_energy_barriers();
    }

    void equilibrate() override {
        // Equilibrate the model, i.e. assume that the model sits at _temperature for an infinite amount of time and
        // update rho with those values.

        // Compute LEMs and energy barriers for current T_.
        _compute_lems_and_energy_barriers();

        // Retrieve a possible equilibrium rho.
        _compute_equilibrium_rho();

        // Use the existing equilibrium rho distribution for the population rho.
        _rho = _rho_eq;
    }

    void update_field(Real H) override {
        update_field(H, _hx, _hy, _hz);
    }

    void update_field(Real H, Real hx, Real hy, Real hz) override {
        // Set the field to the new values.
        _set_field(H, hx, hy, hz);

        // Set temperature dependent constant functions.
        _set_cc1_cc2();

        // Set temperature dependent Stoner-Wohlfarth functions and their 1st/2nd derivatives.
        _set_sw();

        // Compute the LEMs and energy barriers associated with the new values.
        _compute_lems_and_energy_barriers();
    }

    void update_field_and_temperature(Real H, Real T) override {
        update_field_and_temperature(H, _hx, _hy, _hz, T);
    }

    void update_field_and_temperature(Real H, Real hx, Real hy, Real hz, Real T) override {
        // Set the field to the new values.
        _set_field(H, hx, hy, hz);

        // Set the temperature to the new value.
        _temperature = T;

        // Set the temperature dependent constant functions.
        _set_cc1_cc2();

        // Set the temperature dependent Stoner-Wohlfarth functions and their 1st/2nd derivatives.
        _set_sw();

        // Compute the LEMs and energy barriers associated with the new values.
        _compute_lems_and_energy_barriers();
    }

    /**
     * Update the population partition as a function of time spent at the current population temperature.
     * @param dt
     */
    void update_rho(Real dt) override {
        // Update the population magnetization.
        auto P = compute_probability_transition_matrix(dt);

        DEBUG_MSG_SIMPLE_VAR(dt);

        // Save the current rho vector.
        _old_rho = _rho;

        // Update the current rho vector.
        _rho.clear();
        _rho.resize(_lem_states.size(), 0);

        if (_lem_states.size() == 1) {
            // If there is only one LEM state then the occupancy vector _rho is also a single value set to 1 which
            // represent the fact that all grains in a population occupy that state.
            _rho[0] = 1;
        } else {
            // Otherwise, we must update _rho since it
            for (size_t i = 0; i < P.size(); ++i) {
                for (size_t j = 0; j < P[0].size(); ++j) {
                    _rho[i] += _old_rho[j] * P[j][i];
                }
            }
        }

        DEBUG_MSG_STD_VECTOR(_old_rho);
        DEBUG_MSG_STD_VECTOR(_rho);

        // Update the equilibrium magnetization.
        _compute_equilibrium_rho();
    }

    /**
     * Compute the population magnetization vector without projection on to the applied field.
     * @param with_ms include the temperature dependent spontaneous magnetization.
     * @return the population magnetization vector.
     */
    [[nodiscard]] Vector3D population_magnetization_no_project(bool with_ms) const override {
        return _magnetization_no_project(with_ms, POPULATION_MAGNETIZATION);
    }

    /**
     * Compute the population magnetization.
     * @param with_ms include the temperature dependent spontaneous magnetization.
     * @return the scalar projection of the population magnetization on to the applied field.
     */
    [[nodiscard]] Real population_magnetization(bool with_ms) const override {
        return _magnetization(with_ms, POPULATION_MAGNETIZATION);
    }

    /**
     * Compute the equilibrium magnetization (i.e. t -> infity).
     * @param ptype the projection to use, NONE just returns the raw vector, PARALLEL is the magnetization component
     *              parallel to the applied field direction, and PERPENDICULAR is the component perpendicular to the
     *              applied field direction (assuming a right handed rule). By default this is NONE.
     * @param with_ms include the temperature depenent spontaneous magnetization. By default this is false.
     * @return a three dimensional magnetization vector.
     */
    [[nodiscard]] Real equilibrium_magnetization(bool with_ms) const override {
        return _magnetization(with_ms, EQUILIBRIUM_MAGNETIZATION);
    }

    /**
     * Retrieve the vector containing the fraction of grains in one state and the fraction of grains in the other.
     * @return a vector with the proportion of grains in state one and the proportion of grains in state two.
     */
    [[nodiscard]] std::vector<Real> get_rho() const override {
        return _rho;
    }

    /**
     * Retrieve the vector containing the fraction of grains in one state and the fraction of grains in the other at
     * equilibirum.
     * @return a vector with the proportion of grains in state one and the proportion of grains in state two.
     */
    [[nodiscard]] std::vector<Real> get_rho_eq() const override {
        return _rho;
    }

    /**
     * Retrieve energy barriers.
     * @return energy barriers.
     */
    [[nodiscard]] std::vector<std::unordered_map<Index, Real>>
    energy_barriers() const {
        return _energy_barriers;
    }

    /**
     * Retrieve LEM states.
     * @return LEM states.
     */
    [[nodiscard]] std::vector<Real> lem_states() const {
        return _lem_states;
    }

    /**
     * Compute the probability transition matrix based on the current state of the system, i.e. the current set of
     * LEM states and the energy barriers between the current set of LEM states.
     * @param t the length of time (in seconds) that the population has experienced the current field/temperature.
     * @return the probability transition matrix.
     */
    [[nodiscard]] std::vector<std::vector<Real> >
    compute_probability_transition_matrix(Real t) const {
        using namespace std;

        using MPReal = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<500>>;

        using MPMatrix = Eigen::Matrix<MPReal, Eigen::Dynamic, Eigen::Dynamic>;

        using Solver = Eigen::EigenSolver<MPMatrix>;

        DEBUG_MSG_SIMPLE_VAR(t);
        DEBUG_MSG_SIMPLE_VAR(_h_magnitude);
        DEBUG_MSG_SIMPLE_VAR(_phi);
        DEBUG_MSG_SIMPLE_VAR(_v);
        DEBUG_MSG_SIMPLE_VAR(_elong);
        DEBUG_MSG_SIMPLE_VAR(_temperature);

        DEBUG_MSG_STD_VECTOR(_lem_states);
        DEBUG_MSG_STD_VECTOR(_rho);

        DEBUG_MSG_STD_VECTOR_OF_MAPS(_energy_barriers);

        // Check the number of LEM states.
        if (_lem_states.size() == 1) {
            // If there is only one LEM state the probability transition matrix is
            // trivial.
            return {{1}};
        } else {

            // Construct the probability transition matrix.
            MPMatrix Q(_lem_states.size(), _lem_states.size());

            // The rows of Q should sum to zero.
            for (Index i = 0; i < _lem_states.size(); ++i) {
                MPReal off_diagonal_sum(0);
                for (Index j = 0; j < _lem_states.size(); ++j) {
                    if (i != j) {

                        MPReal one_over_tau(1.0 / _tau0);
                        DEBUG_MSG_SIMPLE_VAR(one_over_tau);

                        MPReal delta_E(_energy_barriers[i].at(j));
                        DEBUG_MSG_SIMPLE_VAR(delta_E);

                        MPReal kbT(kb * (_temperature + 273.15));
                        DEBUG_MSG_SIMPLE_VAR(kbT);

                        MPReal delta_t(t);
                        DEBUG_MSG_SIMPLE_VAR(delta_t);

                        MPReal matrix_entry = one_over_tau * boost::multiprecision::exp(-delta_E / kbT) * delta_t;

                        Q(i, j) = matrix_entry;

                        off_diagonal_sum += matrix_entry;
                        DEBUG_MSG_SIMPLE_VAR(off_diagonal_sum);
                    }
                }
                // Set the (i,i)th entry to the negative off-diagonal sum.
                Q(i, i) = -off_diagonal_sum;
            }
            DEBUG_MSG_EIGEN_MATRIX(Q);

            /////////////////////////////////////////////////////////////////////////
            /// Calculate matrix exponential using matrix diagonalization.        ///
            /////////////////////////////////////////////////////////////////////////

            Solver solver;
            solver.compute(Q);

            auto eigen_values = solver.eigenvalues();
            auto eigen_vectors = solver.eigenvectors();

            DEBUG_MSG_EIGEN_MATRIX(eigen_values);
            DEBUG_MSG_EIGEN_MATRIX(eigen_vectors);

            // Verify here that the no. of _states is equal to the number of
            // eigenvalues.
            if (_lem_states.size() != eigen_values.size()) {
                throw std::runtime_error(
                        Formatter() << "FATAL: No. of LEM states and eigenvalues does not match"
                );
            }

            ISize n = _lem_states.size();

            // Calculate the diagonal matrix of eigenvalues.
            MPMatrix expD(n, n);
            for (size_t i = 0; i < n; ++i) {
                expD(i, i) = boost::multiprecision::exp(eigen_values(i).real());
            }
            DEBUG_MSG_EIGEN_MATRIX(expD);

            // Calculate the eigenvector matrix.
            MPMatrix P(n, n);
            for (Index i = 0; i < n; ++i) {
                for (Index j = 0; j < n; ++j) {
                    P(i, j) = eigen_vectors(i, j).real();
                }
            }
            DEBUG_MSG_EIGEN_MATRIX(P);

            auto invP = P.inverse();
            DEBUG_MSG_EIGEN_MATRIX(invP);

            auto expQ = P * expD * invP;
            DEBUG_MSG_EIGEN_MATRIX(expQ);

            std::vector<std::vector<Real> > R(n);
            for (Index i = 0; i < n; ++i) {
                R[i].resize(n);
                for (Index j = 0; j < n; ++j) {
                    R[i][j] = (Real) expQ(i, j);
                }
            }
            DEBUG_MSG_STD_VECTOR_OF_VECTORS(R);

            return R;
        }
    }

private:

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Private member variables.                                                                                   ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Real _esvd{};
    Real _elong{};
    Real _h_magnitude{};
    Real _hx{};
    Real _hy{};
    Real _hz{};
    Real _ux{};
    Real _uy{};
    Real _uz{};
    std::string _str_material;
    Real _initial_temperature{};
    std::vector<Real> _rho;
    std::vector<Real> _old_rho;
    std::vector<Real> _rho_eq;
    Real _tau0;
    Real _eps;
    ISize _n_polish;

    // Derived constant values.
    Real _v;

    // Derived values

    // The angle at which the applied field is directed.
    Real _phi;

    // The material enumeration.
    Material _material;

    // A unit vector pointing in the direction of the applied field.
    Vector3D _h;

    // A unit vector pointing in the direction of the grain's principal axis.
    Vector3D _u;

    // Material dependent functions.
    std::function<Complex(Real)> _c1;
    std::function<Complex(Real)> _c2;

    std::function<Real(Real, Real)> _sw;
    std::function<Real(Real, Real)> _dsw;
    std::function<Real(Real, Real)> _ddsw;

    std::function<Real(Real)> _ms;

    // The temperature.
    Real _temperature;

    // States
    std::vector<Real> _old_lem_states;
    std::vector<Real> _lem_states;

    // Energy barriers.
    std::vector<std::unordered_map<Index, Real> > _old_energy_barriers;
    std::vector<std::unordered_map<Index, Real> > _energy_barriers;

    // Class constants
    const static int EQUILIBRIUM_MAGNETIZATION = 1;
    const static int POPULATION_MAGNETIZATION = 2;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Private member functions.                                                                                   ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void _check_object() const {
        // Epsilon must be non-zero
        if (_eps < 0.0) {
            throw InvalidGrainPopulation("Epsilon is non-zero.");
        }

        // Equivalent spherical volume diameter (ESVD) must be non-zero.
        if (_esvd < _eps) {
            throw InvalidGrainPopulation("ESVD was too small or negative.");
        }

        // Elongation must be non-zero.
        if (_elong < _eps) {
            throw InvalidGrainPopulation("Elongation was too small or negative");
        }
    }

    /**
     * Set the material name string and the material enumeration.
     * @param material the material name.
     */
    void _set_material(const std::string &material) {
        // Set the material string.
        _str_material = material;

        // Set the material enumeration.
        _material = name_to_material(_str_material);

        // Set saturation magnetization function.
        _ms = saturation_magnetization_function(name_to_material(_str_material));
    }

    /**
     * Set the grain populations volume by calculating the volume of a sphere from
     * the equivalent spherical volume diameter (esvd).
     * @param esvd equivalend spherical volume diameter (in meter).
     */
    void _set_volume(Real esvd) {
        _esvd = esvd;
        _v = (pi / 6.0) * _esvd * _esvd * _esvd;
    }

    /**
     * Set the grain population's elongation.
     * @param elong the elongation.
     */
    void _set_elongation(Real elong) {
        _elong = elong;
    }

    /**
     * Set the grain population's orientation.
     * @param ux the x component of the orientation vector.
     * @param uy the y component of the orientation vector.
     * @param uz the z component of the orientation vector.
     */
    void _set_orientation(Real ux, Real uy, Real uz) {
        // Record the values of u.
        _ux = ux;
        _uy = uy;
        _uz = uz;

        auto norm = sqrt(_ux * _ux + _uy * _uy + _uz * _uz);

        // Compute a normalized version of u.
        _u(0) = _ux / norm;
        _u(1) = _uy / norm;
        _u(2) = _uz / norm;
        _u.normalize();
    }

    /**
     * Set the applied field associated variables.
     * @param h_magnitude the applied field strength (A/m).
     * @param hx the x component of the applied field direction vector.
     * @param hy the y component of the applied field direction vector.
     * @param hz the z component of the applied field direction vector.
     */
    void _set_field(Real h_magnitude, Real hx, Real hy, Real hz) {
        // Update applied field strength.
        _h_magnitude = h_magnitude;

        // Record applied field direction components.
        _hx = hx;
        _hy = hy;
        _hz = hz;

        auto norm = sqrt(_hx * _hx + _hy * _hy + _hz * _hz);

        // Compute normalized version of h.
        _h(0) = _hx / norm;
        _h(1) = _hy / norm;
        _h(2) = _hz / norm;

        // Compute planar angle of applied field.
        _phi = acos(_h.dot(_u));
    }

    /**
     * Set the temperature dependent cc1 and cc2 constants.
     */
    void _set_cc1_cc2() {
        _c1 = cc1_function(_elong, _h_magnitude, _phi, _str_material);
        _c2 = cc2_function(_elong, _h_magnitude, _phi, _str_material);
    }

    /**
     * Set the Stoner-Wohlfarth energy functionals.
     */
    void _set_sw() {
        _sw = stoner_wohlfarth_energy_density_function(
                _elong, _h_magnitude, _phi, _str_material
        );
        _dsw = d_stoner_wohlfarth_energy_density_function(
                _elong, _h_magnitude, _phi, _str_material
        );
        _ddsw = dd_stoner_wohlfarth_energy_density_function(
                _elong, _h_magnitude, _phi, _str_material
        );
    }

    /**
     * Set the model's initial temperature.
     * @param initial_temperature the model's initial temperature.
     */
    void _set_initial_temperature(Real initial_temperature) {
        _initial_temperature = initial_temperature;
    }

    /**
     * Set the model's rho - distribution of states.
     * @param rho the rho vector.
     */
    void _set_rho(const std::vector<Real> &rho) {
        // Normalize rho so that each entry sums to 1.
        Real sum = 0;
        for (const auto &v: rho) {
            sum += v;
        }
        for (const auto &v: rho) {
            _rho.push_back(v / sum);
        }
    }

    /**
     * Initialize the model to the initial temperature and compute the LEMs and the energy barriers.
     */
    void _initialize_model() {
        // Set the current time to _initial_temperature.
        _temperature = _initial_temperature;

        // Compute the LEM states and energy barriers.
        // _compute_lems_and_energy_barriers();
    }

    /**
     * Polish a root using `niter` number of iterations of Newton-Raphson root finder.
     * @param theta0 the theta value to polish.
     * @param niter the number of iterations (should be small).
     * @return a polished version of the input `theta`.
     */
    Real _newton_raphson_polisher(Real theta0, ISize niter) {
        Real th = theta0;
        for (Index i = 0; i < niter; ++i) {
            auto nrcorrect = _dsw(th, _temperature) / _ddsw(th, _temperature);
            th = th - nrcorrect;
        }
        return th;
    }

    /**
     * Compute the LEM states (corresponding to angles) using the Stoner-Wolhfath model.
     * This only needs to be done if one of the following changes:
     * 	Applied field magnitude
     * 	Applied field direction
     * 	Temperature
     */
    void _compute_lems_and_energy_barriers() {

        ///////////////////////////////////////////////////////////////////////////
        /// Theta finding routine                                               ///
        /// finds critical points of the S/W energy density functional.         ///
        ///////////////////////////////////////////////////////////////////////////

        std::vector<Real> theta_minima;
        std::vector<Real> theta_maxima;
        std::vector<Real> theta_inflct;

        // Compute the Hessendorf matrix.
        ComplexMatrix4x4 H;
        H(0, 0) = Complex(-1, 0) * _c1(_temperature);
        H(0, 1) = Complex(0.0, 0.0);
        H(0, 2) = _c2(_temperature);
        H(0, 3) = Complex(1.0, 0.0);

        H(1, 0) = Complex(1.0, 0.0);
        H(1, 1) = Complex(0.0, 0.0);
        H(1, 2) = Complex(0.0, 0.0);
        H(1, 3) = Complex(0.0, 0.0);

        H(2, 0) = Complex(0.0, 0.0);
        H(2, 1) = Complex(1.0, 0.0);
        H(2, 2) = Complex(0.0, 0.0);
        H(2, 3) = Complex(0.0, 0.0);

        H(3, 0) = Complex(0.0, 0.0);
        H(3, 1) = Complex(0.0, 0.0);
        H(3, 2) = Complex(1.0, 0.0);
        H(3, 3) = Complex(0.0, 0.0);

        // Compute the eigenvalues of the Hessendorf matrix.
        auto eigen_values = H.eigenvalues();

        for (auto i = 0; i < eigen_values.size(); ++i) {

            // For each eigenvalue of H matrix, apply the transform theta = i * log(x).
            auto theta = Complex(0.0, 1.0) * log(eigen_values(i));

            // If the imaginary part is small enough, then ...
            if (abs(theta.imag()) < _eps) {

                // ... this is a real angle - i.e. a magnetization.
                auto theta_real = theta.real();

                // Try to correct for pi and map -pi to pi if needed.
                if (abs(pi + theta_real) < _eps) {
                    theta_real = pi;
                }

                // If the angle is small enough call it zero.
                if (abs(theta_real) < _eps) {
                    theta_real = 0.0;
                }

                // Polish the theta if required.
                if (_n_polish > 0) {
                    theta_real = _newton_raphson_polisher(theta_real, _n_polish);
                }

                // Check if `theta_real` is a minimum/maximum or inflection.
                if (_ddsw(theta_real, _temperature) < 0.0) {
                    // `theta_real` must be a maximum.
                    theta_maxima.push_back(theta_real);
                } else if (_ddsw(theta_real, _temperature) > 0.0) {
                    // `theta_real` must be a minimum.
                    theta_minima.push_back(theta_real);
                } else {
                    // `theta_real` must be an inflection point.
                    theta_inflct.push_back(theta_real);
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// LEM organisation routine                                                                                ///
        /// for any new minimum theta found, try to match it up with an                                             ///
        /// old/existing theta value to keep track of states. Insert new LEM                                        ///
        /// states as needed.                                                                                       ///
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // If there are no old_lems or new_lems, then we must be at the start of a simulation, so just set
        // _old_lem_states and _lem_states to theta_minima.
        if (_old_lem_states.empty() && _lem_states.empty()) {
            _lem_states = theta_minima;
            _old_lem_states = theta_minima;
        } else {
            // otherwise, we must have theta_minima & theta_maxima already set and we need to calculate them

            // Copy over _lem_states to _old_lem_states.
            _old_lem_states = _lem_states;

            // Copy over _rho to _old_rho
            _old_rho = _rho;

            if (_old_lem_states.size() > theta_minima.size()) {
                auto result = more_old_lems_than_new(
                        _old_lem_states, _old_rho, theta_minima
                );
                _lem_states = result.first;
                _rho = result.second;
            } else if (_old_lem_states.size() < theta_minima.size()) {
                auto result = more_new_lems_than_old(
                        _old_lem_states, _old_rho, theta_minima
                );
                _lem_states = result.first;
                _rho = result.second;
            } else {
                auto result = same_no_of_new_lems_as_old(
                        _old_lem_states, _old_rho, theta_minima
                );
                _lem_states = result.first;
                _rho = result.second;
            }
        }

        ///////////////////////////////////////////////////////////////////////////
        /// Energy barrier calculation routine.                                 ///
        /// this is specific to S/W uniaxial grains - with VIRGIL these are     ///
        /// precomputed.                                                        ///
        ///////////////////////////////////////////////////////////////////////////

        // Clear our current set of energy barriers.
        _energy_barriers.clear();

        Index other_lem_index = 0;
        for (Index i = 0; i < _lem_states.size(); ++i) {
            // For each of the thetas corresponding to energy minima.
            auto lem = _lem_states[i];

            // Energy barriers from the ith LEM to each other LEMs are stored here.
            std::unordered_map<Index, Real> ith_barriers;

            // Begin calculating the smallest barrier value between the current LEM
            // and each other LEM state - this calculation assumes a simple uniaxial
            // model so we can just look for the smallest barrier between all LEMS.
            Real smallest_barrier = 1E20;
            for (auto theta_max : theta_maxima) {
                // For each of the thetas corresponding to energy maxima, compute an
                // energy barrier.
                auto Eb = _v * (_sw(theta_max, _temperature) - _sw(lem, _temperature));
                if (Eb < smallest_barrier) {
                    // If the energy barrier is smaller than the current value for the
                    // energy barrier then use the smaller value as the energy barrier
                    // corresponding with the `theta_min` value.
                    smallest_barrier = Eb;
                }
            }

            // Append the smallest barrer we could find to the barriers associated
            // with the ith LEM. Note: for uniaxial-SD, there will be only one entry.
            other_lem_index = (i + 1) % 2;

            ith_barriers[other_lem_index] = smallest_barrier;

            // Add it to the energy barriers.
            _energy_barriers.push_back(ith_barriers);
        }
    }

    /**
     * Compute equilibrium rho values based on the energy barriers.
     */
    void _compute_equilibrium_rho() {
        using MPReal = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<500>>;

        // If there is only one LEM state ...
        if (_lem_states.size() == 1) {
            // ... then the equilibrium distribution of states is 1
            _rho_eq = {1};
        } else {
            // ... otherwise we compute the Boltzmann partition (eq. 8.11, pg. 213 Dunlop & Ozdemir).

            // First calculate Z, the normalization constant.
            MPReal Z(0.0);
            for (Real lem_state : _lem_states) {
                MPReal energy_density = _sw(lem_state, _temperature);
                MPReal energy = _v * energy_density;
                MPReal kbT = kb * (273.15 + _temperature);
                MPReal q = -1.0 * energy / kbT;
                MPReal exp_q = boost::multiprecision::exp(q);

                Z += exp_q;
            }
            DEBUG_MSG_SIMPLE_VAR(Z);

            // Calculate rho values.
            _rho_eq.clear();
            for (Real lem_state : _lem_states) {
                MPReal energy_density = _sw(lem_state, _temperature);
                MPReal energy = _v * energy_density;
                MPReal kbT = kb * (273.15 + _temperature);
                MPReal q = -1.0 * energy / kbT;
                MPReal exp_q = boost::multiprecision::exp(q);
                MPReal exp_q_over_Z = exp_q / Z;

                DEBUG_MSG_SIMPLE_VAR(exp_q);
                DEBUG_MSG_SIMPLE_VAR(exp_q_over_Z);

                _rho_eq.push_back((Real) (exp_q_over_Z));
            }
        }
    }

    [[nodiscard]] Vector3D _magnetization_no_project(bool with_ms, int magnetization_type) const {

        // Calculate a 3D magnetization vector from the theta values, and the equilibrium distribution.

        Vector3D m = {0.0, 0.0, 0.0};
        for (size_t lem_state_idx = 0; lem_state_idx < _lem_states.size(); ++lem_state_idx) {
            auto lem_state = _lem_states[lem_state_idx];

            Real rho = 0;
            switch (magnetization_type) {
                case EQUILIBRIUM_MAGNETIZATION:
                    rho = _rho_eq[lem_state_idx];
                    break;
                case POPULATION_MAGNETIZATION:
                    rho = _rho[lem_state_idx];
                    break;
                default:
                    throw std::runtime_error("Unknown magnetization type requested.");
                    break;
            }

            // If the lem_state, theta value is ...
            if (std::abs(lem_state) < _eps) {
                // ... equal to zero up to epsilon, then m is the grain axis vector.
                m += rho * _u;
            } else if (std::abs(lem_state - pi) < _eps) {
                // ... equal to pi up to epsilon, then m is the negative grain axis vector.
                m += -1.0 * rho * _u;
            } else if (std::abs(lem_state + pi) < _eps) {
                // ... equal to -pi up to epsilon, then m is the negative grain axis vector.
                m += -1.0 * rho * _u;
            } else {
                // ... theta must be some kind of intermediate value, so

                // 1) calculate the rotation axis
                auto rax = _u.cross(_h).normalized();

                // 2) calculate a rotation matrix
                auto R = rotation_matrix(rax, lem_state);

                // 3) rotate the grains axis vector by the angle of the LEM state and accumulate
                //    the fractional alignment
                m += rho * R * _u;
            }
        }

        // Multiply by Ms if the user has requested this.
        if (with_ms) {
            m = _ms(_temperature) * m;
        }

        return m;
    }

    [[nodiscard]] Real _magnetization(bool with_ms, int magnetization_type) const {
        Vector3D m = _magnetization_no_project(with_ms, magnetization_type);

        // Return the projection of the accumulated m on to the applied field.
        return scalar_proj(_h, m);
    }

};

#endif //SD_COOLING_STONER_WOHLFARTH_HPP
