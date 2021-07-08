//
// Created by L. Nagy on 24/11/2020.
//

#ifndef SD_COOLING_SD_COOLING_HPP
#define SD_COOLING_SD_COOLING_HPP

//#include <ranges>
#include <vector>
#include <algorithm>
#include <numeric>

#include "constants.hpp"
#include "basic_types.hpp"
#include "utilities.hpp"

#include "grain_population.hpp"
#include "stoner_wohlfarth.hpp"

struct GrainPopulationFraction {
    std::shared_ptr<GrainPopulation> grain_population;
    Real fraction;
};

class Assemblage {
public:

    explicit Assemblage(Real Tinit, Real eps = 1E-12, ISize npolish = 0) :
            _initial_temperature(Tinit), _eps(eps), _npolish(npolish) {}

    /**
     * Add a Stoner-Wohlfarth grain population to the assemblage.
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
     * @param rho0 the fraction of grains in state 0 (will be normalized to sum to 1 with rho1)
     * @param rho1 the fraction of grains in state 1 (will be normalized to sum to 1 with rho0)
     * @param fraction the fraction/proportion of the total assemblage that this grain represents
     * @param tau0 the switching frequency in Hz (i.e. per second), see Dunlop & Ozdemir (2001), pp. 202 (eqn. 8.3)
     */
    void add_stoner_wohlfarth_population(Real esvd, Real elong,
                                         Real H, Real hx, Real hy, Real hz,
                                         Real ux, Real uy, Real uz,
                                         const std::string &material,
                                         Real rho0, Real rho1, Real fraction, Real tau0) {
        // Push back a new grain with the give parameters, along with its fraction.
        std::vector<Real> rho = {rho0, rho1};
        _grain_population_fractions.push_back(
                {
                        std::make_shared<StonerWohlfarthGrainPopulation>(
                                esvd, elong,
                                H, hx, hy, hz,
                                ux, uy, uz,
                                material,
                                _initial_temperature,
                                rho, tau0,
                                _eps, _npolish
                        ),
                        fraction
                }
        );
    }

    void update_temperature(Real T) {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [T](GrainPopulationFraction &gp) {
                    gp.grain_population->update_temperature(T);
                }
        );
    }

    void equilibrate() {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [](GrainPopulationFraction &gp) {
                    gp.grain_population->equilibrate();
                }
        );
    }

    void update_field_strength(Real H) {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [H](GrainPopulationFraction &gp) {
                    gp.grain_population->update_field(H);
                }
        );
    }

    void update_field(Real H, Real hx, Real hy, Real hz) {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [H, hx, hy, hz](GrainPopulationFraction &gp) {
                    gp.grain_population->update_field(H, hx, hy, hz);
                }
        );
    }

    void update_rho(Real dt) {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [dt](GrainPopulationFraction &gp) {
                    gp.grain_population->update_rho(dt);
                }
        );
    }

    void update_temperature_and_rho(Real T, Real dt) {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [T, dt](GrainPopulationFraction &gp) {
                    gp.grain_population->update_temperature(T);
                    gp.grain_population->update_rho(dt);
                }
        );
    }

    void update_field_strength_and_rho(Real H, Real dt) {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [H, dt](GrainPopulationFraction &gp) {
                    gp.grain_population->update_field(H);
                    gp.grain_population->update_rho(dt);
                }
        );
    }

    void update_field_and_rho(Real H, Real hx, Real hy, Real hz, Real dt) {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [H, hx, hy, hz, dt](GrainPopulationFraction &gp) {
                    gp.grain_population->update_field(H, hx, hy, hz);
                    gp.grain_population->update_rho(dt);
                }
        );
    }

    void update_field_temperature_and_rho(Real H, Real T, Real dt) {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [H, T, dt](GrainPopulationFraction &gp) {
                    gp.grain_population->update_field_and_temperature(H, T);
                    gp.grain_population->update_rho(dt);
                }
        );
    }

    void update_field_temperature_and_rho(Real H, Real hx, Real hy, Real hz, Real T, Real dt) {
        std::for_each(
                _grain_population_fractions.begin(),
                _grain_population_fractions.end(),
                [H, hx, hy, hz, T, dt](GrainPopulationFraction &gp) {
                    gp.grain_population->update_field_and_temperature(H, hx, hy, hz, T);
                    gp.grain_population->update_rho(dt);
                }
        );
    }

    [[nodiscard]] Real population_magnetization(bool with_ms) const {
        Real result = 0.0;

        for (const auto &gp : _grain_population_fractions) {
            result += gp.fraction * gp.grain_population->population_magnetization(with_ms);
        }

        return result;
    }

    [[nodiscard]] Real equilibrium_magnetization(bool with_ms) const {
        Real result = 0.0;

        for (const auto &gp : _grain_population_fractions) {
            result += gp.fraction * gp.grain_population->equilibrium_magnetization(with_ms);
        }

        return result;
    }

    [[nodiscard]] std::vector<std::vector<Real>> rhos() {
        std::vector<std::vector<Real>> rho_values;

        for (const auto &gp : _grain_population_fractions) {
          auto rho = gp.grain_population->get_rho();
          rho_values.push_back(rho);
        }

        return rho_values;
    }

    [[nodiscard]] std::vector<std::vector<Real>> equilibrium_rhos() {
        std::vector<std::vector<Real>> rho_eq_values;

        for (const auto &gp : _grain_population_fractions) {
            auto rho_eq = gp.grain_population->get_rho_eq();
            rho_eq_values.push_back(rho_eq);
        }

        return rho_eq_values;
    }

private:
    std::vector<GrainPopulationFraction> _grain_population_fractions;

    std::vector< std::array<Real, 9> > _data_cache;

    Real _initial_temperature;
    Real _eps;
    ISize _npolish;
};

#endif //SD_COOLING_SD_COOLING_HPP
