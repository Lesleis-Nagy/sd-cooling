//
// Created by L. Nagy on 25/11/2020.
//

#ifndef SD_COOLING_GRAIN_POPULATION_HPP
#define SD_COOLING_GRAIN_POPULATION_HPP

#include "basic_types.hpp"
#include "utilities.hpp"

class GrainPopulation {
public:
    virtual void update_temperature(Real T) = 0;
    virtual void update_rho(Real t) = 0;
    virtual void update_field(Real H) = 0;
    virtual void update_field(Real H, Real hx, Real hy, Real hz) = 0;
    virtual void update_field_and_temperature(Real H, Real T) = 0;
    virtual void update_field_and_temperature(Real H, Real hx, Real hy, Real hz, Real T) = 0;
    virtual void equilibrate() = 0;

    [[nodiscard]] virtual Vector3D  population_magnetization_no_project(bool with_ms) const = 0;
    [[nodiscard]] virtual Real population_magnetization(bool with_ms) const = 0;
    [[nodiscard]] virtual Real equilibrium_magnetization(bool with_ms) const = 0;


    [[nodiscard]] virtual std::vector<Real> get_rho() const = 0;
    [[nodiscard]] virtual std::vector<Real> get_rho_eq() const = 0;
};

#endif //SD_COOLING_GRAIN_POPULATION_HPP
