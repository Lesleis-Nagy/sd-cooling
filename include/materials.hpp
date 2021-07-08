//
// Created by L. Nagy on 22/10/2020.
//

#ifndef SD_COOLING_MATERIALS_HPP
#define SD_COOLING_MATERIALS_HPP

#include <cmath>
#include <exception>
#include <functional>
#include <locale>
#include <string>

#include "basic_types.hpp"
#include "utilities.hpp"

/**
 * This exception class is thrown if a user specifies an unknown material name.
 */
class UnknownMaterialException: public std::exception
{
    [[nodiscard]] const char* what() const noexcept override {
        return "An unknown material name was requested.";
    }
};

/**
 * Enumeration to hold types of materials.
 */
enum Material {
    MAGNETITE,
    IRON
};


/**
 * Retrieves the Material enumeration value corresponding to an input string.
 * @param name the input string material.
 * @return an enumeration corresponding to the material name input.
 */
Material
name_to_material(const std::string & name) {
    using namespace std;

    string lower_name = to_lower_case(name);

    if (lower_name == "magnetite") {
        return MAGNETITE;
    }

    if (lower_name == "iron") {
        return IRON;
    }

    throw UnknownMaterialException();
}

/**
 * Return a function that will compute the saturation magnetization at temperature (in degrees centigrade) for a given
 * material.
 * @param material the constant associated with the material.
 * @return the saturation magnetization function for the material.
 */
std::function<Real(Real)>
saturation_magnetization_function(Material material)
{
    switch (material) {
        case IRON:
            return [] (Real T) -> Real {
                Real t = 273.0 + T;
                return 1.75221e6 - 1.21716e3 * t + 33.3368 * pow(t, 2.0) - 0.363228 * pow(t, 3.0) +
                    1.96713e-3 * pow(t, 4.0) - 5.98015e-6 * pow(t, 5.0) + 1.06587e-8 * pow(t, 6.0)
                    - 1.1048e-11 * pow(t, 7.0) + 6.16143e-15 * pow(t, 8.0) - 1.42904e-18 * pow(t, 9.0);
            };
        case MAGNETITE:
            return [] (Real T) -> Real {
                Real tmp = 580.0 - T;
                return 737.384 * 51.876 * pow(tmp, 0.4);
            };
        default:
            throw UnknownMaterialException();
    }
}

/**
 * Return a function that will compute the exchange constant at temperature (in degrees centigrade) for a given
 * material.
 * @param material the constant associated with the material.
 * @return the exchange constant function for the material.
 */
std::function<Real(Real)>
exchange_function(Material material)
{
    switch(material) {
        case IRON:
            return [] (Real T) -> Real {
                Real t = 273.0 + T;
                return -1.8952e-12 + 3.0657e-13 * t - 1.599e-15 * pow(t, 2.0) + 4.0151e-18 * pow(t, 3.0)
                    - 5.3728e-21 * pow(t, 4.0) + 3.6501e-24 * pow(t, 5.0) - 9.9515e-28 * pow(t, 6.0);
            };
        case MAGNETITE:
            return [] (Real T) -> Real {
                Real val = 21622.526 + 816.476 * (580.0 - T);
                return (sqrt(val) - 147.046) / 408.238e11;
            };
        default:
            throw UnknownMaterialException();
    }
}

/**
 * Return a function that will compute the magnetocrystalline anisotropy at temperature (in degrees centigrade) for a
 * given material.
 * @param material the constant associated with the material.
 * @return the magnetocrystalline anisotropy function for the material.
 */
std::function<Real(Real)>
anisotropy_function(Material material)
{
    switch(material) {
        case IRON:
            return [] (Real T) -> Real {
                Real t = 273.0 + T;
                Real k1 = 54967.1 + 44.2946 * t - 0.426485 * pow(t, 2.0) + 0.000811152 * pow(t, 3.0)
                        - 1.07579e-6 * pow(t, 4.0) + 8.83207e-10 * pow(t, 5.0) - 2.90947e-13 * pow(t, 6.0);
                Real C = 480.0 / 456.0;
                return C * k1;
            };
        case MAGNETITE:
            return [] (Real T) -> Real {
                Real t = 580.0 - T;
                return -2.13074e-5 * pow(t, 3.2);
            };
        default:
            throw UnknownMaterialException();
    }
}

#endif //SD_COOLING_MATERIALS_HPP
