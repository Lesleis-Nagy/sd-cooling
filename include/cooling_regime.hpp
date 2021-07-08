//
// Created by L.Nagy on 20/05/2021.
//

#ifndef SD_COOLING_COOLING_REGIME_HPP
#define SD_COOLING_COOLING_REGIME_HPP

#include <vector>

#include "nlohmann/json.hpp"

#include "basic_types.hpp"
#include "temperature.hpp"
#include "debug.hpp"

struct TemperatureStep {
    Real time;
    Real temperature;
    Real delta_time;
};

typedef std::vector<TemperatureStep> TemperatureStepList;

TemperatureStepList simple_cooling(
        Real ambient_temperature,
        Real initial_temperature,
        Real temperature_at_t_one,
        Real t_one,
        Real allowable_fractional_drop,
        Real stopping_temperature
) {
    TemperatureStepList regime;

    auto cool_fun = newtonian_temperature_function(
            ambient_temperature,
            initial_temperature,
            temperature_at_t_one,
            t_one
    );

    auto d_cool_fun_dt = d_newtonian_temperature_function_dt(
            ambient_temperature,
            initial_temperature,
            temperature_at_t_one,
            t_one
    );

    Real time = 0;
    Real temperature = initial_temperature;

    while (stopping_temperature < temperature) {
        temperature = cool_fun(time);
        Real delta_time = std::abs(allowable_fractional_drop * temperature / d_cool_fun_dt(time));
        regime.push_back({time, temperature, delta_time});
        time += delta_time;
    }

    return regime;
}

TemperatureStepList simple_cooling_with_A_and_fixed_dt (
        Real ambient_temperature,
        Real initial_temperature,
        Real A,
        Real dt,
        Real stopping_temperature
) {
    TemperatureStepList regime;

    auto cool_fun = newtonian_temperature_function(
            ambient_temperature,
            initial_temperature,
            A
    );

    Real time = 0;
    Real temperature = initial_temperature;

    while (stopping_temperature < temperature) {
        temperature = cool_fun(time);
        Real delta_time = dt;
        regime.push_back({time, temperature, delta_time});
        time += delta_time;
    }

    return regime;
}

TemperatureStepList simple_cooling_from_json(const nlohmann::json &json) {
    TemperatureStepList output;

    if (json.contains("cooling_regime")) {
        auto cooling_regime = json["cooling_regime"];
        Real ambient_temperature = 0;
        Real initial_temperature = 0;
        Real temperature_at_t_one = 0;
        Real t_one = 0;
        Real allowable_fractional_drop = 0;
        Real stopping_temperature = 0;

        if (!cooling_regime.contains("ambient_temperature")) {
            throw std::runtime_error("Program JSON 'cooling_regime' 'ambient_temperature' is missing");
        }
        ambient_temperature = cooling_regime["ambient_temperature"];

        if (!cooling_regime.contains("initial_temperature")) {
            throw std::runtime_error("Program JSON 'cooling_regime' 'initial_temperature' is missing");
        }
        initial_temperature = cooling_regime["initial_temperature"];

        if (!cooling_regime.contains("reference_time")) {
            throw std::runtime_error("Program JSON 'cooling_regime' 'reference_time' is missing");
        }
        t_one = cooling_regime["reference_time"];

        if (!cooling_regime.contains("temperature_at_reference_time")) {
            throw std::runtime_error("Program JSON 'cooling_regime' 'temperature_at_reference_time' is missing");
        }
        temperature_at_t_one = cooling_regime["temperature_at_reference_time"];

        if (cooling_regime.contains("allowable_percentage_drop")) {
            allowable_fractional_drop = (Real)cooling_regime["allowable_percentage_drop"] / 100.0;
        } else {
            allowable_fractional_drop = 0.1/100.0;
        }

        if (cooling_regime.contains("stopping_temperature")) {
            stopping_temperature = cooling_regime["stopping_temperature"];
        } else {
            stopping_temperature = temperature_at_t_one;
        }

        DEBUG_MSG_SIMPLE_VAR(ambient_temperature);
        DEBUG_MSG_SIMPLE_VAR(initial_temperature);
        DEBUG_MSG_SIMPLE_VAR(temperature_at_t_one);
        DEBUG_MSG_SIMPLE_VAR(t_one);
        DEBUG_MSG_SIMPLE_VAR(allowable_fractional_drop);
        DEBUG_MSG_SIMPLE_VAR(stopping_temperature);

        output = simple_cooling(ambient_temperature,
                                initial_temperature,
                                temperature_at_t_one,
                                t_one,
                                allowable_fractional_drop,
                                stopping_temperature);
        return output;
    } else {
        throw std::runtime_error("Program JSON doesn't contain 'cooling_regime'");
    }
}

void save_temperature_regime(const TemperatureStepList &regime, const std::string &file_name) {
    // Write data to file.
    std::ofstream fout;
    fout.open(file_name);
    for (auto &record : regime) {
        fout << std::setw(20) << std::setprecision(15) << record.time << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.temperature << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.delta_time;
        fout << std::endl;
    }
    fout.close();
}

#endif //SD_COOLING_COOLING_REGIME_HPP
