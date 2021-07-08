//
// Created by L.Nagy on 20/05/2021.
//

#ifndef SD_COOLING_EXPERIMENT_REGIME_HPP
#define SD_COOLING_EXPERIMENT_REGIME_HPP

#include <vector>

#include "nlohmann/json.hpp"

#include "basic_types.hpp"
#include "temperature.hpp"
#include "debug.hpp"

struct TemperatureAndFieldStep {
    Real time;
    Real temperature;
    Real delta_time;
    Real field_strength;    // A/m
    Real field_x_direction;
    Real field_y_direction;
    Real field_z_direction;
};

typedef std::vector<TemperatureAndFieldStep> TemperatureAndFieldStepList;

TemperatureAndFieldStepList
izzi_experiment(
        Real lab_ambient_temperature,
        Real lab_cooled_temperature,
        Real lab_heating_time,
        const std::vector<Real>& lab_bake_temperatures,
        Real lab_bake_time,
        Real lab_cooling_time,
        Real lab_field_strength,
        const std::string& field_unit,
        Real lab_field_x_direction,
        Real lab_field_y_direction,
        Real lab_field_z_direction,
        Real allowable_fractional_drop
) {
    std::vector<Real> field_strengths{field_to_amps_per_meter(lab_field_strength, field_unit), 0.0};
    TemperatureAndFieldStepList protocol;
    Real t = 0;
    for (long double bake_temperature : lab_bake_temperatures) {
        for (Real field_strength : field_strengths) {
            // Heating step.
            Real stop_time = t + lab_heating_time;
            auto tfun_heat = linear_temperature_function(t, lab_cooled_temperature, stop_time, bake_temperature);
            auto dtfun_heat = d_linear_temperature_function_dt(t, lab_cooled_temperature, stop_time, bake_temperature);
            while (t < stop_time) {
                Real dt = abs((allowable_fractional_drop * tfun_heat(t)) / dtfun_heat(t));
                protocol.push_back({
                    t, tfun_heat(t), dt, field_strength,
                    lab_field_x_direction, lab_field_y_direction, lab_field_z_direction
                });
                t += dt;
            }
            // Fix up the last entry value.
            auto& last_entry = protocol.back();
            last_entry.delta_time = stop_time - last_entry.time;
            // Add entry to account for stop time.
            protocol.push_back({
                stop_time,
                bake_temperature,
                lab_bake_time,
                field_strength,
                lab_field_x_direction,
                lab_field_y_direction,
                lab_field_z_direction
            });

            // Roll t forward by ex_bake_time.
            t = stop_time + lab_bake_time;

            // Cooling step.
            stop_time = t + lab_cooling_time;
            auto tfun_cool = newtonian_temperature_function(
                    lab_ambient_temperature, t, bake_temperature, stop_time, lab_cooled_temperature);
            auto dtfun_cool = d_newtonian_temperature_function_dt(
                    lab_ambient_temperature, t, bake_temperature, stop_time, lab_cooled_temperature);
            while (t < stop_time) {
                Real dt = std::abs((allowable_fractional_drop * tfun_cool(t)) / dtfun_cool(t));
                protocol.push_back({
                    t, tfun_cool(t), dt, field_strength,
                    lab_field_x_direction, lab_field_y_direction, lab_field_z_direction
                });
                t += dt;
            }
        }
    }
    return protocol;
}

TemperatureAndFieldStepList
izzi_experiment_from_json(const nlohmann::json &lab_demagnetising_regime) {
    Real field_strength = 0;
    Real field_x_direction = 0;
    Real field_y_direction = 0;
    Real field_z_direction = 0;
    std::string field_unit = "";
    Real ambient_temperature = 0;
    Real cooled_temperature = 0;
    Real heating_time = 0;
    std::vector<Real> baking_temperatures;
    Real bake_time;
    Real cooling_time;
    Real allowable_percentage_drop = 0.1/100.0;
    if (lab_demagnetising_regime.contains("applied_field")) {
        auto applied_field = lab_demagnetising_regime["applied_field"];
        if (!applied_field.contains("strength")) {
            throw std::runtime_error("Program JSON, IZZI regime 'applied_field'.'strength' missing");
        }
        if (!applied_field.contains("direction")) {
            throw std::runtime_error("Program JSON, IZZI regime 'applied_field'.'direction' missing");
        }
        if (applied_field["direction"].size() != 3) {
            throw std::runtime_error("Program JSON, IZZI regime 'applied_field'.'direction' has 3 elements");
        }
        if (!applied_field.contains("unit")) {
            throw std::runtime_error("Program JSON, IZZI regime 'applied_field'.'unit' missing");
        }
        field_strength = applied_field["strength"];
        field_x_direction = applied_field["direction"][0];
        field_y_direction = applied_field["direction"][1];
        field_z_direction = applied_field["direction"][2];
        field_unit = applied_field["unit"];
    } else {
        throw std::runtime_error("Program JSON, IZZI regime has no 'applied_field'");
    }
    if (lab_demagnetising_regime.contains("ambient_temperature")) {
        ambient_temperature = lab_demagnetising_regime["ambient_temperature"];
    } else {
        throw std::runtime_error("Program JSON, IZZI regime, 'ambient_temperature' is missing");
    }
    if (lab_demagnetising_regime.contains("cooled_temperature")) {
        cooled_temperature = lab_demagnetising_regime["cooled_temperature"];
    } else {
        throw std::runtime_error("Program JSON, IZZI regime, 'cooled_temperature' is missing");
    }
    if (lab_demagnetising_regime.contains("heating_time")) {
        heating_time = lab_demagnetising_regime["heating_time"];
    } else {
        throw std::runtime_error("Program JSON, IZZI regime, 'heating_time' is missing");
    }
    if (lab_demagnetising_regime.contains("baking_temperatures")) {
        for (Real baking_temperature : lab_demagnetising_regime["baking_temperatures"]) {
            baking_temperatures.push_back(baking_temperature);
        }
    } else {
        throw std::runtime_error("Program JSON, IZZI regime, 'baking_temperatures' is missing");
    }
    if (lab_demagnetising_regime.contains("bake_time")) {
        bake_time = lab_demagnetising_regime["bake_time"];
    } else {
        throw std::runtime_error("Program JSON, IZZI regime, 'bake_time' is missing");
    }
    if (lab_demagnetising_regime.contains("cooling_time")) {
        cooling_time = lab_demagnetising_regime["cooling_time"];
    } else {
        throw std::runtime_error("Program JSON, IZZI regime, 'cooling_time' is missing");
    }
    if (lab_demagnetising_regime.contains("allowable_percentage_drop")) {
        allowable_percentage_drop = (Real)lab_demagnetising_regime["allowable_percentage_drop"]/100.0;
    } else {
        throw std::runtime_error("Program JSON, IZZI regime, 'allowable_percentage_drop' is missing");
    }

    DEBUG_MSG_SIMPLE_VAR(field_strength);
    DEBUG_MSG_SIMPLE_VAR(field_x_direction);
    DEBUG_MSG_SIMPLE_VAR(field_y_direction);
    DEBUG_MSG_SIMPLE_VAR(field_z_direction);
    DEBUG_MSG_SIMPLE_VAR(field_unit);
    DEBUG_MSG_SIMPLE_VAR(ambient_temperature);
    DEBUG_MSG_SIMPLE_VAR(cooled_temperature);
    DEBUG_MSG_SIMPLE_VAR(heating_time);
    DEBUG_MSG_STD_VECTOR(baking_temperatures);
    DEBUG_MSG_SIMPLE_VAR(bake_time);
    DEBUG_MSG_SIMPLE_VAR(cooling_time);
    DEBUG_MSG_SIMPLE_VAR(allowable_percentage_drop);

    return izzi_experiment(
            ambient_temperature,
            cooled_temperature,
            heating_time,
            baking_temperatures,
            bake_time,
            cooling_time,
            field_strength,
            field_unit,
            field_x_direction,
            field_y_direction,
            field_z_direction,
            allowable_percentage_drop);
}

TemperatureAndFieldStepList
experiment_from_json(const nlohmann::json &json) {
    TemperatureAndFieldStepList regime;

    if (json.contains("lab_demagnetising_regime")) {
        auto lab_demagnetising_regime = json["lab_demagnetising_regime"];
        if (!lab_demagnetising_regime.contains("regime")) {
            throw std::runtime_error("Program JSON, 'regime' missing from 'demagnetizing_regime'");
        }
        std::string regime = lab_demagnetising_regime["regime"];
        if (regime == "IZZI") {
            return izzi_experiment_from_json(lab_demagnetising_regime);
        } else {
            throw std::runtime_error("Program JSON, unknown demagnetising regime");
        }
    } else {
        throw std::runtime_error("Program JSON doesn't contain 'lab_demagnetising_regime'");
    }
}

void save_experiment_regime(const TemperatureAndFieldStepList &regime, const std::string &file_name) {
    // Write data to file.
    std::ofstream fout;
    fout.open(file_name);
    for (auto &record : regime) {
        fout << std::setw(20) << std::setprecision(15) << record.time << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.temperature << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.delta_time << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.field_strength << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.field_x_direction << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.field_y_direction << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.field_z_direction;
        fout << std::endl;
    }
    fout.close();
}

#endif //SD_COOLING_EXPERIMENT_REGIME_HPP
