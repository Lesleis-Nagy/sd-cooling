//
// Created by L.Nagy on 20/05/2021.
//

#ifndef SD_COOLING_EXPERIMENT_REGIME_HPP
#define SD_COOLING_EXPERIMENT_REGIME_HPP

#include <vector>

#include <boost/format.hpp>

#include "nlohmann/json.hpp"

#include "basic_types.hpp"
#include "temperature.hpp"
#include "debug.hpp"

#include "models/monodispertion_complex_protocol.hpp"

struct TemperatureAndFieldStep {
    Real time;
    Real temperature;
    Real delta_time;
    Real field_strength;    // A/m
    Real field_x_direction;
    Real field_y_direction;
    Real field_z_direction;
    Real step_bake_temperature;
    int step_type_index;
    std::string notes;
};

typedef std::vector<TemperatureAndFieldStep> TemperatureAndFieldStepList;

std::pair<TemperatureAndFieldStepList, TemperatureAndFieldStepList> simulation_protocol_steps(
    const Protocol &input_protocol
) {
    using boost::format;
    using boost::str;

    TemperatureAndFieldStepList trm_acquisition_protocol;
    TemperatureAndFieldStepList lab_demag_protocol;
    Real t = 0;

    // Add the acquisition part.
    auto cool_fun = newtonian_temperature_function(
            input_protocol.trm_acquisition.ambient_temperature,
            input_protocol.trm_acquisition.initial_temperature,
            input_protocol.trm_acquisition.temperature_at_t1,
            input_protocol.trm_acquisition.t1
    );

    auto d_cool_fun_dt = d_newtonian_temperature_function_dt(
            input_protocol.trm_acquisition.ambient_temperature,
            input_protocol.trm_acquisition.initial_temperature,
            input_protocol.trm_acquisition.temperature_at_t1,
            input_protocol.trm_acquisition.t1
    );

    Real time = 0;
    Real temperature = input_protocol.trm_acquisition.initial_temperature;

    while (input_protocol.trm_acquisition.temperature_at_t1 < temperature) {
        temperature = cool_fun(time);
        Real delta_time = std::abs(input_protocol.trm_acquisition.allowable_fractional_drop * temperature / d_cool_fun_dt(time));
        trm_acquisition_protocol.push_back({time,
                                            temperature,
                                            delta_time,
                                            input_protocol.trm_acquisition.field.strength,
                                            input_protocol.trm_acquisition.field.x_direction,
                                            input_protocol.trm_acquisition.field.y_direction,
                                            input_protocol.trm_acquisition.field.z_direction,
                                            0.0,
                                            0,
                                            "TRM ACQUISITION"});
        time += delta_time;
    }

    // Add the cooling part.
    for (const auto& step : input_protocol.lab_demag_protocol.steps) {
        // Heating step.
        Real stop_time = t + step.heating_time;
        auto tfun_heat = linear_temperature_function(t, step.stop_temperature, stop_time, step.bake_temperature);
        auto dtfun_heat = d_linear_temperature_function_dt(t, step.stop_temperature, stop_time, step.bake_temperature);
        while (t < stop_time) {
            Real dt = abs((input_protocol.lab_demag_protocol.allowable_fractional_drop * tfun_heat(t)) / dtfun_heat(t));
            lab_demag_protocol.push_back({t,
                                          tfun_heat(t),
                                          dt,
                                          step.field.strength,
                                          step.field.x_direction,
                                          step.field.y_direction,
                                          step.field.z_direction,
                                          step.bake_temperature,
                                          step.type_index,
                                          str(format("DEMAG STEP - HEATING - %1%") % step.type)
            });
            t += dt;
        }

        // Fix up the last entry value.
        auto& last_entry = lab_demag_protocol.back();
        last_entry.delta_time = stop_time - last_entry.time;

        // Add entry to account for stop time.
        lab_demag_protocol.push_back({stop_time,
                                      step.bake_temperature,
                                      step.bake_time,
                                      step.field.strength,
                                      step.field.x_direction,
                                      step.field.y_direction,
                                      step.field.z_direction,
                                      step.bake_temperature,
                                      step.type_index,
                                      str(format("DEMAG STEP - BAKING - %1%") % step.type)
        });

        // Roll t forward by ex_bake_time.
        t = stop_time + step.bake_time;

        // Cooling step.
        stop_time = t + step.cool_time;
        auto tfun_cool = newtonian_temperature_function(
                step.ambient_temperature, t, step.bake_temperature, stop_time, step.stop_temperature);
        auto dtfun_cool = d_newtonian_temperature_function_dt(
                step.ambient_temperature, t, step.bake_temperature, stop_time, step.stop_temperature);
        while (t < stop_time) {
            Real dt = std::abs((input_protocol.lab_demag_protocol.allowable_fractional_drop * tfun_cool(t)) / dtfun_cool(t));
            lab_demag_protocol.push_back({t,
                                          tfun_cool(t),
                                          dt,
                                          step.field.strength,
                                          step.field.x_direction,
                                          step.field.y_direction,
                                          step.field.z_direction,
                                          step.bake_temperature,
                                          step.type_index,
                                          str(format("DEMAG STEP - COOLING - %1%") % step.type)
            });
            t += dt;
        }
    }

    return {trm_acquisition_protocol, lab_demag_protocol};
}

void save_experiment_regime(const TemperatureAndFieldStepList &regime, const std::string &file_name) {
    // Write data to file.
    std::ofstream fout;
    fout.open(file_name);

    // Heading.
    fout << "time (s)," << "temperature (C)," << "delta time (s)," << "field strength (A/m),"
         << "field x direction," << "field y direction," << "field z direction,"
         << "bake temperature (C)," << "record type," << "notes"
         << std::endl;

    // Records.
    for (auto &record : regime) {
        fout << std::setw(20) << std::setprecision(15) << record.time << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.temperature << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.delta_time << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.field_strength << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.field_x_direction << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.field_y_direction << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.field_z_direction << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.step_bake_temperature << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.step_type_index << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.notes;
        fout << std::endl;
    }

    fout.close();
}

#endif //SD_COOLING_EXPERIMENT_REGIME_HPP
