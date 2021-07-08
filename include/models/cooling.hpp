//
// Created by L.Nagy on 20/05/2021.
//

#ifndef SD_COOLING_COOLING_HPP
#define SD_COOLING_COOLING_HPP

#include <iostream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include "sd_cooling.hpp"
#include "distribution.hpp"
#include "directions.hpp"
#include "temperature.hpp"
#include "unit.hpp"
#include "cooling_regime.hpp"

struct ModelDataRecord {
    Real time;
    Real temperature;
    Real m;
    Real m_eq;
};

typedef std::vector<ModelDataRecord> ModelDataRecordList;

ModelDataRecordList run_model(
        const std::string &material,
        Real initial_temperature,
        Real applied_field_strength,
        Real applied_field_x_dir,
        Real applied_field_y_dir,
        Real applied_field_z_dir,
        const std::string &applied_field_unit,
        const HistogramBinList &size_distr,
        const std::string &size_unit,
        const HistogramBinList &elong_distr,
        const DirectionBinList &dir_distr,
        const TemperatureStepList &temperature_regime,
        Real tau0,
        Real epsilon,
        unsigned int n_polish
) {
    ModelDataRecordList data;

    // Create the assemblage.
    Assemblage assemblage(initial_temperature, epsilon, n_polish);
    for (const auto &size_bin : size_distr) {
        for (const auto &elong_bin : elong_distr) {
            for (const auto &dir_bin : dir_distr) {
                Real fraction = size_bin.fraction * elong_bin.fraction * dir_bin.fraction;
                assemblage.add_stoner_wohlfarth_population(
                        size_to_meter(size_bin.value, size_unit),
                        elong_bin.value,
                        field_to_amps_per_meter(applied_field_strength, applied_field_unit),
                        applied_field_x_dir,
                        applied_field_y_dir,
                        applied_field_z_dir,
                        dir_bin.dir_x,
                        dir_bin.dir_y,
                        dir_bin.dir_z,
                        material,
                        1.0, 0.0,
                        fraction,
                        tau0
                );
            }
        }
    }

    for (auto & model_step : temperature_regime) {
        assemblage.update_temperature_and_rho(model_step.temperature, model_step.delta_time);
        data.push_back({
            model_step.time,
            model_step.temperature,
            assemblage.population_magnetization(false),
            assemblage.equilibrium_magnetization(false)
        });
    }

    return data;
}

void save_model_output(const std::vector<ModelDataRecord> &output, const std::string &file_name) {
    // Write data to file.
    std::ofstream fout;
    fout.open(file_name);
    for (auto &record : output) {
        fout << std::setw(20) << std::setprecision(15) << record.time << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.temperature << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.m << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.m_eq;
        fout << std::endl;
    }
    fout.close();
}

#endif //SD_COOLING_COOLING_HPP
