//
// Created by L.Nagy on 08/02/2022.
//

#include <fstream>
#define RYML_SINGLE_HDR_DEFINE_NOW
#include "ryml_all.hpp"

#include "models/monodispertion_complex_protocol.hpp"
#include "experiment_regime.hpp"
#include "unit.hpp"
#include "sd_cooling.hpp"

struct ProgramOptions {
    bool real_run = false;
    std::string protocol_file;
};

struct OutputRecord {
    std::string sample_name;
    Real time;
    Real temperature;
    Real m;
    Real mx;
    Real my;
    Real mz;
    Real rho0;
    Real rho1;
    Real bake_temperature;
    int step_type_index;
    std::string note;
};

ProgramOptions get_command_line_parameters(int argc, char *argv[]);

int main(int argc, char *argv[]) {

    ProgramOptions prog_opts;
    prog_opts = get_command_line_parameters(argc, argv);

    std::cout << "Running model with protocol '" << prog_opts.protocol_file << "'" << std::endl;

    std::vector<OutputRecord> output;

    try {
        auto protocol = Protocol::from_file(prog_opts.protocol_file);

        TemperatureAndFieldStepList trm_acquisition_steps;
        TemperatureAndFieldStepList lab_demag_steps;

        std::tie(trm_acquisition_steps, lab_demag_steps) = simulation_protocol_steps(protocol);

        if (prog_opts.real_run) {
            std::cout << "Perform real run" << std::endl;
            Real tau0 = 1E-10;
            Real epsilon = 1E-15;
            unsigned int n_polish = 0;
            Assemblage assemblage(protocol.trm_acquisition.initial_temperature, epsilon, n_polish);
            assemblage.add_stoner_wohlfarth_population(
                    size_to_meter(protocol.geometry.size, protocol.geometry.size_unit),
                    protocol.geometry.elongation,
                    protocol.trm_acquisition.field.strength, // already in A/m
                    protocol.trm_acquisition.field.x_direction,
                    protocol.trm_acquisition.field.y_direction,
                    protocol.trm_acquisition.field.z_direction,
                    1.0,
                    0.0,
                    0.0,
                    protocol.material.name,
                    1.0, 0.0,
                    1.0,
                    tau0);

            // This is the TRM acquisition phase.
            for (const auto &step: trm_acquisition_steps) {
                assemblage.update_temperature_and_rho(step.temperature, step.delta_time);
                Vector3D m = assemblage.population_magnetization_no_project(true);
                std::vector<std::vector<Real>> rho = assemblage.rhos();
                output.push_back({protocol.sample.name,
                                  step.time,
                                  step.temperature,
                                  m.norm(),
                                  m[0], m[1], m[2],
                                  rho[0][0], rho[0][1],
                                  step.step_bake_temperature,
                                  step.step_type_index,
                                  step.notes});
            }

            // This is the lab demagnetization phase.
            for (const auto &step: lab_demag_steps) {
                assemblage.update_field_temperature_and_rho(step.field_strength,
                                                            step.field_x_direction,
                                                            step.field_y_direction,
                                                            step.field_z_direction,
                                                            step.temperature,
                                                            step.delta_time);
                Vector3D m = assemblage.population_magnetization_no_project(true);
                std::vector<std::vector<Real>> rho = assemblage.rhos();
                output.push_back({protocol.sample.name,
                                  step.time,
                                  step.temperature,
                                  m.norm(),
                                  m[0], m[1], m[2],
                                  rho[0][0], rho[0][1],
                                  step.step_bake_temperature,
                                  step.step_type_index,
                                  step.notes});
            }


            // write the output file.
            if (!output.empty()) {
                std::ofstream fout;
                fout.open(protocol.outputs.raw_simulation_data_csv);
                fout << "sample name," << "time," << "temperature," << "m," << "mx," << "my," << "mz," << "rho0,"
                     << "rho1," << "bake temperature," << "step type index," << "note" << std::endl;
                for (const auto &rec: output) {
                    fout << rec.sample_name << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.time << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.temperature << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.m << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.mx << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.my << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.mz << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.rho0 << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.rho1 << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.bake_temperature << ",";
                    fout << std::setw(20) << std::setprecision(15) << rec.step_type_index << ",";
                    fout << rec.note;
                    fout << std::endl;
                }
            }

        }

        if (protocol.outputs.produce_trm_acquisistion_protocol_csv) {
            std::cout << "Writing TRM acquisition protocol to " << protocol.outputs.trm_acquisition_protocol_csv << std::endl;
            save_experiment_regime(trm_acquisition_steps, protocol.outputs.trm_acquisition_protocol_csv);
        }

        if (protocol.outputs.produce_lab_demag_protocol_csv) {
            std::cout << "Writing demag protocol to " << protocol.outputs.lab_demag_protocol_csv << std::endl;
            save_experiment_regime(lab_demag_steps, protocol.outputs.lab_demag_protocol_csv);
        }

    } catch (ProtocolDataInvalidException &exp) {
        std::cout << "An error occurred: " << std::endl;
        std::cout << exp.message() << std::endl;
    }

    return 0;
}

ProgramOptions
get_command_line_parameters(int argc, char *argv[]) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    po::positional_options_description pos_desc;

    ProgramOptions prog_opts;

    try {
        desc.add_options()
                ("help", "produce help message")
                ("protocol-file", po::value<std::string>(), "protocol file")
                ("real-run", po::bool_switch(&prog_opts.real_run), "actually run the model")
                ;

        // Positional.
        pos_desc.add("protocol-file", 1);

        po::variables_map vm;
        po::store(
                po::command_line_parser(argc, argv)
                        .options(desc)
                        .positional(pos_desc)
                        .style(po::command_line_style::unix_style ^ po::command_line_style::allow_short)
                        .run(), vm
        );
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << std::endl;
            exit(EXIT_SUCCESS);
        }

        if (vm.count("protocol-file")) {
            prog_opts.protocol_file = vm["protocol-file"].as<std::string>();
        } else {
            std::cout << "Protocol file is missing." << std::endl;
            exit(EXIT_FAILURE);
        }

        return prog_opts;

    } catch (po::error &e) {
        std::cout << e.what() << std::endl;
        exit(EXIT_FAILURE);
    } catch (...) {
        std::exception_ptr ptr_exc = std::current_exception();
        std::cout << boost::current_exception_diagnostic_information() << std::endl;
        exit(EXIT_FAILURE);
    }
}

