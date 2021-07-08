//
// Created by L.Nagy on 19/05/2021.
//

#include <algorithm>
#include <fstream>

#include "nlohmann/json.hpp"

#include "distribution.hpp"

#include "models/cooling.hpp"

struct ProgramOptions {
    bool real_run = false;
    std::string json_spec_file;
};

ProgramOptions get_command_line_parameters(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    using json = nlohmann::json;

    ProgramOptions prog_opts;
    prog_opts = get_command_line_parameters(argc, argv);

    // Load json
    std::ifstream fin(prog_opts.json_spec_file);
    json json_data;
    fin >> json_data;

    // Model variables

    // Material name
    std::string material;
    if (!json_data.contains("material")) {
        throw std::runtime_error("Program JSON does not contain 'material'");
    }
    material = json_data["material"];

    // Initial temperature
    Real initial_temperature;
    if (json_data.contains("cooling_regime")) {
        auto cooling_regime = json_data["cooling_regime"];
        if (!cooling_regime.contains("initial_temperature")) {
            throw std::runtime_error("Program JSON 'cooling_regime' is missing 'initial_temperature'");
        } else {
            initial_temperature = cooling_regime["initial_temperature"];
        }
    } else {
        throw std::runtime_error("Program JSON does not contain 'cooling_regime'");
    }

    // Size unit
    std::string size_unit = size_unit_from_json(json_data);

    // Applied field.
    std::string applied_field_unit = field_unit_from_json(json_data);
    Real applied_field_strength = 0;
    Real applied_field_x_dir = 0;
    Real applied_field_y_dir = 0;
    Real applied_field_z_dir = 0;
    if (json_data.contains("applied_field")) {
        auto applied_field = json_data["applied_field"];
        if (!applied_field.contains("strength")) {
            throw std::runtime_error("Program JSON, 'applied_field' is given but 'strength' is missing.");
        }
        if (!applied_field.contains("direction")) {
            throw std::runtime_error("Program JSON, 'applied_field' is given but 'direction' is missing.");
        }
        if (applied_field.contains("direction")) {
            if (applied_field["direction"].size() != 3) {
                throw std::runtime_error("Program JSON 'applied_field' 'direction' needs 3 components");
            }
        }
        applied_field_strength = applied_field["strength"];
        applied_field_x_dir = applied_field["direction"][0];
        applied_field_y_dir = applied_field["direction"][1];
        applied_field_z_dir = applied_field["direction"][2];
    }

    // Additional parameters.
    Real tau0 = 0;
    Real epsilon = 0;
    unsigned int n_polish = 0;

    if (json_data.contains("tau0")) {
        tau0 = json_data["tau0"];
    } else {
        // Assign default value.
        tau0 = 1E-10;
    }

    if (json_data.contains("epsilon")) {
        epsilon = json_data["epsilon"];
    } else {
        // Assign default value.
        epsilon = 1E-15;
    }

    if (json_data.contains("n_polish")) {
        n_polish = json_data["n_polish"];
    } else {
        // Assign default value.
        n_polish = 0;
    }

    // Create a size distribution.
    HistogramBinList size_distr = size_distribution_from_json(json_data);
    // Create an elongation distribution.
    HistogramBinList elongation_distr = elongation_distribution_from_json(json_data);
    // Create a direction distribution.
    DirectionBinList direction_distr = directions_from_json(json_data);
    // Create a temperature regime.
    TemperatureStepList temperature_regime = simple_cooling_from_json(json_data);

    // run model.
    if (prog_opts.real_run) {

        // Check that we have a model file.
        std::string model_file;
        if (json_data.contains("outputs")) {
            auto outputs = json_data["outputs"];
            if (outputs.contains("model")) {
                model_file = outputs["model"];
            } else {
                throw std::runtime_error("Program JSON does not contain 'outputs' 'model'");
            }
        } else {
            throw std::runtime_error("Program JSON does not contain 'outputs'");
        }

        ModelDataRecordList model_data = run_model(
                material,
                initial_temperature,
                applied_field_strength,
                applied_field_x_dir,
                applied_field_y_dir,
                applied_field_z_dir,
                applied_field_unit,
                size_distr,
                size_unit,
                elongation_distr,
                direction_distr,
                temperature_regime,
                tau0,
                epsilon,
                n_polish);

        std::cout << "Saving model run to file " << model_file << std::endl;
        save_model_output(model_data, model_file);

    } else {

        Real x = applied_field_x_dir;
        Real y = applied_field_y_dir;
        Real z = applied_field_z_dir;

        std::cout << "The cooling model will run with the following parameters:" << std::endl;
        std::cout << "Material:                              " << material << std::endl;
        std::cout << "Applied field strength:                " << applied_field_strength << std::endl;
        std::cout << "Field direction:                       " << x << ", " << y << ", " << z << std::endl;

        std::cout << std::endl;

        std::cout << "-----------  optional model parameters  ------------" << std::endl;
        std::cout << "tau0:                                  " << tau0 << std::endl;
        std::cout << "Size unit:                             " << size_unit << std::endl;
        std::cout << "Field unit:                            " << applied_field_unit << std::endl;
        std::cout << "Model epsilon value:                   " << epsilon << std::endl;
        std::cout << "No. of Newton-Raphson polishing steps: " << n_polish << std::endl;

        std::cout << std::endl;

    }

    // Deal with outputs.
    if (json_data.contains("outputs")) {
        auto outputs = json_data["outputs"];
        if (outputs.contains("sizes")) {
            auto model_sizes = outputs["sizes"];
            std::cout << "Saving size distribution to file " << model_sizes << std::endl;
            save_distribution(size_distr, model_sizes);
        }
        if (outputs.contains("elongations")) {
            auto model_elongations = outputs["elongations"];
            std::cout << "Saving elongation distribution to file " << model_elongations << std::endl;
            save_distribution(elongation_distr, model_elongations);
        }
        if (outputs.contains("directions")) {
            auto model_directions = outputs["directions"];
            std::cout << "Saving direction distribution to file " << model_directions << std::endl;
            save_direction_distribution(direction_distr, model_directions);
        }
        if (outputs.contains("cooling_regime")) {
            auto model_cooling_regime = outputs["cooling_regime"];
            std::cout << "Saving cooling regime to file " << model_cooling_regime << std::endl;
            save_temperature_regime(temperature_regime, model_cooling_regime);
        }
    }
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
                ("json-spec-file", po::value<std::string>(), "the JSON model spec file")
                ("real-run", po::bool_switch(&prog_opts.real_run), "actually run the model")
                ;

        // Positional.
        pos_desc.add("json-spec-file", 1);

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

        if (vm.count("json-spec-file")) {
            prog_opts.json_spec_file = vm["json-spec-file"].as<std::string>();
        } else {
            std::cout << "Material was not set." << std::endl;
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
