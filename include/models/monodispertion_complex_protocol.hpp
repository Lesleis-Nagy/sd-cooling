//
// Created by lnagy2 on 08/02/2022.
//

#ifndef SD_COOLING_MONODISPERTION_COMPLEX_PROTOCOL_HPP
#define SD_COOLING_MONODISPERTION_COMPLEX_PROTOCOL_HPP

#include <exception>
#include <sys/stat.h>

#include <tuple>
#include <iostream>
#include <iomanip>

#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <utility>

#include "basic_types.hpp"
#include "unit.hpp"
// #include "experiment_regime.hpp"

class ProtocolFileNotFoundException: std::exception {
public:
    /**
     * Retrieve the error string.
     * @return the error string.
     */
    [[nodiscard]] const char *what() const noexcept override {
        return "File could not be found.";
    }
};

class ProtocolDataInvalidException : std::exception {
public:

    /**
     * Constructor creates a new exception with the given message.
     * @param message
     */
    ProtocolDataInvalidException(std::string message) {
        _message = std::move(message);
    }

    /**
     * Retrieve the message associated with this exception.
     * @return exception message.
     */
    [[nodiscard]] const std::string &message() const { return _message; }

    /**
     * Retrieve the error string.
     * @return the error string.
     */
    [[nodiscard]] const char *what() const noexcept override {
        return "File could not be found.";
    }

private:

    std::string _message;
};

struct Sample {
    std::string name;
};

struct InputMaterial {

    std::string name;

};

struct Geometry {
    std::string type;
    double size;
    std::string size_unit;
    double elongation; // As a percentage
};

struct Outputs {

    bool produce_trm_acquisistion_protocol_csv;
    std::string trm_acquisition_protocol_csv;

    bool produce_lab_demag_protocol_csv;
    std::string lab_demag_protocol_csv;

    std::string raw_simulation_data_csv;

};

struct Field {

    double strength;     // *MUST* be in A/m
    double x_direction;
    double y_direction;
    double z_direction;

};

struct TRMAcquisition {

    double allowable_fractional_drop;
    double ambient_temperature;
    double initial_temperature;
    double temperature_at_t1;
    double t1;
    Field field;

};

struct Step {

    std::string type;
    int type_index;
    double bake_temperature;
    double ambient_temperature;
    double stop_temperature;
    double heating_time;
    double bake_time;
    double cool_time;
    Field field;

};

struct LabDemagProtocol {

    double allowable_fractional_drop;
    std::vector<Step> steps;

};

struct Protocol {

    static Protocol from_file(const std::string& file_name);

    Sample sample;

    InputMaterial material;

    Geometry geometry;

    Outputs outputs;

    bool has_trm_acquisition;
    TRMAcquisition trm_acquisition;

    bool has_lab_demag_protocol;
    LabDemagProtocol lab_demag_protocol;

};

std::pair<bool, std::string> check_field(const ryml::NodeRef &field) {
    if (!field.has_child("strength")) {
        return {false, "strength"};
    }

    if (!field.has_child("unit")) {
        return {false, "unit"};
    }

    if (!field.has_child("direction")) {
        return {false, "direction"};
    } else {
        std::vector<double> directions;

        try {
            field["direction"] >> directions;
        } catch (...) {
            return {false, "direction"};
        }

        if (directions.size() != 3) {
            return {false, "direction"};
        }
    }

    return {true, ""};
}

Protocol read_protocol(const ryml::Tree &protocol) {
    using boost::format;
    using boost::str;

    // Initialize the output protocol object.
    Protocol output_protocol{};

    // This part checks that our protocol is valid.
    ryml::NodeRef root = protocol.rootref();

    if (!root.has_child("outputs")) {
        throw ProtocolDataInvalidException("The 'outputs' tag is missing.");
    } else {
        ryml::NodeRef outputs = root["outputs"];

        if (outputs.has_child("trm-acquisition-protocol")) {
            output_protocol.outputs.produce_trm_acquisistion_protocol_csv = true;
            outputs["trm-acquisition-protocol"] >> output_protocol.outputs.trm_acquisition_protocol_csv;
        }

        if (outputs.has_child("lab-demag-protocol")) {
            output_protocol.outputs.produce_lab_demag_protocol_csv = true;
            outputs["lab-demag-protocol"] >> output_protocol.outputs.lab_demag_protocol_csv;
        }

        if (!outputs.has_child("raw-simulation-data")) {
            throw ProtocolDataInvalidException("The 'outputs' tag is missing 'raw-simulation-data'.");
        } else {
            outputs["raw-simulation-data"] >> output_protocol.outputs.raw_simulation_data_csv;
        }
    }

    if (!root.has_child("sample-name")) {
        throw ProtocolDataInvalidException("The 'sample-name' tag is missing.");
    } else {
        root["sample-name"] >> output_protocol.sample.name;
    }

    if (!root.has_child("material")) {
        throw ProtocolDataInvalidException("The 'material' tag is missing.");
    } else {
        root["material"] >> output_protocol.material.name;
    }

    if (!root.has_child("geometry")) {
        throw ProtocolDataInvalidException("The 'geometry' tag is missing.");
    } else {
        ryml::NodeRef geometry = root["geometry"];

        if (!geometry.has_child("type")) {
            throw ProtocolDataInvalidException("The 'geometry' tag is missing 'type'.");
        } else {
            geometry["type"] >> output_protocol.geometry.type;
        }

        if (!geometry.has_child("size")) {
            throw ProtocolDataInvalidException("The 'geometry' tag is missing 'size'.");
        } else {
            geometry["size"] >> output_protocol.geometry.size;
        }

        if (!geometry.has_child("size-unit")) {
            throw ProtocolDataInvalidException("The 'geometry' tag is missing 'size'.");
        } else {
            geometry["size-unit"] >> output_protocol.geometry.size_unit;
        }

        if (!geometry.has_child("elongation")) {
            throw ProtocolDataInvalidException("The 'geometry' tag is missing 'elongation'.");
        } else {
            geometry["elongation"] >> output_protocol.geometry.elongation;
        }

    }

    // The 'trm-acquisition' tag is optional, but if it is present it must contain the child tags.
    if (root.has_child("trm-acquisition")) {
        ryml::NodeRef trm_acquisition = root["trm-acquisition"];
        output_protocol.has_trm_acquisition = true;

        if (!trm_acquisition.has_child("allowable-fractional-drop")) {
            throw ProtocolDataInvalidException("The 'trm-acquisition' tag is missing 'allowable-fractional-drop'.");
        } else {
            trm_acquisition["allowable-fractional-drop"] >> output_protocol.trm_acquisition.allowable_fractional_drop;
        }

        if (!trm_acquisition.has_child("ambient-temperature")) {
            throw ProtocolDataInvalidException("The 'trm-acquisition' tag is missing 'ambient-temperature'.");
        } else {
            trm_acquisition["ambient-temperature"] >> output_protocol.trm_acquisition.ambient_temperature;
        }

        if (!trm_acquisition.has_child("initial-temperature")) {
            throw ProtocolDataInvalidException("The 'trm-acquisition' tag is missing 'initial-temperature'.");
        } else {
            trm_acquisition["initial-temperature"] >> output_protocol.trm_acquisition.initial_temperature;
        }

        if (!trm_acquisition.has_child("temperature-at-t1")) {
            throw ProtocolDataInvalidException("The 'trm-acquisition' tag is missing 'temperature-at-t1'.");
        } else {
            trm_acquisition["temperature-at-t1"] >> output_protocol.trm_acquisition.temperature_at_t1;
        }

        if (!trm_acquisition.has_child("t1")) {
            throw ProtocolDataInvalidException("The 'trm-acquisition' tag is missing 't1'.");
        } else {
            trm_acquisition["t1"] >> output_protocol.trm_acquisition.t1;
        }

        if (!trm_acquisition.has_child("field")) {
            throw ProtocolDataInvalidException("The 'trm-acquisition' tag is missing 'field'.");
        } else {
            ryml::NodeRef field = trm_acquisition["field"];

            bool is_valid;
            std::string tag_name;

            std::tie(is_valid, tag_name) = check_field(field);

            if (!is_valid) {
                throw ProtocolDataInvalidException(str(format("The 'field' tag in 'trm-acquisition' has an invalid/missing '%1%' tag.") % tag_name));
            } else {
                std::string unit;
                double strength;
                std::vector<double> direction;

                field["unit"] >> unit;
                field["strength"] >> strength;
                field["direction"] >> direction;

                double length = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

                output_protocol.trm_acquisition.field.strength = (double)field_to_amps_per_meter((Real)strength, unit);
                output_protocol.trm_acquisition.field.x_direction = direction[0]/length;
                output_protocol.trm_acquisition.field.y_direction = direction[1]/length;
                output_protocol.trm_acquisition.field.z_direction = direction[2]/length;
            }
        }
    }

    // The 'lab-demag-protocol' is optional, but if it is present it must contain the correct tags.
    if (root.has_child("lab-demag-protocol")) {
        ryml::NodeRef lab_demag_protocol = root["lab-demag-protocol"];
        output_protocol.has_lab_demag_protocol = true;

        if (!lab_demag_protocol.has_child("allowable-fractional-drop")) {
            throw ProtocolDataInvalidException("The 'lab-demag-protocol' tag is missing 'allowable-fractional-drop'.");
        } else {
            lab_demag_protocol["allowable-fractional-drop"] >> output_protocol.lab_demag_protocol.allowable_fractional_drop;
        }

        // The 'defaults' tag is optional, but if it is not present, we expect subsequent steps to provide necessary
        // information.

        bool lab_demag_protocol_has_default_ambient_temperature = false;
        double default_ambient_temperature;

        bool lab_demag_protocol_has_default_stop_temperature = false;
        double default_stop_temperature;

        bool lab_demag_protocol_has_default_heating_time = false;
        double default_heating_time;

        bool lab_demag_protocol_has_default_bake_time = false;
        double default_bake_time;

        bool lab_demag_protocol_has_default_cool_time = false;
        double default_cool_time;

        bool lab_demag_protocol_has_default_field = false;
        Field default_lab_field{0.0, 1.0, 1.0, 1.0};

        if (lab_demag_protocol.has_child("defaults")) {
            ryml::NodeRef defaults = lab_demag_protocol["defaults"];

            if (defaults.has_child("ambient-temperature")) {
                lab_demag_protocol_has_default_ambient_temperature = true;
                defaults["ambient-temperature"] >> default_ambient_temperature;
            }

            if (defaults.has_child("stop-temperature")) {
                lab_demag_protocol_has_default_stop_temperature = true;
                defaults["stop-temperature"] >> default_stop_temperature;
            }

            if (defaults.has_child("heating-time")) {
                lab_demag_protocol_has_default_heating_time = true;
                defaults["heating-time"] >> default_heating_time;
            }

            if (defaults.has_child("bake-time")) {
                lab_demag_protocol_has_default_bake_time = true;
                defaults["bake-time"] >> default_bake_time;
            }

            if (defaults.has_child("cool-time")) {
                lab_demag_protocol_has_default_cool_time = true;
                defaults["cool-time"] >> default_cool_time;
            }

            if (defaults.has_child("field")) {

                ryml::NodeRef field = defaults["field"];

                bool is_valid;
                std::string tag_name;

                std::tie(is_valid, tag_name) = check_field(field);

                if (!is_valid) {
                    throw ProtocolDataInvalidException(str(format("The 'field' tag in 'defaults' of 'lab-demag-protocol' has an invalid/missing '%1%' tag.") % tag_name));
                } else {
                    lab_demag_protocol_has_default_field = true;

                    std::string unit;
                    double strength;
                    std::vector<double> direction;

                    field["unit"] >> unit;
                    field["strength"] >> strength;
                    field["direction"] >> direction;

                    double length = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

                    default_lab_field.strength = (double)field_to_amps_per_meter((Real)strength, unit);
                    default_lab_field.x_direction = direction[0]/length;
                    default_lab_field.y_direction = direction[1]/length;
                    default_lab_field.z_direction = direction[2]/length;
                }
            }
        }

        // Demag protocol steps.
        if (!lab_demag_protocol.has_child("steps")) {
            throw ProtocolDataInvalidException("The 'steps' tag is missing from 'lab-demag-protocol'.");
        } else {
            ryml::NodeRef steps = lab_demag_protocol["steps"];
            int step_counter = 0;
            for (auto step : steps) {

                Step protocol_step;

                step_counter += 1;

                if (!step.has_child("type")) {
                    throw ProtocolDataInvalidException(str(format("Step %1% in the 'steps' tag of 'lab-demag-protocol' has no 'type'.")));
                }

                if (!step.has_child("bake-temperature")) {
                    throw ProtocolDataInvalidException(str(format("Step %1% in the 'steps' tag of 'lab-demag-protocol' has no 'bake-temperature'.")));
                }

                ryml::NodeRef type = step["type"];
                if (type.val() == "zero-field") {
                    step["type"] >> protocol_step.type;
                    step["bake-temperature"] >> protocol_step.bake_temperature;
                    protocol_step.type_index = 0;

                    if (step.has_child("ambient-temperature")) {
                        step["ambient-temperature"] >> protocol_step.ambient_temperature;
                    } else if (lab_demag_protocol_has_default_ambient_temperature) {
                        protocol_step.ambient_temperature = default_ambient_temperature;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'ambient-temperature' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("stop-temperature")) {
                        step["stop-temperature"] >> protocol_step.stop_temperature;
                    } else if (lab_demag_protocol_has_default_stop_temperature) {
                        protocol_step.stop_temperature = default_stop_temperature;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'stop-temperature' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("heating-time")) {
                        step["heating-time"] >> protocol_step.heating_time;
                    } else if (lab_demag_protocol_has_default_heating_time) {
                        protocol_step.heating_time = default_heating_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'heating-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("bake-time")) {
                        step["bake-time"] >> protocol_step.bake_time;
                    } else if (lab_demag_protocol_has_default_bake_time) {
                        protocol_step.bake_time = default_bake_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'bake-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("cool-time")) {
                        step["cool-time"] >> protocol_step.cool_time;
                    } else if (lab_demag_protocol_has_default_cool_time) {
                        protocol_step.cool_time = default_cool_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'cool-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    protocol_step.field.strength = 0;
                    protocol_step.field.x_direction = 1.0 / sqrt(3.0);
                    protocol_step.field.y_direction = 1.0 / sqrt(3.0);
                    protocol_step.field.z_direction = 1.0 / sqrt(3.0);

                } else if (type.val() == "in-field") {

                    step["type"] >> protocol_step.type;
                    step["bake-temperature"] >> protocol_step.bake_temperature;
                    protocol_step.type_index = 1;

                    if (step.has_child("ambient-temperature")) {
                        step["ambient-temperature"] >> protocol_step.ambient_temperature;
                    } else if (lab_demag_protocol_has_default_ambient_temperature) {
                        protocol_step.ambient_temperature = default_ambient_temperature;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'ambient-temperature' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("stop-temperature")) {
                        step["stop-temperature"] >> protocol_step.stop_temperature;
                    } else if (lab_demag_protocol_has_default_stop_temperature) {
                        protocol_step.stop_temperature = default_stop_temperature;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'stop-temperature' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("heating-time")) {
                        step["heating-time"] >> protocol_step.heating_time;
                    } else if (lab_demag_protocol_has_default_heating_time) {
                        protocol_step.heating_time = default_heating_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'heating-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("bake-time")) {
                        step["bake-time"] >> protocol_step.bake_time;
                    } else if (lab_demag_protocol_has_default_bake_time) {
                        protocol_step.bake_time = default_bake_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'bake-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("cool-time")) {
                        step["cool-time"] >> protocol_step.cool_time;
                    } else if (lab_demag_protocol_has_default_cool_time) {
                        protocol_step.cool_time = default_cool_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'cool-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("field")) {

                        ryml::NodeRef field = step["field"];

                        bool is_valid;
                        std::string tag_name;

                        std::tie(is_valid, tag_name) = check_field(field);

                        if (is_valid) {
                            std::string unit;
                            double strength;
                            std::vector<double> direction;

                            field["unit"] >> unit;
                            field["strength"] >> strength;
                            field["direction"] >> direction;

                            double length = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

                            protocol_step.field.strength = (double)field_to_amps_per_meter((Real)strength, unit);
                            protocol_step.field.x_direction = direction[0]/length;
                            protocol_step.field.y_direction = direction[1]/length;
                            protocol_step.field.z_direction = direction[2]/length;
                        } else {
                            throw ProtocolDataInvalidException(str(format("Lab demag step %1% is invalid.") % step_counter));
                        }
                    } else if (lab_demag_protocol_has_default_field) {
                        protocol_step.field.strength = default_lab_field.strength;
                        protocol_step.field.x_direction = default_lab_field.x_direction;
                        protocol_step.field.y_direction = default_lab_field.y_direction;
                        protocol_step.field.z_direction = default_lab_field.z_direction;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'field' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                } else if (type.val() == "pTRM-check") {

                    step["type"] >> protocol_step.type;
                    step["bake-temperature"] >> protocol_step.bake_temperature;
                    protocol_step.type_index = 2;

                    if (step.has_child("ambient-temperature")) {
                        step["ambient-temperature"] >> protocol_step.ambient_temperature;
                    } else if (lab_demag_protocol_has_default_ambient_temperature) {
                        protocol_step.ambient_temperature = default_ambient_temperature;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'ambient-temperature' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("stop-temperature")) {
                        step["stop-temperature"] >> protocol_step.stop_temperature;
                    } else if (lab_demag_protocol_has_default_stop_temperature) {
                        protocol_step.stop_temperature = default_stop_temperature;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'stop-temperature' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("heating-time")) {
                        step["heating-time"] >> protocol_step.heating_time;
                    } else if (lab_demag_protocol_has_default_heating_time) {
                        protocol_step.heating_time = default_heating_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'heating-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("bake-time")) {
                        step["bake-time"] >> protocol_step.bake_time;
                    } else if (lab_demag_protocol_has_default_bake_time) {
                        protocol_step.bake_time = default_bake_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'bake-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("cool-time")) {
                        step["cool-time"] >> protocol_step.cool_time;
                    } else if (lab_demag_protocol_has_default_cool_time) {
                        protocol_step.cool_time = default_cool_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'cool-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("field")) {

                        ryml::NodeRef field = step["field"];

                        bool is_valid;
                        std::string tag_name;

                        std::tie(is_valid, tag_name) = check_field(field);

                        if (is_valid) {
                            std::string unit;
                            double strength;
                            std::vector<double> direction;

                            field["unit"] >> unit;
                            field["strength"] >> strength;
                            field["direction"] >> direction;

                            double length = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);

                            protocol_step.field.strength = (double)field_to_amps_per_meter((Real)strength, unit);
                            protocol_step.field.x_direction = direction[0]/length;
                            protocol_step.field.y_direction = direction[1]/length;
                            protocol_step.field.z_direction = direction[2]/length;
                        } else {
                            throw ProtocolDataInvalidException(str(format("Lab demag step %1% is invalid.") % step_counter));
                        }
                    } else if (lab_demag_protocol_has_default_field) {
                        protocol_step.field.strength = default_lab_field.strength;
                        protocol_step.field.x_direction = default_lab_field.x_direction;
                        protocol_step.field.y_direction = default_lab_field.y_direction;
                        protocol_step.field.z_direction = default_lab_field.z_direction;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'field' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                } else if (type.val() == "pTRM-tail-check") {

                    step["type"] >> protocol_step.type;
                    step["bake-temperature"] >> protocol_step.bake_temperature;
                    protocol_step.type_index = 3;

                    if (step.has_child("ambient-temperature")) {
                        step["ambient-temperature"] >> protocol_step.ambient_temperature;
                    } else if (lab_demag_protocol_has_default_ambient_temperature) {
                        protocol_step.ambient_temperature = default_ambient_temperature;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'ambient-temperature' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("stop-temperature")) {
                        step["stop-temperature"] >> protocol_step.stop_temperature;
                    } else if (lab_demag_protocol_has_default_stop_temperature) {
                        protocol_step.stop_temperature = default_stop_temperature;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'stop-temperature' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("heating-time")) {
                        step["heating-time"] >> protocol_step.heating_time;
                    } else if (lab_demag_protocol_has_default_heating_time) {
                        protocol_step.heating_time = default_heating_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'heating-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("bake-time")) {
                        step["bake-time"] >> protocol_step.bake_time;
                    } else if (lab_demag_protocol_has_default_bake_time) {
                        protocol_step.bake_time = default_bake_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'bake-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    if (step.has_child("cool-time")) {
                        step["cool-time"] >> protocol_step.cool_time;
                    } else if (lab_demag_protocol_has_default_cool_time) {
                        protocol_step.cool_time = default_cool_time;
                    } else {
                        throw ProtocolDataInvalidException(str(format("Lab demag step %1% is a '%2%' step but 'cool-time' is not defined (either in the step or in 'defaults').") % step_counter % protocol_step.type));
                    }

                    protocol_step.field.strength = 0;
                    protocol_step.field.x_direction = 1.0 / sqrt(3.0);
                    protocol_step.field.y_direction = 1.0 / sqrt(3.0);
                    protocol_step.field.z_direction = 1.0 / sqrt(3.0);

                } else {

                    throw ProtocolDataInvalidException(str(format("Unknown step %1%") % step.val()));

                }

                output_protocol.lab_demag_protocol.steps.push_back(protocol_step);

            }
        }
    }

    return output_protocol;
}

Protocol Protocol::from_file(const std::string& file_name) {

    struct stat buf{};
    if (stat(file_name.c_str(), &buf) == -1) {
        throw ProtocolFileNotFoundException();
    }

    std::ifstream fin(file_name);
    std::stringstream buffer;
    buffer << fin.rdbuf();

    ryml::Tree protocol = ryml::parse_in_arena(ryml::to_csubstr(buffer.str()));

    return read_protocol(protocol);

}

std::ostream& operator << (std::ostream& out, const Protocol& protocol) {
    out << "+----------------------+" << std::endl;
    out << "| Protocol information |" << std::endl;
    out << "+----------------------+" << std::endl;
    out << "Output files:" << std::endl;

    if (!protocol.outputs.produce_trm_acquisistion_protocol_csv) {
        std::cout << "\tThe TRM acquisition protocol will not be produced." << std::endl;
    } else {
        std::cout << "\tTRM acquisition protocol will be written to: " << protocol.outputs.trm_acquisition_protocol_csv << std::endl;
    }

    if (!protocol.outputs.produce_lab_demag_protocol_csv) {
        std::cout << "\tThe lab demag protocol will not be produced." << std::endl;
    } else {
        std::cout << "\tLab demag protocol will be written to: " << protocol.outputs.lab_demag_protocol_csv << std::endl;
    }

    std::cout << "\tRaw simulation data will be written to: " << protocol.outputs.raw_simulation_data_csv << std::endl;

    out << "Material: " << protocol.material.name << std::endl;

    if (protocol.has_trm_acquisition) {
        out << "TRM acquisition information" << std::endl;
        out << "\tAllowable fractional drop: " << protocol.trm_acquisition.allowable_fractional_drop << std::endl;
        out << "\tAmbient temperature: " << protocol.trm_acquisition.ambient_temperature << std::endl;
        out << "\tInitial temperature: " << protocol.trm_acquisition.initial_temperature << std::endl;
        out << "\tTemperature @ t1: " << protocol.trm_acquisition.temperature_at_t1 << std::endl;
        out << "\tt1: " << protocol.trm_acquisition.t1 << std::endl;
        out << "\tField: <" << protocol.trm_acquisition.field.x_direction << ", "
                            << protocol.trm_acquisition.field.y_direction << ", "
                            << protocol.trm_acquisition.field.z_direction << "> @ "
                            << protocol.trm_acquisition.field.strength << " A/m" << std::endl;
    } else {
        out << "No TRM acquisition information" << std::endl;
    }

    if (protocol.has_lab_demag_protocol) {
        out << "TRM lab demag protocol" << std::endl;
        out << "\tAllowable fractional drop: " << protocol.lab_demag_protocol.allowable_fractional_drop << std::endl;
        for (const auto& step : protocol.lab_demag_protocol.steps) {
            out << "\tStep type: " << step.type << " @ " << step.bake_temperature << " C" << std::endl;
            out << "\t\tStep type index: " << step.type_index << std::endl;
            out << "\t\tAmbient temperature: " << step.ambient_temperature << std::endl;
            out << "\t\tStop temperature: " << step.stop_temperature << std::endl;
            out << "\t\theating time: " << step.heating_time << std::endl;
            out << "\t\tbake time: " << step.bake_time << std::endl;
            out << "\t\tcool time: " << step.cool_time << std::endl;
            out << "\t\tField: <" << step.field.x_direction << ", "
                << step.field.y_direction << ", "
                << step.field.z_direction << "> @ "
                << step.field.strength << " A/m" << std::endl;

        }
    } else {
        out << "No lab demag protocol information" << std::endl;
    }

    return out;
}

#endif //SD_COOLING_MONODISPERTION_COMPLEX_PROTOCOL_HPP
