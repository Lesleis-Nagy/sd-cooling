//
// Created by lnagy2 on 20/05/2021.
//

#ifndef SD_COOLING_FIBONACCI_HPP
#define SD_COOLING_FIBONACCI_HPP

#include <vector>

#include "nlohmann/json.hpp"

#include "basic_types.hpp"
#include "constants.hpp"

struct DirectionBin {
    Real dir_x;
    Real dir_y;
    Real dir_z;
    Real fraction;
};

typedef std::vector<DirectionBin> DirectionBinList;

void normalize_direction_bin_fractions(DirectionBinList &bins) {
    // Calculate the fraction sum;
    Real fraction_sum = 0;
    for (const auto & bin : bins) {
        fraction_sum += bin.fraction;
    }

    for (auto & bin : bins) {
        bin.fraction /= fraction_sum;
    }
}

DirectionBinList
uniform_fibonacci_sphere(size_t m) {
    DirectionBinList points;

    // Edge case, if m == 1
    if (m == 1) {
        points.push_back({1, 0, 0, 1});
        return points;
    }

    // If m > 1,
    for (size_t j = 0; j < m; ++j) {
        Real i = (Real) (j + 1);
        Real n = (Real) m;
        Real theta = (2.0 * pi * j) / PHI;
        Real psi = std::acos(1 - ((2 * i) / n) + (1 / n));

        points.push_back(
                {
                        std::cos(theta) * std::sin(psi),
                        std::sin(theta) * std::sin(psi),
                        std::cos(psi),
                        1.0 / ((Real) m)
                }
        );
    }

    return points;
}

DirectionBinList
directions_from_json(const nlohmann::json &json) {

    DirectionBinList output;

    if (json.contains("directions")) {
        auto directions = json["directions"];
        if (directions.contains("distribution")) {
            auto distribution = directions["distribution"];
            if (distribution.contains("type")) {
                auto type = distribution["type"];
                if (type == "fibonacci") {
                    if (!distribution.contains("nbins")) {
                        throw std::runtime_error("Program JSON is using fibonacci direction dist. but no 'nbins'");
                    }

                    unsigned int nbins = distribution["nbins"];
                    DEBUG_MSG_SIMPLE_VAR(nbins);

                    output = uniform_fibonacci_sphere(nbins);
                }
            }
        } else if (directions.contains("list")) {
            auto list = directions["list"];
            for (const auto &item : list) {
                if (!item.contains("value")) {
                    throw std::runtime_error("Program JSON 'sizes' 'list' entry must have 'value'");
                }
                if (!item.contains("fraction")) {
                    throw std::runtime_error("Program JSON 'sizes' 'list' entry must have 'fraction'");
                }
                if (item["value"].size() != 3) {
                    throw std::runtime_error("Program JSON 'sizes' 'list' 'value' entry must 3 components");
                }
                output.push_back({item["value"][0], item["value"][1], item["value"][2], item["fraction"]});
            }
        } else {
            throw std::runtime_error("Program JSON 'directions' requires 'distribution' or 'list'");
        }
    } else {
        throw std::runtime_error("Program JSON is missing 'directions'");
    }

    normalize_direction_bin_fractions(output);
    return output;
}

void save_direction_distribution(
        const DirectionBinList &direction_distribution_records,
        const std::string &file_name
) {
    // Write data to file.
    std::ofstream fout;
    fout.open(file_name);
    for (auto &record : direction_distribution_records) {
        fout << std::setw(20) << std::setprecision(15) << record.dir_x << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.dir_y << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.dir_z << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.fraction;
        fout << std::endl;
    }
    fout.close();
}

#endif //SD_COOLING_FIBONACCI_HPP
