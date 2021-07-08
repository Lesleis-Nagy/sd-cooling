//
// Created by L.Nagy on 19/05/2021.
//

#ifndef SD_COOLING_DISTRIBUTION_HPP
#define SD_COOLING_DISTRIBUTION_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>
#include <vector>

#include <boost/math/distributions/lognormal.hpp>

#include "nlohmann/json.hpp"

#include "basic_types.hpp"
#include "constants.hpp"

#include "debug.hpp"

struct HistogramBin {
    Real value;
    Real fraction;
};

typedef std::vector<HistogramBin> HistogramBinList;

void normalize_histogram_bin_fractions(HistogramBinList &bins) {
    // Calculate the fraction sum;
    Real fraction_sum = 0;
    for (const auto & bin : bins) {
        fraction_sum += bin.fraction;
    }

    for (auto & bin : bins) {
        bin.fraction /= fraction_sum;
    }
}

/**
 * Return the log-normal probability density function based on input parameters the formula used is
 * given by
 * \f[
 *  p(x|s,l,s_c) = \frac{1}{s (x-l) \sqrt{2\pi}}\exp{\left(-\frac{\log^2{\left(\frac{x-l}{s_c}\right) }}{2s^2}\right)}
 * \f]
 * @param s the shape parameter
 * @param loc the location parameter
 * @param scale is the scaling parameter
 * @return
 */
std::function<Real(Real)>
log_normal_pdf(Real s, Real loc, Real scale) {
    return [s, loc, scale](Real x) -> Real {
        auto y = (x - loc) / scale;

        auto exp_arg_num = std::log(y) * std::log(y);
        auto exp_arg_den = 2 * s * s;

        return ((1 / (s * y * sqrt2pi)) * std::exp(-1.0 * exp_arg_num / exp_arg_den)) / scale;
    };
}

HistogramBinList
log_normal_histogram(Real x_start, Real x_end, Real s, Real loc, Real scale, size_t n) {

    HistogramBinList bins;

    // Edge case, if n == 1 just return the x_start value
    if (n == 1) {
        bins.push_back({
            x_start, 1.0
        });
        return bins;
    }

    // Create a log normal distribution.
    auto pdf = log_normal_pdf(s, loc, scale);

    // Create bins
    auto dx = (x_end - x_start) / (n - 1);
    for (size_t i = 0; i < n - 1; ++i) {
        auto x_i = x_start + i * dx;
        auto x_ip1 = x_start + (i + 1) * dx;
        auto x_mid = (x_i + x_ip1) / 2;

        auto prob_x_mid = pdf(x_mid);

        // Push back the bin's start point as the bin value.
        bins.push_back({x_i, prob_x_mid});
    }

    return bins;
}

HistogramBinList
size_distribution_from_json(const nlohmann::json &json) {

    HistogramBinList output;

    if (json.contains("sizes")) {
        auto sizes = json["sizes"];
        DEBUG_MSG("Found sizes");
        if (sizes.contains("distribution")) {
            auto distribution = sizes["distribution"];
            DEBUG_MSG("Found distribution");
            DEBUG_MSG_SIMPLE_VAR(distribution);
            if (distribution.contains("type")) {
                auto type = distribution["type"];
                DEBUG_MSG_SIMPLE_VAR(type);
                if (type == "lognormal") {
                    if (!distribution.contains("nbins")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'nbins'");
                    }
                    if (!distribution.contains("shape")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'shape'");
                    }
                    if (!distribution.contains("location")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'location'");
                    }
                    if (!distribution.contains("scale")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'scale'");
                    }
                    if (!distribution.contains("start")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'start'");
                    }
                    if (!distribution.contains("end")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'end'");
                    }

                    unsigned int nbins = distribution["nbins"];
                    Real shape = distribution["shape"];
                    Real location = distribution["location"];
                    Real scale = distribution["scale"];
                    Real start = distribution["start"];
                    Real end = distribution["end"];

                    DEBUG_MSG_SIMPLE_VAR(nbins);
                    DEBUG_MSG_SIMPLE_VAR(shape);
                    DEBUG_MSG_SIMPLE_VAR(location);
                    DEBUG_MSG_SIMPLE_VAR(scale);
                    DEBUG_MSG_SIMPLE_VAR(start);
                    DEBUG_MSG_SIMPLE_VAR(end);

                    output = log_normal_histogram(start, end, shape, location, scale, nbins);
                }
            } else {
                throw std::runtime_error("Program JSON, using distribution but no 'type' given");
            }
        } else if (sizes.contains("list")) {
            auto list = sizes["list"];
            for (const auto &item : list) {
                if (!item.contains("value")) {
                    throw std::runtime_error("Program JSON 'sizes' 'list' entry must have 'value'");
                }
                if (!item.contains("fraction")) {
                    throw std::runtime_error("Program JSON 'sizes' 'list' entry must have 'fraction'");
                }
                output.push_back({item["value"], item["fraction"]});
            }
        } else {
            throw std::runtime_error("Program JSON 'sizes' requires 'distribution' or 'list'");
        }
    } else {
        throw std::runtime_error("Program JSON missing 'sizes'");
    }

    normalize_histogram_bin_fractions(output);
    return output;
}

HistogramBinList
elongation_distribution_from_json(const nlohmann::json &json) {

    HistogramBinList output;

    if (json.contains("elongations")) {
        auto elongations = json["elongations"];
        DEBUG_MSG("Found elongations");
        if (elongations.contains("distribution")) {
            auto distribution = elongations["distribution"];
            DEBUG_MSG("Found distribution");
            DEBUG_MSG_SIMPLE_VAR(distribution);
            if (distribution.contains("type")) {
                auto type = distribution["type"];
                DEBUG_MSG_SIMPLE_VAR(type);
                if (type == "lognormal") {
                    if (!distribution.contains("nbins")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'nbins'");
                    }
                    if (!distribution.contains("shape")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'shape'");
                    }
                    if (!distribution.contains("location")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'location'");
                    }
                    if (!distribution.contains("scale")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'scale'");
                    }
                    if (!distribution.contains("start")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'start'");
                    }
                    if (!distribution.contains("end")) {
                        throw std::runtime_error("Program JSON is using log-normal distribution but no 'end'");
                    }

                    unsigned int nbins = distribution["nbins"];
                    Real shape = distribution["shape"];
                    Real location = distribution["location"];
                    Real scale = distribution["scale"];
                    Real start = distribution["start"];
                    Real end = distribution["end"];

                    DEBUG_MSG_SIMPLE_VAR(nbins);
                    DEBUG_MSG_SIMPLE_VAR(shape);
                    DEBUG_MSG_SIMPLE_VAR(location);
                    DEBUG_MSG_SIMPLE_VAR(scale);
                    DEBUG_MSG_SIMPLE_VAR(start);
                    DEBUG_MSG_SIMPLE_VAR(end);

                    output = log_normal_histogram(start, end, shape, location, scale, nbins);
                }
            } else {
                throw std::runtime_error("Program JSON, using distribution but no 'type' given");
            }
        } else if (elongations.contains("list")) {
            auto list = elongations["list"];
            for (const auto &item : list) {
                if (!item.contains("value")) {
                    throw std::runtime_error("Program JSON 'elongations' 'list' entry must have 'value'");
                }
                if (!item.contains("fraction")) {
                    throw std::runtime_error("Program JSON 'elongations' 'list' entry must have 'fraction'");
                }
                output.push_back({item["value"], item["fraction"]});
            }
        } else {
            throw std::runtime_error("Program JSON 'elongations' requires 'distribution' or 'list'");
        }
    } else {
        throw std::runtime_error("Program JSON missing 'elongations'");
    }

    normalize_histogram_bin_fractions(output);
    return output;
}

void
save_distribution(const HistogramBinList &size_distribution_records, const std::string &file_name) {
    // Write data to file.
    std::ofstream fout;
    fout.open(file_name);
    for (auto &record : size_distribution_records) {
        fout << std::setw(20) << std::setprecision(15) << record.value << ", ";
        fout << std::setw(20) << std::setprecision(15) << record.fraction;
        fout << std::endl;
    }
    fout.close();
}



#endif //SD_COOLING_DISTRIBUTION_HPP
