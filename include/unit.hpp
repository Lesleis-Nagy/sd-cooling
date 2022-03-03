//
// Created by L.Nagy on 13/05/2021.
//

#ifndef SD_COOLING_INCLUDE_UNIT_HPP
#define SD_COOLING_INCLUDE_UNIT_HPP

#include <exception>

#include <basic_types.hpp>
#include <constants.hpp>

#include <boost/format.hpp>

#include <nlohmann/json.hpp>
#include <utility>

class UnknownUnitException : public std::exception {
public:
    explicit UnknownUnitException(std::string unit) {
        _message = boost::str(boost::format("Unknown unit '%1%'.") % unit);
    }

    [[nodiscard]] const char *what() const noexcept override {
	    return _message.c_str();
    }

private:
    std::string _message;
};

Real field_to_amps_per_meter(Real h, const std::string &unit) {
  if (unit == "A/m") {
	return h;
  } else if (unit == "T") {
	return h / mu0;
  } else if (unit == "mT") {
	return (h * 1E-3) / mu0;
  } else if (unit == "uT") {
	return (h * 1E-6) / mu0;
  } else if (unit == "nT") {
	return (h * 1E-9) / mu0;
  } else {
	throw UnknownUnitException(unit);
  }
}

Real size_to_meter(Real s, const std::string &unit) {
  if (unit == "m") {
	return s;
  } else if (unit == "cm") {
	return s * 1E-2;
  } else if (unit == "mm") {
	return s * 1E-3;
  } else if (unit == "um") {
	return s * 1E-6;
  } else if (unit == "nm") {
	return s * 1E-9;
  } else {
	throw UnknownUnitException(unit);
  }
}

std::string size_unit_from_json(const nlohmann::json &json) {
    if (json.contains("sizes")) {
        auto sizes = json["sizes"];
        if (sizes.contains("unit")) {
            return sizes["unit"];
        } else {
            return "nm"; // Default size units used are nanometer.
        }
    } else {
        throw std::runtime_error("Program JSON missing 'sizes'");
    }
}

std::string field_unit_from_json(const nlohmann::json &json) {
    if (json.contains("applied_field")) {
        auto applied_field = json["applied_field"];
        if (applied_field.contains("unit")) {
            return applied_field["unit"];
        } else {
            return "uT"; // Default field units used are microtesla.
        }
    } else {
        throw std::runtime_error("Program JSON missing 'applied_field'");
    }
}

#endif //SD_COOLING_INCLUDE_UNIT_HPP
