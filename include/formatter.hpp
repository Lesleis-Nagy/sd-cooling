//
// Created by L.Nagy on 15/05/2021.
//

#ifndef SD_COOLING_INCLUDE_FORMATTER_HPP
#define SD_COOLING_INCLUDE_FORMATTER_HPP

#include <string>
#include <stdexcept>
#include <sstream>

class Formatter {
public:
  Formatter() {}
  ~Formatter() {}

  template<typename Type>
  Formatter &operator<<(const Type &value) {
	stream_ << value;
	return *this;
  }

  std::string str() const { return stream_.str(); }
  operator std::string() const { return stream_.str(); }

  enum ConvertToString {
	to_str
  };
  std::string operator>>(ConvertToString) { return stream_.str(); }

private:
  std::stringstream stream_;

  Formatter(const Formatter &);
  Formatter &operator=(Formatter &);
};

#endif //SD_COOLING_INCLUDE_FORMATTER_HPP
