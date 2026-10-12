#pragma once
// Tiny --key value command-line parser shared by the apps.

#include "NuFitModel.h"

#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

namespace nufit {

class Args {
public:
  Args(int argc, char** argv, const std::map<std::string, std::string>& defaults,
       const std::string& usage)
      : v_(defaults) {
    for (int i = 1; i < argc; ++i) {
      std::string k = argv[i];
      if (k == "-h" || k == "--help") {
        std::cout << usage << "\nOptions (default):\n";
        for (const auto& [key, val] : defaults) std::cout << "  --" << key << " (" << val << ")\n";
        std::exit(0);
      }
      if (k.rfind("--", 0) != 0 || i + 1 >= argc) throw std::invalid_argument("bad argument " + k);
      k = k.substr(2);
      if (!v_.count(k)) throw std::invalid_argument("unknown option --" + k);
      v_[k] = argv[++i];
    }
  }
  std::string str(const std::string& k) const { return v_.at(k); }
  double num(const std::string& k) const { return std::stod(v_.at(k)); }
  long long integer(const std::string& k) const { return static_cast<long long>(std::stod(v_.at(k))); }
  const std::map<std::string, std::string>& all() const { return v_; }

  ModelConfig model() const {
    ModelConfig c;
    c.dataDir = str("datadir");
    c.variant = str("variant");
    const std::string s = str("surrogate");
    if (s != "factorized" && s != "coupled") throw std::invalid_argument("--surrogate factorized|coupled");
    c.coupledT13 = s == "coupled";
    const std::string a = str("angle-prior");
    if (a != "sin2" && a != "angle") throw std::invalid_argument("--angle-prior sin2|angle");
    c.anglePrior = a == "sin2" ? AnglePrior::kSin2 : AnglePrior::kAngle;
    const std::string d = str("dm21-prior");
    if (d != "linear" && d != "log") throw std::invalid_argument("--dm21-prior linear|log");
    c.dm21Prior = d == "linear" ? Dm21Prior::kLinear : Dm21Prior::kLog;
    return c;
  }

  static std::map<std::string, std::string> modelDefaults() {
    return {{"datadir", "data/nufit"},     {"variant", "TBoff"},     {"surrogate", "factorized"},
            {"angle-prior", "sin2"},       {"dm21-prior", "linear"}};
  }

private:
  std::map<std::string, std::string> v_;
};

} // namespace nufit
