#include "NuFitModel.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace nufit {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kLn10 = std::log(10.0);

std::string trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r");
  if (b == std::string::npos) return "";
  return s.substr(b, s.find_last_not_of(" \t\r") - b + 1);
}

// Check that the delta_CP = -180 and +180 slices of a table agree (they are the same point).
void checkPeriodic(const GridTable& t, int dcpAxis, const std::string& what) {
  const auto& a = t.axis(dcpAxis);
  if (a.front() != -180 || a.back() != 180)
    throw std::runtime_error(what + ": delta_CP axis does not span [-180,180]");
  const size_t last = (a.size() - 1) * t.stride(dcpAxis);
  for (size_t k = 0; k < t.size(); ++k) {
    if ((k / t.stride(dcpAxis)) % a.size() != 0) continue;
    if (std::abs(t.values()[k] - t.values()[k + last]) > 1e-6)
      throw std::runtime_error(what + ": chi2 at delta_CP=-180 and +180 differ");
  }
}
} // namespace

const char* parName(int p) {
  static const char* names[kNPar] = {"s12", "s13", "s23", "dcp", "dm3l", "dm21"};
  return names[p];
}

NuFitFile NuFitFile::read(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  NuFitFile f;
  std::string line, tag;
  std::vector<std::string> cols;
  std::vector<std::vector<double>> rows;
  auto flush = [&]() {
    if (!tag.empty()) f.tables.emplace(tag, GridTable::fromRows(cols, rows));
    rows.clear();
  };
  while (std::getline(in, line)) {
    const std::string s = trim(line);
    if (s.empty()) continue;
    if (s.rfind("###", 0) == 0) {
      f.header = trim(s.substr(3));
      continue;
    }
    if (s[0] == '#') {
      const auto pos = s.find(" projection:");
      if (pos == std::string::npos) continue;
      flush();
      tag = trim(s.substr(1, pos - 1));
      std::istringstream ss(s.substr(pos + 12));
      cols.clear();
      for (std::string c; ss >> c;) cols.push_back(c);
      cols.pop_back(); // Delta_chi^2 column
      continue;
    }
    std::vector<double> r;
    r.reserve(4);
    const char* p = s.c_str();
    char* end = nullptr;
    for (double v = std::strtod(p, &end); end != p; v = std::strtod(p, &end)) {
      r.push_back(v);
      p = end;
    }
    rows.push_back(std::move(r));
  }
  flush();
  return f;
}

const GridTable& NuFitFile::at(const std::string& tag) const {
  auto it = tables.find(tag);
  if (it == tables.end()) throw std::runtime_error("NuFIT file has no section " + tag);
  return it->second;
}

std::string ModelConfig::describe() const {
  std::ostringstream ss;
  ss << "variant=" << variant << " surrogate=" << (coupledT13 ? "coupled" : "factorized")
     << " anglePrior=" << (anglePrior == AnglePrior::kSin2 ? "sin2" : "angle")
     << " dm21Prior=" << (dm21Prior == Dm21Prior::kLinear ? "linear" : "log");
  return ss.str();
}

NuFitModel::NuFitModel(const ModelConfig& cfg) : cfg_(cfg) {
  if (cfg.variant != "TBoff" && cfg.variant != "TByes")
    throw std::invalid_argument("variant must be TBoff or TByes");
  for (int o = 0; o < 2; ++o) {
    const std::string path =
        cfg.dataDir + "/v60.release-" + cfg.variant + (o == kNO ? "-NO.txt" : "-IO.txt");
    NuFitFile f = NuFitFile::read(path);
    t3_[o] = f.at("T23/DMA/DCP");
    t2_[o] = f.at("T12/DMS");
    t1_[o] = f.at("T13");
    t1323_[o] = f.at("T13/T23");
    t23_[o] = f.at("T23");
    min_[o] = t1_[o].minValue();
    for (const GridTable* t : {&t3_[o], &t2_[o], &t1323_[o], &t23_[o]})
      if (std::abs(t->minValue() - min_[o]) > 1e-6)
        throw std::runtime_error(path + ": projections do not share a common minimum");
    checkPeriodic(t3_[o], 2, path + " T23/DMA/DCP");
  }
  if (t3_[kNO].axis(1).front() <= 0 || t3_[kIO].axis(1).back() >= 0)
    throw std::runtime_error("unexpected sign of the Delta m^2_3l axes");

  lo_[kS12] = t2_[kNO].axis(0).front();
  hi_[kS12] = t2_[kNO].axis(0).back();
  lo_[kS13] = t1_[kNO].axis(0).front();
  hi_[kS13] = t1_[kNO].axis(0).back();
  lo_[kS23] = t3_[kNO].axis(0).front();
  hi_[kS23] = t3_[kNO].axis(0).back();
  lo_[kDCP] = -180;
  hi_[kDCP] = 180;
  lo_[kDM3L] = t3_[kNO].axis(1).front();
  hi_[kDM3L] = t3_[kNO].axis(1).back();
  if (-t3_[kIO].axis(1).back() != lo_[kDM3L] || -t3_[kIO].axis(1).front() != hi_[kDM3L])
    throw std::runtime_error("NO and IO Delta m^2_3l windows are not mirror images");
  lo_[kDM21] = std::pow(10.0, t2_[kNO].axis(1).front()) / 1e-5;
  hi_[kDM21] = std::pow(10.0, t2_[kNO].axis(1).back()) / 1e-5;
}

double NuFitModel::wrapDeg(double dcp) {
  double w = std::fmod(dcp + 180.0, 360.0);
  if (w < 0) w += 360.0;
  return w - 180.0;
}

bool NuFitModel::inSupport(const double* p) const {
  for (int k : {kS12, kS13, kS23, kDM21})
    if (!(p[k] >= lo_[k] && p[k] <= hi_[k])) return false;
  const double a = std::abs(p[kDM3L]);
  return a >= lo_[kDM3L] && a <= hi_[kDM3L] && std::isfinite(p[kDCP]);
}

double NuFitModel::deltaChi2(const double* p) const {
  if (!inSupport(p)) return kNaN;
  const Ordering o = ordering(p);
  const double x3[3] = {p[kS23], p[kDM3L], wrapDeg(p[kDCP])};
  const double x2[2] = {p[kS12], std::log10(p[kDM21] * 1e-5)};
  double chi2 = t3_[o].eval(x3) + (t2_[o].eval(x2) - min_[o]);
  if (cfg_.coupledT13) {
    const double x13[2] = {p[kS13], p[kS23]};
    chi2 += t1323_[o].eval(x13) - t23_[o].eval1(p[kS23]);
  } else {
    chi2 += t1_[o].eval1(p[kS13]) - min_[o];
  }
  return chi2;
}

double NuFitModel::logPriorAngle(double s) const {
  if (cfg_.anglePrior == AnglePrior::kSin2) return 0;
  // theta = asin(sqrt(s)); d theta / ds = 1 / (2 sqrt(s (1-s)))
  return -0.5 * std::log(s * (1 - s)) - std::log(2.0);
}

double NuFitModel::logPriorLog10Dm21(double l) const {
  // density in l = log10(dm21/eV^2): flat in l (log prior) or d dm21/dl ~ 10^l (linear prior)
  return cfg_.dm21Prior == Dm21Prior::kLog ? 0 : l * kLn10;
}

double NuFitModel::logPrior(const double* p) const {
  if (!inSupport(p)) return -std::numeric_limits<double>::infinity();
  double lp = logPriorAngle(p[kS12]) + logPriorAngle(p[kS13]) + logPriorAngle(p[kS23]);
  if (cfg_.dm21Prior == Dm21Prior::kLog) lp -= std::log(p[kDM21]);
  return lp;
}

} // namespace nufit
