#include "Output.h"

#include "NuFitModel.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>

namespace nufit {

void writeNpy(const std::string& path, const std::vector<double>& data, size_t cols) {
  if (cols == 0 || data.size() % cols != 0) throw std::invalid_argument("writeNpy: bad shape");
  std::ostringstream h;
  h << "{'descr': '<f8', 'fortran_order': False, 'shape': (" << data.size() / cols << ", " << cols
    << "), }";
  std::string header = h.str();
  const size_t pre = 10; // magic(6) + version(2) + header length(2)
  header.append(64 - (pre + header.size() + 1) % 64, ' ');
  header.push_back('\n');
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot write " + path);
  out.write("\x93NUMPY\x01\x00", 8);
  const uint16_t len = static_cast<uint16_t>(header.size());
  const char lenBytes[2] = {static_cast<char>(len & 0xff), static_cast<char>(len >> 8)};
  out.write(lenBytes, 2);
  out << header;
  out.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
}

namespace {
std::string quote(const std::string& s) {
  std::string q = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') q.push_back('\\');
    q.push_back(c);
  }
  return q + "\"";
}
std::string number(double v) {
  if (!std::isfinite(v)) return "null";
  std::ostringstream s;
  s << std::setprecision(17) << v;
  return s.str();
}
} // namespace

Json& Json::num(const std::string& k, double v) {
  items_.emplace_back(k, number(v));
  return *this;
}
Json& Json::str(const std::string& k, const std::string& v) {
  items_.emplace_back(k, quote(v));
  return *this;
}
Json& Json::arr(const std::string& k, const std::vector<double>& v) {
  std::string s = "[";
  for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + number(v[i]);
  items_.emplace_back(k, s + "]");
  return *this;
}
Json& Json::strs(const std::string& k, const std::vector<std::string>& v) {
  std::string s = "[";
  for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + quote(v[i]);
  items_.emplace_back(k, s + "]");
  return *this;
}
Json& Json::obj(const std::string& k, const Json& v) {
  items_.emplace_back(k, v.dump());
  return *this;
}
Json& Json::objs(const std::string& k, const std::vector<Json>& v) {
  std::string s = "[";
  for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + v[i].dump();
  items_.emplace_back(k, s + "]");
  return *this;
}
std::string Json::dump() const {
  std::string s = "{";
  for (size_t i = 0; i < items_.size(); ++i)
    s += (i ? ", " : "") + quote(items_[i].first) + ": " + items_[i].second;
  return s + "}";
}
void Json::write(const std::string& path) const {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << dump() << "\n";
}

std::vector<std::string> sampleColumns() {
  std::vector<std::string> c;
  for (int p = 0; p < kNPar; ++p) c.push_back(parName(p));
  c.push_back("logL");
  c.push_back("logPrior");
  return c;
}

void writeRootTree(const std::string& path, const std::string& treeName,
                   const std::vector<std::string>& cols, const std::vector<double>& data) {
  const size_t nc = cols.size();
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "RECREATE"));
  if (!f || f->IsZombie()) throw std::runtime_error("cannot write " + path);
  auto* tree = new TTree(treeName.c_str(), "NuFIT surrogate posterior samples");
  std::vector<double> buf(nc);
  for (size_t c = 0; c < nc; ++c) tree->Branch(cols[c].c_str(), &buf[c]);

  auto find = [&](const char* name) {
    auto it = std::find(cols.begin(), cols.end(), name);
    return it == cols.end() ? -1 : static_cast<int>(it - cols.begin());
  };
  const int i3l = find("dm3l"), i21 = find("dm21");
  const bool derived = i3l >= 0 && i21 >= 0;
  double dm31 = 0, dm32 = 0;
  int ordering = 0;
  if (derived) {
    tree->Branch("dm31_eV2", &dm31);
    tree->Branch("dm32_eV2", &dm32);
    tree->Branch("ordering", &ordering);
  }
  for (size_t r = 0; r < data.size() / nc; ++r) {
    std::copy(data.begin() + r * nc, data.begin() + (r + 1) * nc, buf.begin());
    if (derived) {
      const double d3l = buf[i3l] * 1e-3, d21 = buf[i21] * 1e-5;
      ordering = d3l > 0 ? 0 : 1;
      dm31 = ordering == 0 ? d3l : d3l + d21; // IO: dm3l is Dm32, Dm31 = Dm32 + Dm21
      dm32 = ordering == 0 ? d3l - d21 : d3l;
    }
    tree->Fill();
  }
  tree->Write();
  f->Close();
}

} // namespace nufit
