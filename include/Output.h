#pragma once
// Output helpers: .npy arrays (read by the Python validation suite), a small JSON
// metadata writer, and ROOT trees of samples for downstream (pseudo-experiment) use.

#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace nufit {

// Write a row-major (rows x cols) float64 array in NumPy .npy format.
void writeNpy(const std::string& path, const std::vector<double>& data, size_t cols);

// Minimal ordered JSON object builder (numbers, strings, arrays, nested objects).
class Json {
public:
  Json& num(const std::string& k, double v);
  Json& str(const std::string& k, const std::string& v);
  Json& arr(const std::string& k, const std::vector<double>& v);
  Json& strs(const std::string& k, const std::vector<std::string>& v);
  Json& obj(const std::string& k, const Json& v);
  Json& objs(const std::string& k, const std::vector<Json>& v);
  std::string dump() const;
  void write(const std::string& path) const;

private:
  std::vector<std::pair<std::string, std::string>> items_;
};

// Write samples (rows of the given columns) as a ROOT TTree with one double branch per
// column, adding derived dm31/dm32 [eV^2] and an ordering flag (0 = NO, 1 = IO) when the
// standard parameter columns are present.
void writeRootTree(const std::string& path, const std::string& treeName,
                   const std::vector<std::string>& cols, const std::vector<double>& data);

std::vector<std::string> sampleColumns(); // s12 ... dm21, logL, logPrior

} // namespace nufit
