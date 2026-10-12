#pragma once
// Dense Delta chi^2 table on a complete rectilinear grid (1-3 dimensions),
// evaluated by N-linear interpolation.
//
// N-linear interpolation is exact at the grid nodes, never over/undershoots the
// node values, and attains its extrema over any cell at the cell corners. The
// last property is what makes the exact rejection sampler in DirectSampler.h work.

#include <cstddef>
#include <string>
#include <vector>

namespace nufit {

class GridTable {
public:
  GridTable() = default;
  // values are row-major with the last axis varying fastest.
  GridTable(std::vector<std::string> names, std::vector<std::vector<double>> axes,
            std::vector<double> values);

  // Build from scattered rows (x_0 ... x_{n-1}, value); throws unless the rows
  // form a complete rectilinear grid with each node given exactly once.
  static GridTable fromRows(std::vector<std::string> names,
                            const std::vector<std::vector<double>>& rows);

  int ndim() const { return static_cast<int>(axes_.size()); }
  const std::vector<double>& axis(int d) const { return axes_[d]; }
  const std::string& name(int d) const { return names_[d]; }
  const std::vector<double>& values() const { return values_; }
  size_t size() const { return values_.size(); }
  size_t stride(int d) const { return strides_[d]; }
  double node(const int* idx) const { return values_[flatIndex(idx)]; }
  size_t flatIndex(const int* idx) const;
  double minValue() const { return minValue_; }

  // Interpolated value; NaN outside the grid (no extrapolation).
  double eval(const double* x) const;
  double eval1(double x) const { return eval(&x); }

  // Lower node index i and fractional position t in [0,1] of x along axis d.
  // Returns false outside [axis.front(), axis.back()].
  bool locate(int d, double x, int& i, double& t) const;

private:
  std::vector<std::string> names_;
  std::vector<std::vector<double>> axes_;
  std::vector<size_t> strides_;
  std::vector<double> values_;
  double minValue_ = 0;
};

} // namespace nufit
