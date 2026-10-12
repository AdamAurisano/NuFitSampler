#include "GridTable.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nufit {

GridTable::GridTable(std::vector<std::string> names, std::vector<std::vector<double>> axes,
                     std::vector<double> values)
    : names_(std::move(names)), axes_(std::move(axes)), values_(std::move(values)) {
  if (axes_.empty() || axes_.size() > 3) throw std::invalid_argument("GridTable: 1-3 axes supported");
  strides_.assign(axes_.size(), 1);
  size_t n = 1;
  for (int d = ndim() - 1; d >= 0; --d) {
    if (axes_[d].size() < 2) throw std::invalid_argument("GridTable: axis needs >= 2 nodes");
    if (!std::is_sorted(axes_[d].begin(), axes_[d].end()) ||
        std::adjacent_find(axes_[d].begin(), axes_[d].end()) != axes_[d].end())
      throw std::invalid_argument("GridTable: axis must be strictly increasing");
    strides_[d] = n;
    n *= axes_[d].size();
  }
  if (n != values_.size()) throw std::invalid_argument("GridTable: value count does not match axes");
  minValue_ = *std::min_element(values_.begin(), values_.end());
}

GridTable GridTable::fromRows(std::vector<std::string> names,
                              const std::vector<std::vector<double>>& rows) {
  if (rows.empty()) throw std::invalid_argument("GridTable::fromRows: no rows");
  const int nd = static_cast<int>(rows[0].size()) - 1;
  std::vector<std::vector<double>> axes(nd);
  for (int d = 0; d < nd; ++d) {
    auto& a = axes[d];
    a.reserve(rows.size());
    for (const auto& r : rows) a.push_back(r[d]);
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
  }
  size_t n = 1;
  for (const auto& a : axes) n *= a.size();
  if (n != rows.size())
    throw std::runtime_error("GridTable::fromRows: rows do not form a complete grid");

  std::vector<double> values(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<size_t> strides(nd, 1);
  for (int d = nd - 2; d >= 0; --d) strides[d] = strides[d + 1] * axes[d + 1].size();
  for (const auto& r : rows) {
    if (static_cast<int>(r.size()) != nd + 1) throw std::runtime_error("GridTable::fromRows: ragged rows");
    size_t k = 0;
    for (int d = 0; d < nd; ++d) {
      auto it = std::lower_bound(axes[d].begin(), axes[d].end(), r[d]);
      k += static_cast<size_t>(it - axes[d].begin()) * strides[d];
    }
    if (!std::isnan(values[k])) throw std::runtime_error("GridTable::fromRows: duplicate node");
    values[k] = r[nd];
  }
  return GridTable(std::move(names), std::move(axes), std::move(values));
}

size_t GridTable::flatIndex(const int* idx) const {
  size_t k = 0;
  for (int d = 0; d < ndim(); ++d) k += static_cast<size_t>(idx[d]) * strides_[d];
  return k;
}

bool GridTable::locate(int d, double x, int& i, double& t) const {
  const auto& a = axes_[d];
  if (!(x >= a.front() && x <= a.back())) return false; // also rejects NaN
  i = static_cast<int>(std::upper_bound(a.begin(), a.end(), x) - a.begin()) - 1;
  i = std::min(i, static_cast<int>(a.size()) - 2);
  t = (x - a[i]) / (a[i + 1] - a[i]);
  return true;
}

double GridTable::eval(const double* x) const {
  const int nd = ndim();
  int i[3];
  double t[3];
  for (int d = 0; d < nd; ++d)
    if (!locate(d, x[d], i[d], t[d])) return std::numeric_limits<double>::quiet_NaN();

  const size_t base = flatIndex(i);
  double sum = 0;
  for (int c = 0; c < (1 << nd); ++c) {
    double w = 1;
    size_t k = base;
    for (int d = 0; d < nd; ++d) {
      if (c >> d & 1) {
        w *= t[d];
        k += strides_[d];
      } else {
        w *= 1 - t[d];
      }
    }
    if (w != 0) sum += w * values_[k];
  }
  return sum;
}

} // namespace nufit
