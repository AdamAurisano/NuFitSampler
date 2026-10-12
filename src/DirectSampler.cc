#include "DirectSampler.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nufit {

namespace {
// Cells whose best corner is this far above the minimum carry no mass (exp(-700)).
constexpr double kChi2Cut = 1400;
} // namespace

GridBlock::GridBlock(const GridTable* table, std::vector<bool> sampled, double offset,
                     const GridTable* offsetTable,
                     std::vector<std::function<double(double)>> logPrior)
    : table_(table), sampled_(std::move(sampled)), offset_(offset), offsetTable_(offsetTable),
      logPrior_(std::move(logPrior)) {
  const int nd = table_->ndim();
  if (static_cast<int>(sampled_.size()) != nd || static_cast<int>(logPrior_.size()) != nd)
    throw std::invalid_argument("GridBlock: per-axis arguments do not match table dimension");
  for (int d = 0; d < nd; ++d) {
    if (sampled_[d]) {
      sAxes_.push_back(d);
      cellShape_.push_back(static_cast<int>(table_->axis(d).size()) - 1);
    } else {
      if (condAxis_ >= 0) throw std::invalid_argument("GridBlock: at most one conditioning axis");
      condAxis_ = d;
    }
  }
  if (offsetTable_) {
    if (condAxis_ < 0 || offsetTable_->ndim() != 1 ||
        offsetTable_->axis(0) != table_->axis(condAxis_))
      throw std::invalid_argument("GridBlock: offset table must live on the conditioning axis");
  }
  const int ns = static_cast<int>(sAxes_.size());
  const int nCond = condAxis_ >= 0 ? static_cast<int>(table_->axis(condAxis_).size()) : 1;

  size_t nCells = 1;
  for (int n : cellShape_) nCells *= n;
  cumMass_.resize(nCells);
  logEnv_.resize(nCells);

  std::vector<int> ci(ns, 0);
  int idx[3] = {0, 0, 0};
  double cum = 0;
  for (size_t c = 0; c < nCells; ++c) {
    // decode cell index (last sampled axis fastest)
    size_t r = c;
    for (int k = ns - 1; k >= 0; --k) {
      ci[k] = static_cast<int>(r % cellShape_[k]);
      r /= cellShape_[k];
    }
    double minChi = std::numeric_limits<double>::infinity();
    for (int j = 0; j < nCond; ++j) {
      const double off = offset_ + (offsetTable_ ? offsetTable_->values()[j] : 0.0);
      if (condAxis_ >= 0) idx[condAxis_] = j;
      for (int corner = 0; corner < (1 << ns); ++corner) {
        for (int k = 0; k < ns; ++k) idx[sAxes_[k]] = ci[k] + (corner >> k & 1);
        minChi = std::min(minChi, table_->node(idx) - off);
      }
    }
    double logPriorMax = 0, vol = 1;
    for (int k = 0; k < ns; ++k) {
      const auto& a = table_->axis(sAxes_[k]);
      const double x0 = a[ci[k]], x1 = a[ci[k] + 1];
      logPriorMax += std::max(logPrior_[sAxes_[k]](x0), logPrior_[sAxes_[k]](x1));
      vol *= x1 - x0;
    }
    double mass = 0;
    logEnv_[c] = std::numeric_limits<double>::infinity();
    if (minChi < kChi2Cut) {
      logEnv_[c] = -0.5 * minChi + logPriorMax;
      if (!std::isfinite(logEnv_[c]))
        throw std::runtime_error("GridBlock: unbounded prior in a cell with non-negligible likelihood");
      mass = std::exp(logEnv_[c]) * vol;
    }
    cum += mass;
    cumMass_[c] = cum;
  }
  if (!(cum > 0)) throw std::runtime_error("GridBlock: envelope has zero mass");
}

double GridBlock::logDensity(const double* x) const {
  double off = offset_;
  if (offsetTable_) off += offsetTable_->eval1(x[condAxis_]);
  double lf = -0.5 * (table_->eval(x) - off);
  for (int d : sAxes_) lf += logPrior_[d](x[d]);
  return lf;
}

double GridBlock::propose(std::mt19937_64& rng, double* x) const {
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  const double target = u01(rng) * cumMass_.back();
  size_t c = static_cast<size_t>(std::upper_bound(cumMass_.begin(), cumMass_.end(), target) -
                                 cumMass_.begin());
  if (c >= cumMass_.size()) c = cumMass_.size() - 1;
  size_t r = c;
  for (int k = static_cast<int>(sAxes_.size()) - 1; k >= 0; --k) {
    const int i = static_cast<int>(r % cellShape_[k]);
    r /= cellShape_[k];
    const auto& a = table_->axis(sAxes_[k]);
    x[sAxes_[k]] = a[i] + u01(rng) * (a[i + 1] - a[i]);
  }
  const double lf = logDensity(x);
  return std::isfinite(lf) ? std::exp(lf - logEnv_[c]) : 0.0;
}

DirectSampler::DirectSampler(const NuFitModel& model) : model_(model) {
  const auto angle = [&model](double s) { return model.logPriorAngle(s); };
  const auto flat = [](double) { return 0.0; };
  const auto l21 = [&model](double l) { return model.logPriorLog10Dm21(l); };
  double mass[2];
  for (int oi = 0; oi < 2; ++oi) {
    const Ordering o = static_cast<Ordering>(oi);
    const double m = model.chi2Min(o);
    b3_.emplace_back(&model.t3(o), std::vector<bool>{true, true, true}, 0.0, nullptr,
                     std::vector<std::function<double(double)>>{angle, flat, flat});
    b2_.emplace_back(&model.t2(o), std::vector<bool>{true, true}, m, nullptr,
                     std::vector<std::function<double(double)>>{angle, l21});
    if (model.config().coupledT13)
      b1_.emplace_back(&model.t1323(o), std::vector<bool>{true, false}, 0.0, &model.t23(o),
                       std::vector<std::function<double(double)>>{angle, flat});
    else
      b1_.emplace_back(&model.t1(o), std::vector<bool>{true}, m, nullptr,
                       std::vector<std::function<double(double)>>{angle});
    mass[oi] = b3_[oi].totalMass() * b2_[oi].totalMass() * b1_[oi].totalMass();
  }
  orderingCut_ = mass[kNO] / (mass[kNO] + mass[kIO]);
}

void DirectSampler::sample(std::mt19937_64& rng, double* p) {
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  const bool coupled = model_.config().coupledT13;
  for (;;) {
    ++nProposed_;
    const int o = u01(rng) < orderingCut_ ? kNO : kIO;
    double x3[3], x2[2], x1[2];
    double a = b3_[o].propose(rng, x3);
    a *= b2_[o].propose(rng, x2);
    if (coupled) x1[1] = x3[0]; // condition s13 on the sampled s23
    a *= b1_[o].propose(rng, x1);
    if (u01(rng) >= a) continue;
    ++nAccepted_;
    p[kS23] = x3[0];
    p[kDM3L] = x3[1];
    p[kDCP] = x3[2];
    p[kS12] = x2[0];
    p[kDM21] = std::pow(10.0, x2[1]) / 1e-5;
    p[kS13] = x1[0];
    return;
  }
}

} // namespace nufit
