// Unit tests that need no NuFIT data: interpolation, exact rejection sampler, PT on a toy.

#include "DirectSampler.h"
#include "GridTable.h"
#include "PTSampler.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>

using namespace nufit;

static int failures = 0;
#define CHECK(cond)                                                                          \
  do {                                                                                       \
    if (!(cond)) {                                                                           \
      std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond << "\n";          \
      ++failures;                                                                            \
    }                                                                                        \
  } while (0)
#define CHECK_NEAR(a, b, tol)                                                                \
  do {                                                                                       \
    const double va = (a), vb = (b);                                                         \
    if (!(std::abs(va - vb) <= (tol))) {                                                     \
      std::cerr << __FILE__ << ":" << __LINE__ << ": " #a " = " << va << " vs " #b " = " << vb \
                << " (tol " << (tol) << ")\n";                                               \
      ++failures;                                                                            \
    }                                                                                        \
  } while (0)

static std::vector<double> irregular(double lo, double hi, int n, double power) {
  std::vector<double> a(n);
  for (int i = 0; i < n; ++i) a[i] = lo + (hi - lo) * std::pow(double(i) / (n - 1), power);
  return a;
}

static void testInterpolation() {
  // N-linear interpolation reproduces any function that is affine in each coordinate.
  auto f = [](double x, double y, double z) { return 1 + 2 * x - 3 * y + 0.5 * z + 0.7 * x * y * z; };
  const auto ax = irregular(0, 1, 7, 1.5), ay = irregular(-2, 3, 11, 0.7), az = irregular(-180, 180, 9, 1);
  std::vector<std::vector<double>> rows;
  for (double x : ax)
    for (double y : ay)
      for (double z : az) rows.push_back({x, y, z, f(x, y, z)});
  std::shuffle(rows.begin(), rows.end(), std::mt19937_64(3));
  const GridTable t = GridTable::fromRows({"x", "y", "z"}, rows);
  CHECK(t.ndim() == 3 && t.axis(1).size() == ay.size());
  std::mt19937_64 rng(1);
  std::uniform_real_distribution<double> u(0, 1);
  for (int i = 0; i < 1000; ++i) {
    const double p[3] = {u(rng), -2 + 5 * u(rng), -180 + 360 * u(rng)};
    CHECK_NEAR(t.eval(p), f(p[0], p[1], p[2]), 1e-9);
  }
  for (double x : ax) // exact at nodes, including the upper edges
    for (double z : az) {
      const double p[3] = {x, ay.back(), z};
      CHECK_NEAR(t.eval(p), f(x, ay.back(), z), 1e-12);
    }
  const double out1[3] = {1.0001, 0, 0}, out2[3] = {0.5, NAN, 0};
  CHECK(std::isnan(t.eval(out1)));
  CHECK(std::isnan(t.eval(out2)));

  bool threw = false;
  rows.pop_back();
  try {
    GridTable::fromRows({"x", "y", "z"}, rows);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);
}

// Integral of exp(-chi2/2) * exp(lp(x)) over [a,b] with chi2 interpolated linearly.
template <class F>
static double quad(const GridTable& t, F lp, double a, double b, int n, double power) {
  double s = 0;
  const double h = (b - a) / n;
  for (int i = 0; i < n; ++i) {
    const double x = a + (i + 0.5) * h;
    s += std::pow(x, power) * std::exp(-0.5 * t.eval1(x) + lp(x)) * h;
  }
  return s;
}

static void testDirectBlock() {
  // 1D: chi2 = ((x-0.3)/0.1)^2 on an irregular coarse grid, with a monotone prior factor.
  const auto ax = irregular(0, 1, 15, 1.3);
  std::vector<double> v;
  for (double x : ax) v.push_back(std::pow((x - 0.3) / 0.1, 2));
  const GridTable t({"x"}, {ax}, v);
  auto lp = [](double x) { return 2.0 * x; };
  GridBlock b(&t, {true}, 0.0, nullptr, {lp});
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> u01(0, 1);
  const int N = 400000;
  double m1 = 0, m2 = 0;
  int acc = 0;
  for (int i = 0; i < 4 * N && acc < N; ++i) {
    double x;
    if (u01(rng) < b.propose(rng, &x)) {
      m1 += x;
      m2 += x * x;
      ++acc;
    }
  }
  m1 /= acc;
  m2 /= acc;
  const double z = quad(t, lp, 0, 1, 200000, 0);
  const double e1 = quad(t, lp, 0, 1, 200000, 1) / z, e2 = quad(t, lp, 0, 1, 200000, 2) / z;
  const double sd = std::sqrt(e2 - e1 * e1);
  CHECK_NEAR(m1, e1, 5 * sd / std::sqrt(acc));
  CHECK_NEAR(m2 - m1 * m1, e2 - e1 * e1, 5 * sd * sd * std::sqrt(2.0 / acc));
  // envelope mass bounds the integral
  CHECK(b.totalMass() >= z);

  // 2D conditional block: chi2(x, c) = ((x - 0.2 - 0.4 c)/0.1)^2 sampled given c, with
  // offset table min_x chi2 = 0. At fixed c the conditional must be the 1D law.
  const auto cx = irregular(0, 1, 41, 1), cc = irregular(0, 1, 5, 1);
  std::vector<double> v2, off(cc.size(), 0.0);
  for (double x : cx)
    for (double c : cc) v2.push_back(std::pow((x - 0.2 - 0.4 * c) / 0.1, 2));
  const GridTable t2({"x", "c"}, {cx, cc}, v2);
  const GridTable to({"c"}, {cc}, off);
  auto flat = [](double) { return 0.0; };
  GridBlock bc(&t2, {true, false}, 0.0, &to, {flat, flat});
  const double c = 0.5; // node: conditional mean 0.4
  double s = 0;
  acc = 0;
  for (int i = 0; i < 4 * N && acc < N / 4; ++i) {
    double x[2] = {0, c};
    if (u01(rng) < bc.propose(rng, x)) {
      s += x[0];
      ++acc;
    }
  }
  CHECK_NEAR(s / acc, 0.4, 5 * 0.1 / std::sqrt(acc) + 1e-3);
}

// Two well-separated 2D Gaussians with weights 0.3 / 0.7: a random walk cannot cross.
class ToyMixture : public Target {
public:
  int dim() const override { return 2; }
  std::vector<std::string> names() const override { return {"x", "y"}; }
  void evaluate(const double* x, double& logL, double& logPrior) const override {
    if (std::abs(x[0]) > 10 || std::abs(x[1]) > 10) {
      logPrior = -INFINITY;
      logL = 0;
      return;
    }
    logPrior = 0;
    auto g = [&](double mx) { return std::exp(-0.5 * ((x[0] - mx) * (x[0] - mx) + x[1] * x[1]) / (0.3 * 0.3)); };
    logL = std::log(0.3 * g(-4) + 0.7 * g(4) + 1e-300);
  }
  void initialScale(double* s) const override { s[0] = s[1] = 0.3; }
  void initialPoint(std::mt19937_64& rng, double* x) const override {
    std::uniform_real_distribution<double> u(-5, 5);
    x[0] = u(rng);
    x[1] = u(rng);
  }
};

static void testPT() {
  ToyMixture toy;
  PTConfig cfg;
  cfg.nReplicas = 4;
  cfg.nRungs = 12;
  cfg.tMax = 300;
  cfg.nBurn = 20000;
  cfg.nSteps = 200000;
  cfg.thin = 10;
  cfg.useJumps = false;
  const auto res = PTSampler(toy, cfg).run();
  double fracAll = 0;
  std::vector<double> fr;
  for (const auto& r : res) {
    const size_t rows = r.cold.size() / 4;
    double nRight = 0;
    for (size_t i = 0; i < rows; ++i) nRight += r.cold[i * 4] > 0;
    fr.push_back(nRight / rows);
    fracAll += nRight / rows / res.size();
    CHECK(r.roundTrips > 10);
  }
  double var = 0;
  for (double f : fr) var += (f - fracAll) * (f - fracAll) / (fr.size() - 1);
  const double mcse = std::sqrt(var / fr.size());
  std::printf("PT toy: P(right mode) = %.4f +- %.4f (truth 0.7), round trips %llu\n", fracAll, mcse,
              static_cast<unsigned long long>(res[0].roundTrips));
  CHECK_NEAR(fracAll, 0.7, std::max(4 * mcse, 0.01));
}

int main() {
  testInterpolation();
  testDirectBlock();
  testPT();
  if (failures) {
    std::cerr << failures << " failure(s)\n";
    return 1;
  }
  std::cout << "all core tests passed\n";
  return 0;
}
