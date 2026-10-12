// Exact i.i.d. sampling of the NuFIT surrogate posterior.
//   sample_direct --variant TBoff --n 1000000 --out results/direct_TBoff
// writes <out>_samples.npy, <out>_meta.json and <out>.root (TTree "samples").

#include "Args.h"
#include "DirectSampler.h"
#include "NuFitModel.h"
#include "Output.h"

#include <chrono>
#include <iostream>

using namespace nufit;

int main(int argc, char** argv) try {
  auto defaults = Args::modelDefaults();
  defaults.insert({{"n", "1000000"}, {"seed", "1"}, {"out", "results/direct"}, {"root", "1"}});
  const Args args(argc, argv, defaults, "Exact sampler for the NuFIT 6.0 surrogate posterior");
  const ModelConfig cfg = args.model();
  const auto t0 = std::chrono::steady_clock::now();
  const NuFitModel model(cfg);
  DirectSampler sampler(model);
  const double tBuild = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  const long long n = args.integer("n");
  std::seed_seq seq{static_cast<uint64_t>(args.integer("seed")), uint64_t{0xd1ec7}};
  std::mt19937_64 rng(seq);
  const auto cols = sampleColumns();
  std::vector<double> data;
  data.reserve(n * cols.size());
  double p[kNPar];
  long long nIO = 0;
  for (long long i = 0; i < n; ++i) {
    sampler.sample(rng, p);
    data.insert(data.end(), p, p + kNPar);
    data.push_back(-0.5 * model.deltaChi2(p));
    data.push_back(model.logPrior(p));
    nIO += p[kDM3L] < 0;
  }
  const double tTot = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  const std::string out = args.str("out");
  writeNpy(out + "_samples.npy", data, cols.size());
  if (args.integer("root")) writeRootTree(out + ".root", "samples", cols, data);
  Json opts;
  for (const auto& [k, v] : args.all()) opts.str(k, v);
  Json()
      .str("sampler", "direct")
      .str("model", cfg.describe())
      .obj("options", opts)
      .strs("columns", cols)
      .num("n", n)
      .num("acceptance", double(sampler.nAccepted()) / sampler.nProposed())
      .num("fracIO", double(nIO) / n)
      .num("chi2MinNO", model.chi2Min(kNO))
      .num("chi2MinIO", model.chi2Min(kIO))
      .num("secondsBuild", tBuild)
      .num("secondsTotal", tTot)
      .write(out + "_meta.json");
  std::cout << cfg.describe() << "\n"
            << n << " samples, acceptance " << double(sampler.nAccepted()) / sampler.nProposed()
            << ", P(IO) = " << double(nIO) / n << ", " << tTot << " s -> " << out << "_*\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << "error: " << e.what() << "\n";
  return 1;
}
