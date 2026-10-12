// Parallel-tempering MCMC of the NuFIT surrogate posterior.
//   sample_pt --variant TBoff --replicas 4 --rungs 8 --burn 50000 --steps 500000 --out results/pt_TBoff
// writes <out>_cold.npy (replica, iter, params, logL, logPrior), <out>_rungs.npy,
// <out>_labels.npy, <out>_ladder.npy, <out>_meta.json and <out>.root (TTree "samples").

#include "Args.h"
#include "NuFitModel.h"
#include "NuFitTarget.h"
#include "Output.h"
#include "PTSampler.h"

#include <iostream>

using namespace nufit;

int main(int argc, char** argv) try {
  auto defaults = Args::modelDefaults();
  defaults.insert({{"replicas", "4"},
                   {"rungs", "8"},
                   {"tmax", "50"},
                   {"burn", "50000"},
                   {"steps", "500000"},
                   {"thin", "10"},
                   {"pjump", "0.1"},
                   {"jumps", "1"},
                   {"adapt-ladder", "1"},
                   {"seed", "12345"},
                   {"out", "results/pt"},
                   {"root", "1"}});
  const Args args(argc, argv, defaults, "Parallel-tempering sampler for the NuFIT 6.0 surrogate posterior");
  const ModelConfig mcfg = args.model();
  const NuFitModel model(mcfg);
  const NuFitTarget target(model);

  PTConfig cfg;
  cfg.nReplicas = static_cast<int>(args.integer("replicas"));
  cfg.nRungs = static_cast<int>(args.integer("rungs"));
  cfg.tMax = args.num("tmax");
  cfg.nBurn = args.integer("burn");
  cfg.nSteps = args.integer("steps");
  cfg.thin = static_cast<int>(args.integer("thin"));
  cfg.pJump = args.num("pjump");
  cfg.useJumps = args.integer("jumps") != 0;
  cfg.adaptLadder = args.integer("adapt-ladder") != 0;
  cfg.seed = args.integer("seed");

  const auto res = PTSampler(target, cfg).run();

  const int K = cfg.nRungs;
  std::vector<std::string> coldCols = {"replica", "iter"};
  for (const auto& c : sampleColumns()) coldCols.push_back(c);
  const size_t rowLen = sampleColumns().size();
  std::vector<double> cold, rungs, labels, ladder, coldRoot;
  for (int r = 0; r < cfg.nReplicas; ++r) {
    const auto& R = res[r];
    const size_t rows = R.cold.size() / rowLen;
    for (size_t i = 0; i < rows; ++i) {
      const double it = cfg.nBurn + double(i) * cfg.thin;
      cold.push_back(r);
      cold.push_back(it);
      cold.insert(cold.end(), R.cold.begin() + i * rowLen, R.cold.begin() + (i + 1) * rowLen);
      rungs.push_back(r);
      rungs.push_back(it);
      rungs.insert(rungs.end(), R.rungLogL.begin() + i * K, R.rungLogL.begin() + (i + 1) * K);
      labels.push_back(r);
      labels.push_back(it);
      labels.insert(labels.end(), R.labelRung.begin() + i * K, R.labelRung.begin() + (i + 1) * K);
    }
    for (size_t i = 0; i < R.ladderHistory.size() / (K + 1); ++i) {
      ladder.push_back(r);
      ladder.insert(ladder.end(), R.ladderHistory.begin() + i * (K + 1), R.ladderHistory.begin() + (i + 1) * (K + 1));
    }
  }
  const std::string out = args.str("out");
  writeNpy(out + "_cold.npy", cold, coldCols.size());
  writeNpy(out + "_rungs.npy", rungs, K + 2);
  writeNpy(out + "_labels.npy", labels, K + 2);
  if (!ladder.empty()) writeNpy(out + "_ladder.npy", ladder, K + 2);
  if (args.integer("root")) writeRootTree(out + ".root", "samples", coldCols, cold);

  Json opts;
  for (const auto& [k, v] : args.all()) opts.str(k, v);
  std::vector<std::string> jumpNames;
  for (int j = 0; j < target.nJumps(); ++j) jumpNames.push_back(target.jumpName(j));
  std::vector<Json> reps;
  for (const auto& R : res) {
    Json j;
    j.arr("betas", R.betas).arr("rwAccept", R.rwAccept).arr("rwScale", R.rwScale).arr("swapAccept", R.swapAccept);
    std::vector<Json> ja;
    for (const auto& v : R.jumpAccept) ja.push_back(Json().arr("accept", v));
    j.objs("jumpAccept", ja).num("roundTrips", R.roundTrips).num("seconds", R.seconds);
    reps.push_back(j);
  }
  Json()
      .str("sampler", "pt")
      .str("model", mcfg.describe())
      .obj("options", opts)
      .strs("columns", coldCols)
      .strs("jumps", jumpNames)
      .num("nRungs", K)
      .num("chi2MinNO", model.chi2Min(kNO))
      .num("chi2MinIO", model.chi2Min(kIO))
      .objs("replicas", reps)
      .write(out + "_meta.json");

  std::cout << mcfg.describe() << "\n";
  for (int r = 0; r < cfg.nReplicas; ++r) {
    std::cout << "replica " << r << ": " << res[r].seconds << " s, round trips " << res[r].roundTrips
              << ", cold RW acc " << res[r].rwAccept[0] << ", swap acc";
    for (double a : res[r].swapAccept) std::cout << " " << a;
    std::cout << "\n";
  }
  std::cout << "-> " << out << "_*\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << "error: " << e.what() << "\n";
  return 1;
}
