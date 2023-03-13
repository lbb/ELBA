/* Created by Giulia Guidi on 4/14/2021 */

#pragma once

#include "../../ipuma-lib/src/driver.hpp"
#include "../../ipuma-lib/src/swatlib/vector.hpp"
#include "../AlignmentInfo.hpp"
#include "PairwiseFunction.hpp"
#include <string_view>

typedef std::tuple<int, int> SeedPair;

// template <typename TSequenceValue, typename TSpec>
class IPUAligner : public PairwiseFunction {
 public:
  IPUAligner(ScoringScheme scoring_scheme,
             ushort seed_length, int xdrop, int seed_count);

  void
  PostAlignDecision(const LoganAlignmentInfo &ai, bool &passed, float &ratioScoreOverlap,
                    int &dir, int &dirT, int &sfx, int &sfxT, uint32_t &overlap, const bool noAlign, std::vector<int64_t> &ContainedSeqMyThread);

  void
  apply(uint64_t l_col_idx, uint64_t g_col_idx,
        uint64_t l_row_idx, uint64_t g_row_idx,
        seqan::Dna5String *seq_h, seqan::Dna5String *seq_v,
        ushort k,
        elba::CommonKmers &cks, std::stringstream &ss) override;

  void
  apply_batch(seqan::StringSet<seqan::Dna5String> &seqsh,
              seqan::StringSet<seqan::Dna5String> &seqsv,
              uint64_t *lids,
              uint64_t col_offset,
              uint64_t row_offset,
              PSpMat<elba::CommonKmers>::ref_tuples *mattuples,
              std::ofstream &lfs,
              const bool noAlign,
              ushort k,
              uint64_t nreads,
              std::vector<int64_t> &ContainedSeqPerThread,
              float ratioScoreOverlap = 0.99,  // GGGG: Precomputed for error rate = 15% and default scoring matrix (1,-1,-1) (0.445 for CLR, 0.99 for CCS)
              int debugThr = 50) override;     // GGGG: Fixed threshold, this is convenient only for debugging
  void runIPUAlign(
      const std::vector<std::string> &seqHs, 
      const std::vector<std::string> &seqVs, 
      std::vector<std::tuple<SeedPair, SeedPair>> seeds, 
      std::vector<int> &xscores, 
      int xdrop, 
      int seed_length);

 private:
  ipu::batchaffine::SWAlgorithm *driver_algo = nullptr;
  ScoringScheme scoring_scheme;
  ushort seed_length;
  int xdrop;
  int seed_count;
};

