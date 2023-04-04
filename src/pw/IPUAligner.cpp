/* Created by Saliya Ekanayake on 2019-07-05 and modified by Giulia Guidi on 4/14/2021. */

#include "../../include/pw/IPUAligner.hpp"

#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Log.h>

#include <nlohmann/json.hpp>
#include <string_view>

using json = nlohmann::json;

uint _minOverlapLenL = 5000;

// for (int32_t i = 0; i < npairs; ++i) {
//   sequences[2 * i] = seqVs[i];
//   sequences[2 * i + 1] = seqHs[i];
//   int offsetH1, offsetV1;
//   int offsetH2, offsetV2;
//   SeedPair s1, s2;
//   std::tie(s1, s2) = seeds[i];

//   std::tie(offsetH1, offsetV1) = s1;
//   std::tie(offsetH2, offsetV2) = s2;
//   // if (!(seqVs[i].length() >= offsetV)) {
//   //   std::cerr << "REMAP (V) len:" << seqVs[i].length() << ", offset" << offsetV << std::endl;
//   // }
//   // if (!(seqHs[i].length() >= offsetH)) {
//   //   std::cerr << "REMAP (H) len:" << seqHs[i].length() << ", offset" << offsetH << std::endl;
//   // }
//   assert(seqVs[i].length() >= offsetV && "VSeq offset");
//   assert(seqHs[i].length() >= offsetH && "HSeq offset");
//   comparisons[i] = {
//       i,
//       2 * i,
//       sequences[2 * i].size(),
//       2 * i + 1,
//       sequences[2 * i + 1].size(),
//       {{{offsetV1, offsetH1}, {offsetV2, offsetH2}}}};
// }
// {
//   ofstream myfile;
//   myfile.open ("seqs_seed_H.txt");
//   for (auto &&i : seqHs) {
//     myfile << i << '\n';
//   }
//   myfile.close();
// }
// {
//   ofstream myfile;
//   myfile.open ("seqs_seed_V.txt");
//   for (auto &&i : seqVs) {
//     myfile << i << '\n';
//   }
//   myfile.close();
// }
// {
//   ofstream myfile;
//   myfile.open ("seeds_H1.txt");
//   for (int32_t i = 0; i < npairs; ++i) {
//     int offsetH1, offsetV1;
//     int offsetH2, offsetV2;
//     SeedPair s1, s2;
//     std::tie(s1, s2) = seeds[i];

//     std::tie(offsetH1, offsetV1) = s1;
//     std::tie(offsetH2, offsetV2) = s2;
//     myfile << offsetH1 << '\n';
//   }
//   myfile.close();
// }
// {
//   ofstream myfile;
//   myfile.open ("seeds_V1.txt");
//   for (int32_t i = 0; i < npairs; ++i) {
//     int offsetH1, offsetV1;
//     int offsetH2, offsetV2;
//     SeedPair s1, s2;
//     std::tie(s1, s2) = seeds[i];

//     std::tie(offsetH1, offsetV1) = s1;
//     std::tie(offsetH2, offsetV2) = s2;
//     myfile << offsetV1 << '\n';
//   }
//   myfile.close();
// }
// {
//   ofstream myfile;
//   myfile.open ("seeds_H2.txt");
//   for (int32_t i = 0; i < npairs; ++i) {
//     int offsetH1, offsetV1;
//     int offsetH2, offsetV2;
//     SeedPair s1, s2;
//     std::tie(s1, s2) = seeds[i];

//     std::tie(offsetH1, offsetV1) = s1;
//     std::tie(offsetH2, offsetV2) = s2;
//     myfile << offsetH2 << '\n';
//   }
//   myfile.close();
// }
// {
//   ofstream myfile;
//   myfile.open ("seeds_V2.txt");
//   for (int32_t i = 0; i < npairs; ++i) {
//     int offsetH1, offsetV1;
//     int offsetH2, offsetV2;
//     SeedPair s1, s2;
//     std::tie(s1, s2) = seeds[i];

//     std::tie(offsetH1, offsetV1) = s1;
//     std::tie(offsetH2, offsetV2) = s2;
//     myfile << offsetV2 << '\n';
//   }
//   myfile.close();
// }

char _complementbase(char n) {
  switch (n) {
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
  }
}

std::string
_reversecomplement(const std::string &seq) {
  std::string cpyseq = seq;
  std::reverse(cpyseq.begin(), cpyseq.end());

  std::transform(
      std::begin(cpyseq),
      std::end(cpyseq),
      std::begin(cpyseq),
      _complementbase);

  return cpyseq;
}

int ipu_total_cmps = 0;

void IPUAligner::PostAlignDecision(const LoganAlignmentInfo &ai, bool &passed, float &ratioScoreOverlap,
                                   int &dir, int &dirT, int &sfx, int &sfxT, uint32_t &overlap, const bool noAlign, std::vector<int64_t> &ContainedSeqMyThread) {
  // {begin/end}Position{V/H}: Returns the begin/end position of the seed in the seqVs (vertical/horizonral direction)
  // these four return seqan:Tposition objects
  int begpV = ai.begSeedV;
  int endpV = ai.endSeedV;
  int begpH = ai.begSeedH;
  int endpH = ai.endSeedH;

  unsigned short int overlapLenH = ai.seq_h_seed_length;
  unsigned short int overlapLenV = ai.seq_v_seed_length;

  unsigned short int rlenH = ai.seq_h_length;
  unsigned short int rlenV = ai.seq_v_length;

  unsigned short int minLeft = min(begpV, begpH);
  unsigned short int minRight = min(rlenV - endpV, rlenH - endpH);

  int64_t seqV = ai.seq_v_g_idx;
  int64_t seqH = ai.seq_h_g_idx;

  overlap = minLeft + minRight + (overlapLenV + overlapLenH) / 2;

  if ((seqV == 5 && seqH == 99) ||
      (seqV == 5 && seqH == 100) ||
      (seqV == 5 && seqH == 184) ||
      (seqV == 24 && seqH == 40)) {
    std::cout << seqV + 1 << "\t" << seqH + 1 << "\t" << ai.rc << "\t" << begpV
              << "\t" << endpV << "\t" << begpH << "\t" << endpH << "\t" << ai.xscore << "\t" << overlap << std::endl;
  }

#ifndef FIXEDTHR
  float myThr = (1 - DELTACHERNOFF) * (ratioScoreOverlap * (float)overlap);

  // Contained overlaps removed for now, reintroduce them later
  // @GGGG-TODO: identify chimeric sequences
  bool contained = false;
  bool chimeric = false;

  if (begpV <= begpH && (rlenV - endpV) <= (rlenH - endpH)) {
    ContainedSeqMyThread.push_back(seqV);
    contained = true;
  } else if (begpV >= begpH && (rlenV - endpV) >= (rlenH - endpH)) {
    ContainedSeqMyThread.push_back(seqH);
    contained = true;
  } else if (!noAlign) {
    passed = ((float)ai.xscore >= myThr && overlap >= _minOverlapLenL);

    if (passed) {
      if (begpV > begpH) {
        dir = ai.rc ? 0 : 1;
        dirT = ai.rc ? 0 : 2;
        sfx = ((rlenH - endpH) - (rlenV - endpV));
        sfxT = begpV - begpH;
      } else {
        dir = ai.rc ? 3 : 2;
        dirT = ai.rc ? 3 : 1;
        sfx = begpH - begpV;
        sfxT = ((rlenV - endpV) - (rlenH - endpH));
      }
    }
  }

#else
  if (ai.xscore >= FIXEDTHR)
    passed = true;
#endif
}

  std::string getEnvVar( std::string const & key ) {
    char * val = getenv( key.c_str() );
    return val == NULL ? std::string("") : std::string(val);
}

IPUAligner::IPUAligner(
    ScoringScheme scoring_scheme,
    ushort seed_length, int xdrop, int seed_count) : PairwiseFunction(),
                                                     scoring_scheme(scoring_scheme),
                                                     seed_length(seed_length),
                                                     xdrop(xdrop),
                                                     seed_count(seed_count) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);
  std::cout << "Seed Length: " << seed_length << std::endl;
  SW_CONFIGURATION = {
      .gapInit = -1,
      .gapExtend = -1,
      .matchValue = 1,
      .mismatchValue = -1,
      .ambiguityValue = -1,
      .similarity = swatlib::Similarity::nucleicAcid,
      .datatype = swatlib::DataType::nucleicAcid,
      .seedLength = seed_length,
      .xDrop = xdrop,
  };

  ALGOCONFIG = {
      .numVertices = 1472,
      .maxSequenceLength = 19295,
      .maxComparisonsPerVertex = 200,
      .vertexBufferSize = 160000,
      .vtype = ipu::VertexType::xdroprestrictedseedextend,
      .fillAlgo = ipu::Algorithm::greedy,
      .complexityAlgo = ipu::Complexity::xdrop,
      .partitionadd = ipu::PartitionAdd::alternating,
      .partitioningSortComparisons = true,
      .forwardOnly = false,
      .ioTiles = 0,
      .bandPercentageXDrop = 0.45,
  };
  PLOGW << "ALGOCONFIG" << json{ALGOCONFIG}.dump();
  PLOGW << "SWCONFIG" << json{SW_CONFIGURATION}.dump();
  auto ipus = 1;
  auto env = getEnvVar("NIPUS");
  if (env != "") {
    PLOGI << "Evironment NIPUS = " << env;
    ipus = stoi(env);
  }

  this->driver_algo = new ipu::batchaffine::SWAlgorithm(SW_CONFIGURATION, ALGOCONFIG, 0, ipus);
}

void IPUAligner::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Dna5String *seqH, seqan::Dna5String *seqV, ushort k,
    elba::CommonKmers &cks, std::stringstream &ss) {
  // ...
}

double total_time = 0;
void IPUAligner::runIPUAlign(
    ColPacks &columnpacks,
    const std::vector<int> &seqHs_idx,
    const std::vector<int> &seqVs_idx,
    const std::vector<std::string> &sequences,
    const std::vector<std::string> &seqHs,
    const std::vector<std::string> &seqVs,
    std::vector<std::tuple<SeedPair, SeedPair>> seeds, std::vector<int> &xscores, int xdrop, int seed_length) {
  auto maxbatch_size = 100;
  std::vector<int> hist(maxbatch_size, 0);
  std::vector<ipu::MultiComparison> comparisons;

  swatlib::TickTock multibatch_timer;
  multibatch_timer.tick();
  auto cnt_comparisons = 0;
  for (int cpi = 0; cpi < columnpacks.size(); cpi++) {
    auto colinfo = columnpacks[cpi];
    // Twin / plain
    for (int tp = 0; tp < 2; tp++) {
      auto seqstore = 0;

      std::vector<ipu::Comparison> tmpcmps;
      auto coldirection = colinfo[tp];
      for (int colt = 0; colt < coldirection.size(); colt++) {
        auto [id, seqH_idx, seqV_idx] = coldirection[colt];
        int offsetH1, offsetV1;
        int offsetH2, offsetV2;
        SeedPair s1, s2;
        std::tie(s1, s2) = seeds[id];

        std::tie(offsetH1, offsetV1) = s1;
        std::tie(offsetH2, offsetV2) = s2;

        if (sequences[seqH_idx].size() + seqstore + sequences[seqV_idx].size() >= ALGOCONFIG.vertexBufferSize || tmpcmps.size() >= maxbatch_size) {
          comparisons.emplace_back(tmpcmps, seed_length);
          hist[tmpcmps.size() - 1] += 1;
          tmpcmps.resize(0);
          seqstore = 0;
        }

        seqstore += sequences[seqV_idx].size();
        cnt_comparisons += 2;
        tmpcmps.push_back({
            id,
            seqV_idx,
            sequences[seqV_idx].size(),
            seqH_idx,
            sequences[seqH_idx].size(),
            {{{offsetV1, offsetH1}, {offsetV2, offsetH2}}},
            0,
        });
      }
      // if (tmpcmps.size() == 7) {
      //     // Reduce to 6
      //     comparisons.emplace_back(std::vector<ipu::Comparison>{tmpcmps[tmpcmps.size()  - 1]}, seed_length);
      //     hist[1 - 1] += 1;
      //     tmpcmps.resize(tmpcmps.size() - 1);
      // }
      // if (tmpcmps.size() == 6) {
      //     // Reduce to 4
      //     comparisons.emplace_back(std::vector<ipu::Comparison>{tmpcmps[tmpcmps.size()  - 1], tmpcmps[tmpcmps.size()  - 2]}, seed_length);
      //     hist[2 - 1] += 1;
      //     tmpcmps.resize(tmpcmps.size() - 2);
      // }
      // if (tmpcmps.size() == 5) {
      //     // Reduce to 4
      //     comparisons.emplace_back(std::vector<ipu::Comparison>{tmpcmps[tmpcmps.size()  - 1]}, seed_length);
      //     hist[1 - 1] += 1;
      //     tmpcmps.resize(tmpcmps.size() - 1);
      // }
      // if (tmpcmps.size() == 4) {
      //     // Reduce to 2
      //     comparisons.emplace_back(std::vector<ipu::Comparison>{tmpcmps[tmpcmps.size()  - 1], tmpcmps[tmpcmps.size()  - 2]}, seed_length);
      //     hist[2 - 1] += 1;
      //     tmpcmps.resize(tmpcmps.size() - 2);
      // }
      // if (tmpcmps.size() == 3) {
      //     // Reduce to 2
      //     comparisons.emplace_back(std::vector<ipu::Comparison>{tmpcmps[tmpcmps.size()  - 1]}, seed_length);
      //     hist[1 - 1] += 1;
      //     tmpcmps.resize(tmpcmps.size() - 1);
      // }
      // if (tmpcmps.size() == 2) {
      //     // Reduce to 1
      //     comparisons.emplace_back(std::vector<ipu::Comparison>{tmpcmps[tmpcmps.size()  - 1], tmpcmps[tmpcmps.size()  - 2]}, seed_length);
      //     hist[2 - 1] += 1;
      //     tmpcmps.resize(tmpcmps.size() - 2);
      // }
      // if (tmpcmps.size() != 0) {
      //   auto s = tmpcmps.size();
      //   for (size_t i = 0; i < s; i++) {
      //     comparisons.emplace_back(std::vector<ipu::Comparison>{tmpcmps[tmpcmps.size()  - 1]}, seed_length);
      //     tmpcmps.resize(tmpcmps.size() - 1);
      //     hist[1 - 1] += 1;
      //   }
      // }
      if (tmpcmps.size() != 0) {
        comparisons.emplace_back(tmpcmps, seed_length);
        hist[tmpcmps.size() - 1] += 1;
        tmpcmps.resize(0);
      }
    }
  }
  multibatch_timer.tock();
  PLOGI << "MultiBatching took [ms]: " << multibatch_timer.duration();
  PLOGI << "Sigular comparison count: " << cnt_comparisons;

  PLOGD << "We got MultiComparison: " << comparisons.size();
  for (size_t i = 0; i < hist.size(); i++) {
    PLOGD.printf("Bucket %5d has %7d entries.", i + 1, hist[i]);
  }

  // std::cerr << "maxV: " << maxV << std::endl;
  // std::cerr << "maxSeedV: " << maxSeedV << std::endl;
  std::vector<std::string_view> view_seqs(sequences.size());
#pragma omp parallel for
  for (int i = 0; i < sequences.size(); i++) {
    view_seqs[i] = sequences[i];
  }

  // {
  //   ofstream myfile;
  //   myfile.open ("seqs.txt");
  //   json s = sequences;
  //   myfile << s << '\n';
  //   // for (auto i = 0; i < sequences.size(); ++i) {
  //   //   myfile << sequences[i] << '\n';
  //   // }
  //   myfile.close();
  // }
  // {
  //   ofstream myfile;
  //   myfile.open ("cmps.txt");
  //   json c = comparisons;
  //   myfile << c << '\n';
  //   // for (auto i = 0; i < comparisons.size(); ++i) {
  //   //   json j{};
  //   //   ipu::to_json(j, comparisons[i]);
  //   //   myfile << j << '\n';
  //   // }
  //   myfile.close();
  // }

  std::vector<ipu::partition::BatchMapping> mappings;
  auto &cmps = comparisons;
  auto &seqs = view_seqs;

  // std::vector<ipu::Batch> batches;
  swatlib::TickTock mapping_timer;
  mapping_timer.tick();
  PLOGI << "Use multicomparisons for creating batches";
  mappings = ipu::partition::mapBatches(ALGOCONFIG, seqs, cmps);
  // batches = ipu::create_batches(seqs, mcmps, config.ipuconfig, config.swconfig);
  mapping_timer.tock();
  PLOGI << "Mapping took [ms]: " << mapping_timer.duration();

  std::vector<int32_t> scores(xscores.size());
  std::vector<ipu::batchaffine::Job *> jobs(mappings.size());
  std::vector<ipu::Batch> batches(mappings.size());
  bool lazyGeneration = true;

  if (mappings.size() > 750) {
    lazyGeneration = true;
    PLOGE << "Generate batches lazy, as batches " << mappings.size();
  }

  if (!lazyGeneration) {
    swatlib::TickTock batch_timer;
    batch_timer.tick();
#pragma omp parallel for
    for (int i = 0; i < mappings.size(); ++i) {
      batches[i] = ipu::create_batch(mappings[i], seqs, ALGOCONFIG, SW_CONFIGURATION);
    }
    batch_timer.tock();
    PLOGI << "Batch copies took [ms]: " << batch_timer.duration();
  }

  auto runt = std::chrono::system_clock::now();
  int64_t gcells = 0;
  int progress = 0;
#pragma omp parallel for
  for (int i = 0; i < mappings.size(); ++i) {
    ipu::Batch batch;
    if (lazyGeneration) {
      batch = ipu::create_batch(mappings[i], seqs, ALGOCONFIG, SW_CONFIGURATION);
    } else {
      batch = std::move(batches[i]);
    }
    jobs[i] = driver_algo->async_submit(&batch);
    driver_algo->blocking_join(*jobs[i]);

#pragma omp critical
    {
      progress++;
      gcells += batch.cellCount / 1e9;
      PLOGI << "Received batch " << progress << " / " << jobs.size();
    }

    auto result = batch.get_result();
    delete jobs[i];
    for (int ii = 0; ii < batch.origin_comparison_index.size(); ++ii) {
      auto [orig_i, orig_seed] = ipu::unpackOriginIndex(batch.origin_comparison_index[ii]);
      if (orig_i >= 0) {
        int lpartScore = result.a_range_result[ii];
        int rpartScore = result.b_range_result[ii];
        // PLOGF << orig_i << ":" << batch.origin_comparison_index[i] << " " << lpartScore << " " << rpartScore;
        scores[orig_i] = std::max(lpartScore + rpartScore + SW_CONFIGURATION.seedLength, scores[orig_i]);
      }
    }
  }

  PLOGI << "Cell count: " << gcells;

  // std::vector<ipu::Batch> batches = ipu::create_batches(view_seqs, comparisons, ALGOCONFIG, SW_CONFIGURATION);
  //// std::vector<ipu::Batch> batches = ipu::create_batches(sequences_old, real_cmps, ALGOCONFIG, SW_CONFIGURATION);

  // auto runt = std::chrono::system_clock::now();
  // std::vector<ipu::batchaffine::Job *> jobs(batches.size());
  // for (size_t i = 0; i < batches.size(); i++) {
  //   jobs[i] = driver_algo->async_submit(&batches[i]);
  // }

  // std::vector<ipu::BlockAlignmentResults> results(batches.size());
  //// #pragma omp parallel for
  //// }

  // for (size_t i = 0; i < batches.size(); i++) {
  //   driver_algo->blocking_join(*jobs[i]);
  //   int j = i;
  //   auto &batch = batches[j];
  //   ipu::BlockAlignmentResults res = batch.get_result();
  //   for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
  //     auto [orig_i, seed_j] = ipu::unpackOriginIndex(batch.origin_comparison_index[i]);
  //     if (orig_i >= 0 && seed_j == 0) {
  //       xscores[orig_i] = std::max( (int32_t)res.a_range_result[i]+(int32_t)res.b_range_result[i], (int32_t)res.a_range_result[i+1]+(int32_t)res.b_range_result[i+1]);
  //       // PLOGW <<  std::max( (int32_t)res.a_range_result[i]+(int32_t)res.b_range_result[i], (int32_t)res.a_range_result[i+1]+(int32_t)res.b_range_result[i+1]);
  //     }
  //   }
  // }
  auto end_runt = std::chrono::system_clock::now();
  std::cout << "Local TIMEms:::::::::::::::::  " << (ms_t(end_runt - runt)).count() << std::endl;
  total_time += (ms_t(end_runt - runt)).count();
  std::cout << "TIMEms:::::::::::::::::  " << total_time << std::endl;
}

void IPUAligner::apply_batch(
    seqan::StringSet<seqan::Dna5String> &seqsh,
    seqan::StringSet<seqan::Dna5String> &seqsv,
    uint64_t *lids,
    uint64_t col_offset,
    uint64_t row_offset,
    PSpMat<elba::CommonKmers>::ref_tuples *mattuples,
    std::ofstream &lfs,
    const bool noAlign,
    ushort k,
    uint64_t nreads,
    std::vector<int64_t> &ContainedSeqPerBatch,
    float ratioScoreOverlap,  // GGGG: this is my ratioScoreOverlap variable change name later
    int debugThr) {
  throw std::runtime_error("This should never be called");
}

// @NOTE This is hard-coded to the number of seeds being <= 2
void IPUAligner::apply_batch(
    std::vector<int> &seqsh_idx,
    std::vector<int> &seqsv_idx,
    std::vector<seqan::Dna5String *> &seqs_db_h,
    std::vector<seqan::Dna5String *> &seqs_db_v,
    seqan::StringSet<seqan::Dna5String> &dep_seqsh,
    seqan::StringSet<seqan::Dna5String> &dep_seqsv,
    uint64_t *lids,
    uint64_t col_offset,
    uint64_t row_offset,
    PSpMat<elba::CommonKmers>::ref_tuples *mattuples,
    std::ofstream &lfs,
    const bool noAlign,
    ushort k,
    uint64_t nreads,
    std::vector<int64_t> &ContainedSeqPerBatch,
    float ratioScoreOverlap,  // GGGG: this is my ratioScoreOverlap variable change name later
    int debugThr) {
  PLOGD << "Apply batch";
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int numThreads = 1;
#ifdef THREADED
#pragma omp parallel
  {
    numThreads = omp_get_num_threads();
  }
#endif

  uint64_t npairs = seqsh_idx.size();

  lfs << "processing batch of size " << npairs << " with " << numThreads << " threads " << std::endl;

  // for multiple seeds we store the seed with the highest identity
  std::vector<IPumaAlignmentInfo> ai(npairs);

  // bool *strands = new bool[npairs];
  // int  *xscores = new int[npairs];
  // TSeed  *seeds = new TSeed[npairs];

  /* GGGG: seed_count is hardcoded here (2) */
  int count_errors = 0;

  // TODO: We dont need this anymore/
  std::vector<string> seqHs(npairs);
  std::vector<string> seqVs(npairs);

  auto totalUniqueSeqs = 0;
  for (int i = 0; i < npairs; i++) {
    totalUniqueSeqs = std::max(totalUniqueSeqs, seqsh_idx[i]);
    totalUniqueSeqs = std::max(totalUniqueSeqs, seqsv_idx[i]);
  }
  totalUniqueSeqs += 1;
  PLOGD << "totalUniqueSeqs " << totalUniqueSeqs;
  auto getHid = [&](int i) { return seqsh_idx[i]; };
  auto getVid = [&](int i) { return totalUniqueSeqs + seqsv_idx[i]; };
  auto getHTid = [&](int i) { return totalUniqueSeqs * 2 + seqsh_idx[i]; };

  swatlib::TickTock copy_time;
  copy_time.tick();
  std::vector<string> seqs(totalUniqueSeqs * 3);
// TODO: Here, we can copy over sequences....
#pragma omp parallel for
  for (auto i = 0; i < totalUniqueSeqs; i++) {
    std::string seqH;
    seqan::assign(seqH, seqan::Dna5String(*(seqs_db_h[i])));

    std::string seqV;
    seqan::assign(seqV, seqan::Dna5String(*(seqs_db_v[i])));

    std::string twinseqH;
    twinseqH.resize(seqH.size()); 
    // std::reverse(std::begin(twinseqH), std::end(twinseqH));
    // std::transform(std::begin(twinseqH), std::end(twinseqH), std::begin(twinseqH), _complementbase);
    std::transform(std::rbegin(seqH), std::rend(seqH), std::begin(twinseqH), _complementbase);


    seqs[i] = seqH;
    seqs[totalUniqueSeqs + i] = seqV;
    seqs[totalUniqueSeqs * 2 + i] = twinseqH;
  }
  copy_time.tock();
  PLOGI << "Copy time in [ms]: " << copy_time.duration();

  // We will have _totalUniqueSeqs_ dimensions
  ColPacks columnpacks;
  std::vector<uint64_t> colstartrange;
  auto lastcol = -1;
  PLOGI << "npairs = " << npairs;
  for (auto i = 0; i < npairs; i++) {
    if (seqsh_idx[i] != lastcol) {
      columnpacks.push_back({{}});
      colstartrange.push_back(i);
      lastcol = seqsh_idx[i];
    }
  }

  swatlib::TickTock colprep_time;
  colprep_time.tick();
  std::vector<bool> rcs(npairs);
  std::vector<int> xscores(npairs);
  std::vector<std::tuple<SeedPair, SeedPair>> seeds(npairs);
  auto start_time = std::chrono::system_clock::now();
  int skipps = 0;

  #pragma omp parallel for
  for (uint64_t j = 0; j < columnpacks.size(); j++) {
    auto current_column = seqsh_idx[colstartrange[j]];
    for (int i = colstartrange[j]; seqsh_idx[i] == current_column && i < npairs; i++) {
      elba::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);
      bool rc = false;
      std::array<SeedPair, 2> tmpSeeds;
      for (int count = 0; count < seed_count; ++count) {
        // argA (see KmerIntersectSR.hpp) == row == seqV
        ushort LocalSeedVOffset = (count == 0) ? cks->first.first : cks->second.first;
        // argB (see KmerIntersectSR.hpp) == col == seqH
        ushort LocalSeedHOffset = (count == 0) ? cks->first.second : cks->second.second;

        // Get sequences
        std::string seqH;
        std::string seqV;

        seqH = seqs[getHid(i)];
        // seqan::assign(seqH, dep_seqsh[i]);
        seqV = seqs[getVid(i)];
        // seqan::assign(seqV, dep_seqsv[i]);

        uint lenH = seqH.length();
        uint lenV = seqV.length();

        std::string seedH = seqH.substr(LocalSeedHOffset, seed_length);
        std::string seedV = seqV.substr(LocalSeedVOffset, seed_length);

        std::string twinH = _reversecomplement(seedH);

        if (twinH == seedV) {
          // std::string twinseqH(seqH);
          // std::reverse(std::begin(twinseqH), std::end(twinseqH));
          // std::transform(std::begin(twinseqH), std::end(twinseqH), std::begin(twinseqH), _complementbase);

          std::string twinseqH = seqs[getHTid(i)];

          LocalSeedHOffset = twinseqH.length() - LocalSeedHOffset - seed_length;
          assert(LocalSeedHOffset >= 0);
          assert(LocalSeedHOffset <= twinseqH.length());
          assert(LocalSeedVOffset <= seqV.length());
          if (count == 0) {
            rc = true;
            // seqs[getVid(i)] = seqV;
            // seqVs[i] = seqV;
            // seqs[getHTid(i)] = twinseqH;
            // seqHs[i] = twinseqH;
            columnpacks[j][0].push_back({i, getHTid(i), getVid(i)});
            rcs[i] = rc;
            tmpSeeds[count] = {LocalSeedHOffset, LocalSeedVOffset};
          } else {
            if (rc != true) {
              tmpSeeds[count] = {-1, -1};
              count_errors++;
            }
            tmpSeeds[count] = {LocalSeedHOffset, LocalSeedVOffset};
          }
        } else if (seedH == seedV) {
          assert(LocalSeedHOffset <= seqH.length());
          assert(LocalSeedVOffset <= seqV.length());
          if (count == 0) {
            rc = false;
            // seqVs[i] = seqV;
            // seqs[getVid(i)] = seqV;
            // seqHs[i] = seqH;
            // seqs[getHid(i)] = seqH;
            rcs[i] = rc;
            tmpSeeds[count] = {LocalSeedHOffset, LocalSeedVOffset};
            columnpacks[j][1].push_back({i, getHid(i), getVid(i)});
          } else {
            if (rc != false) {
              tmpSeeds[count] = {-1, -1};
              count_errors++;
            }
            tmpSeeds[count] = {LocalSeedHOffset, LocalSeedVOffset};
          }
        } else {
          tmpSeeds[count] = {-1, -1};
          #pragma omp critical
          {
            skipps++;
          }
        }
      }
      seeds[i] = {tmpSeeds[0], tmpSeeds[1]};
    }
  }

  colprep_time.tock();
  PLOGI << "Colprep time in [ms]: " << colprep_time.duration();
  auto end_time = std::chrono::system_clock::now();
  add_time("XA:LoganPreprocess", (ms_t(end_time - start_time)).count());

  start_time = std::chrono::system_clock::now();

  std::cout << "Skipped: " << skipps << std::endl;
  std::cout << "FLIP ERRORS: " << count_errors << std::endl;

  PLOGD.printf("We got %d columpacks", columnpacks.size());
  if (!noAlign) {
    this->runIPUAlign(columnpacks, seqsh_idx, seqsv_idx, seqs, seqHs, seqVs, seeds, xscores, xdrop, seed_length);
  } else {
    std::cout << "What are we doing?" << std::endl;
    exit(1);
  }

  // std::cout << "IPUAligner.cpp HERE 3" << std::endl;

  end_time = std::chrono::system_clock::now();
  add_time("XA:LoganAlign", (ms_t(end_time - start_time)).count());

  start_time = std::chrono::system_clock::now();

  // Compute stats
  /*
  if (count == 0)  //overwrite in the first seed
  {
    for (uint64_t i = 0; i < npairs; ++i) {
      ai[i].xscore = xscores[i];
      ai[i].rc = rcs[i];

      ai[i].begSeedH = xscores[i].begSeedH;
      ai[i].endSeedH = xscores[i].endSeedH;
      ai[i].begSeedV = xscores[i].begSeedV;
      ai[i].endSeedV = xscores[i].endSeedV;

      ai[i].seq_h_length = seqan::length(seqsh[i]);
      ai[i].seq_v_length = seqan::length(seqsv[i]);

      this is a bit redundant since we can extract it from seed
      ai[i].seq_h_seed_length = ai[i].endSeedH - ai[i].begSeedH;
      ai[i].seq_v_seed_length = ai[i].endSeedV - ai[i].begSeedV;

      GGGG: global idx over here to use in the FullDistVect for removing contained vertices/seqs
      ai[i].seq_h_g_idx = col_offset + std::get<1>(mattuples[lids[i]]);
      ai[i].seq_v_g_idx = row_offset + std::get<0>(mattuples[lids[i]]);
    }
  } else {
      for (uint64_t i = 0; i < npairs; ++i)
      {
              if (xscores[i].score > ai[i].xscore)
              {
                      std::cout << "Does this happen?" << std::endl;

                      ai[i].xscore = xscores[i].score;
                      ai[i].rc     = xscores[i].rc;

              ai[i].begSeedH = xscores[i].begSeedH;
              ai[i].endSeedH = xscores[i].endSeedH;
              ai[i].begSeedV = xscores[i].begSeedV;
              ai[i].endSeedV = xscores[i].endSeedV;

                      // @GGGG: this is a bit redundant since we can extract it from seed
                      ai[i].seq_h_seed_length = ai[i].endSeedH - ai[i].begSeedH;
                      ai[i].seq_v_seed_length = ai[i].endSeedV - ai[i].begSeedV;
              }
      }
  }



  */

  end_time = std::chrono::system_clock::now();
  add_time("XA:ComputeStats", (ms_t(end_time - start_time)).count());

  start_time = std::chrono::system_clock::now();
  std::vector<std::vector<int64_t>> ContainedSeqPerThread(numThreads);

  // Dump alignment info
  {
    // std::cout << "IPUAligner.cpp HERE 7" << std::endl;
    for (uint64_t i = 0; i < npairs; ++i) {
      // Only keep alignments that meet BELLA criteria
      bool passed = false;
      int tid = omp_get_thread_num();

      // GGGG: ai stores global idx to to store in ContainedSeqPerBatch
      // GGGG: in PostAlignDecision() we can mark as contained sequences as removable in ContainedSeqPerBatch and their local contained edges
      // GGGG: ContainedSeqPerBatch global indexes of contained sequences

      elba::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);
      // PostAlignDecision(ai[i], passed, ratioScoreOverlap, cks->dir, cks->dirT, cks->sfx, cks->sfxT, cks->overlap, noAlign, ContainedSeqPerThread[tid]);

      if (passed) {
        // GGGG: store updated seed start/end position in the CommonKmers pairs (the semantics of these pairs change wrt the original semantics but that's okay)
        cks->first.first = ai[i].begSeedV;   // start on ver sequence
        cks->second.first = ai[i].begSeedH;  // start on hor sequence

        cks->first.second = ai[i].endSeedV;   // end on ver sequence
        cks->second.second = ai[i].endSeedH;  // end on hor sequence

        cks->lenv = ai[i].seq_v_length;
        cks->lenh = ai[i].seq_h_length;

        cks->score = ai[i].xscore;
        cks->passed = passed;  // keep this
      }
    }
    // std::cout << "IPUAligner.cpp HERE 8" << std::endl;
  }

  int readcount = 0;
  for (int t = 0; t < numThreads; ++t) {
    readcount += ContainedSeqPerThread[t].size();
  }

  unsigned int readssofar = 0;
  ContainedSeqPerBatch.resize(readcount);

  // std::cout << "IPUAligner.cpp HERE 9" << std::endl;
  // Concatenate per-thread result
  for (int t = 0; t < numThreads; ++t) {
    copy(ContainedSeqPerThread[t].begin(), ContainedSeqPerThread[t].end(), ContainedSeqPerBatch.begin() + readssofar);
    readssofar += ContainedSeqPerThread[t].size();
  }

  end_time = std::chrono::system_clock::now();
  add_time("XA:StringOp", (ms_t(end_time - start_time)).count());
  std::cout << "total_cmps=" << ipu_total_cmps << std::endl;

  std::cout << "IPUAligner.cpp EXIT" << std::endl;
  return;
}
