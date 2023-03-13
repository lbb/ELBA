/* Created by Saliya Ekanayake on 2019-07-05 and modified by Giulia Guidi on 4/14/2021. */

#include "../../include/pw/IPUAligner.hpp"

#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Log.h>
#include <string_view>

uint _minOverlapLenL = 5000;


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
  assert(false);
  return ' ';
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
  const ipu::SWConfig SW_CONFIGURATION = {
      .gapInit = -1,
      .gapExtend = -1,
      .matchValue = 1,
      .mismatchValue = -1,
      .ambiguityValue = -1,
      .similarity = swatlib::Similarity::nucleicAcid,
      .datatype = swatlib::DataType::nucleicAcid,
  };

  const ipu::IPUAlgoConfig ALGOCONFIG = {
      .numVertices = 1472,
      .maxSequenceLength = 19295,
      .maxComparisonsPerVertex = 20,
      .vertexBufferSize = 160000,
      .vtype = ipu::VertexType::xdroprestrictedseedextend,
      .fillAlgo = ipu::Algorithm::greedy,
      .forwardOnly = false,
      .ioTiles = 0,
      .xDrop = xdrop,
      .bandPercentageXDrop = 0.45,
      .seedLength = seed_length,
  };
  const auto ipus = 1;
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
void IPUAligner::runIPUAlign(const std::vector<std::string> &seqHs, const std::vector<std::string> &seqVs, std::vector<std::tuple<SeedPair, SeedPair>> seeds, std::vector<int> &xscores, int xdrop, int seed_length) {
  ipu_total_cmps +=  seqVs.size() ;
  
  auto npairs = seqVs.size();
  std::vector<ipu::Comparison> comparisons(npairs);
  std::vector<std::string_view> sequences(npairs * 2);

  // #pragma omp parallel for
  for (int32_t i = 0; i < npairs; ++i) {
    sequences[2 * i] = seqVs[i];
    sequences[2 * i + 1] = seqHs[i];
    int offsetH1, offsetV1;
    int offsetH2, offsetV2;
    SeedPair s1, s2;
    std::tie(s1, s2) = seeds[i];

    std::tie(offsetH1, offsetV1) = s1;
    std::tie(offsetH2, offsetV2) = s2;
    // if (!(seqVs[i].length() >= offsetV)) {
    //   std::cerr << "REMAP (V) len:" << seqVs[i].length() << ", offset" << offsetV << std::endl;
    // }
    // if (!(seqHs[i].length() >= offsetH)) {
    //   std::cerr << "REMAP (H) len:" << seqHs[i].length() << ", offset" << offsetH << std::endl;
    // }
    assert(seqVs[i].length() >= offsetV && "VSeq offset");
    assert(seqHs[i].length() >= offsetH && "HSeq offset");
    comparisons[i] = {
        2 * i, 2 * i + 1,
        {{{offsetV1, offsetH1}, {offsetV2, offsetH2}}}};
  }

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

  // std::cerr << "maxV: " << maxV << std::endl;
  // std::cerr << "maxSeedV: " << maxSeedV << std::endl;

  std::vector<ipu::Batch> batches = driver_algo->create_batches(sequences, comparisons);

  auto runt = std::chrono::system_clock::now();
  std::vector<ipu::batchaffine::Job *> jobs(batches.size());
  for (size_t i = 0; i < batches.size(); i++) {
    jobs[i] = driver_algo->async_submit(&batches[i]);
  }

  std::vector<ipu::BlockAlignmentResults> results(batches.size());
  // #pragma omp parallel for
  for (size_t i = 0; i < batches.size(); i++) {
    driver_algo->blocking_join(*jobs[i]);
  }
  auto end_runt = std::chrono::system_clock::now();
  std::cout << "Local TIMEms:::::::::::::::::  " << (ms_t(end_runt - runt)).count() << std::endl;
  total_time += (ms_t(end_runt - runt)).count();
  std::cout << "TIMEms:::::::::::::::::  " << total_time << std::endl;

  for (size_t i = 0; i < batches.size(); i++) {
    int j = i;
    auto &batch = batches[j];
    ipu::BlockAlignmentResults res = batch.get_result();
    // #pragma omp for
    for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
      auto orig_i = batch.origin_comparison_index[i];
      if (orig_i >= 0) {
        xscores[batch.origin_comparison_index[i]] = res.scores[i][0] > res.scores[i][1] ? res.scores[i][0] : res.scores[i][1];
      }
    }
  }
}

// @NOTE This is hard-coded to the number of seeds being <= 2
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
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int numThreads = 1;
#ifdef THREADED
#pragma omp parallel
  {
    numThreads = omp_get_num_threads();
  }
#endif

  uint64_t npairs = seqan::length(seqsh);

  lfs << "processing batch of size " << npairs << " with " << numThreads << " threads " << std::endl;

  // for multiple seeds we store the seed with the highest identity
  std::vector<IPumaAlignmentInfo> ai(npairs);


  // bool *strands = new bool[npairs];
  // int  *xscores = new int[npairs];
  // TSeed  *seeds = new TSeed[npairs];

  /* GGGG: seed_count is hardcoded here (2) */
    int count_errors = 0;
    std::vector<string> seqHs(npairs);
    std::vector<string> seqVs(npairs);
    std::vector<bool> rcs(npairs);
    std::vector<int> xscores(npairs);
    std::vector<std::tuple<SeedPair, SeedPair>> seeds(npairs);
    auto start_time = std::chrono::system_clock::now();
    int skipps = 0;

    for (uint64_t i = 0; i < npairs; ++i) {
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

      seqan::assign(seqH, seqsh[i]);
      seqan::assign(seqV, seqsv[i]);

      uint lenH = seqH.length();
      uint lenV = seqV.length();

      // Get seed string
      assert(LocalSeedHOffset <= seqH.length());
      assert(LocalSeedHOffset+seed_length <= seqH.length());
      assert(LocalSeedVOffset <= seqV.length());
      assert(LocalSeedVOffset+seed_length <= seqV.length());
      std::string seedH = seqH.substr(LocalSeedHOffset, seed_length);
      std::string seedV = seqV.substr(LocalSeedVOffset, seed_length);

      std::string twinH = _reversecomplement(seedH);

      if (twinH == seedV) {
        std::string twinseqH(seqH);
        std::reverse(std::begin(twinseqH), std::end(twinseqH));
        std::transform(std::begin(twinseqH), std::end(twinseqH), std::begin(twinseqH), _complementbase);
        LocalSeedHOffset = twinseqH.length() - LocalSeedHOffset - seed_length;
        assert(LocalSeedHOffset >= 0);
        assert(LocalSeedHOffset <= twinseqH.length());
        assert(LocalSeedVOffset <= seqV.length());
        if (count == 0) {
          rc = true;
          seqVs[i] = seqV;
          seqHs[i] = twinseqH;
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
          seqVs[i] = seqV;
          seqHs[i] = seqH;
          rcs[i] = rc;
          tmpSeeds[count] = {LocalSeedHOffset, LocalSeedVOffset};
        } else {
          if (rc != false) {
            tmpSeeds[count] = {-1, -1};
            count_errors++;
          }
          tmpSeeds[count] = {LocalSeedHOffset, LocalSeedVOffset};
        }
      } else {
        tmpSeeds[count] = {-1, -1};
        skipps++;
      }
    }
    seeds[i] = {tmpSeeds[0], tmpSeeds[1]};
  }

    auto end_time = std::chrono::system_clock::now();
    add_time("XA:LoganPreprocess", (ms_t(end_time - start_time)).count());

    start_time = std::chrono::system_clock::now();

    std::cout << "Skipped: " <<  skipps << std::endl;
    std::cout << "FLIP ERRORS: " << count_errors << std::endl;

    if (!noAlign) {
      this->runIPUAlign(seqHs, seqVs, seeds, xscores, xdrop, seed_length);
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

  // std::cout << "IPUAligner.cpp EXIT" << std::endl;
  return;
}
