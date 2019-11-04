
#include "all_all_dist.h"
#include "all_all_dist.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>

#define DEFAULT_ALIGNMENT_BATCH 500

using namespace std;

void AllAllDist::EnqueueAlignment(const WorkItem& item) {
  if (cur_num_alignments_ == 0) {
    req_ = MarshalledRequest();
    req_.CreateAlignmentRequest();
  }

  // append new request
  int abs1 = get<0>(item)->ID();
  int abs2 = get<1>(item)->ID();
  req_.AddAlignment(abs1, abs2);
  cur_num_alignments_++;
  total_alignments_++;

  // if cur_num == batch size, ship it off

  if (cur_num_alignments_ == DEFAULT_ALIGNMENT_BATCH) {
    request_queue_->push(std::move(req_));
    if (req_.buf.data() != nullptr) {
      cout << "req buf was not null after move??\n";
    }
    cur_num_alignments_ = 0;
    outstanding_++;
  }
}

void AllAllDist::ProcessResult(const char* result_buf) {
  cout << "processing result\n";
  const DistMatchResult* matches;
  size_t num_matches;
  AlignmentResults::DeserializeResults(&matches, &num_matches, result_buf);

  for (size_t i = 0; i < num_matches; i++) {
    const auto& match = matches[i];
    auto& seq1 = sequences_[match.abs_seq_1];
    auto& seq2 = sequences_[match.abs_seq_2];
    auto genome_pair = std::make_pair(seq1.Seq(), seq2.Seq());

    mu_.Lock();
    if (file_map.find(genome_pair) == file_map.end()) {
      // create the file
      string path = absl::StrCat(output_dir_, "/", genome_pair.first);
      struct stat info;
      if (stat(path.c_str(), &info) != 0) {
        // doesnt exist, create
        // cout << "creating dir " << path << "\n";
        int e = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (e != 0) {
          cout << "could not create output dir " << path << ", exiting ...\n";
          exit(0);
        }
      } else if (!(info.st_mode & S_IFDIR)) {
        // exists but not dir
        cout << "output dir exists but is not dir, exiting ...\n";
        exit(0);
      }  // else, dir exists,

      absl::StrAppend(&path, "/", genome_pair.second);
      cout << "opening file " << path << "\n";

      LockedStream ls;
      ls.out_stream = std::ofstream(path);
      file_map[genome_pair] = std::move(ls);

      string line = absl::StrCat("# AllAll of ", genome_pair.first, " vs ",
                                 genome_pair.second, ";\nRefinedMatches(\n[");
      file_map[genome_pair].out_stream << line;
    }

    auto& out_file = file_map[genome_pair];
    mu_.Unlock();

    // write the match to the file now that we have a file handle

    std::ostringstream ss;
    ss << "[" << seq1.GenomeIndex() + 1 << ", " << seq2.GenomeIndex() + 1
       << ", ";
    ss << std::fixed;
    ss.precision(7);
    auto& the_match = match.m;
    ss << round(the_match.score * 10000000.0f) / 10000000.0f << ", ";
    if (the_match.distance >= 45.0f)
      ss << int(the_match.distance);
    else if (the_match.distance > 0.1f) {
      ss.precision(4);
      ss << round(the_match.distance * 10000.0f) / 10000.0f;
    } else {
      ss.precision(8);
      ss << round(the_match.distance * 100000000.0f) / 100000000.0f;
    }
    ss.precision(8);

    ss << ", " << the_match.seq1_min + 1 << ".." << the_match.seq1_max + 1
       << ", " << the_match.seq2_min + 1 << ".." << the_match.seq2_max + 1
       << ", " << round(the_match.variance * 100000000.0f) / 100000000.0f
       << ", " << the_match.cluster_size << "],\n";

    {
      // is this lock necessary? are streams thread safe?
      absl::MutexLock l(&out_file.mu);
      out_file.out_stream << ss.str();
    }
    total_matches_++;
  }

  outstanding_--;
}

void AllAllDist::Finish() {
  // wait until all outstanding alignments are completed and written to the
  // files

  if (cur_num_alignments_ > 0) {
    request_queue_->push(std::move(req_));
    if (req_.buf.data() != nullptr) {
      cout << "req buf was not null after move??\n";
    }
    cur_num_alignments_ = 0;
    outstanding_++;
  }

  while (outstanding_.load() > 0)
    ;

  // finalize output files
  for (auto& v : file_map) {
    auto& ls = v.second;

    std::ofstream& of = ls.out_stream;
    long pos = of.tellp();
    of.seekp(pos - 2);
    of << " ]):";
  }

  cout << "Total matches: " << total_matches_.load() << "\n";
  cout << "Total All-All full alignments: " << total_alignments_.load() << "\n";
}