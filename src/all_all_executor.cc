
#include "all_all_executor.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include <iostream>
#include "aligner.h"

using std::cout;
using std::get;
using std::make_pair;
using std::string;

void AllAllExecutor::EnqueueAlignment(const WorkItem& item) {
  if (!work_queue_->push(item)) {
    cout << "failed to push work item to queue.\n";
  }
}

void AllAllExecutor::FinishAndOutput(absl::string_view output_folder) {
  while (!work_queue_->empty()) {
    ;
    ;
  }
  cout << "Queue emptied, unblocking...\n";
  run_ = false;
  work_queue_->unblock();
  for (auto& f : threads_) {
    f.join();
  }
  cout << "All threads finished.\n";

  std::unordered_map<GenomePair, std::ofstream, PairHash> file_map;
  string output_dir("output_matches");
  size_t total_candidates = 0;

  for (const auto& matches_map : matches_per_thread_) {
    for (const auto& matches_kv : matches_map) {
      auto& genome_pair = matches_kv.first;
      auto& candidate_map = matches_kv.second;

      if (file_map.find(genome_pair) == file_map.end()) {
        // create the file
        string path = output_dir + "/" + genome_pair.first;
        struct stat info;
        if (stat(path.c_str(), &info) != 0) {
          // doesnt exist, create
          cout << "creating dir " << path << "\n";
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
        file_map[genome_pair] = std::ofstream(path);

        file_map[genome_pair]
            << absl::StrCat("# AllAll of ", genome_pair.first, " vs ",
                            genome_pair.second, ";\nRefinedMatches(\n[");
      }

      auto& out_file = file_map[genome_pair];

      for (auto& match : candidate_map) {
        total_candidates++;
        std::ostringstream ss;
        ss << "[" << match.first.first + 1 << ", " << match.first.second + 1
           << ", ";
        ss << std::fixed;
        ss.precision(7);
        auto& the_match = match.second;
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
           << "],\n";

        /*if (std::distance(it, candidate_map.end()) == 1)
          ss << "] ):";
        else
          ss << ",\n";*/

        out_file << ss.str();
      }
    }
  }

  for (auto& of_kv : file_map) {
    std::ofstream& of = of_kv.second;
    long pos = of.tellp();
    of.seekp(pos - 2);
    of << " ):";
  }
}

AllAllExecutor::AllAllExecutor(size_t num_threads, size_t capacity,
                               AlignmentEnvironments* envs, Parameters* params)
    : envs_(envs), params_(params), num_threads_(num_threads) {
  work_queue_.reset(new ConcurrentQueue<WorkItem>(capacity));
  matches_per_thread_.resize(num_threads);

  cout << "Start executor, id is " << id_.load() << "\n";

  num_active_threads_ = num_threads;
}

void AllAllExecutor::Initialize() {
  threads_.reserve(num_threads_);
  for (size_t i = 0; i < num_threads_; i++) {
    threads_.push_back(std::thread(&AllAllExecutor::Worker, this));
  }
}