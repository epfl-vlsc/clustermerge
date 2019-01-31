
#include "all_all_executor.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <ctime>
#include "json.hpp"
#include "absl/strings/str_cat.h"
#include "aligner.h"

using std::cout;
using std::get;
using std::make_pair;
using std::string;
using namespace std::literals::chrono_literals;

void AllAllExecutor::EnqueueAlignment(const WorkItem& item) {
  if (!work_queue_->push(item)) {
    cout << "failed to push work item to queue.\n";
  }
}

void AllAllExecutor::FinishAndOutput(const string& output_dir) {
  cout << "waiting for work queue to empty\n";
  while (!work_queue_->empty()) {
    std::this_thread::sleep_for(1s);
    //cout << "work queue has " << work_queue_->size() << " entries left\n";
  }
  //cout << "Queue emptied, unblocking...\n";
  run_ = false;
  work_queue_->unblock();
  for (auto& f : threads_) {
    f.join();
  }
  queue_measure_thread_.join();
  cout << "All threads finished.\n";

  absl::flat_hash_map<GenomePair, std::ofstream> file_map;
  struct stat info;
  if (stat(output_dir.c_str(), &info) != 0) {
    // doesnt exist, create
    cout << "creating dir " << output_dir << "\n";
    int e = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (e != 0) {
      cout << "could not create output dir " << output_dir << ", exiting ...\n";
      exit(0);
    }
  } else if (!(info.st_mode & S_IFDIR)) {
    // exists but not dir
    cout << "output dir exists but is not dir, exiting ...\n";
    exit(0);
  } else {
    // dir exists, nuke
    // im too lazy to do this the proper way
    string cmd = absl::StrCat("rm -rf ", output_dir, "/*");
    cout << "dir " << output_dir << " exists, nuking ...\n";
    int nuke_result = system(cmd.c_str());
    if (nuke_result != 0) {
      cout << "Could not nuke dir " << output_dir << "\n";
      exit(0);
    }
  }

  size_t total_candidates = 0;

  for (const auto& matches_map : matches_per_thread_) {
    for (const auto& matches_kv : matches_map) {
      auto& genome_pair = matches_kv.first;
      auto& candidate_map = matches_kv.second;

      if (file_map.find(genome_pair) == file_map.end()) {
        // create the file
        string path = absl::StrCat(output_dir, "/", genome_pair.first);
        struct stat info;
        if (stat(path.c_str(), &info) != 0) {
          // doesnt exist, create
          //cout << "creating dir " << path << "\n";
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
        //cout << "opening file " << path << "\n";
        file_map[genome_pair] = std::ofstream(path);

        string line = absl::StrCat("# AllAll of ", genome_pair.first, " vs ",
                                   genome_pair.second, ";\nRefinedMatches(\n[");
        file_map[genome_pair] << line;
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
           << ", " << the_match.cluster_size
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
  cout << "Total matches: " << total_candidates << "\n";
  cout << "Total All-All full alignments: " << num_full_alignments_.load()
       << "\n";
  cout << "Total threshold alignments: " << num_pass_threshold_.load() << "\n";

  // queue size stats
  std::vector<std::pair<size_t, size_t>> values;
  for (size_t i = 0; i < timestamps_.size(); i++) {
    values.push_back(std::make_pair(timestamps_[i], queue_sizes_[i]));
  }

  nlohmann::json j(values);

  //std::cout << "dumping queue sizes ...\n";
  std::ofstream o("queue.json");

  o << std::setw(2) << j << std::endl;

}

AllAllExecutor::AllAllExecutor(size_t num_threads, size_t capacity,
                               AlignmentEnvironments* envs, const Parameters* params)
    : envs_(envs), params_(params), num_threads_(num_threads) {
  work_queue_.reset(new ConcurrentQueue<WorkItem>(capacity));
  matches_per_thread_.resize(num_threads);

  // cout << "Start executor, id is " << id_.load() << "\n";

  num_active_threads_ = num_threads;
  timestamps_.reserve(100000);
  queue_sizes_.reserve(100000);

  queue_measure_thread_ = std::thread([this](){
      // get timestamp, queue size
      cout << "queue measure thread starting ...\n";
      while(run_) {
        time_t result = std::time(nullptr);
        timestamps_.push_back(static_cast<long int>(result));
        queue_sizes_.push_back(work_queue_->size());
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
      }
      cout << "queue measure thread finished\n";
    }
  );


}

void AllAllExecutor::Initialize() {
  threads_.reserve(num_threads_);
  for (size_t i = 0; i < num_threads_; i++) {
    threads_.push_back(std::thread(&AllAllExecutor::Worker, this));
  }
}
