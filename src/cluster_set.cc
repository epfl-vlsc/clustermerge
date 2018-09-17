
#include "cluster_set.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include "agd/json.hpp"
#include "aligner.h"
#include "debug.h"

using std::vector;

ClusterSet ClusterSet::MergeClusters(ClusterSet& other,
                                     ProteinAligner* aligner) {
  // this is the money method

  // merge clusters, clusters can "disappear" from either
  // set, so we just create a new one and resize its internal
  // cluster vector for a single alloc

  ClusterSet new_cluster_set(clusters_.size() + other.clusters_.size());

  ProteinAligner::Alignment alignment;
  agd::Status s;
  for (auto& c : clusters_) {
    for (auto& c_other : other.clusters_) {
      // s = c.AlignReps(c_other, alignment);
      // analyze 4 different merge cases
      // construct new_cluster_set

      if (!c_other.IsFullyMerged() && c.PassesThreshold(c_other, aligner)) {
        // std::cout << "passed threshold, aligning ...\n";
        s = c.AlignReps(c_other, &alignment, aligner);
        // does c contain c_other fully
        if (alignment.seq1_max - alignment.seq1_min == c.Rep().Seq().size()) {
          // c rep is fully covered, probably contained in c_other rep
          // std::cout << "C rep is fully covered\n";
          // add all seqs in c into c_other, mark c fully merged, add c_other to
          // new cluster set

          for (const auto& seq : c.Sequences()) {
            c_other.AddSequence(seq);
          }
          c.SetFullyMerged();
          break;  // c is gone, move to next

        } else if (alignment.seq2_max - alignment.seq2_min ==
                   c_other.Rep().Seq().size()) {
          // c_other rep is fully covered, probably contained in c rep
          // std::cout << "C other rep is fully covered\n";
          // add all seqs in c_other into c, mark c_other fully merged, add c to
          // new cluster set
          for (const auto& seq : c_other.Sequences()) {
            c.AddSequence(seq);
          }
          c_other.SetFullyMerged();

        } else {
          // situation is :
          // |-------------------|
          //            |-------------------|
          // or opposite. If the coverage of one is within X
          // of total residues, merge completely. Otherwise, we just
          // add matching seqs from one to the other
          // std::cout << "reps are partially overlapped\n";

          auto c_num_uncovered =
              c.Rep().Seq().size() - (alignment.seq1_max - alignment.seq1_min);
          auto c_other_num_uncovered =
              c_other.Rep().Seq().size() -
              (alignment.seq2_max - alignment.seq2_min);

          if (c_num_uncovered < aligner->Params()->max_n_aa_not_covered) {
            // they are _almost_ overlapped, merge completely
            // std::cout << "Nearly complete overlap, merging c into c_other\n";
            for (const auto& seq : c.Sequences()) {
              c_other.AddSequence(seq);
            }
            c.SetFullyMerged();
            break;

          } else if (c_other_num_uncovered <
                     aligner->Params()->max_n_aa_not_covered) {
            // std::cout << "Nearly complete overlap, merging c_other into c\n";
            for (const auto& seq : c_other.Sequences()) {
              c.AddSequence(seq);
            }
            c_other.SetFullyMerged();
          } else {
            // add c_other_rep into c
            // for each sequence in c_other, add if it matches c rep
            // keep both clusters
            // std::cout << "merging and keeping both clusters\n";
            c.Merge(c_other, aligner);
          }
        }
      }  // if passes threshold
    }
    if (!c.IsFullyMerged()) {
      new_cluster_set.clusters_.push_back(std::move(c));
    }
  }

  for (auto& c_other : other.clusters_) {
    if (!c_other.IsFullyMerged()) {
      // push any not fully merged cluster into the new set and we are done
      new_cluster_set.clusters_.push_back(std::move(c_other));
    }
  }
  // std::cout << "new cluster set is \n";
  // new_cluster_set.DebugDump();

  return new_cluster_set;
}

void ClusterSet::DebugDump() const {
  std::cout << "Dumping " << clusters_.size() << " clusters in set... \n";
  for (const auto& cluster : clusters_) {
    std::cout << "\tCluster seqs:\n";
    for (const auto& seq : cluster.Sequences()) {
      std::cout << "\t\tGenome: " << seq.Genome() << ", sequence: "
                << PrintNormalizedProtein(seq.Seq().data(), seq.Seq().length())
                << "\n\n";
    }
  }
}

void ClusterSet::ScheduleAlignments(AllAllExecutor* executor) {
  for (const auto& cluster : clusters_) {
    for (auto it = cluster.Sequences().begin(); it != cluster.Sequences().end();
         it++) {
      for (auto itt = next(it); itt != cluster.Sequences().end(); itt++) {
        // std::cout << "Scheduling alignment...\n";
        AllAllExecutor::WorkItem item = std::make_tuple(&(*it), &(*itt));
        executor->EnqueueAlignment(item);
      }
    }
  }
}

void ClusterSet::DumpJson() const {
  vector<vector<size_t>> cluster_seqs;
  for (const auto& c : clusters_) {
    vector<size_t> seq_ids;
    for (const auto& s : c.Sequences()) {
      seq_ids.push_back(s.ID());
    }

    cluster_seqs.push_back(seq_ids);
  }

  nlohmann::json j(cluster_seqs);

  std::cout << "dumping clusters ...\n";
  std::ofstream o("clusters.json");

  o << std::setw(2) << j << std::endl;
}
