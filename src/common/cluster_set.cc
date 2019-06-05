
#include "cluster_set.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "absl/container/flat_hash_set.h"
#include "aligner.h"
#include "candidate_map.h"
#include "debug.h"
#include "json.hpp"
#include "merge_executor.h"

using std::make_tuple;
using std::vector;

void free_func(void* data, void* hint) {
  delete reinterpret_cast<char*>(data);
}

void ClusterSet::BuildMarshalledResponse(int id, MarshalledResponse* response) {
  agd::Buffer buf;
  buf.reserve(256);
  buf.AppendBuffer(reinterpret_cast<char*>(&id), sizeof(int));

  ClusterSetHeader h;
  h.num_clusters = 0;  // set after once we know the value
  buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(ClusterSetHeader));

  for (const auto& c : clusters_) { // keep the fully merged for controller to remove
    ClusterHeader ch;
    ch.fully_merged = c.IsFullyMerged();
    ch.num_seqs = c.Sequences().size();
    buf.AppendBuffer(reinterpret_cast<char*>(&ch), sizeof(ClusterHeader));
    for (auto s : c.Sequences()) {
      uint32_t i = s.ID();
      buf.AppendBuffer(reinterpret_cast<char*>(&i), sizeof(int));
    }
  }

  char* data = buf.mutable_data() + sizeof(int);
  ClusterSetHeader* hp = reinterpret_cast<ClusterSetHeader*>(data);
  hp->num_clusters = clusters_.size();
  // hand the buf pointer to the message
  response->msg = zmq::message_t(buf.release_raw(), buf.size(), free_func, NULL);
}

ClusterSet::ClusterSet(MarshalledClusterSet& marshalled_set,
                       const std::vector<Sequence>& sequences) {
  //std::vector<Cluster> clusters(marshalled_set.NumClusters());
  MarshalledClusterView cluster;
  while (marshalled_set.NextCluster(&cluster)) {
  //for (size_t cs_i = 0; cs_i < set_proto.clusters_size(); cs_i++) {
    //const auto& cluster_proto = set_proto.clusters(cs_i);
    // std::cout << "cluster has " << cluster_proto.indexes_size() << " seqs\n";
    //std::cout << "adding cluster with rep idx " << cluster.SeqIndex(0) << " and num seqs = " << cluster.NumSeqs() << "\n";
    Cluster c(sequences[cluster.SeqIndex(0)]);
    uint32_t num_seqs = cluster.NumSeqs();
    for (size_t seq_i = 1; seq_i < num_seqs; seq_i++) {
      c.AddSequence(sequences[cluster.SeqIndex(seq_i)]);
    }
    if (cluster.IsFullyMerged()) {
      c.SetFullyMerged();
    }
    clusters_.push_back(std::move(c));
  }
}

ClusterSet::ClusterSet(MarshalledClusterSetView& marshalled_set,
                       const std::vector<Sequence>& sequences) {
  // yeah its copied from above idc
  //std::vector<Cluster> clusters(marshalled_set.NumClusters());
  MarshalledClusterView cluster;
  while (marshalled_set.NextCluster(&cluster)) {
  //for (size_t cs_i = 0; cs_i < set_proto.clusters_size(); cs_i++) {
    //const auto& cluster_proto = set_proto.clusters(cs_i);
    // std::cout << "cluster has " << cluster_proto.indexes_size() << " seqs\n";
    //std::cout << "adding cluster with rep idx " << cluster.SeqIndex(0) << " and num seqs = " << cluster.NumSeqs() << "\n";
    Cluster c(sequences[cluster.SeqIndex(0)]);
    uint32_t num_seqs = cluster.NumSeqs();
    for (size_t seq_i = 1; seq_i < num_seqs; seq_i++) {
      c.AddSequence(sequences[cluster.SeqIndex(seq_i)]);
    }
    if (cluster.IsFullyMerged()) {
      c.SetFullyMerged();
    }
    clusters_.push_back(std::move(c));
  }
}

/*void ClusterSet::ConstructProto(cmproto::ClusterSet* set_proto) {
  for (const auto& c : clusters_) {
    auto* c_proto = set_proto->add_clusters();
    for (const auto& s : c.Sequences()) {
      c_proto->add_indexes(s.ID());
    }
    c_proto->set_fully_merged(c.IsFullyMerged());
  }
}*/

ClusterSet ClusterSet::MergeClustersParallel(ClusterSet& other,
                                             MergeExecutor* executor) {
  ClusterSet new_cluster_set(clusters_.size() + other.clusters_.size());

  MultiNotification n;
  for (auto& c : clusters_) {
    // enqueue each comparison between c and all cluster in other
    MergeExecutor::WorkItem item = make_tuple(&c, &other, &n);
    executor->EnqueueMerge(item);
  }

  n.SetMinNotifies(clusters_.size());
  n.WaitForNotification();

  for (auto& c_other : other.clusters_) {
    if (!c_other.IsFullyMerged()) {
      // push any not fully merged cluster into the new set and we are done
      new_cluster_set.clusters_.push_back(std::move(c_other));
    }
  }
  for (auto& c : clusters_) {
    if (!c.IsFullyMerged()) {
      // push any not fully merged cluster into the new set and we are done
      new_cluster_set.clusters_.push_back(std::move(c));
    }
  }

  // sort so that larger rep clusters come first, leading to
  // better scheduling overlap of work
  std::sort(new_cluster_set.clusters_.begin(), new_cluster_set.clusters_.end(),
            [](Cluster& a, Cluster& b) {
              return a.Rep().Seq().size() > b.Rep().Seq().size();
            });

  // new_cluster_set.RemoveDuplicates();
  return new_cluster_set;
}

void ClusterSet::MergeClusterLocked(Cluster* cluster, ProteinAligner* aligner) {
  // this func called from multiple threads, but we are guaranteed that
  // `cluster` is accessed exclusively clusters of `this` must be locked before
  // modifying, because other threads may be accessing
  agd::Status s;
  ProteinAligner::Alignment alignment;
  for (auto& c_other : clusters_) {
    if (!c_other.IsFullyMerged() &&
        cluster->PassesThreshold(c_other, aligner)) {
      // std::cout << "passed threshold, aligning ...\n";
      s = cluster->AlignReps(c_other, &alignment, aligner);

      // situation is :
      // |-------------------|
      //            |-------------------|
      // or opposite. If the coverage of one is within X
      // of total residues, merge completely. Otherwise, we just
      // add matching seqs from one to the other
      // std::cout << "reps are partially overlapped\n";

      auto c_num_uncovered = cluster->Rep().Seq().size() -
                             (alignment.seq1_max - alignment.seq1_min);
      auto c_other_num_uncovered = c_other.Rep().Seq().size() -
                                   (alignment.seq2_max - alignment.seq2_min);

      if (c_num_uncovered < aligner->Params()->max_n_aa_not_covered &&
          alignment.score > aligner->Params()->min_full_merge_score) {
        // they are _almost_ overlapped, merge completely
        // std::cout << "Nearly complete overlap, merging c into c_other,
        // score is " << alignment.score << "\n";

        // TODO move this code into class Cluster
        c_other.Lock();
        if (c_other.IsFullyMerged()) {
          c_other.Unlock();
          continue;
        }
        for (const auto& seq : cluster->Sequences()) {
          c_other.AddSequence(seq);
        }
        cluster->SetFullyMerged();
        c_other.Unlock();
        break;

      } else if (c_other_num_uncovered <
                     aligner->Params()->max_n_aa_not_covered &&
                 alignment.score > aligner->Params()->min_full_merge_score) {
        // std::cout << "Nearly complete overlap, merging c_other into c,
        // score is " << alignment.score << "\n";
        // TODO move this code into class Cluster
        c_other.Lock();
        if (c_other.IsFullyMerged()) {
          c_other.Unlock();
          continue;
        }
        for (const auto& seq : c_other.Sequences()) {
          cluster->AddSequence(seq);
        }
        c_other.SetFullyMerged();
        c_other.Unlock();
      } else {
        // add c_other_rep into c
        // for each sequence in c_other, add if it matches c rep
        // keep both clusters
        // std::cout << "merging and keeping both clusters\n";
        c_other.Lock();
        if (c_other.IsFullyMerged()) {
          c_other.Unlock();
          continue;
        }
        cluster->Merge(&c_other, aligner);
        c_other.Unlock();
      }
    }  // if passes threshold
  }    // for c_other in clusters
}

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
      if (!c_other.IsFullyMerged() && c.PassesThreshold(c_other, aligner)) {
        // std::cout << "passed threshold, aligning ...\n";
        s = c.AlignReps(c_other, &alignment, aligner);

        // situation is :
        // |-------------------|
        //            |-------------------|
        // or opposite. If the coverage of one is within X
        // of total residues, merge completely. Otherwise, we just
        // add matching seqs from one to the other
        // std::cout << "reps are partially overlapped\n";

        auto c_num_uncovered =
            c.Rep().Seq().size() - (alignment.seq1_max - alignment.seq1_min);
        auto c_other_num_uncovered = c_other.Rep().Seq().size() -
                                     (alignment.seq2_max - alignment.seq2_min);

        if (c_num_uncovered < aligner->Params()->max_n_aa_not_covered &&
            alignment.score > aligner->Params()->min_full_merge_score) {
          // they are _almost_ overlapped, merge completely
          // std::cout << "Nearly complete overlap, merging c into c_other,
          // score is " << alignment.score << "\n";
          for (const auto& seq : c.Sequences()) {
            c_other.AddSequence(seq);
          }
          c.SetFullyMerged();
          break;

        } else if (c_other_num_uncovered <
                       aligner->Params()->max_n_aa_not_covered &&
                   alignment.score > aligner->Params()->min_full_merge_score) {
          // std::cout << "Nearly complete overlap, merging c_other into c,
          // score is " << alignment.score << "\n";
          for (const auto& seq : c_other.Sequences()) {
            c.AddSequence(seq);
          }
          c_other.SetFullyMerged();
        } else {
          // add c_other_rep into c
          // for each sequence in c_other, add if it matches c rep
          // keep both clusters
          // std::cout << "merging and keeping both clusters\n";
          c.Merge(&c_other, aligner);
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

ClusterSet ClusterSet::MergeCluster(Cluster& c_other, ProteinAligner* aligner, bool& worker_signal_) {
  // for the dist version, we keep fully merged clusters around,
  // because this is a partial merge of two large sets,
  // the results of which need to be merged by the controller
  ClusterSet new_cluster_set(clusters_.size() + 1);

  ProteinAligner::Alignment alignment;
  agd::Status s;
  for (auto& c : clusters_) {
    if(worker_signal_)  {
      std::cout << "Breaking partial merge.\n";
      break;
    }
    if (!c_other.IsFullyMerged() && c.PassesThreshold(c_other, aligner)) {
      // std::cout << "passed threshold, aligning ...\n";
      s = c.AlignReps(c_other, &alignment, aligner);

      // situation is :
      // |-------------------|
      //            |-------------------|
      // or opposite. If the coverage of one is within X
      // of total residues, merge completely. Otherwise, we just
      // add matching seqs from one to the other
      // std::cout << "reps are partially overlapped\n";

      auto c_num_uncovered =
          c.Rep().Seq().size() - (alignment.seq1_max - alignment.seq1_min);
      auto c_other_num_uncovered = c_other.Rep().Seq().size() -
                                   (alignment.seq2_max - alignment.seq2_min);

      if (c_num_uncovered < aligner->Params()->max_n_aa_not_covered &&
          alignment.score > aligner->Params()->min_full_merge_score) {
        // they are _almost_ overlapped, merge completely
        // std::cout << "Nearly complete overlap, merging c into c_other,
        // score is " << alignment.score << "\n";
        for (const auto& seq : c.Sequences()) {
          c_other.AddSequence(seq);
        }
        c.SetFullyMerged();
      } else if (c_other_num_uncovered <
                     aligner->Params()->max_n_aa_not_covered &&
                 alignment.score > aligner->Params()->min_full_merge_score) {
        // std::cout << "Nearly complete overlap, merging c_other into c,
        // score is " << alignment.score << "\n";
        for (const auto& seq : c_other.Sequences()) {
          c.AddSequence(seq);
        }
        c_other.SetFullyMerged();
      } else {
        // add c_other_rep into c
        // for each sequence in c_other, add if it matches c rep
        // keep both clusters
        // std::cout << "merging and keeping both clusters\n";
        c.Merge(&c_other, aligner);
      }
    }  // if passes threshold

    new_cluster_set.clusters_.push_back(std::move(c));
  }

  // we can leave out without confusing the controller
  if (!c_other.IsFullyMerged())
    new_cluster_set.clusters_.push_back(std::move(c_other));

  return new_cluster_set;
}

void ClusterSet::ScheduleAlignments(AllAllExecutor* executor) {
  // removing duplicate clusters (clusters with same sequences)
  // for some reason, absl::InlinedVector doesnt work here
  absl::flat_hash_set<std::vector<size_t>> set_map;

  size_t num_dups_found = 0;
  for (auto& c : clusters_) {
    std::vector<size_t> cluster_set;
    for (const auto& s : c.Sequences()) {
      cluster_set.push_back(s.ID());
    }
    std::sort(cluster_set.begin(), cluster_set.end());

    auto result = set_map.insert(std::move(cluster_set));
    if (!result.second) {
      c.SetDuplicate();
      num_dups_found++;
    }
  }

  std::cout << "Found " << num_dups_found << " duplicate clusters.\n";
  // sort by residue total first
  // to schedule the heaviest computations first
  std::cout << "sorting clusters ...\n";
  std::sort(clusters_.begin(), clusters_.end(), [](Cluster& a, Cluster& b) {
    return a.Sequences().size() > b.Sequences().size();
  });
  std::cout << "done sorting clusters.\n";

  CandidateMap candidate_map(10000000);  // only a few MB
  int num_avoided = 0;

  for (const auto& cluster : clusters_) {
    // std::cout << "Cluster has " << cluster.Sequences().size() << " seqs\n";
    if (cluster.IsDuplicate()) {
      continue;
    }
    for (auto it = cluster.Sequences().begin(); it != cluster.Sequences().end();
         it++) {
      for (auto itt = next(it); itt != cluster.Sequences().end(); itt++) {
        auto* seq1 = &(*it);
        auto* seq2 = &(*itt);

        if (seq1->Genome() == seq2->Genome() &&
            seq1->GenomeIndex() == seq2->GenomeIndex()) {
          // not sure if this can actually happen yet, but no need to align
          // against self
          continue;
        }

        if (seq1->GenomeSize() > seq2->GenomeSize() ||
            ((seq1->GenomeSize() == seq2->GenomeSize()) &&
             seq1->Genome() > seq2->Genome())) {
          std::swap(seq1, seq2);
        }

        if (seq1->Genome() == seq2->Genome() &&
            seq1->GenomeIndex() > seq2->GenomeIndex()) {
          std::swap(seq1, seq2);
        }

        auto abs_seq_pair = std::make_pair(seq1->ID(), seq2->ID());
        if (!candidate_map.ExistsOrInsert(abs_seq_pair)) {
          AllAllExecutor::WorkItem item =
              std::make_tuple(seq1, seq2, cluster.Sequences().size());
          executor->EnqueueAlignment(item);
        } else {
          num_avoided++;
        }
      }
    }
  }
  std::cout << "Avoided " << num_avoided << " alignments.\n";
}

void ClusterSet::DumpJson(const std::string& filename) const {
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
  std::ofstream o(filename);

  o << std::setw(2) << j << std::endl;
}

void ClusterSet::RemoveDuplicates() {
  absl::flat_hash_set<std::vector<size_t>> set_map;

  auto cluster_it = clusters_.begin();
  while (cluster_it != clusters_.end()) {
    std::vector<size_t> cluster_set;

    for (const auto& s : cluster_it->Sequences()) {
      cluster_set.push_back(s.ID());
    }
    // std::sort(cluster_set.begin(), cluster_set.end());

    auto result = set_map.insert(std::move(cluster_set));
    if (!result.second) {
      cluster_it = clusters_.erase(cluster_it);
    } else {
      cluster_it++;
    }
  }
}
