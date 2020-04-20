
#include "cluster_set.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "absl/container/flat_hash_set.h"
#include "aligner.h"
#include "candidate_map.h"
#include "debug.h"
#include "json.hpp"
#include "merge_executor.h"
#include "src/common/buffered_compressor.h"

using std::make_tuple;
using std::vector;

void cm_free_func(void* data, void* hint) {
  delete[] reinterpret_cast<char*>(data);
}

// type is implicit in the call, retained for future improvments
void ClusterSet::BuildMarshalledResponse(int id, RequestType type,
                                         MarshalledResponse* response) {
  // calculate buffer size for a single alloc
  size_t buf_size = sizeof(ResponseHeader) + sizeof(ClusterSetHeader);
  for (const auto& c : clusters_) {
    buf_size += sizeof(ClusterHeader) + c.Sequences().size() * sizeof(int);
  }
  agd::Buffer buf(buf_size);

  ResponseHeader rh;
  rh.id = id;
  rh.type = type;
  buf.AppendBuffer(reinterpret_cast<char*>(&rh), sizeof(ResponseHeader));

  ClusterSetHeader h;
  h.num_clusters = 0;  // set after once we know the value
  buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(ClusterSetHeader));

  for (const auto& c : clusters_) {
    // keep the fully merged for controller to remove
    ClusterHeader ch;
    ch.fully_merged = c.IsFullyMerged();
    ch.num_seqs = c.Sequences().size();
    buf.AppendBuffer(reinterpret_cast<char*>(&ch), sizeof(ClusterHeader));
    for (const auto& s : c.Sequences()) {
      uint32_t i = s;
      buf.AppendBuffer(reinterpret_cast<char*>(&i), sizeof(int));
    }
  }

  char* data = buf.mutable_data() + sizeof(ResponseHeader);
  ClusterSetHeader* hp = reinterpret_cast<ClusterSetHeader*>(data);
  hp->num_clusters = clusters_.size();
  // hand the buf pointer to the message
  response->msg =
      zmq::message_t(buf.release_raw(), buf.size(), cm_free_func, NULL);
}

// RequestType is implicit in the call
void ClusterSet::BuildMarshalledResponse(int id, int start_index, int end_index,
                                         int cluster_index,
                                         MarshalledResponse* response) {
  // calculate buffer size for a single alloc
  size_t buf_size =
      sizeof(ResponseHeader) + 3 * sizeof(int) + sizeof(ClusterSetHeader);
  for (const auto& c : clusters_) {
    buf_size += sizeof(ClusterHeader) + c.Sequences().size() * sizeof(int);
  }
  agd::Buffer buf(buf_size);
  ResponseHeader rh;
  rh.id = id;
  rh.type = RequestType::Partial;

  // buf contains response header - three integers - cluster set
  buf.AppendBuffer(reinterpret_cast<char*>(&rh), sizeof(ResponseHeader));
  buf.AppendBuffer(reinterpret_cast<char*>(&start_index), sizeof(int));
  buf.AppendBuffer(reinterpret_cast<char*>(&end_index), sizeof(int));
  buf.AppendBuffer(reinterpret_cast<char*>(&cluster_index), sizeof(int));

  ClusterSetHeader h;
  h.num_clusters = 0;  // set after once we know the value
  buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(ClusterSetHeader));

  for (const auto& c : clusters_) {
    // keep the fully merged for controller to remove
    ClusterHeader ch;
    ch.fully_merged = c.IsFullyMerged();
    ch.num_seqs = c.Sequences().size();
    buf.AppendBuffer(reinterpret_cast<char*>(&ch), sizeof(ClusterHeader));
    for (const auto& s : c.Sequences()) {
      uint32_t i = s;
      buf.AppendBuffer(reinterpret_cast<char*>(&i), sizeof(int));
    }
  }

  char* data = buf.mutable_data() + sizeof(ResponseHeader) + 3 * sizeof(int);
  ClusterSetHeader* hp = reinterpret_cast<ClusterSetHeader*>(data);
  hp->num_clusters = clusters_.size();
  // hand the buf pointer to the message
  response->msg =
      zmq::message_t(buf.release_raw(), buf.size(), cm_free_func, NULL);
}

ClusterSet::ClusterSet(MarshalledClusterSet& marshalled_set,
                       const std::vector<Sequence>& sequences) {
  // std::vector<Cluster> clusters(marshalled_set.NumClusters());
  MarshalledClusterView cluster;
  while (marshalled_set.NextCluster(&cluster)) {
    // for (size_t cs_i = 0; cs_i < set_proto.clusters_size(); cs_i++) {
    // const auto& cluster_proto = set_proto.clusters(cs_i);
    // std::cout << "cluster has " << cluster_proto.indexes_size() << " seqs\n";
    // std::cout << "adding cluster with rep idx " << cluster.SeqIndex(0) << "
    // and num seqs = " << cluster.NumSeqs() << "\n";
    Cluster c(cluster.SeqIndex(0), sequences);
    uint32_t num_seqs = cluster.NumSeqs();
    for (size_t seq_i = 1; seq_i < num_seqs; seq_i++) {
      c.AddSequence(cluster.SeqIndex(seq_i));
    }
    if (cluster.IsFullyMerged()) {
      c.SetFullyMerged();
    }
    clusters_.push_back(std::move(c));
  }
}

ClusterSet::ClusterSet(MarshalledClusterSetView& marshalled_set,
                       const vector<size_t>& set_offsets, int start_index,
                       int end_index, const std::vector<Sequence>& sequences) {
  MarshalledClusterView cluster;
  clusters_.reserve(end_index - start_index);
  for (int i = start_index; i <= end_index; i++) {
    marshalled_set.ClusterAtOffset(&cluster, set_offsets[i]);
    Cluster c(cluster.SeqIndex(0), sequences);
    uint32_t num_seqs = cluster.NumSeqs();
    c.Reserve(num_seqs);
    for (size_t seq_i = 1; seq_i < num_seqs; seq_i++) {
      c.AddSequence(cluster.SeqIndex(seq_i));
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
  // std::vector<Cluster> clusters(marshalled_set.NumClusters());
  MarshalledClusterView cluster;
  while (marshalled_set.NextCluster(&cluster)) {
    // for (size_t cs_i = 0; cs_i < set_proto.clusters_size(); cs_i++) {
    // const auto& cluster_proto = set_proto.clusters(cs_i);
    // std::cout << "cluster has " << cluster_proto.indexes_size() << " seqs\n";
    // std::cout << "adding cluster with rep idx " << cluster.SeqIndex(0) << "
    // and num seqs = " << cluster.NumSeqs() << "\n";
    Cluster c(cluster.SeqIndex(0), sequences);
    uint32_t num_seqs = cluster.NumSeqs();
    for (size_t seq_i = 1; seq_i < num_seqs; seq_i++) {
      c.AddSequence(cluster.SeqIndex(seq_i));
    }
    if (cluster.IsFullyMerged()) {
      c.SetFullyMerged();
    }
    clusters_.push_back(std::move(c));
  }
}

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
              return a.SeqRep().Seq().size() > b.SeqRep().Seq().size();
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

      auto c_num_uncovered = cluster->SeqRep().Seq().size() -
                             (alignment.seq1_max - alignment.seq1_min);
      auto c_other_num_uncovered = c_other.SeqRep().Seq().size() -
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
        break;
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
            c.SeqRep().Seq().size() - (alignment.seq1_max - alignment.seq1_min);
        auto c_other_num_uncovered = c_other.SeqRep().Seq().size() -
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

void ClusterSet::DebugDump(const std::vector<Sequence>& sequences) const {
  std::cout << "Dumping " << clusters_.size() << " clusters in set... \n";
  for (const auto& cluster : clusters_) {
    std::cout << "\tCluster seqs:\n";
    for (const auto& seq : cluster.Sequences()) {
      std::cout << "\t\tGenome: " << sequences[seq].Genome() << ", sequence: "
                << PrintNormalizedProtein(sequences[seq].Seq().data(),
                                          sequences[seq].Seq().length())
                << "\n\n";
    }
  }
}

ClusterSet ClusterSet::MergeCluster(Cluster& c_other, ProteinAligner* aligner,
                                    bool& worker_signal_) {
  // for the dist version, we keep fully merged clusters around,
  // because this is a partial merge of two large sets,
  // the results of which need to be merged by the controller
  ClusterSet new_cluster_set(clusters_.size() + 1);

  ProteinAligner::Alignment alignment;
  agd::Status s;
  bool fully_merged = false;
  int c_other_old_nseqs = c_other.Sequences().size();
  for (auto& c : clusters_) {
    if (worker_signal_) {
      std::cout << "Breaking partial merge.\n";
      break;
    }
    Cluster c_standin(
        c_other.AllSequences());  // stand in cluster which holds diffs
    if (!fully_merged && c.PassesThreshold(c_other, aligner)) {
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
          c.SeqRep().Seq().size() - (alignment.seq1_max - alignment.seq1_min);
      auto c_other_num_uncovered = c_other.SeqRep().Seq().size() -
                                   (alignment.seq2_max - alignment.seq2_min);
      if (c_num_uncovered < aligner->Params()->max_n_aa_not_covered &&
          alignment.score > aligner->Params()->min_full_merge_score) {
        // they are _almost_ overlapped, merge completely
        // std::cout << "Nearly complete overlap, merging c into c_other,
        // score is " << alignment.score << "\n";
        for (const auto& seq : c.Sequences()) {
          c_other.AddSequence(seq);
        }
        c_standin.SetFullyMerged();
        fully_merged = true;
      } else if (c_other_num_uncovered <
                     aligner->Params()->max_n_aa_not_covered &&
                 alignment.score > aligner->Params()->min_full_merge_score) {
        // std::cout << "Nearly complete overlap, merging c_other into c,
        // score is " << alignment.score << "\n";
        for (const auto& seq : c_other.Sequences()) {
          c_standin.AddSequence(seq);
        }
        c_other.SetFullyMerged();
        fully_merged = true;
      } else {
        // add c_other_rep into c
        // for each sequence in c_other, add if it matches c rep
        // keep both clusters
        // std::cout << "merging and keeping both clusters\n";
        int num_old_seqs = c.Sequences().size();
        c.Merge(&c_other, aligner);
        auto seqs = c.Sequences();
        auto it = seqs.begin();
        // push only newly added sequences
        std::advance(it, num_old_seqs);
        while (it != seqs.end()) {
          c_standin.AddSequence(*it);
          it++;
        }
      }
    }  // if passes threshold

    new_cluster_set.clusters_.push_back(std::move(c_standin));
  }

  // we can leave out without confusing the controller
  // adding only diffs
  if (!c_other.IsFullyMerged()) {
    Cluster c_standin(c_other.AllSequences());
    auto seqs = c_other.Sequences();
    auto it = seqs.begin();
    // push only newly added sequences
    std::advance(it, c_other_old_nseqs);
    while (it != seqs.end()) {
      c_standin.AddSequence(*it);
      it++;
    }
    new_cluster_set.clusters_.push_back(std::move(c_standin));
  }

  return new_cluster_set;
}

void ClusterSet::ScheduleAlignments(AllAllBase* executor,
                                    std::vector<Sequence>& sequences) {
  // removing duplicate clusters (clusters with same sequences)
  // for some reason, absl::InlinedVector doesnt work here
  absl::flat_hash_set<std::vector<size_t>> set_map;

  size_t num_dups_found = 0;
  for (auto& c : clusters_) {
    std::vector<size_t> cluster_set;
    for (const auto& s : c.Sequences()) {
      cluster_set.push_back(s);
    }
    std::sort(cluster_set.begin(), cluster_set.end());

    auto result = set_map.insert(std::move(cluster_set));
    if (!result.second) {
      c.SetDuplicate();
      num_dups_found++;
    }
  }

  std::cout << "Found " << num_dups_found << " duplicate clusters."
            << std::endl;
  // sort by residue total first
  // to schedule the heaviest computations first
  std::cout << "sorting clusters ...\n";
  std::sort(clusters_.begin(), clusters_.end(), [](Cluster& a, Cluster& b) {
    return a.Sequences().size() > b.Sequences().size();
  });
  std::cout << "done sorting clusters." << std::endl;

  CandidateMap candidate_map(40000000);  // only a few MB
  int num_avoided = 0, num_scheduled = 0;
  auto time_last_candmap_status = std::chrono::high_resolution_clock::now();

  for (const auto& cluster : clusters_) {
    // std::cout << "Cluster has " << cluster.Sequences().size() << " seqs\n";
    if (cluster.IsDuplicate()) {
      continue;
    }
    for (auto it = cluster.Sequences().begin(); it != cluster.Sequences().end();
         it++) {
      for (auto itt = next(it); itt != cluster.Sequences().end(); itt++) {
        auto seq1 = *it;
        auto seq2 = *itt;

        if (sequences[seq1].Genome() == sequences[seq2].Genome() &&
            sequences[seq1].GenomeIndex() == sequences[seq2].GenomeIndex()) {
          // not sure if this can actually happen yet, but no need to align
          // against self
          continue;
        }

        if (sequences[seq1].GenomeSize() > sequences[seq2].GenomeSize() ||
            ((sequences[seq1].GenomeSize() == sequences[seq2].GenomeSize()) &&
             sequences[seq1].Genome() > sequences[seq2].Genome())) {
          std::swap(seq1, seq2);
        }

        if (sequences[seq1].Genome() == sequences[seq2].Genome() &&
            sequences[seq1].GenomeIndex() > sequences[seq2].GenomeIndex()) {
          std::swap(seq1, seq2);
        }

        auto abs_seq_pair = std::make_pair(seq1, seq2);
        if (!candidate_map.ExistsOrInsert(abs_seq_pair)) {
          AllAllExecutor::WorkItem item = std::make_tuple(
              &sequences[seq1], &sequences[seq2], cluster.Sequences().size());
          //std::cout << "enqueueing alignment between " << sequences[seq1].ID() << 
          // " and " << sequences[seq2].ID() << std::endl;
          executor->EnqueueAlignment(item);
          num_scheduled++;
        } else {
          num_avoided++;
        }
        auto cur_time = std::chrono::high_resolution_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - time_last_candmap_status).count() > 60000){
          time_t now_time = std::time(0);
          std::cout << "[" << std::put_time(std::localtime(&now_time), "%F %T") << "]" 
                    << " current candidate map size: " << candidate_map.size() 
                    << " scheduled/avoided pairs: " << num_scheduled << "/" 
                    << num_avoided << std::endl;
          time_last_candmap_status = cur_time;
        }
      }
    }
  }
  std::cout << "Avoided " << num_avoided << " / Scheduled " << num_scheduled << " alignments." << std::endl;
}

#define RETURN_ON_ERROR(...) \
  {int e = (__VA_ARGS__);     \
  if (e) return e;}

int ClusterSet::DumpJson(const std::string& filename,
                         std::vector<std::string>& dataset_file_names) const {
  std::cout << "dumping clusters ...\n";

  // total amount of clusters can be large. So we buffer and compress data to
  // file as we go, with default 2 MB buffers
  BufferedCompressor compressor(2 * 1024 * 1024);
  compressor.Init(filename);

  std::stringstream ss;

  std::string dsets("{\n\"datasets\": [\n");
  RETURN_ON_ERROR(compressor.Write(dsets.data(), dsets.size()));

  for (auto n_it = dataset_file_names.begin(); n_it != dataset_file_names.end();
       n_it++) {
    RETURN_ON_ERROR(compressor.Write("\t\"", 2));
    RETURN_ON_ERROR(compressor.Write(n_it->data(), n_it->size()));
    if (std::next(n_it) == dataset_file_names.end()) {
      RETURN_ON_ERROR(compressor.Write("\"\n],\n", 5));
    } else {
      RETURN_ON_ERROR(compressor.Write("\",\n", 3));
    }
  }

  std::string clstrs("\"clusters\": [\n");
  RETURN_ON_ERROR(compressor.Write(clstrs.data(), clstrs.size()));

  // im using tabs/space here so it's at least somewhat human readable when
  // decompressed
  for (auto c_it = clusters_.begin(); c_it != clusters_.end(); c_it++) {
    RETURN_ON_ERROR(compressor.Write("[\n", 2));

    for (auto s_it = c_it->Sequences().begin(); s_it != c_it->Sequences().end();
         s_it++) {
      const auto& s = *s_it;

      ss << "{\n\t\"Genome\": \"" << c_it->SeqAt(s).Genome()
         << "\",\n\t\"Index\": " << c_it->SeqAt(s).GenomeIndex()
         << ",\n\t\"AbsoluteIndex\": " << s << "\n}";
      if (std::next(s_it) != c_it->Sequences().end()) {
        ss << ",";
      }
      ss << "\n";
      std::string out = ss.str();
      RETURN_ON_ERROR(compressor.Write(out.data(), out.size()));
      ss.str(std::string());  // clear the stream buffer
    }

    if (std::next(c_it) == clusters_.end()) {
      RETURN_ON_ERROR(compressor.Write("]\n", 2));
    } else {
      RETURN_ON_ERROR(compressor.Write("],\n", 3));
    }
  }
  RETURN_ON_ERROR(compressor.Write("]}\n", 3));
  return 0;
}

void ClusterSet::RemoveDuplicates() {
  absl::flat_hash_set<std::vector<size_t>> set_map;

  auto cluster_it = clusters_.begin();
  while (cluster_it != clusters_.end()) {
    std::vector<size_t> cluster_set;

    for (const auto& s : cluster_it->Sequences()) {
      cluster_set.push_back(s);
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
