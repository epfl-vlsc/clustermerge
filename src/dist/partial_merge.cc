
#include "src/dist/partial_merge.h"

void PartialMergeSet::Init(MarshalledClusterSet& set1, MarshalledClusterSet& set2) {
  MarshalledClusterView cluster;
  while (set1.NextCluster(&cluster)) {
    IndexedCluster ic(cluster);
    clusters_set1_.push_back(std::move(ic));
  }
  while (set2.NextCluster(&cluster)) {
    IndexedCluster ic(cluster);
    clusters_set2_.push_back(std::move(ic));
  }
}

void PartialMergeSet::BuildMarshalledSet(MarshalledClusterSet* set) {
  ClusterSetHeader h;
  h.num_clusters = 0;  // set after once we know the value
  set->buf.reset();
  set->buf.reserve(512);
  set->buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(ClusterSetHeader));

  uint32_t total_not_merged = 0;
  absl::flat_hash_set<uint32_t> seq_set;
  
  for (const auto& c : clusters_set1_) {
    if (c.IsFullyMerged()) {
      continue;
    }
    total_not_merged++;

    seq_set.clear();
    for (auto s : c.SeqIndexes()) {
      seq_set.insert(s);
    }
    c.AddNewSeqs(&seq_set);

    ClusterHeader ch;
    ch.fully_merged = false;
    ch.num_seqs = seq_set.size();
    set->buf.AppendBuffer(reinterpret_cast<char*>(&ch), sizeof(ClusterHeader));
    auto r = c.Representative();
    set->buf.AppendBuffer(reinterpret_cast<char*>(&r), sizeof(uint32_t));
    // remove duplicates in new seqs using the set
    int i = 1;
    for (auto s : seq_set) {
      if (s != r) {
        i++;
        set->buf.AppendBuffer(reinterpret_cast<char*>(&s), sizeof(uint32_t));
      }
    }
    assert(i == seq_set.size());
  }

  for (const auto& c : clusters_set2_) {
    if (c.IsFullyMerged()) {
      continue;
    }
    total_not_merged++;

    seq_set.clear();
    for (auto s : c.SeqIndexes()) {
      seq_set.insert(s);
    }
    c.AddNewSeqs(&seq_set);

    ClusterHeader ch;
    ch.fully_merged = false;
    ch.num_seqs = seq_set.size();
    set->buf.AppendBuffer(reinterpret_cast<char*>(&ch), sizeof(ClusterHeader));
    auto r = c.Representative();
    set->buf.AppendBuffer(reinterpret_cast<char*>(&r), sizeof(uint32_t));
    // remove duplicates in new seqs using the set
    int i = 1;
    for (auto s : seq_set) {
      if (s != r) {
        i++;
        set->buf.AppendBuffer(reinterpret_cast<char*>(&s), sizeof(uint32_t));
      }
    }
    assert(i == seq_set.size());
  }

  char* data = set->buf.mutable_data();
  ClusterSetHeader* hp = reinterpret_cast<ClusterSetHeader*>(data);
  hp->num_clusters = total_not_merged;
}

void PartialMergeSet::MergeClusterSet(MarshalledClusterSetView set, int start_index, 
  int end_index, int cluster_index) {
  MarshalledClusterView cluster;
  uint32_t clusters_in_chunk = end_index - start_index + 1;
  assert(set.NumClusters() >= clusters_in_chunk);
  for (int i = start_index; i <= end_index; i++) {
    // auto cluster = set.Cluster(i);
    if (!set.NextCluster(&cluster)) {
      std::cout << "error next cluster returned false? on index " << i
                << " with clusters size " << set.NumClusters()
                << " and orig clusters size " << clusters_in_chunk << "\n";
      exit(0);
    }

    if (cluster.IsFullyMerged()) {
      clusters_set2_[i].SetFullyMerged();
    } else {
      // cluster has only diffs, push them directly
      for (uint32_t x = 0; x < cluster.NumSeqs(); x++) {
        clusters_set2_[i].Insert(cluster.SeqIndex(x));
      }
    }
  }

  uint32_t num_clusters = set.NumClusters();
  // insert the new cluster, if there is one
  if (num_clusters > clusters_in_chunk) {
    assert(clusters_in_chunk == num_clusters - 1);
    // IndexedCluster c;
    // const auto& cluster = set.Cluster(num_clusters - 1);
    if (!set.NextCluster(&cluster)) {
      std::cout << "error next cluster returned false? second\n";
      exit(0);
    }

    for (uint32_t x = 0; x < cluster.NumSeqs(); x++) {
      clusters_set1_[cluster_index].Insert(cluster.SeqIndex(x));
    }
  } else {
      clusters_set1_[cluster_index].SetFullyMerged();
  }
}
