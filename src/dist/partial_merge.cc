
#include "src/dist/partial_merge.h"

/*void PartialMergeSet::RemoveFullyMerged() {
  auto cluster_it = clusters_.begin();
  while (cluster_it != clusters_.end()) {
    if (cluster_it->IsFullyMerged()) {
      cluster_it = clusters_.erase(cluster_it);
      continue;
    }
    cluster_it++;
  }
  
  auto new_cluster_it = new_clusters_.begin();
  while (new_cluster_it != new_clusters_.end()) {
    if (new_cluster_it->fully_merged()) {
      new_cluster_it = new_clusters_.erase(new_cluster_it);
      continue;
    }
    new_cluster_it++;
  }

}*/

void PartialMergeSet::Init(MarshalledClusterSet& set) {
  MarshalledClusterView cluster;
  while (set.NextCluster(&cluster)) {
    IndexedCluster ic(cluster);
    clusters_.push_back(std::move(ic));
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
    for (const auto& c : clusters_) {
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
      for (auto s : seq_set) {
        if (s != r) {
          set->buf.AppendBuffer(reinterpret_cast<char*>(&s), sizeof(uint32_t));
        }
      }
    }

    // dont forget the new clusters, they are already marshalled
    for (auto& c : new_clusters_) {
      set->buf.AppendBuffer(c.buf.data(), c.buf.size());
      total_not_merged++;
    }

    char* data = set->buf.mutable_data();
    ClusterSetHeader* hp = reinterpret_cast<ClusterSetHeader*>(data);
    hp->num_clusters = total_not_merged;
}

void PartialMergeSet::MergeClusterSet(MarshalledClusterSetView set) {
  MarshalledClusterView cluster;
  assert(set.NumClusters() >= clusters_.size());
  for (uint32_t i = 0; i < clusters_.size(); i++) {
    //auto cluster = set.Cluster(i);
    if (!set.NextCluster(&cluster)) {
      std::cout << "error next cluster returned false? on index " << i << " with clusters size " << set.NumClusters() 
        << " and orig clusters size " << clusters_.size() << "\n";
      exit(0);
    }

    if (cluster.IsFullyMerged()) {
      clusters_[i].SetFullyMerged();
    }
    auto num_seqs = cluster.NumSeqs();
    // only do new seqs, orig seqs are already present and can be skipped
    auto orig_seqs = clusters_[i].NumOrigSeqs();
    if (orig_seqs > num_seqs) { 
      std::cout << "error num seqs " << num_seqs << " < " << orig_seqs << " orig seqs, cluster " << i << "\n"; 
      exit(0); 
    }
    for (auto x = orig_seqs; x < num_seqs; x++) {
      clusters_[i].Insert(cluster.SeqIndex(x));
    }
  }

  uint32_t num_clusters = set.NumClusters();
  // insert the new cluster, if there is one
  if (num_clusters > clusters_.size()) {
    assert(clusters_.size() == num_clusters - 1);
    //IndexedCluster c;
    //const auto& cluster = set.Cluster(num_clusters - 1);
    if (!set.NextCluster(&cluster)) {
      std::cout << "error next cluster returned false? second\n";
      exit(0);
    }
    // copy out the marshalled cluster into the new clusters array
    MarshalledCluster mc(cluster.data, cluster.TotalSize());

    absl::MutexLock l(&mu_);
    new_clusters_.push_back(std::move(mc));
  }

}
