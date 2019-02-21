
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

void PartialMergeSet::Init(const MarshalledClusterSet& set) {
  int num_clusters = set.NumClusters();
  for (int i = 0; i < num_clusters; i++) {
    auto c = set.Cluster(i);
    IndexedCluster ic(c);
    clusters_.push_back(std::move(ic));
  }
}

void PartialMergeSet::BuildMarshalledSet(MarshalledClusterSet* set) {
    ClusterSetHeader h;
    h.num_clusters = 0;  // set after once we know the value
    set->buf.reset();
    set->buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(ClusterSetHeader));

    uint32_t total_not_merged = 0;
    for (const auto& c : clusters_) {
      if (c.IsFullyMerged()) {
        continue;
      }
      total_not_merged++;
      ClusterHeader ch;
      ch.fully_merged = false;
      ch.num_seqs = c.SeqIndexes().size();
      set->buf.AppendBuffer(reinterpret_cast<char*>(&ch), sizeof(ClusterHeader));
      int r = c.Representative();
      set->buf.AppendBuffer(reinterpret_cast<char*>(&r), sizeof(uint32_t));
      for (auto s : c.SeqIndexes()) {
        if (s != r) {
          set->buf.AppendBuffer(reinterpret_cast<char*>(&s), sizeof(decltype(s)));
        }
      }
    }

    char* data = set->buf.mutable_data();
    ClusterSetHeader* hp = reinterpret_cast<ClusterSetHeader*>(data);
    hp->num_clusters = total_not_merged;
}

void PartialMergeSet::MergeClusterSet(MarshalledClusterSetView set) {
  MarshalledClusterView cluster;
  for (int i = 0; i < clusters_.size(); i++) {
    //auto cluster = set.Cluster(i);
    if (!set.NextCluster(&cluster)) {
      std::cout << "error next cluster returned false?\n";
      exit(0);
    }

    if (cluster.IsFullyMerged()) {
      clusters_[i].SetFullyMerged();
    }
    auto num_seqs = cluster.NumSeqs();
    for (auto x = 0; x < num_seqs; x++) {
      clusters_[i].Insert(cluster.SeqIndex(x));
    }
  }

  uint32_t num_clusters = set.NumClusters();
  // insert the new cluster, if there is one
  if (num_clusters > clusters_.size()) {
    assert(clusters_.size() == num_clusters - 1);
    IndexedCluster c;
    //const auto& cluster = set.Cluster(num_clusters - 1);
    if (!set.NextCluster(&cluster)) {
      std::cout << "error next cluster returned false?\n";
      exit(0);
    }
    // copy out the marshalled cluster into the new clusters array
    MarshalledCluster mc(cluster.data, cluster.TotalSize());

    absl::MutexLock l(&mu_);
    new_clusters_.push_back(std::move(mc));
  }

}
