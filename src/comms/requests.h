#pragma once

#include "src/agd/buffer.h"
#include "zmq.hpp"
#include <iostream>
#include <tuple>

enum RequestType { Batch = 0, Partial, LargePartial, SubLargePartial };

struct __attribute__((packed)) BatchRequestHeader {
  int id;
  RequestType type;
  int num_cluster_sets;
};

struct __attribute__((packed)) PartialRequestHeader {
  int id;
  RequestType type;
  // cluster
  // set
};

struct __attribute__((packed)) ClusterSetHeader {
  uint32_t num_clusters;
};

struct __attribute__((packed)) ClusterHeader {
  char fully_merged;
  uint32_t num_seqs;
};

struct __attribute__((packed)) ResponseHeader {
  int id;
  RequestType type;
};

struct MarshalledClusterView {
  MarshalledClusterView() = default;
  MarshalledClusterView(const char* d) : data(d) {}
  const char* data;
  uint32_t NumSeqs() const {
    const ClusterHeader* h = reinterpret_cast<const ClusterHeader*>(data);
    return h->num_seqs;
  }
  size_t TotalSize() const {
    uint32_t seqs = NumSeqs();
    return sizeof(ClusterHeader) + sizeof(uint32_t) * seqs;
  }
  // get sequence index at idx
  uint32_t SeqIndex(uint32_t idx) const {
    const char* p = data + sizeof(ClusterHeader) + sizeof(uint32_t) * idx;
    return *reinterpret_cast<const uint32_t*>(p);
  }
  bool IsFullyMerged() const {
    const ClusterHeader* h = reinterpret_cast<const ClusterHeader*>(data);
    return h->fully_merged;
  }
};

// for sets to merge
struct MarshalledClusterSetView {
  MarshalledClusterSetView() = default;
  MarshalledClusterSetView(const char* d) : data(d) {}
  const char* data;
  const char* cur_cluster_ptr = nullptr;
  uint32_t cur_cluster_idx = 0;

  void Reset() { cur_cluster_ptr = nullptr; }

  uint32_t NumClusters() const {
    const ClusterSetHeader* h = reinterpret_cast<const ClusterSetHeader*>(data);
    return h->num_clusters;
  }
  bool NextCluster(MarshalledClusterView* cluster) {
    if (cur_cluster_ptr == nullptr) {
      cur_cluster_ptr = data + sizeof(ClusterSetHeader);
      *cluster = MarshalledClusterView(cur_cluster_ptr);
      cur_cluster_idx = 0;
      return true;
    }
    // advance to next cluster and return
    uint32_t cluster_size =
        reinterpret_cast<const ClusterHeader*>(cur_cluster_ptr)->num_seqs;
    cur_cluster_ptr += sizeof(ClusterHeader) + sizeof(uint32_t) * cluster_size;
    cur_cluster_idx++;
    if (cur_cluster_idx >= NumClusters()) {
      return false;
    } else {
      *cluster = MarshalledClusterView(cur_cluster_ptr);
      return true;
    }
  }

  void Offsets(std::vector<size_t>& offsets)  {
    // for cluster at index zero
    size_t offset = sizeof(ClusterSetHeader);
    offsets.push_back(offset);
    const char* cluster_ptr = data + offset;
    
    uint32_t cluster_index = 1;
    while(cluster_index < NumClusters())  {
      uint32_t cluster_size = 
        reinterpret_cast<const ClusterHeader*>(cluster_ptr)->num_seqs;
      offset += sizeof(ClusterHeader) + sizeof(uint32_t) * cluster_size;   
      offsets.push_back(offset);
      cluster_ptr = data + offset;
      cluster_index++;
    }
  } 

  // no check on offset, assumed correct
  bool ClusterAtOffset(MarshalledClusterView* cluster, size_t offset) {
    const char* cluster_ptr = data + offset;
    *cluster = MarshalledClusterView(cluster_ptr);
    return true;
  }
};

// for response queue, keep msg to prevent copies
struct MarshalledResponse {
  // contains marshalled id and cluster set
  zmq::message_t msg;
  MarshalledResponse() = default;
  MarshalledResponse(MarshalledResponse&&) = default;
  MarshalledResponse& operator=(MarshalledResponse&&) = default;
  int ID() { 
    const ResponseHeader* rh = reinterpret_cast<const ResponseHeader*>(msg.data());
    return rh->id;
  }
  RequestType Type()  {
    const ResponseHeader* rh = reinterpret_cast<const ResponseHeader*>(msg.data());
    return rh->type;
  }

  MarshalledClusterSetView Set() {
    const ResponseHeader* rh = reinterpret_cast<const ResponseHeader*>(msg.data());
    const char* data;
    if(rh->type == RequestType::SubLargePartial)  {
      data = reinterpret_cast<const char*>(msg.data()) + sizeof(ResponseHeader) + 3*sizeof(int);
    } else {
      data = reinterpret_cast<const char*>(msg.data()) + sizeof(ResponseHeader);
    }
    return MarshalledClusterSetView(data);
  }
  
  std::tuple<int, int, int> Indexes() {
    std::tuple<int, int, int> indexes;
    if(this->Type() != RequestType::SubLargePartial)  {
      std::cout << "Function called on incorrect reponse type.\n";
      exit(0);
    }  else {
      const char* data = reinterpret_cast<const char*>(msg.data()) + sizeof(ResponseHeader);
      int start_index = *reinterpret_cast<const int*>(data);
      int end_index = *reinterpret_cast<const int*>(data + sizeof(int));
      int cluster_index = *reinterpret_cast<const int*>(data + 2*sizeof(int));
      indexes = std::make_tuple(start_index, end_index, cluster_index);
    }
    return indexes;
  }
};

struct MarshalledClusterSet {
  MarshalledClusterSet() : buf(agd::Buffer(256, 128)) {}
  MarshalledClusterSet(uint32_t idx)
      : buf(agd::Buffer(
            sizeof(ClusterSetHeader) + sizeof(ClusterHeader) + sizeof(uint32_t),
            512)) {  // create from single sequence

    ClusterSetHeader h;
    h.num_clusters = 1;
    buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(ClusterSetHeader));

    ClusterHeader ch;
    ch.fully_merged = false;
    ch.num_seqs = 1;
    buf.AppendBuffer(reinterpret_cast<char*>(&ch), sizeof(ClusterHeader));

    buf.AppendBuffer(reinterpret_cast<char*>(&idx), sizeof(uint32_t));
  }

  //only called for batch requests, assumes user has checked response type
  MarshalledClusterSet(MarshalledResponse& response)
      : buf(agd::Buffer(response.msg.size() - sizeof(ResponseHeader), 128)) {
    // a response buffer is an int id and a marshalled clusterset
    const char* data = reinterpret_cast<const char*>(response.msg.data());
    data += sizeof(ResponseHeader);
    auto data_size = response.msg.size() - sizeof(ResponseHeader);
    // buf.reserve(data_size);
    buf.AppendBuffer(data, data_size);
  }

  void SortSet()  {
    std::vector<std::pair<uint32_t, size_t>> clusters;
    uint32_t num_clusters = 0, total_clusters = NumClusters();
    size_t offset = sizeof(ClusterSetHeader);
    const char* cur_cluster_ptr = buf.data() + offset;
    
    while(num_clusters < total_clusters) {
      uint32_t cluster_size = 
        reinterpret_cast<const ClusterHeader*>(cur_cluster_ptr)->num_seqs;
      clusters.push_back({cluster_size, offset});
      offset += sizeof(ClusterHeader) + sizeof(uint32_t) * cluster_size;
      cur_cluster_ptr = buf.data() + offset;
      num_clusters++;  
    }

    // std::sort(clusters.begin(), clusters.end(), 
    //   [](const std::pair<uint32_t, size_t>& p1, const std::pair<uint32_t, size_t>& p2) {
    //     return p1.first >= p2.first;
    // });

    std::sort(clusters.begin(), clusters.end(), std::greater<std::pair<uint32_t, size_t>>());
    
    agd::Buffer new_buf(buf.size());
    const char* data = buf.data();
    new_buf.AppendBuffer(data, sizeof(ClusterSetHeader));
    for(const auto& pr: clusters)  {
      data = buf.data() + pr.second;
      new_buf.AppendBuffer(data, sizeof(ClusterHeader) + pr.first * sizeof(uint32_t));  
    }
    assert(new_buf.size() == buf.size());
    buf = std::move(new_buf);
    Reset();  
  }

  uint32_t NumClusters() const {
    const char* data = buf.data();
    const ClusterSetHeader* h = reinterpret_cast<const ClusterSetHeader*>(data);
    return h->num_clusters;
  }

  bool NextCluster(MarshalledClusterView* cluster) {
    if (cur_cluster_ptr == nullptr) {
      cur_cluster_ptr = buf.data() + sizeof(ClusterSetHeader);
      *cluster = MarshalledClusterView(cur_cluster_ptr);
      return true;
    }
    // advance to next cluster and return
    uint32_t cluster_size =
        reinterpret_cast<const ClusterHeader*>(cur_cluster_ptr)->num_seqs;
    cur_cluster_ptr += sizeof(ClusterHeader) + sizeof(uint32_t) * cluster_size;
    if (cur_cluster_ptr >= buf.data() + buf.size()) {
      return false;
    } else {
      *cluster = MarshalledClusterView(cur_cluster_ptr);
      return true;
    }
  }

  size_t TotalSize() const { return buf.size(); }

  void Reset() { cur_cluster_ptr = nullptr; }

  const char* cur_cluster_ptr = nullptr;
  agd::Buffer buf;
};

struct MarshalledCluster {
  MarshalledCluster(const char* data, uint32_t size)
      : buf(agd::Buffer(size, 128)) {
    buf.AppendBuffer(data, size);
  }
  agd::Buffer buf;
};

// for request queue
struct MarshalledRequest {
  MarshalledRequest() : buf(agd::Buffer(256, 128)) {}

  void CreateBatchRequest(int id) {
    BatchRequestHeader h;
    h.id = id;
    h.type = RequestType::Batch;
    h.num_cluster_sets = 0;
    buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(BatchRequestHeader));
  }

  void AddSetToBatch(const MarshalledClusterSet& set) {
    buf.AppendBuffer(set.buf.data(), set.buf.size());
    BatchRequestHeader* h =
        reinterpret_cast<BatchRequestHeader*>(buf.mutable_data());
    h->num_cluster_sets++;
  }

  void CreatePartialRequest(int id, MarshalledClusterView cluster,
                            const MarshalledClusterSet& set) {
    size_t total = cluster.TotalSize() + set.TotalSize();
    buf.reserve(total + sizeof(PartialRequestHeader));
    PartialRequestHeader h;
    h.id = id;
    h.type = RequestType::Partial;
    buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(PartialRequestHeader));
    buf.AppendBuffer(cluster.data, cluster.TotalSize());
    buf.AppendBuffer(set.buf.data(), set.buf.size());
    assert(buf.size() == total + sizeof(PartialRequestHeader));
  }

  // Does not have a clusterset, but is still large :P
  void CreateLargePartialRequest(int id, MarshalledClusterView cluster) {
    buf.reserve(cluster.TotalSize() + sizeof(PartialRequestHeader));
    PartialRequestHeader h;
    h.id = id;
    h.type = RequestType::LargePartial;
    buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(PartialRequestHeader));
    buf.AppendBuffer(cluster.data, cluster.TotalSize());
    assert(buf.size() == cluster.TotalSize() + sizeof(PartialRequestHeader));
  }

  void CreateSubLargePartialRequest(int id, MarshalledClusterView cluster, 
    int start_index, int end_index, int cluster_index)  {
    PartialRequestHeader h;
    h.id = id;
    h.type = RequestType::SubLargePartial;
    buf.reserve(sizeof(PartialRequestHeader) + sizeof(int)*3 + cluster.TotalSize());
    buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(PartialRequestHeader));
    buf.AppendBuffer(reinterpret_cast<char*>(&start_index), sizeof(int));
    buf.AppendBuffer(reinterpret_cast<char*>(&end_index), sizeof(int));
    buf.AppendBuffer(reinterpret_cast<char*>(&cluster_index), sizeof(int));
    buf.AppendBuffer(cluster.data, cluster.TotalSize());
    assert(buf.size() == cluster.TotalSize() + sizeof(PartialRequestHeader) + 3*sizeof(int));  
  }

  agd::Buffer buf;
};

struct MarshalledRequestView {
  MarshalledRequestView(const char* d, size_t sz) : data(d), size(sz) {}
  const char* data;
  size_t size;
  const char* cur_clusterset_ptr = nullptr;

  void Reset() { cur_clusterset_ptr = nullptr; }

  RequestType Type() {
    return reinterpret_cast<const PartialRequestHeader*>(data)->type;
  }

  int ID() { return reinterpret_cast<const PartialRequestHeader*>(data)->id; }

  uint32_t NumClusterSets() {
    // relies on user to call Type and ensure this is correct
    return reinterpret_cast<const BatchRequestHeader*>(data)->num_cluster_sets;
  }

  // for partial requests
  void ClusterAndSet(MarshalledClusterSetView* cluster_set,
                     MarshalledClusterView* cluster) {
    const char* cluster_ptr = data + sizeof(PartialRequestHeader);
    const ClusterHeader* h =
        reinterpret_cast<const ClusterHeader*>(cluster_ptr);
    uint32_t sz = h->num_seqs * sizeof(uint32_t) + sizeof(ClusterHeader);
    *cluster = MarshalledClusterView(cluster_ptr);

    const char* cluster_set_ptr = cluster_ptr + sz;
    *cluster_set = MarshalledClusterSetView(cluster_set_ptr);
  }

  // for Large Partial requests
  void Cluster(MarshalledClusterView* cluster) {
    const char* cluster_ptr = data + sizeof(PartialRequestHeader);
    *cluster = MarshalledClusterView(cluster_ptr);
  }

  // relies on user to call Type and ensure this is correct
  void IndexesAndCluster(int* start_index, int* end_index, int* cluster_index, MarshalledClusterView* cluster) {
    const char* ptr = data + sizeof(PartialRequestHeader);
    *start_index = *reinterpret_cast<const int*>(ptr);
    *end_index = *reinterpret_cast<const int*>(ptr + sizeof(int));
    *cluster_index = *reinterpret_cast<const int*>(ptr + 2*sizeof(int));
    const char* cluster_ptr = ptr + 3*sizeof(int);
    *cluster = MarshalledClusterView(cluster_ptr);
  }

  // for batch
  bool NextClusterSet(MarshalledClusterSetView* cluster_set) {
    if (cur_clusterset_ptr == nullptr) {
      cur_clusterset_ptr = data + sizeof(BatchRequestHeader);
      *cluster_set = MarshalledClusterSetView(cur_clusterset_ptr);
      return true;
    } else {
      // find the start of next cluster set
      const char* p = cur_clusterset_ptr;
      uint32_t num_clusters =
          reinterpret_cast<const ClusterSetHeader*>(p)->num_clusters;
      p += sizeof(ClusterSetHeader);
      for (uint32_t i = 0; i < num_clusters; i++) {
        uint32_t num_seqs = reinterpret_cast<const ClusterHeader*>(p)->num_seqs;
        p += sizeof(ClusterHeader) + num_seqs * sizeof(uint32_t);
      }
      if (p >= data + size) {
        return false;
      } else {
        cur_clusterset_ptr = p;
        *cluster_set = MarshalledClusterSetView(p);
        return true;
      }
    }
  }
};
