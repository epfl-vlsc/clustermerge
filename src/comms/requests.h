#pragma once

#include "src/agd/buffer.h"
#include "zmq.hpp"

enum RequestType { Batch = 0, Partial, LargePartial };

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
};

// for response queue, keep msg to prevent copies
struct MarshalledResponse {
  // contains marshalled id and cluster set
  zmq::message_t msg;
  MarshalledResponse() = default;
  MarshalledResponse(MarshalledResponse&&) = default;
  MarshalledResponse& operator=(MarshalledResponse&&) = default;
  int ID() { return *reinterpret_cast<int*>(msg.data()); }
  MarshalledClusterSetView Set() {
    const char* data = reinterpret_cast<const char*>(msg.data()) + sizeof(int);
    return MarshalledClusterSetView(data);
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
  
  MarshalledClusterSet(MarshalledResponse& response)
      : buf(agd::Buffer(response.msg.size() - sizeof(uint32_t), 128)) {
    // a response buffer is an int id and a marshalled clusterset
    char* data = reinterpret_cast<char*>(response.msg.data());
    data += sizeof(uint32_t);
    auto data_size = response.msg.size() - sizeof(uint32_t);
    // buf.reserve(data_size);
    buf.AppendBuffer(data, data_size);
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

  //Does not have a clusterset, but is still large :P
  void CreateLargePartialRequest(int id, MarshalledClusterView cluster) {
    buf.reserve(cluster.TotalSize() + sizeof(PartialRequestHeader));
    PartialRequestHeader h;
    h.id = id;
    h.type = RequestType::LargePartial;
    buf.AppendBuffer(reinterpret_cast<char*>(&h), sizeof(PartialRequestHeader));
    buf.AppendBuffer(cluster.data, cluster.TotalSize());
    assert(buf.size() == cluster.TotalSize() + sizeof(PartialRequestHeader));  
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
  void Cluster(MarshalledClusterView* cluster)  {
    const char* cluster_ptr = data + sizeof(PartialRequestHeader);
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
