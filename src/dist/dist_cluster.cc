
// Entry point for the master server for the distributed
// clustering implementation

#include <iostream>
#include "src/dist/proto/requests.pb.h"
#include "zmq.hpp"

using namespace std;

/*
Server cluster format example
{
  "controller": "<ip/addr>",
  "push_port": <port num>,
  "pull_port": <port num>,
  "servers": [
    "<ip/addr>", ...
  ]
}
*/

int main(int argc, char* argv[]) {

  cmproto::Cluster c;
  c.add_indexes(1);

  cmproto::MergeBatch batch;
  auto* set = batch.add_sets();
  auto* cluster = set->add_clusters();
  *cluster = c;

  cmproto::MergeRequest r;
  /*auto* b = r.mutable_batch();
  *b = batch;*/

  if (r.has_batch()) {
    cout << "its a batch\n";
  } else if (r.has_partial()) {
    cout << "has partial\n";
  } else {
    cout << "has none\n";
  }

  // parse servers json, determine role
  // load datasets

  // if controller
  // connect to zmq queues
  // push initial batch mergers to ingress queue
  // loop (possibly with threads), pulling from egress queue
  //   construct new batches
  //   if received cluster set is large (some threshold?), construct partial merge requests
  
  // else if worker server
  return 0;

}