
cc_binary(
  name = "dist_cluster",
  srcs = glob(["*.cc", "*.h"]),
  linkopts = ["-L/usr/local/lib", "-L/cluster/apps/zmq/4.3.1/x86_64/lib", "-lzmq"],
  copts = ["-I/usr/local/include", "-L/usr/local/lib", "-L/cluster/apps/zmq/4.3.1/x86_64/lib"],
  deps = [
          #"//src/proto:clustermerge_cc_proto",
          "@com_google_absl//absl/container:node_hash_map",
          "@com_google_absl//absl/container:flat_hash_set",
          "@zmq//:zmq-cpp",
          "@args//:args",
          "//src/common:common",
          "//src/comms:comms"
         ]
)
