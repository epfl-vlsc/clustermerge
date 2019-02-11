
cc_binary(
  name = "dist_cluster",
  srcs = glob(["*.cc", "*.h"]),
  linkopts = ["-L/usr/local/lib", "-lzmq"],
  copts = ["-I/usr/local/include", "-L/usr/local/lib"],
  deps = [
          "//src/proto:clustermerge_cc_proto",
          "@com_google_absl//absl/container:node_hash_map",
          "@com_google_absl//absl/container:flat_hash_set",
          "@zmq//:zmq-cpp",
          "@args//:args",
          "//src/common:common"
         ]
)