
cc_binary(
  name = "dist_cluster",
  srcs = glob(["*.cc", "*.h"]),
  linkopts = ["-L/usr/local/lib", "-lzmq"],
  copts = ["-I/usr/local/include", "-L/usr/local/lib"],
  deps = [
          "//src/proto:clustermerge_cc_proto",
          "@zmq//:zmq-cpp",
          "@args//:args",
          "//src/common:common"
         ]
)