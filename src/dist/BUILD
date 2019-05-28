
cc_library(
  name = "dist_lib",
  srcs = glob(["*.cc", "*.h"], exclude=["dist_cluster.cc"]),
  linkopts = ["-L/usr/local/lib", "-lzmq"],
  copts = ["-I/usr/local/include", "-L/usr/local/lib"],
  visibility = ["//visibility:public"],
  deps = [
          "@com_google_absl//absl/container:node_hash_map",
          "@com_google_absl//absl/container:flat_hash_set",
          "@com_google_absl//absl/strings",
          "@zmq//:zmq-cpp",
          "//src/common:common",
          "//src/comms:comms"
         ]
)

cc_binary(
  name = "dist_cluster",
  srcs = glob(["dist_cluster.cc"]),
  linkopts = ["-L/usr/local/lib", "-lzmq"],
  copts = ["-I/usr/local/include", "-L/usr/local/lib"],
  deps = [
          "@args//:args",
          ":dist_lib"
         ]
)