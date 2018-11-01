
cc_binary(
  name = "controller",
  srcs = glob(["*.cc"]),
  linkopts = ["-L/usr/local/lib", "-lzmq"],
  copts = ["-I/usr/local/include", "-L/usr/local/lib"],
  deps = ["//src/dist/proto:clustermerge_cc_proto",
          "@zmq//:zmq-cpp"
  ]
)