
cc_library(
  name = "zmq-cpp",
  hdrs = ["zmq.hpp"],
  copts = ["-std=c++14", "-I/cluster/apps/zmq/4.3.1/x86_64/include", "-L/cluster/apps/zmq/4.3.1/x86_64/lib"],
  linkopts = ["-L/cluster/apps/zmq/4.3.1/x86_64/lib", "-lzmq"],
  visibility = ["//visibility:public"]
)

