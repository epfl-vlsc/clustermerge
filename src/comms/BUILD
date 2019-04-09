
cc_library(
    name = "comms",
    srcs = glob(["*.h"]),
    hdrs = glob(["*.h"]),
    # bazel currently doesn't look in /usr/local for some reason
    # so we need to include it for ZMQ
    #copts = ["-std=c++14"],
    linkopts = ["-lzmq"],
    copts = ["-I/usr/local/include", "-L/usr/local/lib"],
    visibility = ["//visibility:public"],
    deps = [
      "//src/agd:agd",
      "@zmq//:zmq-cpp"
    ]
)