
# compile swps3 into a lib with just the things we need (SSE, scalar)
cc_library(
    name = "swps3_lib",
    srcs = [
            "swps3/DynProgr_sse_short.c",
            "swps3/DynProgr_sse_byte.c",
            "swps3/DynProgr_sse_double.c",
            "swps3/DynProgr_scalar.c",
            "swps3/matrix.c",
            "swps3/extras.c",
            "swps3/fasta.c",
            "swps3/debug.c",
            "swps3/EstimatePam.c",
            "swps3/Page_size.c"
            ],
    hdrs = ["swps3/swps3.h",
            "swps3/DynProgr_sse_short.h",
            "swps3/DynProgr_sse_byte.h",
            "swps3/DynProgr_sse_double.h",
            "swps3/DynProgr_scalar.h",
            "swps3/extras.h",
            "swps3/matrix.h",
            "swps3/debug.h",
            "swps3/fasta.h",
            "swps3/EstimatePam.h",
            "swps3/Page_size.h"
            ],
    linkopts = [
        "-msse3"
    ],
    linkstatic = 1,
    alwayslink = 1,
    visibility = ["//visibility:public"],
    copts = ["-msse3"]
)

cc_library(
    name = "common",
    srcs = glob(["*.cc", "*.h"]),
    # bazel currently doesn't look in /usr/local for some reason
    # so we need to include it for ZMQ
    linkopts = ["-lpthread"],
    copts = ["-std=c++14", "-I/usr/local/include", "-L/usr/local/lib"],
    visibility = ["//visibility:public"],
    deps = ["@com_google_absl//absl/strings",
            "@com_google_absl//absl/container:flat_hash_map",
            "@com_google_absl//absl/container:flat_hash_set",
            "@json//:json-cpp",
            "@args//:args",
            #"//src/proto:clustermerge_cc_proto",
            "//src/agd:agd",
            "//src/dataset:dataset",
            "//src/comms:comms",
            ":swps3_lib",
           ]
)
