
cc_library(
    name = "agd",
    srcs = glob(["*.cc"]),
    hdrs = glob(["*.h", "*.hpp"]),
    #copts = ["-std=c++17"],
    linkopts = ["-lz"],
    visibility = ["//visibility:public"],
    deps = ["@com_google_absl//absl/strings",
        "@com_google_absl//absl/synchronization",
    ]
)