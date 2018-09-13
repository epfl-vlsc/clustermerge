
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

# actually generate debug symbols ffs
# only applies to MacOS
genrule(
    name = "main_dsym",
    srcs = [":clustermerge"],
    outs = ["clustermerge.dSYM"],
    output_to_bindir = True,
    cmd = "dsymutil $(location :clustermerge) -o $@"
)

cc_binary(
    name = "clustermerge",
    srcs = glob(["*.cc", "*.h"]),
    linkopts = ["-lpthread"],
    copts = ["-std=c++14"],
    deps = ["@com_google_absl//absl/strings",
            ":swps3_lib",
            "//src/agd:agd"]
)


