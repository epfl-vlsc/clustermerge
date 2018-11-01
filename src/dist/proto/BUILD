
cc_proto_library(
    name = "clustermerge_cc_proto",
    deps = [":clustermerge_proto"],
    visibility = ["//visibility:public"]
)
proto_library(
    name = "clustermerge_proto",
    srcs = ["requests.proto", "responses.proto"],
    deps = [":cluster_proto"]
)

proto_library(
    name = "cluster_proto",
    srcs = ["cluster.proto"]
)
