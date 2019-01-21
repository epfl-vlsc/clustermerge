
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
  srcs = ["main.cc"],
  copts = ["-std=c++14"],
  deps = ["@com_google_absl//absl/strings",
          "@args//:args",
          "//src/common:common"
  ],
)