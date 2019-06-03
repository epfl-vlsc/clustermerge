
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

#local_repository(
  # Name of the Abseil repository. This name is defined within Abseil's
  # WORKSPACE file, in its `workspace()` metadata
  #name = "com_google_absl",

  # NOTE: Bazel paths must be absolute paths. E.g., you can't use ~/Source
  #path = "third_party/abseil-cpp",
#)

local_repository(
  name = "zmq",

  path = "third_party/zmq",
)

local_repository(
  name = "json",

  path = "third_party/json",
)

local_repository(
  name = "args",

  path = "third_party/args",
)

http_archive(
    name = "com_google_protobuf",
    sha256 = "9510dd2afc29e7245e9e884336f848c8a6600a14ae726adb6befdb4f786f0be2",
    strip_prefix = "protobuf-3.6.1.3",
    urls = ["https://github.com/google/protobuf/archive/v3.6.1.3.zip"],
)

git_repository(
  name="com_google_absl",
  remote="https://github.com/abseil/abseil-cpp",
  branch="lts_2018_12_18"
)
