
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

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
    name = "gtest",
    url = "https://github.com/google/googletest/archive/release-1.7.0.zip",
    sha256 = "b58cb7547a28b2c718d1e38aee18a3659c9e3ff52440297e965f5edffe34b6d0",
    build_file = "@//src:gtest.BUILD",
    strip_prefix = "googletest-release-1.7.0",
)

git_repository(
  name="com_google_absl",
  remote="https://github.com/abseil/abseil-cpp",
  branch="lts_2018_12_18"
)
