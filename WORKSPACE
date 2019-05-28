
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

local_repository(
  # Name of the Abseil repository. This name is defined within Abseil's
  # WORKSPACE file, in its `workspace()` metadata
  name = "com_google_absl",

  # NOTE: Bazel paths must be absolute paths. E.g., you can't use ~/Source
  path = "third_party/abseil-cpp",
)

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
