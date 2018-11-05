# ClusterMerge

Bottom up merging to cluster similar protein sequences. 

Install dependencies:

```sh
$ sudo apt-get install bazel
```

Bazel will use the default C++ compiler on your system. 

## Building

Linux (dbg):
`bazel build -c dbg src:clustermerge`

Note: I'm not sure how to increase the default stack size limit on Mac, so currently clustermerge will only work on Linux.
Mac (dbg):
`bazel build -c dbg --spawn_strategy=standalone src:main_dsym`

Mac/Linux (opt)
`bazel build -c opt src:clustermerge`

## Running 

First, set the thread stack limit to 64MB:

```sh
$ ulimit -s 65532000
```

This is to accomodate some alignment routines that stack allocate a lot of data. 

Clustermerge usage:

```sh
$ ./bazel-bin/src/clustermerge file1.fa file2.fa

$ ./bazel-bin/src/clustermerge -h  # view help
```

By default clustermerge will use as many threads as are available on your system.
You can optionally provide a list of your datasets to cluster to avoid typing so much crap on the command line:

```json
# data.json
[
  "file1.fa",
  "file2.fa",
  "fileN.fa"
]
```

and run with

```sh
$ ./bazel-bin/src/clustermerge -i data.json
```

PAM matrices for optimal alignments using SWPS3 are provided in the repo under `data/matrices/json`,
if you run the tool elsewhere you will need to use the `-d / --data_dir` option to provide the path.


If you have complete all vs. all results you want to compare against, `tools` has a compare script:

```sh
python2 tools/compare_results.py path/to/AllAll/ path/to/clustermerge/output_matches/
```


