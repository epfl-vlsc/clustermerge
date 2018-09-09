# clustermerge
test bottom up cluster merging

Linux (dbg):
`bazel build -c dbg src:clustermerge`

Mac (dbg):
`bazel build -c dbg --spawn_strategy=standalone src:main_dsym`

Mac/Linux (opt)
`bazel build -c opt src:clustermerge`