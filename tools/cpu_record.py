
import subprocess

proc = subprocess.Popen(["./bazel-bin/src/clustermerge", "-c", "48", "-m", "48", "-t", "48",
  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chll2/chll2_metadata.json",
  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chlch/chlch_metadata.json", 
  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chll7/chll7_metadata.json", 
  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chlp8/chlp8_metadata.json"])

print("the pid is {}".format(proc.pid))

print("process returned with {}".format(proc.wait()))