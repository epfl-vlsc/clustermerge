
import subprocess
from subprocess import PIPE

pidstat_cmd = ["pidstat", "-hrdu", "-p", "pid_list", "1", "|", "sed", "'1d;/^[#]/{{4,$d}};/^[#]/s/^[#][ ]*//;/^$/d;s/^[ ]*//;s/[ ]\+/,/g'"]

proc = subprocess.Popen(["./bazel-bin/src/clustermerge", "-c", "48", "-m", "48", "-t", "48",
  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chll2/chll2_metadata.json",
  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chlch/chlch_metadata.json", 
  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chll7/chll7_metadata.json", 
  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chlp8/chlp8_metadata.json"])

print("the pid is {}".format(proc.pid))

#cmd = pidstat_cmd.format(pid_list=proc.pid))
pidstat_cmd[3] = proc.pid

pid_proc = subprocess.Popen(pidstat_cmd, stdout=PIPE)

print("process returned with {}".format(proc.wait()))

output = pid_proc.stdout.read()

print("the output is \n{}".format(output))