
import subprocess
from subprocess import PIPE

pidstat_cmd = ["pidstat", "-hrdu", "-p", "pid_list", "1", ]

sed_cmd = ["sed", r'1d;/^[#]/{{4,$d}};/^[#]/s/^[#][ ]*//;/^$/d;s/^[ ]*//;s/[ ]\+/,/g']

#proc = subprocess.Popen(["./bazel-bin/src/clustermerge", "-c", "48", "-m", "48", "-t", "48",
#  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chll2/chll2_metadata.json",
#  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chlch/chlch_metadata.json",
#  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chll7/chll7_metadata.json", 
#  "/scratch/proteomes/agd/bacteria_chlorobiaceae/chlp8/chlp8_metadata.json"])

proc = subprocess.Popen(["./bazel-bin/src/clustermerge", "-c", "48", "-m", "48", "-t", "48", "-i", "bacteria_datasets.json"]) 

print("the pid is {}".format(proc.pid))

#cmd = pidstat_cmd.format(pid_list=proc.pid))
pidstat_cmd[3] = str(proc.pid)

print("pid stat command is: {}".format(pidstat_cmd))
print("sed command is: {}".format(sed_cmd))

pid_output = open("pidout.txt", "w+")
pid_proc = subprocess.Popen(pidstat_cmd, stdout=pid_output)
print("process returned with {}".format(proc.wait()))
pid_proc.wait()
pid_output.close()

#pidstat_output = pid_proc.stdout.read()

#print("the pidstat output is \n{}".format(pidstat_output))

sed_in = open("pidout.txt", "r")
sed_out = open("cpu.csv", "w+")
sed_proc = subprocess.Popen(sed_cmd, stdin=sed_in, stdout=sed_out)


sed_proc.wait()

    
