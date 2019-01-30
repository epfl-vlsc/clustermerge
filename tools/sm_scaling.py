import subprocess
from subprocess import PIPE


#proc = subprocess.Popen(["./bazel-bin/src/clustermerge", "-x", "-c", "48", "-m", "48", "-t", "48", "-i", "bacteria_datasets.json"], stdout=subprocess.PIPE) 
#proc = subprocess.Popen(["./bazel-bin/src/clustermerge", "-x", "-c", "48", "-m", "48", "-t", "48", "-i", "four_datasets.json"], stdout=subprocess.PIPE) 


#print("the output is: {}".format(out))

threads = [1, 2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 47, 48]
#threads = [48, 47]

outfile = open("sm_scaling.csv", "w")

for t in threads:

  total_time = 0
  for i in range(1):

    proc = subprocess.Popen(["./bazel-bin/src/clustermerge", "-x", "-c", "48", "-m", str(t), "-t", "48", "-i", "bacteria_medium.json"], stdout=subprocess.PIPE) 
    print("the pid is {}".format(proc.pid))

    out, err = proc.communicate()

    lines = out.split('\n')
    for l in lines:
      if (l.startswith('Clustering')):
        s = l.split(' ')
        time = s[-1]
        print("the time is {}".format(time))
        total_time += int(time)

  print("average time is {}".format(total_time)) 
  outfile.write("{}, {}\n".format(t, total_time))

outfile.close()
