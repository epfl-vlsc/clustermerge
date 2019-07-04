import subprocess
import psutil
import time
import csv
import pathlib
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__file__)
log.setLevel(level=logging.DEBUG)

proc = subprocess.Popen(
    [
        "./bazel-bin/src/dist/dist_cluster",
        "-C",
        "-t",
        "46",
        "-b",
        "40",
        "-x",
        "-i",
        "bacteria_fa_full.json",
    ]
)
"""
proc = subprocess.Popen(
    [
        "./bazel-bin/src/clustermerge",
        "-x",
        "-i",
        "bacteria_datasets.json",
    ]
)
"""

prior_time = None
outpath = pathlib.Path("./net_stat.csv")
prior_counters = None
iface_name = "eno1"
net_stats = []
record_interval = 1.0

while proc.poll() is None:

    current_time = time.time()
    current_counters = psutil.net_io_counters(pernic=True)
    if prior_counters is not None:
        prior_nic = prior_counters[iface_name]
        current_nic = current_counters[iface_name]
        elapsed = current_time - prior_time
        assert elapsed > 0, "Time diff was not >0. Got {}".format(elapsed)
        elapsed = float(elapsed)
        rx_bytes = current_nic.bytes_recv - prior_nic.bytes_recv
        assert rx_bytes >= 0, "Got a negative rx byte diff {}".format(rx_bytes)
        tx_bytes = current_nic.bytes_sent - prior_nic.bytes_sent
        assert tx_bytes >= 0, "Got a negative tx byte diff {}".format(tx_bytes)
        rx_packets = current_nic.packets_recv - prior_nic.packets_recv
        assert rx_packets >= 0, "Got a negative rx byte diff {}".format(rx_packets)
        tx_packets = current_nic.packets_sent - prior_nic.packets_sent
        assert tx_packets >= 0, "Got a negative tx byte diff {}".format(tx_packets)
        net_stats.append(
            {
                "time": current_time,
                "rx_bytes/sec": rx_bytes / elapsed,
                "tx_bytes/sec": tx_bytes / elapsed,
                "rx_packets/sec": rx_packets / elapsed,
                "tx_packets/sec": tx_packets / elapsed,
                "rx_bytes": current_nic.bytes_recv,
                "tx_bytes": current_nic.bytes_sent,
                "rx_packets": current_nic.packets_recv,
                "tx_packets": current_nic.packets_sent,
            }
        )

    prior_counters = current_counters
    prior_time = current_time
    time.sleep(record_interval)

print("process returned {}".format(proc.returncode))
if len(net_stats) > 0:
    with outpath.open("w") as f:
        writer = csv.DictWriter(f=f, fieldnames=net_stats[0].keys())
        writer.writeheader()
        writer.writerows(net_stats)
else:
    log.warning("Net stats didn't record anything")

