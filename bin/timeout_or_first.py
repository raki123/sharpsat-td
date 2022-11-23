import subprocess
import sys
import time

start = time.time()
p = subprocess.Popen(["./flow_cutter_pace17"], stdin=sys.stdin, stdout=subprocess.PIPE, close_fds = True)
timeout = float(sys.argv[1])
first = False
while p.poll() is None and not first:
    line = p.stdout.readline().decode()
    first = line.startswith("c status")
time_left = timeout - (time.time() - start)
try:
    p.wait(time_left)
except subprocess.TimeoutExpired:
    p.terminate()

print(p.stdout.read().decode()[:-1])

