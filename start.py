__author__ = 'Felipe'

import subprocess
import time

for i in range(10):
    print "Running Process " + str(i)
    p = subprocess.Popen("python Main.py", shell=True)
    p.communicate()
    time.sleep(10)

print "Finish"