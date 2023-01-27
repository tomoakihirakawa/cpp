from os.path import expanduser
import os
import math
import json
from math import pi
import platform

home = expanduser("~")

if platform.system() == "Linux":
    program_home = home + "/research"
else:
    program_home = home + "/Dropbox/markdown"

rho = 1000.
g = 9.81
