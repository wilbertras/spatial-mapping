from ppadb.client import Client as AdbClient
import matplotlib.pyplot as plt
from datetime import datetime
import functions as f
from copy import copy
import numpy as np
import itertools
import pygame
import time
import json
import sys
import os

try:
    os.startfile("scrcpy-win64-v211\scrcpy.exe")
    client = AdbClient(host="127.0.0.1", port=5037) # Default is "127.0.0.1" and 5037
except:
    print('No phone connected')
devices = client.devices()
if len(devices) == 0:
    print('No devices')
    quit()
device = devices[0]
print(f'Connected to {device}')

output = device.shell("dumpsys battery")
for line in output.splitlines():
    if "level" in line:
        battery_level = int(line.split(":")[1].strip())
        print(f"Battery level: {battery_level}%")
        break