from ppadb.client import Client as AdbClient
import os

try:
    # os.startfile("scrcpy-win64-v211\scrcpy.exe")
    client = AdbClient(host="127.0.0.1", port=5037) # Default is "127.0.0.1" and 5037
except:
    print('No phone connected')
devices = client.devices()
if len(devices) == 0:
    print('No devices')
    quit()
else:
    print('devices: ', devices)
device = devices[0]
device.shell("reboot -p")