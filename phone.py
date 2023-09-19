import time
from ppadb.client import Client as AdbClient
import pygame

client = AdbClient(host="127.0.0.1", port=5037) # Default is "127.0.0.1" and 5037
devices = client.devices()

if len(devices) == 0:
    print('No devices')
    quit()

device = devices[0]
screen = pygame.display.set_mode((500, 500))
print(f'Connected to {device}')
# Prompt the user for the keycode to send
running = 1
while running:
    event = pygame.event.poll()
    if event.type == pygame.QUIT:
        running = 0 

    if event.type == pygame.KEYDOWN:
        if event.key == pygame.K_DOWN:
            device.shell("input keyevent KEYCODE_DPAD_DOWN")
            time.sleep(1)

        if event.key == pygame.K_RIGHT:
            device.shell("input keyevent KEYCODE_DPAD_RIGHT")
            time.sleep(1)

        if event.key == pygame.K_SPACE:
            device.shell("input keyevent KEYCODE_SPACE")
            time.sleep(0)