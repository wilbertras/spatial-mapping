import time
from ppadb.client import Client as AdbClient
import pygame
# import vna
import pandas as pd
from datetime import datetime
import os


def timestamp():
    year = datetime.now().year
    month = datetime.now().month
    day = datetime.now().day
    hour = datetime.now().hour
    minute = datetime.now().minute
    return '%dh%d_%d-%d-%d' % (hour, minute, day, month, year)

try:
    client = AdbClient(host="127.0.0.1", port=5037) # Default is "127.0.0.1" and 5037
    os.startfile("scrcpy-win64-v211\scrcpy.exe")
except:
    print('No phone connected')

devices = client.devices()
cols = ['nr', 'x', 'y', 'freqs', 's21']
df = pd.DataFrame(columns=cols)

if len(devices) == 0:
    print('No devices')
    quit()

device = devices[0]
screen = pygame.display.set_mode((500, 500))
print(f'Connected to {device}')
# text_color = 0, 0, 0
# rect_exc = pygame.draw.rect(100, 100, 30, 50)
# text_esc = 'ESC: close and save'
# base_font = pygame.font.Font(None, 32)
# text_surface = base_font.render(text_esc, True, text_color)

running = 1
nr_s21 = 0
x = 0
y = 0
while running:
    # screen.blit(text_surface, rect_exc)
    event = pygame.event.poll()
    if event.type == pygame.QUIT:
        running = 0 

    if event.type == pygame.KEYDOWN:
        if event.key == pygame.K_ESCAPE:
            running = 0
            date = datetime.today()
            name = 'S21s/Scan_'+ timestamp()
            df.to_pickle(name)
            print('Saved: %s' % name)
        if event.key == pygame.K_BACKSPACE:
            device.shell("input keyevent 67")
            x = 0
            y = 0
            nr_s21 = 0
        if event.key == pygame.K_DOWN:
            device.shell("input keyevent KEYCODE_DPAD_DOWN")
            y += 1
        if event.key == pygame.K_UP:
            device.shell("input keyevent KEYCODE_DPAD_UP")
            if y == 0:
                y = 0
            else:
                y -= 1
        if event.key == pygame.K_RIGHT:
            device.shell("input keyevent KEYCODE_DPAD_RIGHT")
            x += 1
        if event.key == pygame.K_LEFT:
            device.shell("input keyevent KEYCODE_DPAD_LEFT")
            if x == 0:
                x = 0
            else:
                x -= 1
        if event.key == pygame.K_e:
            device.shell("input keyevent KEYCODE_E")
        if event.key == pygame.K_q:
            device.shell("input keyevent KEYCODE_Q")
        if event.key == pygame.K_c:
            device.shell("input keyevent KEYCODE_C")
        if event.key == pygame.K_z:
            device.shell("input keyevent KEYCODE_Z")
        if event.key == pygame.K_x:
            x = 0
            print('x = 0')
        if event.key == pygame.K_u:
            device.shell("input keyevent KEYCODE_U")
        if event.key == pygame.K_t:
            device.shell("input keyevent KEYCODE_T")
        if event.key == pygame.K_y:
            device.shell("input keyevent KEYCODE_T")
            y = 0
            print('y = 0')
        if event.key == pygame.K_i:
            device.shell("input keyevent KEYCODE_I")
        if event.key == pygame.K_RETURN:
            if (nr_s21 == 0) & (x != 0 | y != 0):
                print('Please set x and y to zero for the first scan')
            elif (x != 0) & (y != 0):
                print('Please set either x or y to zero for a scan')
            else:
                # freqs, s21 = vna.get_s21(4, 8, 101, 1000)
                freqs, s21 = 0
                pygame.time.wait(0)
                nr_s21 += 1
                print('Scanned: ', x, y)
                df.loc[len(df.index)] = [nr_s21, x, y, freqs, s21]
