from ppadb.client import Client as AdbClient
import pygame
# import vna
import pandas as pd
from datetime import datetime
import os
import sys


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

print(f'Connected to {device}')

pygame.init()
# Constants
WIDTH, HEIGHT = 500, 500
BACKGROUND_COLOR = (0, 0, 0)
TEXT_COLOR = (255, 255, 255)
FONT_SIZE = 20
FONT = pygame.font.Font(None, FONT_SIZE)
COLUMN_WIDTH = WIDTH // 2
BOX_HEIGHT = HEIGHT // 5

screen = pygame.display.set_mode((WIDTH, HEIGHT))


def draw_text_boxes(text_list):
    for i, text in enumerate(text_list):
        row = i // 2
        col = i % 2
        x = col * COLUMN_WIDTH
        y = row * BOX_HEIGHT
        pygame.draw.rect(screen, BACKGROUND_COLOR, (x, y, COLUMN_WIDTH, BOX_HEIGHT), 0)
        
        # Split the text into two lines
        lines = text.split('\n')
        
        # Render each line of text
        for j, line in enumerate(lines):
            text_surface = FONT.render(line, True, TEXT_COLOR)
            text_rect = text_surface.get_rect(center=(x + COLUMN_WIDTH // 2, y + (j * FONT_SIZE) + (BOX_HEIGHT // 4)))
            screen.blit(text_surface, text_rect)

running = 1
nr_s21 = 0
x = 0
y = 0
scan_x = 0
scan_y = 0
linewidth = 1
stepsize_x = 1
stepsize_y = 1
nr_scans = 0
restart = 1
while running:
    if restart:
        device.shell("input keyevent 67")
        restart = 0

    text_list = [
    '(X, Y) = (' + str(x) + ', ' + str(y) + ')',
    '(scan X, scan Y) = (' + str(scan_x) + ', ' + str(scan_y) + ') \n Set scan X to 0 = X, Set scan Y to 0 = Y',
    'stepsize X = ' + str(stepsize_x) + '\n -1 = Z; +1 = C',
    'stepsize Y = ' + str(stepsize_y) + '\n -1 = T; +1 = U',
    'linewidth = ' + str(linewidth) + '\n -1 = Q, +1 = E',
    'nr of scans = ' + str(nr_s21) + '\n Make scan = Enter',
    'invert screen = I',
    'green/white line = G',
    'reset lines = Backspace',
    'quit and save = Esc',
    ]
    screen.fill(BACKGROUND_COLOR)
    draw_text_boxes(text_list)

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
            linewidth = 1
            stepsize_x = 1
            stepsize_y = 1
            nr_s21 = 0
        if event.key == pygame.K_DOWN:
            device.shell("input keyevent KEYCODE_DPAD_DOWN")
            y += stepsize_y
            scan_y += stepsize_y
        if event.key == pygame.K_UP:
            device.shell("input keyevent KEYCODE_DPAD_UP")
            if y - stepsize_y <= 0:
                scan_y -= y
                y = 0
            else:
                y -= stepsize_y
                scan_y -= stepsize_y
        if event.key == pygame.K_RIGHT:
            device.shell("input keyevent KEYCODE_DPAD_RIGHT")
            x += stepsize_x
            scan_x += stepsize_x
        if event.key == pygame.K_LEFT:
            device.shell("input keyevent KEYCODE_DPAD_LEFT")
            if x - stepsize_x <= 0:
                scan_x -= x
                x = 0
            else:
                x -= stepsize_x
                scan_x -= stepsize_x
        if event.key == pygame.K_e:
            device.shell("input keyevent KEYCODE_E")
            linewidth += 1
        if event.key == pygame.K_q:
            device.shell("input keyevent KEYCODE_Q")
            linewidth -= 1
        if event.key == pygame.K_c:
            device.shell("input keyevent KEYCODE_C")
            stepsize_x += 1
        if event.key == pygame.K_z:
            device.shell("input keyevent KEYCODE_Z")
            if stepsize_x == 1:
                stepsize_x = 1
            else:
                stepsize_x -= 1
        if event.key == pygame.K_u:
            device.shell("input keyevent KEYCODE_U")
            stepsize_y += 1
        if event.key == pygame.K_t:
            device.shell("input keyevent KEYCODE_T")
            if stepsize_y == 1:
                stepsize_y = 1
            else:
                stepsize_y -= 1
        if event.key == pygame.K_x:
            scan_x = 0
            print('scan_x = 0')
        if event.key == pygame.K_y:
            device.shell("input keyevent KEYCODE_T")
            scan_y = 0
            print('scan_y = 0')
        if event.key == pygame.K_i:
            device.shell("input keyevent KEYCODE_I")
        if event.key == pygame.K_RETURN:
            if (nr_s21 == 0) & (scan_x != 0 | scan_y != 0):
                print('Please set scan_x and scan_y to zero for the first scan')
            elif (scan_x != 0) & (scan_y != 0):
                print('Please set either scan_x or scan_y to zero for a scan')
            else:
                # freqs, s21 = vna.get_s21(4, 8, 101, 1000)
                freqs, s21 = 0, 0
                pygame.time.wait(0)
                nr_s21 += 1
                print('Scanned: ', scan_x, scan_y)
                df.loc[len(df.index)] = [nr_s21, x, y, scan_x, scan_y, freqs, s21]
    pygame.display.flip()

# Quit Pygame
pygame.quit()
sys.exit()

