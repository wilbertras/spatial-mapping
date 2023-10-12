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
    # os.startfile("scrcpy-win64-v211\scrcpy.exe")
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
TEXT_COLOR = (0, 255, 0)
COLUMN_WIDTH = WIDTH // 2
SMALL_BOX_HEIGHT = HEIGHT // 8
BIG_BOX_HEIGHT = HEIGHT // 2
SMALL_FONT_SIZE = 20
BIG_FONT_SIZE = 30

screen = pygame.display.set_mode((WIDTH, HEIGHT))


def draw_text_boxes(text_list, BACKGROUND_COLOR, TEXT_COLOR):
    y=0
    for i, text in enumerate(text_list):
        row = i // 2
        col = i % 2
        if row==0:
            BOX_HEIGHT = BIG_BOX_HEIGHT
            FONT_SIZE = BIG_FONT_SIZE
            y = 0
        else:
            y = BIG_BOX_HEIGHT + (row-1)*SMALL_BOX_HEIGHT
            BOX_HEIGHT = SMALL_BOX_HEIGHT
            FONT_SIZE = SMALL_FONT_SIZE
        x = col * COLUMN_WIDTH
        pygame.draw.rect(screen, BACKGROUND_COLOR, (x, y, COLUMN_WIDTH, BOX_HEIGHT), 0)
        pygame.draw.rect(screen, TEXT_COLOR, (x, y, COLUMN_WIDTH, BOX_HEIGHT), 1)
        # Split the text into two lines
        lines = text.split('\n')
        
        # Render each line of text
        for j, line in enumerate(lines):
            if j==0:
                FONT_SIZE = BIG_FONT_SIZE
            else:
                FONT_SIZE = SMALL_FONT_SIZE
            FONT = pygame.font.Font(None, FONT_SIZE)
            text_surface = FONT.render(line, True, TEXT_COLOR)
            text_rect = text_surface.get_rect(center=(x + COLUMN_WIDTH // 2, y + (j * FONT_SIZE) + (BOX_HEIGHT // 2)))
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
inverted = 0
green = 1
while running:
    if restart:
        device.shell("input keyevent 67")
        restart = 0

    text_list = [
    'X, Y = ' + str(x) + ', ' + str(y) + '\n Move = up,down,left,right',
    'X_scan, Y_scan = ' + str(scan_x) + ', ' + str(scan_y) + '\n Set2Zero = X, Y',
    'dX = ' + str(stepsize_x) + '\n -/+ = Z/C',
    'dY = ' + str(stepsize_y) + '\n -/+ = T/U',
    'w = ' + str(linewidth) + '\n -/+ = Q/E',
    '# scans = ' + str(nr_s21) + '\n Make scan = Enter',
    'invert screen = I',
    'green/white line = G',
    'reset lines = Backspace',
    'quit and save = Esc',
    ]
    screen.fill(BACKGROUND_COLOR)
    draw_text_boxes(text_list, BACKGROUND_COLOR, TEXT_COLOR)

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
            if y - stepsize_y <= 0:
                scan_y -= y
                y = 0
            else:
                y -= stepsize_y
                scan_y -= stepsize_y
        if event.key == pygame.K_UP:
            device.shell("input keyevent KEYCODE_DPAD_UP")
            y += stepsize_y
            scan_y += stepsize_y
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
            if linewidth == 1:
                linewidth = 1
            else:
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
            if inverted:
                BACKGROUND_COLOR = 0, 0, 0
                TEXT_COLOR = 0, 255, 0
                inverted = 0
            else:
                BACKGROUND_COLOR = 255, 255, 255
                TEXT_COLOR = 0, 0, 0
                inverted = 1
        if event.key == pygame.K_g:
            device.shell("input keyevent KEYCODE_G")
            if inverted:
                pass
            else:
                if green:
                    TEXT_COLOR = 255, 255, 255
                    green = 0
                else:
                    TEXT_COLOR = 0, 255, 0
                    green = 1
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

