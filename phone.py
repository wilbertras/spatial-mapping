from ppadb.client import Client as AdbClient
import pygame
import functions as f
import numpy as np
from datetime import datetime
import os
import matplotlib.pyplot as plt
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

if len(devices) == 0:
    print('No devices')
    quit()

device = devices[0]

print(f'Connected to {device}')

pygame.init()
# Constants
red = 255, 0, 0
green = 0, 255, 0
blue = 0, 0, 255
white = 255, 255, 255
black = 0, 0, 0
colors = [green, white, red, blue]
color_cycle = 0
WIDTH, HEIGHT = 500, 500
BACKGROUND_COLOR = black
TEXT_COLOR = green
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

## Initiate scanning paramaters
running = 1
nr_s21 = 0
x = 0
y = 0
nr_x_scanned = 0
nr_y_scanned = 0
linewidth = 1
stepsize_x = 1
stepsize_y = 1
nr_scans = 0
restart = 1
inverted = 0
green = 1
nr_x_scans = 0
nr_y_scans = 0
measure = 0
wait = 100

## Input S21 parameters
fstop = 6  # GHz
fstart = 4  # GHz
totscanbw = fstop - fstart
num_points = 3201
subscanbw = 0.1  # GHz
num_subscans = int(np.ceil(totscanbw / subscanbw))
realfstart = fstart
realfstop = fstart + num_subscans * subscanbw
len_s21 = int(num_subscans * num_points)
kidpower = -110 # dBm
ifbw = 10000  # Hz
freqs = np.linspace(realfstart, realfstop, num_points*num_subscans)

## Connect to Virtual Intstruments
vna = f.connect2vi("GPIB0::16::INSTR", timeout=3000000)
weinschell = f.connect2vi("GPIB0::10::INSTR", timeout=300000)

while running:
    if restart:
        device.shell("input keyevent 67")
        pygame.time.wait(wait)
        restart = 0

    text_list = [
    'X, Y = ' + str(x) + ', ' + str(y) + '\n Move = up,down,left,right',
    '# X scans, # Y scans = ' + str(nr_x_scanned) + ', ' + str(nr_y_scanned),
    'dX = ' + str(stepsize_x) + '\n -1/0/+1 = Z/X/C',
    'dY = ' + str(stepsize_y) + '\n -1/0/+1 = T/Y/U',
    'w = ' + str(linewidth) + '\n -1/0/+1 = Q/W/E',
    '# scans = ' + str(nr_s21) + '\n Make scan = Enter',
    'invert screen = I',
    'Change color = G \nCycles White,Red,Blue,Green',
    'reset lines = Backspace',
    'quit and save = Esc',
    ]
    screen.fill(BACKGROUND_COLOR)
    draw_text_boxes(text_list, BACKGROUND_COLOR, TEXT_COLOR)
    
    if measure:
        if nr_x_scanned < nr_x_scans:
            plt.close()
            s21 = np.zeros((1, 1, len_s21))
            # freqs, s21 = f.get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw)
            s21s[nr_x_scanned, 0, :] = s21
            nr_x_scanned += 1
            if nr_x_scanned < nr_x_scans:
                device.shell("input keyevent KEYCODE_DPAD_RIGHT")
                pygame.time.wait(wait)
                x += 1
            fig, ax = plt.subplots()
            ax.plot(s21[0, 0, :])
            plt.draw()
        if (nr_x_scanned == nr_x_scans) & (nr_y_scanned < nr_y_scans):
            plt.close()
            s21 = np.zeros((1, 1, len_s21))
            # freqs, s21 = f.get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw)
            s21s[0, nr_y_scanned, :] = s21
            nr_y_scanned += 1
            if nr_y_scanned < nr_y_scans:
                device.shell("input keyevent KEYCODE_DPAD_UP")
                pygame.time.wait(wait)
                y += 1
            fig, ax = plt.subplots()
            ax.plot(s21[0, 0, :])
            plt.draw()
        if nr_y_scanned == nr_y_scans:
            measure = 0
            date = datetime.today()
            name = 'S21s/Scan_'+ timestamp() + '.npy'
            np.save(name, s21s)
            print('Saved: %s' % name)
    
    event = pygame.event.poll()
    if event.type == pygame.QUIT:
        running = 0 

    if event.type == pygame.KEYDOWN:
        if event.key == pygame.K_ESCAPE:
            running = 0
        if event.key == pygame.K_BACKSPACE:
            device.shell("input keyevent 67")
            pygame.time.wait(wait)
            x = 0
            y = 0
            linewidth = 1
            stepsize_x = 1
            stepsize_y = 1
            nr_s21 = 0
        if event.key == pygame.K_DOWN:
            device.shell("input keyevent KEYCODE_DPAD_DOWN")
            pygame.time.wait(wait)
            if y - stepsize_y <= 0:
                y = 0
            else:
                y -= stepsize_y
        if event.key == pygame.K_UP:
            device.shell("input keyevent KEYCODE_DPAD_UP")
            pygame.time.wait(wait)
            y += stepsize_y
        if event.key == pygame.K_RIGHT:
            device.shell("input keyevent KEYCODE_DPAD_RIGHT")
            pygame.time.wait(wait)
            x += stepsize_x
        if event.key == pygame.K_LEFT:
            device.shell("input keyevent KEYCODE_DPAD_LEFT")
            pygame.time.wait(wait)
            if x - stepsize_x <= 0:
                x = 0
            else:
                x -= stepsize_x
        if event.key == pygame.K_e:
            device.shell("input keyevent KEYCODE_E")
            pygame.time.wait(wait)
            linewidth += 1
        if event.key == pygame.K_q:
            device.shell("input keyevent KEYCODE_Q")
            pygame.time.wait(wait)
            if linewidth == 1:
                linewidth = 1
            else:
                linewidth -= 1
        if event.key == pygame.K_c:
            device.shell("input keyevent KEYCODE_C")
            pygame.time.wait(wait)
            stepsize_x += 1
        if event.key == pygame.K_z:
            device.shell("input keyevent KEYCODE_Z")
            pygame.time.wait(wait)
            if stepsize_x == 1:
                stepsize_x = 1
            else:
                stepsize_x -= 1
        if event.key == pygame.K_u:
            device.shell("input keyevent KEYCODE_U")
            pygame.time.wait(wait)
            stepsize_y += 1
        if event.key == pygame.K_t:
            device.shell("input keyevent KEYCODE_T")
            pygame.time.wait(wait)
            if stepsize_y == 1:
                stepsize_y = 1
            else:
                stepsize_y -= 1
        if event.key == pygame.K_x:
            device.shell("input keyevent KEYCODE_X")
            pygame.time.wait(wait)
            stepsize_x = 1
            print('dX = 1')
        if event.key == pygame.K_y:
            device.shell("input keyevent KEYCODE_Y")
            pygame.time.wait(wait)
            stepsize_y = 1
            print('dY = 1')
        if event.key == pygame.K_w:
            device.shell("input keyevent KEYCODE_W")
            pygame.time.wait(wait)
            linewidth = 1
            print('w = 1')
        if event.key == pygame.K_i:
            device.shell("input keyevent KEYCODE_I")
            pygame.time.wait(wait)
            if inverted:
                BACKGROUND_COLOR = black
                TEXT_COLOR = green
                inverted = 0
            else:
                BACKGROUND_COLOR = white
                TEXT_COLOR = black
                inverted = 1
        if event.key == pygame.K_g:
            device.shell("input keyevent KEYCODE_G")
            pygame.time.wait(wait)
            if inverted:
                pass
            else:
                color_cycle += 1
                TEXT_COLOR = colors[color_cycle % 4]
        if event.key == pygame.K_RETURN:
                if linewidth != 1:
                    print('Putting w to 1')
                    device.shell("input keyevent KEYCODE_W")
                    pygame.time.wait(wait)
                elif stepsize_x != 1:
                    print('Putting dX to 1')
                    device.shell("input keyevent KEYCODE_X")
                    pygame.time.wait(wait)
                elif stepsize_y != 1:
                    print('Putting dY to 1')
                    device.shell("input keyevent KEYCODE_Y")
                    pygame.time.wait(wait) 
                 
                nr_x_scans = int(input('Please input the number of scans in x: '))
                nr_y_scans = int(input('Please input the number of scans in y: '))
                s21s = np.zeros((nr_x_scans, nr_y_scans, len_s21))
                measure = 1
                
    pygame.display.flip()

# Quit Pygame
pygame.quit()
sys.exit()

