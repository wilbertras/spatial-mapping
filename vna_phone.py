from ppadb.client import Client as AdbClient
import pygame
import functions as f
import numpy as np
from datetime import datetime
import os
import matplotlib.pyplot as plt
import sys
import json
from copy import copy


def timestamp():
    year = datetime.now().year
    month = datetime.now().month
    day = datetime.now().day
    hour = datetime.now().hour
    minute = datetime.now().minute
    return '%dh%d_%d-%d-%d' % (hour, minute, day, month, year)

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

pygame.init()
# Constants
red = 255, 0, 0
green = 0, 255, 0
blue = 0, 0, 255
white = 255, 255, 255
black = 0, 0, 0
colors = [white, blue, green, red]
nr_colors = len(colors)
color_cycler = 0
width, height = 600, 500
bgcolor = black
linecolor = white
column_width = width // 2
small_box_height = height // 8
big_box_height = height // 2
small_font_size = 20
big_font_size = 30

screen = pygame.display.set_mode((width, height))

def draw_text_boxes(text_list, bgcolor, linecolor):
    y=0
    for i, text in enumerate(text_list):
        row = i // 2
        col = i % 2
        if row==0:
            BOX_height = big_box_height
            font_size = big_font_size
            y = 0
        else:
            y = big_box_height + (row-1)*small_box_height
            BOX_height = small_box_height
            font_size = small_font_size
        x = col * column_width
        pygame.draw.rect(screen, bgcolor, (x, y, column_width, BOX_height), 0)
        pygame.draw.rect(screen, linecolor, (x, y, column_width, BOX_height), 1)
        # Split the text into two lines
        lines = text.split('\n')
        
        # Render each line of text
        for j, line in enumerate(lines):
            if j==0:
                font_size = big_font_size
            else:
                font_size = small_font_size
            FONT = pygame.font.Font(None, font_size)
            text_surface = FONT.render(line, True, linecolor)
            text_rect = text_surface.get_rect(center=(x + column_width // 2, y + (j * font_size) + (BOX_height // 2)))
            screen.blit(text_surface, text_rect)


screen_width = 1080
screen_height = 2240
x_centre = int(screen_width / 2)
y_centre = int(screen_height / 2)
length_line_mm = 50
pxl_mm = 0.06
length_line_pxl = int(length_line_mm / pxl_mm)
half_length = int(length_line_pxl / 2)
y_min = y_centre - half_length
y_max = y_centre + half_length
x_min = x_centre - half_length
x_max = x_centre + half_length
x_init = x_centre
y_init = y_centre
dx_init = 10
dx_min = 1
dy_min = 1
dy_init = 10
w_init = 5
w_min = 1
step = 1

## Initiate scanning paramaters
running = 1
nr_s21 = 0
nr_x_scanned = 0
nr_y_scanned = 0
w = w_init
x = x_init
y = y_init
dx = dx_init
dy = dy_init
nr_scans = 0
restart = 1
inverted = 0
green = 1
nr_x_scans = 0
nr_y_scans = 0
measure = 0
wait = 50
dark = 0

## Input S21 parameters
fstart = 4  # GHz
fstop = 6  # GHz
totscanbw = fstop - fstart
num_points = 6401
subscanbw = 0.1  # GHz
num_subscans = int(np.ceil(totscanbw / subscanbw))
realfstart = fstart
realfstop = fstart + num_subscans * subscanbw
len_s21 = int(num_subscans * num_points)
kidpower = -110 # dBm
ifbw = 10000  # Hz
freqs = np.linspace(realfstart, realfstop, num_points*num_subscans)
date = datetime.today()

## Test connection to Virtual Intstruments
vna = f.connect2vi("GPIB0::16::INSTR", timeout=3000000)
weinschell = f.connect2vi("GPIB0::10::INSTR", timeout=300000)

while running:
    if restart:
        device.shell("input keyevent 67")
        pygame.time.wait(wait)
        restart = 0

    text_list = [
    'X, Y = ' + str(x) + ', ' + str(y) + '\n Move = up,down,left,right',
    '#X_scans, #Y_scans = ' + str(nr_x_scanned) + ', ' + str(nr_y_scanned),
    'dX = ' + str(dx) + '\n -1/=1/+1 = Z/X/C',
    'dY = ' + str(dy) + '\n -1/=1/+1 = T/Y/U',
    'w = ' + str(w) + '\n -1/=1/+1 = Q/W/E',
    '# scans = ' + str(nr_s21) + '\n Make scan = Enter',
    'invert screen = I',
    'Change color = G \nCycles White-Blue-Green-Red-Off',
    'reset lines = Backspace',
    'quit and save = Esc',
    ]
    screen.fill(bgcolor)
    draw_text_boxes(text_list, bgcolor, linecolor)
    
    if measure:
        if nr_x_scanned < nr_scans:
            if nr_x_scanned == 0:
                device.shell("input keyevent KEYCODE_B")
                pygame.time.wait(wait)
            _, s21 = f.get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw)
            s21s[nr_x_scanned, 0, :] = s21
            nr_x_scanned += 1
            if nr_x_scanned < nr_scans:
                device.shell("input keyevent KEYCODE_DPAD_RIGHT")
                pygame.time.wait(wait)
                x += 1
        if (nr_x_scanned == nr_scans) & (nr_y_scanned < nr_scans):
            if nr_y_scanned == 0:
                device.shell("input keyevent KEYCODE_B")
                pygame.time.wait(wait)
            _, s21 = f.get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw)
            s21s[nr_y_scanned, 1, :] = s21
            nr_y_scanned += 1
            if nr_y_scanned < nr_scans:
                device.shell("input keyevent KEYCODE_DPAD_UP")
                pygame.time.wait(wait)
                y += 1
        if nr_y_scanned == nr_scans:
            measure = 0
            np.save(name, s21s)
            print('Saved: %s' % name)
            device.shell("input keyevent KEYCODE_B")
            pygame.time.wait(wait)
    
    event = pygame.event.poll()
    if event.type == pygame.QUIT:
        running = 0 

    if event.type == pygame.KEYDOWN:
        if event.key == pygame.K_ESCAPE:
            running = 0
        if event.key == pygame.K_BACKSPACE:
            device.shell("input keyevent 67")
            pygame.time.wait(wait)
            x = x_init
            y = y_init
            w = w_init
            dx = w_init
            dy = w_init
            nr_s21 = 0
        if event.key == pygame.K_DOWN:
            device.shell("input keyevent KEYCODE_DPAD_DOWN")
            pygame.time.wait(wait)
            y -= dy
            if y < y_min:
                y = y_min
        if event.key == pygame.K_UP:
            device.shell("input keyevent KEYCODE_DPAD_UP")
            pygame.time.wait(wait)
            y += dy
            if y > y_max:
                y = y_max
        if event.key == pygame.K_RIGHT:
            device.shell("input keyevent KEYCODE_DPAD_RIGHT")
            pygame.time.wait(wait)
            x += dx
            if x > x_max:
                x = x_max
        if event.key == pygame.K_LEFT:
            device.shell("input keyevent KEYCODE_DPAD_LEFT")
            pygame.time.wait(wait)
            x -= dx
            if x < x_min:
                x = x_min
        if event.key == pygame.K_e:
            device.shell("input keyevent KEYCODE_E")
            pygame.time.wait(wait)
            w += step
        if event.key == pygame.K_q:
            device.shell("input keyevent KEYCODE_Q")
            pygame.time.wait(wait)
            w -= step
            if w < w_min:
                w = w_min
        if event.key == pygame.K_c:
            device.shell("input keyevent KEYCODE_C")
            pygame.time.wait(wait)
            dx += step
        if event.key == pygame.K_z:
            device.shell("input keyevent KEYCODE_Z")
            pygame.time.wait(wait)
            dx -= step
            if dx < dx_min:
                dx = dx_min
        if event.key == pygame.K_u:
            device.shell("input keyevent KEYCODE_U")
            pygame.time.wait(wait)
            dy += step
        if event.key == pygame.K_t:
            device.shell("input keyevent KEYCODE_T")
            pygame.time.wait(wait)
            dy -= step
            if dy < dy_min:
                dy = dy_min
        if event.key == pygame.K_x:
            device.shell("input keyevent KEYCODE_X")
            pygame.time.wait(wait)
            dx = 1
        if event.key == pygame.K_y:
            device.shell("input keyevent KEYCODE_Y")
            pygame.time.wait(wait)
            dy = 1
        if event.key == pygame.K_w:
            device.shell("input keyevent KEYCODE_W")
            pygame.time.wait(wait)
            w = 1
        if event.key == pygame.K_i:
            device.shell("input keyevent KEYCODE_I")
            pygame.time.wait(wait)
            if not dark:
                current_linecolor = copy(linecolor)
                current_bgcolor = copy(bgcolor)
                bgcolor = current_linecolor
                linecolor = current_bgcolor
                axcolor = current_bgcolor
            if inverted:
                inverted = 0
            else:
                inverted = 1
        if event.key == pygame.K_g:
            device.shell("input keyevent KEYCODE_G")
            if not dark:
                if inverted:
                    color_cycler += 1
                    bgcolor = colors[color_cycler % nr_colors]
                    linecolor = black
                    axcolor = linecolor
                else:
                    color_cycler += 1
                    bgcolor = black
                    linecolor = colors[color_cycler % nr_colors]
                    axcolor = linecolor
        if event.key == pygame.K_b:
            device.shell("input keyevent KEYCODE_B")
            pygame.time.wait(wait)
            if not inverted:
                if dark:
                    linecolor = current_linecolor
                    axcolor = current_axcolor
                    dark = 0
                else:
                    current_linecolor = copy(linecolor)
                    current_axcolor = copy(axcolor)
                    bgcolor = black
                    linecolor = black
                    dark = 1
        if event.key == pygame.K_RETURN:
                # Set screen and stepsizes to correct values
                if inverted:
                    print('Putting inverted screen back to normal')
                    device.shell("input keyevent KEYCODE_I")
                    pygame.time.wait(wait)
                    inverted = 0
                if dx != 1:
                    print('Putting dX to 1')
                    device.shell("input keyevent KEYCODE_X")
                    pygame.time.wait(wait)
                    dx = 1
                if dy != 1:
                    print('Putting dY to 1')
                    device.shell("input keyevent KEYCODE_Y")
                    pygame.time.wait(wait)
                    dy = 1

                # Ask input on number of scans
                nr_scans = int(input('Please input the number of scans in x and y: '))
                
                # Initiate array
                dir = 'S21s/17-10-23/'
                s21s = np.zeros((nr_scans, 2, len_s21))
                date = timestamp()
                name = '%sS21_w%d_%s.npy' % (dir, w, date)
                freqsname = '%sS21_w%d_%s_freqs.npy' % (dir, w, date)
                darkname = '%sS21_w%d_%s_dark.npy' % (dir, w, date)
                settingsname = '%sS21_w%d_%s_settings.txt' % (dir, w, date)

                dict = {'color':colors[color_cycler % nr_colors], 
                        'fstart':realfstart, 'fstop':realfstop, 'subscanbw':subscanbw, 
                        'kidpower':kidpower, 'ifbw':ifbw}
                with open(settingsname, 'w') as file:
                    json.dump(dict, file)

                # Make dark scan and save it
                device.shell("input keyevent KEYCODE_B")
                pygame.time.wait(wait)
                freqs, dark_s21 = f.get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw)
                np.save(darkname, dark_s21)
                np.save(freqsname, freqs)
                print('Saved: %s' % darkname)
                fig, ax = plt.subplots()
                ax.plot(freqs, dark_s21)
                plt.show()
                measure = 1
    pygame.display.flip()
pygame.quit()
sys.exit()

