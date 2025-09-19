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


## Input S21 parameters
fstart = 4 # GHz
fstop = 8.3  # GHz
totscanbw = fstop - fstart
num_points = 6401
subscanbw = 100  # MHz
num_subscans = int(np.ceil(totscanbw / subscanbw))
realfstart = fstart
realfstop = fstart + num_subscans * (subscanbw*1e-3)
len_s21 = int(num_subscans * num_points)
kidpower = -110 # dBm
ifbw = 1000  # Hz
freqs = np.linspace(realfstart, realfstop, num_points*num_subscans)
date = datetime.today()
xsteps = []   
ysteps = []  


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


# Constants
pygame.init()
colors = {
    'white': (255, 255, 255),
    'blue': (0, 0, 255),
    'green': (0, 255, 0),
    'red': (255, 0, 0),
    'black': (0, 0, 0)
}
colorkeys = ['white', 'blue', 'green', 'red']
colorkey_cycler = itertools.cycler(colorkeys)
linecolor = next(colorkey_cycler)
bgcolor = 'black'

width, height = 600, 500
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
running = True
nr_y2scan = len(ysteps)
nr_x2scan = len(xsteps)
nr_xscanned = 0
nr_yscanned = 0
w = w_init
x = x_init
y = y_init
dx = dx_init
dy = dy_init
restart = True
false_true = [False, True]
inverted_cycler = itertools.cycler(false_true)
inverted = next(inverted_cycler)
measure = False
wait = 1
linetypes = ['both', 'none', 'x', 'y']
linectype_cycler = itertools.cycle(linetypes)
linetype = next(linectype_cycler)
square_cycler = itertools.cycler(false_true)
square = next(square_cycler)
scanline = False
datadir = None

# ## Test connection to Virtual Intstruments
try:
    vna = f.connect2vi("GPIB0::16::INSTR", timeout=3000000)
    weinschell = f.connect2vi("GPIB0::10::INSTR", timeout=300000)
except:
    print('No VNA connected')

while running:
    if restart:
        device.shell("input keyevent 67")
        pygame.time.wait(wait)
        restart = 0

    text_list = [
    'X, Y = ' + str(x) + ', ' + str(y) + '\n Move = up,down,left,right',
    '# scans = ' + str(nr_yscanned + nr_xscanned),
    'dX = ' + str(dx) + '\n -1/=1/+1 = Z/X/C',
    'dY = ' + str(dy) + '\n -1/=1/+1 = T/Y/U',
    'w = ' + str(w) + '\n -1/=1/+1 = Q/W/E',
    '# steps = ' + str(len(xsteps)) + ', ' + str(len(ysteps)) + '\n add step = L',
    'invert screen = I',
    'Change color = G \nCycles White-Blue-Green-Red-Off',
    'reset lines = Backspace',
    'quit and save = Esc',
    ]
    screen.fill(colors[bgcolor])
    draw_text_boxes(text_list, colors[bgcolor], colors[linecolor])
    
    if measure:
        if nr_xscanned < nr_x2scan:
            if nr_xscanned == 0:
                device.shell("input keyevent KEYCODE_B")
                pygame.time.wait(wait)
            xstep = xsteps[nr_xscanned]
            print('Scan %d/%d, stepping %d' % (nr_xscanned+1, nr_x2scan, xstep))
            for i in range(xstep):
                device.shell("input keyevent KEYCODE_DPAD_RIGHT")
                pygame.time.wait(wait)
                x += 1
            freqs, s21 = f.get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw)
            name = '%s/S21_x%02d.npy' % (datadir, nr_xscanned)
            np.save(name, np.stack((freqs, s21), axis=-1).T)
            nr_xscanned += 1
        elif nr_yscanned < nr_y2scan and nr_xscanned == nr_x2scan:
            if nr_yscanned == 0:
                device.shell("input keyevent KEYCODE_B")
                pygame.time.wait(wait)
            ystep = ysteps[nr_yscanned]
            print('Scan %d/%d, stepping %d' % (nr_yscanned+1, nr_y2scan, ystep))
            for i in range(ystep):
                device.shell("input keyevent KEYCODE_DPAD_UP")
                pygame.time.wait(wait)
                y += 1
            freqs, s21 = f.get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw)
            name = '%s/S21_y%02d.npy' % (datadir, nr_yscanned)
            np.save(name, np.stack((freqs, s21), axis=-1).T)
            nr_yscanned += 1
        elif nr_xscanned == nr_x2scan and nr_yscanned == nr_y2scan:
            measure = 0
            print('Helemaal f*cking klaar met de meting')
            device.shell("input keyevent KEYCODE_B")
            pygame.time.wait(wait)
        else:
            measure = 0
            device.shell("input keyevent KEYCODE_B")
            print('WARNING: Measurement stopped but might not be complete')

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
            colorkey_cycler = itertools.cycler(colorkeys)
            bgcolor = 'black'
            linecolor = next(colorkey_cycler)
            linetype_cycler = itertools.cycle(linetypes)
            linetype = next(linetype_cycler)
            inverted_cycler = itertools.cycler(false_true)
            inverted = next(inverted_cycler)
            square_cycler = itertools.cycler(false_true)
            square = next(square_cycler)

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
            current_linecolor = copy(linecolor)
            current_bgcolor = copy(bgcolor)
            bgcolor = current_linecolor
            linecolor = current_bgcolor
            inverted = next(inverted_cycler)
        if event.key == pygame.K_g:
            device.shell("input keyevent KEYCODE_G")
            pygame.time.wait(wait)
            if inverted:
                bgcolor = next(colorkey_cycler)
                linecolor = 'black'
            else:
                bgcolor = 'black'
                linecolor = next(colorkey_cycler)
        if event.key == pygame.K_b:
            if not square:
                device.shell("input keyevent KEYCODE_B")
                pygame.time.wait(wait)
                linetype = next(linectype_cycler)
        if event.key == pygame.K_l:
            if linetype == 'x':
                if not len(xsteps):
                    xsteps.append(0)
                    xstart = copy(int(x))
                else:
                    xsteps.append(int(x - xstart))
                    xstart = copy(x)
            elif linetype == 'y':
                if not len(ysteps):
                    ysteps.append(0)
                    ystart = copy(int(y))
                else:
                    ysteps.append(int(y - ystart))
                    ystart = copy(y)
                print('Stpes in X: ', xsteps)
                print('steps in Y: ', ysteps)
        if event.key == pygame.K_s:
            device.shell("input keyevent KEYCODE_S")
            pygame.time.wait(wait)
            square = next(square_cycler)
        if event.key == pygame.K_m:
            datadir = f.select_directory(initial_dir=datadir)
            if datadir:
                print(f"Selected directory: {datadir}")
            freqs, s21 = f.get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw)
            date = f.timestamp()
            name = '%s/S21_x%dy%d_w%d.npy' % (datadir, x, y, w)
            np.save(name, np.stack((freqs, s21), axis=-1).T)
            print('Saved: %s' % (name))
        if event.key == pygame.K_RETURN:
                measure = 1

                print('# steps in X: ', nr_x2scan)
                print('# steps in Y: ', nr_y2scan)
                print('linewidth: ', w)
                print('color: ', linecolor)
                maindir = f.select_directory()
                if maindir:
                    print(f"Selected directory: {maindir}")
                date = f.timestamp()
                datadir = maindir + '/S21s' + date
                f.create_directory(datadir)
                
                # Initiate array
                freqsname = '%s/freqs.npy' % (datadir)
                darkname = '%s/S21_dark.npy' % (datadir)
                settingsname = '%s/settings.txt' % (datadir)
                dict = {'color': linecolor, 'width':w,
                        'fstart':realfstart, 'fstop':realfstop, 'subscanbw':subscanbw, 
                        'kidpower':kidpower, 'ifbw':ifbw, 'nr points':num_points}
                with open(settingsname, 'w') as file:
                    json.dump(dict, file)

                # Make dark scan and save it
                device.shell("input keyevent KEYCODE_B")
                pygame.time.wait(wait)
                st = time.time()
                freqs, dark_s21 = f.get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw)
                et = time.time()
                scan_time = et - st
                print('Time 1 scan = %d seconds' % scan_time)
                print('Expected duration measurement: %d minutes' % ((nr_x2scan + nr_y2scan) * scan_time / 60))
                np.save(darkname, np.stack((freqs, dark_s21), axis=-1).T)
                print('Saved: %s' % darkname)
                
                fig, ax = plt.subplots()
                ax.plot(freqs, dark_s21)
                plt.show()
    pygame.display.flip()
pygame.quit()
sys.exit()

