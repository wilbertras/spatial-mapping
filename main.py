import pygame
import numpy as np


def scan():
    w = 1
    change = 1
    dx = 1
    width = 1080
    height = 2340
    y_start = int(height / 2) - int(width / 2)
    y_stop = y_start + width - 1
    x_start = 0
    x_stop = x_start + width - 1
    dir = 1
    t = 0
    running = 1
    screen = pygame.display.set_mode((width, height))
    linecolor = 0, 255, 0
    bgcolor = 0, 0, 0
    while running:
        event = pygame.event.poll()
        if event.type == pygame.QUIT:
            running = 0

        x_range = np.arange(x_start, x_stop, dx)
        y_range = np.arange(y_start, y_stop, dx)

        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_ESCAPE:
                running = 0
            if event.key == pygame.K_DOWN:
                for y in y_range:
                    screen.fill(bgcolor)
                    pygame.draw.line(screen, linecolor, (0, y), (width-1, y), width=w)
                    pygame.display.flip()
                    pygame.time.wait(t)
                    new_event = pygame.event.poll()
                    if new_event.type == pygame.KEYDOWN:
                        if new_event.key == pygame.K_SPACE:
                            break


            if event.key == pygame.K_RIGHT:
                for x in x_range:
                    screen.fill(bgcolor)
                    pygame.draw.line(screen, linecolor, (x, y_start), (x, y_stop-1), width=w)
                    pygame.display.flip()
                    pygame.time.wait(t)
                    new_event = pygame.event.poll()
                    if new_event.type == pygame.KEYDOWN:
                        if new_event.key == pygame.K_SPACE:
                            break

            if event.key == pygame.K_e:
                w += change
            if event.key == pygame.K_q:
                w -= change
                if w <= 1:
                    w = 1
            if event.key == pygame.K_y:
                t += change
            if event.key == pygame.K_r:
                t -= change
                if t <= 0:
                    t = 0
            if event.key == pygame.K_c:
                dx += change
            if event.key == pygame.K_z:
                dx -= change
                if dx <= 1:
                    dx = 1
        screen.fill(bgcolor)
        pygame.display.flip()

scan()