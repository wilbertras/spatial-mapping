import pygame
import numpy as np

y = 0
x = 0
w = 5
dx = 1
dy = 1
width = 800
height = 600
x_range = np.arange(x, width, dx)
y_range = np.arange(y, height, dy)

dir = 1
dT = 1
running = 1
screen = pygame.display.set_mode((width, height))
linecolor = 0, 255, 0
bgcolor = 0, 0, 0
while running:
    event = pygame.event.poll()
    if event.type == pygame.QUIT:
        running = 0

    if event.type == pygame.KEYDOWN:
        if event.key == pygame.K_DOWN:
            for y in y_range:
                screen.fill(bgcolor)
                pygame.draw.line(screen, linecolor, (0, y), (width-1, y), width=w)
                pygame.display.flip()
                pygame.time.wait(dT)

        if event.key == pygame.K_RIGHT:
            for x in x_range:
                screen.fill(bgcolor)
                pygame.draw.line(screen, linecolor, (x, 0), (x, height-1), width=w)
                pygame.display.flip()
                pygame.time.wait(dT)

    screen.fill(bgcolor)
    pygame.display.flip()