from __future__ import division, print_function
from physics_helper import PhysicsWorld, INITHEIGHT
import pygame as pg
from pygame.constants import QUIT
import os
import json
import numpy as np
from geometry import find_concave_outline


def read_shapes_dict(d):
    world = PhysicsWorld(0, 0.05)
    sa = np.sqrt(2)
    for i in range(len(d['Right'])):
        stack = d['Right'][i]
        xpos = -(i+1)
        hght = INITHEIGHT
        for block in stack:
            hght = world.new_shapes_object(block[0], xpos, hght, block[1]*sa,
                                           block[2], 0., 'W', False)
    for i in range(len(d['Left'])):
        stack = d['Left'][i]
        xpos = i+1
        hght = INITHEIGHT
        for block in stack:
            hght = world.new_shapes_object(block[0], xpos, hght, block[1]*sa,
                                           block[2], 0., 'W', False)
    return world



#world = PhysicsWorld(0, 0.05)
#h1 = world.new_i(-2)
#print (world.run())
flpth = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                     '..','..','..','Scenes','TrJSON_Shapes','CB_1_BC_N.json')
world = read_shapes_dict(json.load(open(flpth, 'rU')))
for b in world.sp.bodies:
    if b.position.x < -3:
        vlist = [[np.array(b.local_to_world(v)) for v in s.get_vertices()] for s in b.shapes]
        print (find_concave_outline(vlist))
#world.new_box(-1,INITHEIGHT)

pg.init()
d = pg.display.set_mode((700,300))
world.demonstrate(d)
'''
d.blit(s, (0,0))
pg.display.flip()

running = True
while running:
    for e in pg.event.get():
        if e.type == QUIT:
            running = False

pg.quit()
'''
