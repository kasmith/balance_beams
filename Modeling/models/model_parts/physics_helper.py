from __future__ import division, print_function
import pymunk as pm
import numpy as np
import pygame as pg
from pygame.constants import QUIT
from geometry import rescale_and_recenter_system, convex_area, rotate_system, \
    mirror_system, find_concave_outline, recenter_system
from numpy import sqrt

COLTYPE_BEAM = 101
COLTYPE_LEFT = 102
COLTYPE_RIGHT = 103

DEF_FRIC = 0.5
DEF_AREA = 0.32

INITHEIGHT = 0.711
SPACING = 0.01

OLINE_WID = 2

DENSDICT = {"W": 1.,
            "B": 2.5,
            "I": 9.}

COLORDICT = {"W": (182, 155, 76, 255),
             "B": (150, 22, 11, 255),
             "I": (67, 75, 77, 255)}

OUTLINE_COLORDICT = {"W": (130, 82, 1, 255),
                     "B": (97, 0, 0, 255),
                     "I": (55, 56, 56, 255)}


def _position_standard_block(lower_right, edge_len=1):
    lr = np.array(lower_right)
    standard = [np.array([0, 0]), np.array([0, edge_len]),
                np.array([edge_len, edge_len]), np.array([edge_len, 0])]
    return [s + lr for s in standard]


def _get_pm_vertices(poly):
    return [np.array(v.rotated(poly.body.angle) +
                     poly.body.position) for v in poly.get_vertices()]

def _just_rescale(system, tot_area):
    new_sys, center = rescale_and_recenter_system(system, tot_area)
    return [[v + center for v in shape] for shape in new_sys]

def make_block(scale=1):
    tot_area = scale * DEF_AREA
    system = [_position_standard_block([-1, 0]), _position_standard_block([0, 0])]
    return _just_rescale(system, tot_area)


def make_box(scale=1):
    tot_area = scale * DEF_AREA
    system = [_position_standard_block([-1., -1.], 2.)]
    return _just_rescale(system, tot_area)


def make_u(scale=1):
    tot_area = scale * DEF_AREA
    system = [_position_standard_block([-3., -2.], 2.),
              _position_standard_block([-3., 0.], 2.),
              _position_standard_block([-1., -2.], 2.),
              _position_standard_block([1., -2.], 2.),
              _position_standard_block([1., 0.], 2.)]
    return _just_rescale(system, tot_area)


def make_l(scale=1):
    tot_area = scale * DEF_AREA
    system = [_position_standard_block([-3.5, -1.5], 2.),
              _position_standard_block([-1.5, -1.5], 2.),
              _position_standard_block([0.5, -1.5], 2.),
              _position_standard_block([0.5, 0.5], 2.)]
    return _just_rescale(system, tot_area)


def make_t(scale=1):
    tot_area = scale * DEF_AREA
    system = [_position_standard_block([-3., -5./3.], 2.),
              _position_standard_block([-1., -5./3.], 2.),
              _position_standard_block([1., -5./3.], 2.),
              _position_standard_block([-1., 1./3.], 2.)]
    return _just_rescale(system, tot_area)


def make_i(scale=1):
    tot_area = scale * DEF_AREA
    system = [_position_standard_block([-3., -1.], 2.),
              _position_standard_block([-1., -1.], 2.),
              _position_standard_block([1., -1.], 2.)]
    return _just_rescale(system, tot_area)


def make_thumb(scale=1):
    tot_area = scale * DEF_AREA
    system = [_position_standard_block([-3, -20./11.], 2.),
              _position_standard_block([-1, -20./11.], 2.),
              _position_standard_block([1., 20./11.], 2.),
              _position_standard_block([-3., 2./11.], 2.),
              _position_standard_block([-1., 2./11.], 2.)]
    return _just_rescale(system, tot_area)


def make_tri(scale=1):
    tot_area = scale * DEF_AREA
    system = [[[-2, -0.5], [0, 1.5], [2, -0.5]]]
    return _just_rescale(system, tot_area)


def _get_system_heights(system):
    min_y = np.min([np.min([p[1] for p in vs]) for vs in system])
    max_y = np.max([np.max([p[1] for p in vs]) for vs in system])
    return -min_y, (max_y - min_y)


class PhysicsWorld(object):

    def __init__(self, com_pos, com_width, beam_mass=None, timestep=.02):
        self.sp = pm.Space()
        self.sp.gravity = 0, -9.81

        self.ts = timestep
        self.objects = []

        # Set up the strut & beam
        hw = com_width
        self.strut_body = pm.Body(1, 1, pm.Body.STATIC)
        self.strut_shape = pm.Poly(self.strut_body, [(-hw, .6), (hw, .6),
                                                     (hw, 0), (-hw, 0)])
        self.strut_shape.color = COLORDICT["W"]
        self.strut_shape.outline = OUTLINE_COLORDICT["W"]
        self.strut_shape.friction = .5
        self.strut_body.position = (com_pos, 0)
        self.sp.add(self.strut_shape)

        if beam_mass is None:
            beam_mass = 1.024 * DENSDICT['W']
        beam_moi = pm.moment_for_box(beam_mass, (10., .1))
        self.beam_body = pm.Body(beam_mass, beam_moi)
        self.beam_shape = pm.Poly(self.beam_body, [(-5.12, .05), (5.12, .05),
                                                   (5.12, -.05), (-5.12, -.05)])
        self.beam_body.position = 0, .661
        self.beam_shape.collision_type = COLTYPE_BEAM
        self.beam_shape.color = COLORDICT["W"]
        self.beam_shape.outline = OUTLINE_COLORDICT["W"]
        self.beam_shape.friction = .5
        self.sp.add(self.beam_body)
        self.sp.add(self.beam_shape)

        self.ground_shape = pm.Poly(self.sp.static_body,
                                    [[-10., 0.], [10., 0.],
                                     [10, -.5], [-10., -.5]])
        self.ground_shape.friction = 0.5
        self.sp.add(self.ground_shape)

        # Set up the sensor shapes
        self.sensor_left = pm.Poly(self.sp.static_body, [(1, .01), (10, .01), (10, 0.), (1, 0)])
        self.sensor_left.sensor = True
        self.sensor_left.collision_type = COLTYPE_LEFT
        self.sp.add(self.sensor_left)
        self.flag_left = False

        self.sensor_right = pm.Poly(self.sp.static_body, [(-10, .01), (-1, .01), (-1, 0), (-10, 0)])
        self.sensor_right.sensor = True
        self.sensor_right.collision_type = COLTYPE_RIGHT
        self.sp.add(self.sensor_right)
        self.flag_right = False

        def l_h(space, arbiter, x):
            self.flag_left = True
            return 0

        def r_h(space, arbiter, x):
            self.flag_right = True
            return 0

        self.sp.add_collision_handler(COLTYPE_BEAM, COLTYPE_LEFT).begin = l_h
        self.sp.add_collision_handler(COLTYPE_BEAM, COLTYPE_RIGHT).begin = r_h

        # Holder for perception through physics model
        self.perc = None

    def run(self, tlim=3.):
        t = 0.
        while t < tlim:
            t += self.ts
            self.sp.step(self.ts)
            if self.flag_left:
                return 'L'
            if self.flag_right:
                return 'R'
        return 'B'

    def draw(self, dims=(700, 300), exist_surf=None):
        if exist_surf is None:
            surf = pg.Surface(dims)
        else:
            surf = exist_surf
            dims = surf.get_size()
        surf.fill(pg.Color('white'))

        def rescale_point(p):
            x = (-p[0] + 7.) / 14. * dims[0]
            y = (p[1] + 0.5) / 6. * dims[1]
            y = dims[1] - y  # Invert for drawing
            return x, y

        ground_pt = rescale_point((7, 0))
        ground_r = pg.Rect(0, ground_pt[1], dims[0], dims[1])
        pg.draw.rect(surf, pg.Color('black'), ground_r)

        strut_verts = _get_pm_vertices(self.strut_shape)
        strut_draw = [rescale_point(v) for v in strut_verts]
        pg.draw.polygon(surf, self.strut_shape.color, strut_draw)
        pg.draw.polygon(surf, self.strut_shape.outline, strut_draw, OLINE_WID)
        beam_verts = _get_pm_vertices(self.beam_shape)
        beam_draw = [rescale_point(v) for v in beam_verts]
        pg.draw.polygon(surf, self.beam_shape.color, beam_draw)
        pg.draw.polygon(surf, self.beam_shape.outline, beam_draw, OLINE_WID)

        for o in self.objects:
            shapes = o[1]
            shape_verts = []
            o_col = shapes[0].outline
            for s in o[1]:
                s_verts = _get_pm_vertices(s)
                s_draw = [rescale_point(v) for v in s_verts]
                shape_verts.append(s_draw)
                pg.draw.polygon(surf, s.color, s_draw)
            outline = find_concave_outline(shape_verts)
            pg.draw.polygon(surf, o_col, outline, OLINE_WID)

        return surf

    def demonstrate(self, pgscreen=None):
        hz = 1 / self.ts
        self.draw(exist_surf=pgscreen)
        pg.display.flip()
        running = True
        clk = pg.time.Clock()
        clk.tick(hz)
        while running:
            self.sp.step(self.ts)
            self.draw(exist_surf=pgscreen)
            pg.display.flip()
            for e in pg.event.get():
                if e.type == QUIT:
                    running = False
            clk.tick(hz)
        pg.quit()

    def _make_pm_object(self, vertexlist, position, scale, material):
        tot_mass = scale * DENSDICT[material]
        vertexlist = recenter_system(vertexlist)[0]
        shape_areas = np.array([convex_area(s) for s in vertexlist])
        shape_areas /= sum(shape_areas)
        # Get the moments out of the shapes
        moments = [pm.moment_for_poly(a * tot_mass, s) for s, a in
                   zip(vertexlist, shape_areas)]
        pm_body = pm.Body(tot_mass, sum(moments))
        # Create the shapes
        pm_shapes = []
        for s in vertexlist:
            new_s = pm.Poly(pm_body, s)
            new_s.color = COLORDICT[material]
            new_s.outline = OUTLINE_COLORDICT[material]
            new_s.friction = 0.5
            pm_shapes.append(new_s)
        # Move the body
        pm_body.position = position
        self.sp.add(pm_body)
        self.sp.add(pm_shapes)
        return pm_body, pm_shapes

    def _do_rotation(self, shapes, rottype):
        ymirror = rottype >= 4
        xmirror = (rottype % 4) == 2
        shapes = mirror_system(shapes, [xmirror, ymirror])
        return shapes

    def new_block(self, xloc, zbottom=INITHEIGHT, scale=1., rottype=0, nexttx=0.,
                  material='W', rev=False):
        if rev:
            rottype = (rottype + 4) % 8
        shapes = make_block(scale)
        shapes = self._do_rotation(shapes, rottype)
        c_height, obj_height = _get_system_heights(shapes)
        pos = np.array([xloc, zbottom + c_height])
        self.objects.append(self._make_pm_object(shapes, pos, scale, material))
        return zbottom + obj_height

    def new_box(self, xloc, zbottom=INITHEIGHT, scale=1., rottype=0, nexttx=0.,
                material='W', rev=False):
        shapes = make_box(scale)
        shapes = self._do_rotation(shapes, rottype)
        c_height, obj_height = _get_system_heights(shapes)
        pos = np.array([xloc, zbottom + c_height])
        self.objects.append(self._make_pm_object(shapes, pos, scale, material))
        return zbottom + obj_height

    def new_u(self, xloc, zbottom=INITHEIGHT, scale=1., rottype=0, nexttx=0.,
              material='W', rev=False):
        if rev:
            rottype = (rottype + 4) % 8
        shapes = make_u(scale)
        shapes = self._do_rotation(shapes, rottype)
        c_height, obj_height = _get_system_heights(shapes)
        pos = np.array([xloc, zbottom + c_height])
        self.objects.append(self._make_pm_object(shapes, pos, scale, material))
        return zbottom + obj_height

    def new_l(self, xloc, zbottom=INITHEIGHT, scale=1., rottype=0, nextx=0.,
              material='W', rev=False):
        shapes = make_l(scale)
        if rev:
            if rottype % 4 > 1:
                rottype -= 2
            else:
                rottype += 2
        inv = rottype >= 4
        if rottype >= 4:
            rottype -= 4
            shapes = mirror_system(shapes, [False, True])
        shapes = self._do_rotation(shapes, rottype)
        c_height, obj_height = _get_system_heights(shapes)
        if inv:
            if nextx <= -.5 and rottype == 0 or nextx >= .5 and rottype == 2:
                newz = zbottom + obj_height / 3. + SPACING
            else:
                newz = zbottom + obj_height + SPACING
        else:
            if nextx <= -.5 and rottype == 0 or nextx >= .5 and rottype == 2:
                newz = zbottom + obj_height + SPACING
            else:
                newz = zbottom + obj_height / 2. + SPACING
        pos = np.array([xloc, zbottom + c_height])
        print(pos)
        self.objects.append(self._make_pm_object(shapes, pos, scale, material))
        return newz

    def new_t(self, xloc, zbottom=INITHEIGHT, scale=1., rottype=0, nexttx=0.,
              material='W', rev=False):
        shapes = make_t(scale)
        shapes = self._do_rotation(shapes, rottype)
        if rev:
            rotate_system(shapes, [True, False])
        c_height, obj_height = _get_system_heights(shapes)
        pos = np.array([xloc, zbottom + c_height])
        self.objects.append(self._make_pm_object(shapes, pos, scale, material))
        return zbottom + obj_height

    def new_i(self, xloc, zbottom=INITHEIGHT, scale=1., rottype=0, nexttx=0.,
              material='W', rev=False):
        shapes = make_i(scale)
        if rottype >= 4:
            rottype -= 4
            shapes = rotate_system(shapes, -np.pi / 2.)
        c_height, obj_height = _get_system_heights(shapes)
        pos = np.array([xloc, zbottom + c_height])
        self.objects.append(self._make_pm_object(shapes, pos, scale, material))
        return zbottom + obj_height

    def new_thumb(self, xloc, zbottom=INITHEIGHT, scale=1., rottype=0, nexttx=0.,
              material='W', rev=False):
        if rev:
            if rottype % 4 > 1:
                rottype -= 2
            else:
                rottype += 2
        shapes = make_thumb(scale)
        if rottype >= 4:
            rottype -= 4
            shapes = rotate_system(shapes, np.pi / 2.)
        if (rottype == 2):
            shapes = mirror_system(shapes, [True, False])
        c_height, obj_height = _get_system_heights(shapes)
        pos = np.array([xloc, zbottom + c_height])
        self.objects.append(self._make_pm_object(shapes, pos, scale, material))
        return zbottom + obj_height

    def new_hat(self, xloc, zbottom=INITHEIGHT, scale=1., rottype=0, nexttx=0.,
              material='W', rev=False):
        if rev:
            if rottype % 4 > 1:
                rottype -= 2
            else:
                rottype += 2
        if rottype >= 8:
            rotang = -np.pi*.75
            rottype -= 8
        elif rottype >= 4:
            rotang = np.pi*.75
            rottype -= 4
        else:
            rotang = 0
        shapes = make_tri(scale)
        shapes = rotate_system(shapes, rotang)
        if rottype == 2:
            shapes = mirror_system(shapes, [True, False])
        c_height, obj_height = _get_system_heights(shapes)
        pos = np.array([xloc, zbottom + c_height])
        self.objects.append(self._make_pm_object(shapes, pos, scale, material))
        return zbottom + obj_height

    def comb_2a(self, xloc, scale=1., h=INITHEIGHT, mat="W", rev=False):
        h1 = self.new_box(xloc, h, .75*scale, 6, 0, mat)
        self.new_hat(xloc, h1, .25*scale, 10, 0, mat)

    def comb_2b(self, xloc, scale=1., h=INITHEIGHT, mat="W", rev=False):
        h1 = self.new_u(xloc, h, .5*scale, 2, -1, mat)
        self.new_thumb(xloc, h1, .5*scale, 2, 0, mat, rev)

    def comb_2c(self, xloc, scale=1., h=INITHEIGHT, mat="W", rev=False):
        o1 = -.05
        o2 = .15
        if not rev:
            o1 *= -1
            o2 *= -1
        sscl = scale / sqrt(2)
        h1 = self.new_l(xloc+o1*sqrt(sscl),h,.75*scale,0,.2,mat,rev=rev)
        self.new_hat(xloc+o2*sqrt(sscl),h1,.25*scale,8,0,mat,rev=rev)

    def comb_2d(self, xloc, scale=1., h=INITHEIGHT, mat="W", rev=False):
        h1 = self.new_box(xloc,h,.5*scale,6,0,mat)
        self.new_i(xloc,h1,.5*scale,6,0,mat)

    def comb_2e(self, xloc, scale=1., h=INITHEIGHT, mat="W", rev=False):
        xL = .252/12. + 1./8. + .05
        xT = -1./8. - .05
        if not rev:
            xL = -xL
            xT = -xT
        sscl = scale / sqrt(2)
        h1 = self.new_l(xloc+xL*sqrt(sscl),h,.5*scale,2,0,mat,rev=rev)
        self.new_t(xloc+xT*sqrt(sscl),h1,.5*scale,6,0,mat,rev=rev)

    def comb_2f(self, xloc, scale=1., h=INITHEIGHT, mat="W", rev=False):
        h1 = self.new_box(xloc,h,2.*scale/3.,6,0,mat)
        self.new_hat(xloc,h1,scale/3.,0,0,mat)

    def comb_2g(self, xloc, scale=1., h=INITHEIGHT, mat="W", rev=False):
        h1 = self.new_i(xloc,h,scale/2.,4,0,mat)
        self.new_i(xloc,h1,scale/2.,6,0,mat)

    def new_shapes_object(self, otype, xloc, zbottom=INITHEIGHT, scale=1.,
                          rottype=0, nexttx=0., material='W', rev=False):
        if otype == 'B':
            return self.new_box(xloc, zbottom, scale, rottype, nexttx, material, rev)
        elif otype == 'I':
            return self.new_i(xloc, zbottom, scale, rottype, nexttx, material, rev)
        elif otype == 'H':
            return self.new_hat(xloc, zbottom, scale, rottype, nexttx, material, rev)
        elif otype == 'U':
            return self.new_u(xloc, zbottom, scale, rottype, nexttx, material, rev)
        elif otype == 'Th':
            return self.new_thumb(xloc, zbottom, scale, rottype, nexttx, material, rev)
        elif otype == 'T':
            return self.new_t(xloc, zbottom, scale, rottype, nexttx, material, rev)
        elif otype == 'L':
            return self.new_l(xloc, zbottom, scale, rottype, nexttx, material, rev)
        elif otype == 'C':
            selfcall = getattr(self, "comb_"+rottype)
            selfcall(xloc, scale, zbottom, material, rev)

    def readjust_strut(self, com_u=None):
        if com_u:
            cuadj = np.random.normal(0, com_u)
            self.sp.remove(self.strut_shape)

            old_vs = self.strut_shape.get_vertices()
            new_vs = map(lambda p: (p[0] + cuadj, p[1]), old_vs)

            self.strut_shape = pm.Poly(self.strut_body, new_vs)
            self.strut_shape.color = (0, 0, 0)
            self.strut_shape.friction = .5
            self.sp.add(self.strut_shape)

            return cuadj
        return None


if __name__ == '__main__':
    print(make_box(2))
