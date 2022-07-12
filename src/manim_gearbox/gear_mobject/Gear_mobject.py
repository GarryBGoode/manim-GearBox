import numpy as np
from manim import *
from typing import Optional, Sequence, Union
from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.optimize import fmin
from scipy.optimize import fmin_powell

__all__ = [
    "involute_func",
    "involute_deriv_func",
    "involute_height_func",
    "involute_point_gen",
    "Gear",
    "Rack"
]



def involute_func(t, r, a=0, rad_offs=0,tan_offs=0):
    '''
    Returns the x-y-z values of the involute function.
    t: input angle
    r: base circle radius
    a: offset angle
    '''
    def involute_val(val):
        offs_v = rotate_vector(np.array([rad_offs,tan_offs,0]),val)
        # x = r * (np.cos(val) + (val-a) * np.sin(val-a)) + offs_v[0]
        # y = r * (np.sin(val) - (val-a) * np.cos(val-a)) + offs_v[1]
        x = r * (np.cos(val) + (val - a) * np.sin(val - a)) + \
            rad_offs * np.cos(val) - tan_offs*np.sin(val)
        y = r * (np.sin(val) - (val - a) * np.cos(val - a)) + \
            rad_offs * np.sin(val) + tan_offs*np.cos(val)
        z = 0
        return np.array((x,y,z))
    if hasattr(t, '__iter__'):
        ret = np.empty((0,3))
        for u in t:
            point = involute_val(u)
            point = np.reshape(point,(1,3))
            ret = np.concatenate((ret,point),0)
        return ret
    else:
        return involute_val(t)

def involute_deriv_func(t,r,a=0,rad_offs=0,tan_offs=0):
    def diff_val(val):
        x = r * (-np.sin(val) + (val - a) * np.cos(val - a) + np.sin(val - a)) - \
            rad_offs * np.sin(val) - tan_offs * np.cos(val)
        y = r * (np.cos(val) + (val - a) * np.sin(val - a) - np.cos(val - a)) + \
            rad_offs * np.cos(val) - tan_offs * np.sin(val)
        z = 0
        return np.array((x,y,z))
    if hasattr(t, '__iter__'):
        ret = np.empty((0, 3))
        for u in t:
            point = diff_val(u)
            point = np.reshape(point,(1,3))
            ret = np.concatenate((ret,point),0)
        return ret
    else:
        return diff_val(t)


def involute_height_func(k, r, **kwargs):
    '''
    Returns the radial height of the involute compared to the base circle.
    '''
    return np.linalg.norm(involute_func(k, r, **kwargs)) - r


def involute_point_gen(t,r,**kwargs):
    '''
    Returns a list of points to be for cubic bezier approximation of the involute curve.
    Output is compatible with Mobject.points.
    Input t is a list where the involute shall be evaluated, it can be unevenly spaced.
    Anchors are added automatically.
    '''
    end_points = involute_func(t,r,**kwargs)
    diff_points = involute_deriv_func(t,r,**kwargs)
    out_points = np.empty((0,3))
    for i in range(len(t)-1):
        t_ratio =  (t[i+1]-t[i]) / 3
        point1 = end_points[i,:]
        point2 = end_points[i+1,:]
        anchor_1 = point1 + diff_points[i,:] * t_ratio
        anchor_2 = point2 - diff_points[i+1,:] * t_ratio
        out_points = np.append(out_points,[end_points[i,:],anchor_1,anchor_2, end_points[i+1,:]],axis=0)

    return out_points


class Gear(VMobject):
    def __init__(self,
                 num_of_teeth,
                 module=0.2,
                 alpha=20,
                 h_a=1,
                 h_f=1.17,
                 inner_teeth=False,
                 profile_shift=0,
                 cutout_teeth_num=0,
                 nppc=5,
                 **kwargs):
        '''
        Basic involute gear. 2 gears need to have the same module and alpha parameters to mesh properly.
        h_a and h_f may be slightly different but should be close.

        Parameters
        ----------
        num_of_teeth: number of gear teeth.
        module: standard size scaling parameter. Diameter = module * num_of_teeth.
        alpha: pressure angle in degrees, affects tooth curvature. Suggested values between 10-30
        h_a: addendum / module coefficient (tooth height above pitch circle)
        h_f: dedendum / module coefficient (tooth height below pitch circle)
        inner_teeth: generate gear where the teeth point inward (for example for planetary gear setup)
        profile_shift: profile shift coefficient x. Changes shape and diameter slightly, reduces undercut.
        cutout_teeth_num: number of teeth not realized
        nppc: number of points per curve. 1 tooth is constructed from 4~6 curve pieces depending on undercut.

        Examples
        --------
        class gear_example(Scene):
            def construct(self):
                # small gear
                gear1=Gear(15, stroke_opacity=0, fill_color=WHITE,fill_opacity=1)
                # larger gear
                gear2=Gear(25,  stroke_opacity=0, fill_color=RED, fill_opacity=1)
                # shifting gear1 away from center
                gear1.shift(-gear1.rp * 1.5 * RIGHT)
                # position gear2 next to gear1 so that they mesh together
                gear2.mesh_to(gear1)

                self.add(gear1, gear2)
                self.play(Rotate(gear1, gear1.pitch_angle, rate_func=linear),
                          Rotate(gear2, - gear2.pitch_angle, rate_func=linear),
                          run_time=4)
        '''
        self.z = num_of_teeth
        self.z_cut = cutout_teeth_num
        self.m = module

        # rp = pitch circle
        # when 2 gears mesh, their pitch circles need to be tangent
        self.rp = module*self.z/2
        # pressure angle
        self.alpha = alpha
        # tooth height
        self.h = (h_a+h_f)*self.m
        # addendum and dedendum coefficients
        self.h_a = h_a
        self.h_f = h_f
        # arc length of a tooth-period
        self.pitch = self.m * PI
        # base circle of involute function
        self.rb = self.rp * np.cos(self.alpha*DEGREES)
        self.X = profile_shift * module

        # for inner teeth, the top / bottom extensions are reversed
        if inner_teeth:
            # ra : outer radius (top of teeth)
            self.ra = self.rp + self.m * (h_f + profile_shift)
            # rf: inner radius (bottom of teeth)
            self.rf = self.rp - self.m * (h_a - profile_shift)
        else:
            self.ra = self.rp + self.m * (h_a + profile_shift)
            self.rf = self.rp - self.m * (h_f - profile_shift)
        self.inner_teeth = inner_teeth

        # angle_ofs: to be used with the construction of involutes
        self.angle_ofs = 0
        # angular period of teeth
        self.pitch_angle = self.pitch / self.rp
        # number of points per involute curve and per arc
        self.nppc = nppc

        # note: points are created by the 'generate_points' function, which is called by some of the supers upon init
        super().__init__(**kwargs)

        # these submobjects are a bit of a hack
        # they are used to track the center and angular position of the gear
        self.submobjects.append(VMobject(stroke_opacity=0, fill_opacity=0))
        self.submobjects[0].points=ORIGIN
        self.submobjects.append(VMobject(stroke_opacity=0, fill_opacity=0))
        self.submobjects[1].points = RIGHT

    def get_center(self):
        return self.submobjects[0].points.copy()

    def get_angle_vector(self):
        return self.submobjects[1].points-self.submobjects[0].points

    def get_angle(self):
        v = self.get_angle_vector()
        return np.arctan2(v[1], v[0])

    def generate_points(self):

        # involute starts at 0 angle at rb, but it should be at 0 on rp, so need an offset angle
        angle_base = fsolve(lambda u: involute_height_func(u, self.rb) - (self.rp - self.rb), self.alpha * DEGREES,
                            xtol=1e-10)
        self.angle_ofs = angle_base[0] - self.alpha * DEGREES

        # from tec-science article
        # https://www.tec-science.com/mechanical-power-transmission/involute-gear/profile-shift/
        # thicknes of the tooth on the pitch circle
        s0 = self.pitch / 2 + 2 * self.X * np.tan(self.alpha * DEGREES)
        ds = s0-self.pitch/2
        da = ds/2 / self.rp

        self.angle_ofs = angle_base[0] - self.alpha * DEGREES + da

        # find t-range for the involute that lies inside the rf-ra range

        def invo_cross_diff(t):
            # rotate involute into a position where the tooth would be symmetrically on the x axis
            p1 = rotate_vector(involute_func(t[0], self.rb),self.pitch_angle/4+self.angle_ofs)
            # when y coordinate is 0, the 2 involutes of the tooth would intersect because of the symmetry
            return p1[1]
        # find max height
        t_hmax = fsolve(invo_cross_diff,angle_base[0]*2)
        hmax=involute_height_func(t_hmax[0],self.rb)

        undercut = False
        if self.ra > self.rb+hmax:
            self.ra = self.rb+hmax
        res = fsolve(lambda u: involute_height_func(u,self.rb)-(self.ra-self.rb) , self.alpha * DEGREES,xtol=1e-9)
        tmax = res[0]
        if(self.rf>self.rb):
            res = fsolve(lambda u: involute_height_func(u,self.rb)-(self.rf-self.rb) , self.alpha * DEGREES,xtol=1e-9)
            tmin = res[0]
        else:
            tmin=0

        ucut_amount = (self.rf / np.cos(self.alpha * DEGREES) - self.rb)
        v_loc = (self.rb + ucut_amount) * RIGHT
        v_loc_2 = rotate_vector(v_loc, -self.alpha * DEGREES)
        ofs_vector = -self.rp * RIGHT + v_loc_2
        rad_ucut = ofs_vector[0]
        tan_ucut = ofs_vector[1]

        def undercut_func(t):
            return involute_func(t, self.rp, rad_offs=rad_ucut, tan_offs=tan_ucut)

        # undercut happening according to standard criteria
        if self.z < 2 / (np.sin(self.alpha * DEGREES) ** 2):
            undercut = True

            def diff_val_func(t):

                if t[1]>0:
                    invo_val = rotate_vector(involute_func(-t[1], self.rb), - self.alpha * DEGREES)
                else:
                    invo_val = rotate_vector(involute_func(t[1], self.rb), - self.alpha * DEGREES)
                ucut_val = undercut_func(t[0])
                diff = ucut_val - invo_val
                return diff[0:2]

            tres_ucut = fsolve(diff_val_func,np.array([0.01,0.05]))
            tmin = tres_ucut[1]
            if tmin<0:
                tmin = -tmin
            tmax_ucut = tres_ucut[0]

            # find where the undercut goes down to the root
            [tmin_ucut] = fsolve(lambda t: np.linalg.norm(undercut_func(t))-self.rf,0)
            t_range_ucut = np.linspace(tmin_ucut,tmax_ucut,self.nppc)
            undercut_curve = VMobject()
            undercut_curve.points = involute_point_gen(t_range_ucut,self.rp,rad_offs=rad_ucut, tan_offs=tan_ucut)

        # if the root circle is smaller than the base, I'm using the undercut curve to smooth out the transition between
        # base and root, simply because it provides a nice tangent curve.
        elif self.rf < self.rb :
            [tmin_ucut] = fsolve(lambda t: np.linalg.norm(undercut_func(t)) - self.rf, 0)
            [tmax_ucut] = fsolve(lambda t: np.linalg.norm(undercut_func(t)) - self.rb, 0)
            t_range_ucut = np.linspace(tmin_ucut, tmax_ucut, self.nppc)
            undercut_curve = VMobject()
            undercut_curve.points = involute_point_gen(t_range_ucut, self.rp, rad_offs=rad_ucut, tan_offs=tan_ucut)
            undercut = True

        trange_invo = np.linspace(-tmax,-tmin,self.nppc)
        involute_curve = VMobject()
        involute_curve.points = involute_point_gen(trange_invo,self.rb)
        involute_curve.rotate_about_origin(-self.alpha*DEGREES)

        if undercut:
            undercut_curve.reverse_direction()
            mid_point = (undercut_curve.points[1,:] + involute_curve.points[-2,:])/2
            undercut_curve.points[0, :] = mid_point
            involute_curve.points[-1, :] = mid_point
            involute_curve.append_points(undercut_curve.points)

        # rotate to construction position
        involute_curve.rotate(angle=self.pitch_angle/4 + self.angle_ofs + self.alpha*DEGREES,
                              about_point=ORIGIN)

        involute_curve2 = involute_curve.copy().flip(axis=RIGHT, about_point=ORIGIN)

        angle_bot_point = involute_curve.points[-1]
        angle_bot = np.arctan2(angle_bot_point[1],angle_bot_point[0])
        arc_bot_1 = Arc(radius=self.rf,
                        start_angle=angle_bot,
                        angle= self.pitch_angle/2-angle_bot, num_components=self.nppc//2+1)
        arc_bot_2 = arc_bot_1.copy().flip(axis=RIGHT, about_point=ORIGIN)
        arc_bot_1.reverse_points()
        arc_top = ArcBetweenPoints(radius=self.rf,
                                   start=involute_curve2.points[0],
                                   end=involute_curve.points[0], num_components=self.nppc)
        arc_top.reverse_points()

        involute_curve.reverse_direction()
        def smooth_curve_joint(curve1: VMobject, curve2: VMobject):
            mid_point = (curve2.points[1, :] + curve1.points[-2, :]) / 2
            curve2.points[0, :] = mid_point
            curve1.points[-1, :] = mid_point

        smooth_curve_joint(arc_bot_1,involute_curve)
        smooth_curve_joint(involute_curve, arc_top)
        smooth_curve_joint(arc_top,involute_curve2)
        smooth_curve_joint(involute_curve2,arc_bot_2)

        tooth_curve_points = np.concatenate((
            arc_bot_1.points,
            involute_curve.points,
            arc_top.points,
            involute_curve2.points,
            arc_bot_2.points
        ))

        self.points = np.empty((0,3))
        # self.points=tooth_curve_points
        for k in range(self.z-self.z_cut):
            self.points = np.concatenate((self.points,tooth_curve_points),0)
            self.rotate(self.pitch / self.rp, about_point=ORIGIN)

        if self.z_cut!=0:
            self.rotate(-self.pitch_angle*(self.z-self.z_cut+1)/2,about_point=ORIGIN)
            arc_patch = Arc(start_angle=np.arctan2(self.points[-1,1],self.points[-1,0]),
                            angle=-self.z_cut*self.pitch_angle,
                            radius=self.rf,
                            arc_center=ORIGIN)

            self.append_points(arc_patch.points)

        if self.inner_teeth:
            Outer_ring = Circle(radius=self.ra*1.1)
            # self.reverse_direction()
            self.append_points(Outer_ring.points)

    def mesh_to(self, gear2: 'Gear', offset: float = 0, positive_bias: bool = True):
        ''' This will position and rotate the gear (self) next to the input gear2 so that they mesh properly.

        Parameters
        ----------
        gear2: the other gear this gear (self) will mesh to. gear2 will not move due to meshing, only the 'self'.
        offset: axial distance offset coefficient. The gears will be offset*module further apart than default.
        positive_bias: When offset is used, there will play between gears. If positive_bias= True,
            this function meshes 'self' gear to gear2 as if there was a positive rotation torque on 'self'.'''

        # get the basic distance vector
        # remember: diff vect points towards self
        diff_vect = self.get_center() - gear2.get_center()
        distance = np.linalg.norm(diff_vect)
        if distance != 0:
            # making it unit vector
            diff_vect = diff_vect / distance
        else:
            diff_vect = RIGHT

        # calculate necessary axial distance. Inside-gears complicate things, as usual.
        # The pitch point is not in the middle between pitch circles. The calculation is based on triangle relations.
        if gear2.inner_teeth:
            pitch_dist = ( gear2.rp - self.rp - offset * self.m - self.X + gear2.X)
            rp1 = self.rb * pitch_dist / (-self.rb + gear2.rb)
            rp2 = gear2.rb * pitch_dist / (-self.rb + gear2.rb)
        elif self.inner_teeth:
            pitch_dist = (-gear2.rp + self.rp - offset * self.m + self.X - gear2.X)
            rp1 = self.rb * pitch_dist / (self.rb - gear2.rb)
            rp2 = gear2.rb * pitch_dist / (self.rb - gear2.rb)
        else:
            pitch_dist = (self.rp+gear2.rp + offset * self.m + self.X + gear2.X)
            rp1 = self.rb * pitch_dist / (self.rb+gear2.rb)
            rp2 = gear2.rb * pitch_dist / (self.rb+gear2.rb)

        self.shift(diff_vect * (-distance + pitch_dist))

        if offset != 0 or gear2.X != 0 or self.X != 0:
            # find the invo roll-angle where the curve goes as high (out) as the pitch point
            if self.inner_teeth:
                invo_offset_1 = fsolve(lambda t: involute_height_func(t,self.rb) - (rp1 - self.rb),
                                       self.angle_ofs + self.alpha*DEGREES)
            else:
                invo_offset_1 = fsolve(lambda t: involute_height_func(t, self.rb) - (rp1 - self.rb),
                                       self.angle_ofs + self.alpha * DEGREES)
            if gear2.inner_teeth:
                invo_offset_2 = fsolve(lambda t: involute_height_func(t, gear2.rb) - (rp2 - gear2.rb),
                                       gear2.angle_ofs + gear2.alpha * DEGREES)
            else:
                invo_offset_2 = fsolve(lambda t: involute_height_func(t, gear2.rb) - (rp2 - gear2.rb),
                                       gear2.angle_ofs + gear2.alpha*DEGREES)
            invo_point_1 = involute_func(invo_offset_1[0], self.rb)
            invo_point_2 = involute_func(invo_offset_2[0], gear2.rb)
            angle_offset_1 = (np.arctan2(invo_point_1[1], invo_point_1[0]) - self.angle_ofs)
            angle_offset_2 = (np.arctan2(invo_point_2[1], invo_point_2[0]) - gear2.angle_ofs)
            if not positive_bias:
                angle_offset_2 = -angle_offset_2
                angle_offset_1 = -angle_offset_1
        else:
            angle_offset_1 = 0
            angle_offset_2 = 0

        # angle of the diff vector
        diff_angle = np.arctan2(diff_vect[1], diff_vect[0])

        # get the 'tooth-phases'
        # these mods represent how much a gear has turned within the repeating pattern of the teeth
        # diff_angle needs to be considered for turning due to movement. Think planetary gear movement.
        # In some places -PI is involved due to diff vector pointing towards or away from the pitch point
        # (and sometimes added PI due to experimentation)
        if self.inner_teeth:
            mod1 = (self.get_angle() - diff_angle - PI - angle_offset_1) % self.pitch_angle / self.pitch_angle
            mod2 = (gear2.get_angle() - diff_angle - PI - angle_offset_2) % gear2.pitch_angle / gear2.pitch_angle
        elif gear2.inner_teeth:
            mod1 = (self.get_angle() - diff_angle - angle_offset_1) % self.pitch_angle / self.pitch_angle
            mod2 = (gear2.get_angle() - diff_angle - angle_offset_2) % gear2.pitch_angle / gear2.pitch_angle
        else:
            mod1 = (self.get_angle() - diff_angle - PI - angle_offset_1) % self.pitch_angle / self.pitch_angle
            mod2 = (gear2.get_angle() - diff_angle - angle_offset_2) % gear2.pitch_angle / gear2.pitch_angle

        # with inside gears, the tooth goes into a tooth-hole, they overlap
        if self.inner_teeth or gear2.inner_teeth:
            self.rotate((+mod2 - mod1) * self.pitch_angle)
        # with outside gears, the tooth pattern needs to shift half-cycle, and the rotation is reversed
        else:
            self.rotate((-mod2 - mod1 + 0.5) * self.pitch_angle)

    def rotate(
            self,
            angle: float,
            axis: np.ndarray = OUT,
            about_point: Optional[Sequence[float]] = None,
            **kwargs,
            ):
        '''Override of original rotate function so that if about_point is not specified, use the gear center'''
        if about_point is None:
            ret = super().rotate(angle, axis, about_point=self.get_center(), **kwargs)
        else:
            ret = super().rotate(angle, axis, about_point=about_point, **kwargs)
        return ret


class Rack(VMobject):
    def __init__(self, num_of_teeth, module=0.2, alpha=20, h_a=1, h_f=1.17, **kwargs):
        '''
        Basic rack for involute gears (pinion-rack connection).
        Must have the same module and alpha as the gear for proper meshing.

        Parameters
        ----------
        num_of_teeth: number of gear teeth.
        module: standard size scaling parameter. Diameter = module * num_of_teeth.
        alpha: pressure angle in degrees, affects tooth curvature. Suggested values between 10-30
        h_a: addendum / module coefficient (tooth height above pitch circle)
        h_f: dedendum / module coefficient (tooth height below pitch circle)

        Examples
        --------
        class test_Rack(Scene):
            def construct(self):
                rack1 = Rack(10,module=1, h_a=1.17, stroke_opacity=0,fill_opacity=1,fill_color=RED)
                gear1 = Gear(10, module=1,stroke_opacity=0,fill_opacity=1,fill_color=WHITE)
                gear1.shift(RIGHT*gear1.rp)
                rack1.shift(UP*rack1.pitch*0.5)
                self.add(gear1,rack1)
                self.play(Rotate(gear1, gear1.pitch_angle * 2),
                          rack1.animate.shift(DOWN*rack1.pitch * 2),
                          rate_func=linear, run_time=10)
        '''
        self.z = num_of_teeth
        self.m = module

        # pressure angle
        self.alpha = alpha
        # tooth height
        self.h = (h_a + h_f) * self.m
        # addendum and dedendum coefficients
        self.h_a = h_a
        self.h_f = h_f
        # arc length of a tooth-period
        self.pitch = self.m * PI

        # note: points are created by the 'generate_points' function, which is called by some of the supers upon init
        super().__init__(**kwargs)

        # these submobjects are a bit of a hack
        # they are used to track the center and angular position of the gear
        self.submobjects.append(VMobject(stroke_opacity=0, fill_opacity=0))
        self.submobjects[0].points = ORIGIN
        self.submobjects.append(VMobject(stroke_opacity=0, fill_opacity=0))
        self.submobjects[1].points = UP

    def generate_points(self):

        h_amax = self.pitch / 4 / np.tan(self.alpha*DEGREES)
        da = self.pitch/4*(h_amax-self.h_a*self.m)/h_amax
        h_fmax = h_amax
        df = self.pitch / 4 * (h_fmax - self.h_f*self.m) / h_amax

        tooth_points = [UP*(self.pitch/2)+LEFT*self.h,
                        UP*(self.pitch/2-df)+LEFT*self.h,
                        UP*da, ORIGIN, DOWN*da,
                        DOWN * (self.pitch / 2 - df) + LEFT * self.h,
                        DOWN * (self.pitch / 2) + LEFT * self.h
                        ]

        self.set_points_as_corners(tooth_points)
        for k in range(self.z-1):
            self.shift(UP * self.pitch)
            self.add_points_as_corners(tooth_points)

        self.shift(DOWN*self.pitch*self.z/2+RIGHT*self.h_a*self.m)
        point2 = LEFT * self.h/2 + self.points[-1, :]
        self.add_line_to(point2)
        self.add_line_to(self.points[0, :]+LEFT*self.h/2)
        self.add_line_to(self.points[0, :])

    def get_center(self):
        return self.submobjects[0].points.copy()

    def get_angle_vector(self):
        return self.submobjects[1].points - self.submobjects[0].points

    def get_angle(self):
        v = self.get_angle_vector()
        return np.arctan2(v[1], v[0])

    def rotate(
            self,
            angle: float,
            axis: np.ndarray = OUT,
            about_point: Optional[Sequence[float]] = None,
            **kwargs,
            ):
        '''Override of original rotate function so that if about_point is not specified, use the gear center'''
        if about_point is None:
            ret = super().rotate(angle, axis, about_point=self.get_center(), **kwargs)
        else:
            ret = super().rotate(angle, axis, about_point=about_point, **kwargs)
        return ret
