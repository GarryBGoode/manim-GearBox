import numpy as np
from manim import *
from scipy.optimize import fsolve
from scipy.optimize import fmin
from scipy.optimize import fmin_powell

__all__ = [
    "involute_func",
    "involute_height_func",
    "Gear"
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
        x = r * (np.cos(val) + (val-a) * np.sin(val-a)) + offs_v[0]
        y = r * (np.sin(val) - (val-a) * np.cos(val-a)) + offs_v[1]
        z = 0
        return np.array((x,y,z))
    if hasattr(t, '__iter__'):
        return np.array([involute_val(u) for u in t])
    else:
        return involute_val(t)

def involute_height_func(k, r, **kwargs):
    '''
    Returns the radial height of the involute compared to the base circle.
    '''
    return np.linalg.norm(involute_func(k, r, **kwargs)) - r





class Gear(VMobject):
    def __init__(self,  num_of_teeth, module=0.2, alpha=20, h_a=1,h_f=1.17, inner_teeth = False, **kwargs):
        '''
        Basic involue gear.

        Parameters
        ----------
        num_of_teeth: number of gear teeth. Use more than ~17, otherwise it might look bad.
        module: standard size parameter. Diameter = module * num_of_teeth. 2 gears must have same module in order to mesh.
        alpha: pressure angle in degrees, affects tooth curvature. Suggested values between 10-30
        h_a: addendum / module coefficient (tooth height above pitch circle)
        h_f: dedendum / module coefficient (tooth height below pitch circle)
        inner_teeth: generate gear where the teeth point inward (for example for planetary gear setup)

        Examples
        --------
        class gear_test(Scene):
            def construct(self):
                # 20 tooth gear
                gear1=Gear(20, stroke_opacity=0, fill_color=WHITE,fill_opacity=1)
                # 40 tooth gear
                gear2=Gear(40, stroke_opacity=0, fill_color=RED, fill_opacity=1)

                # shifting gears away from center
                gear2.shift(gear2.rp * RIGHT)
                gear1.shift(-gear1.rp * RIGHT)

                self.add(gear1,gear2)
                self.play(Rotate(gear1, PI / 10, rate_func=linear),
                          Rotate(gear2, - PI / 10 / 2, rate_func=linear),
                          run_time=4)
        '''
        self.z = num_of_teeth
        self.m = module
        self.rp = module*self.z/2
        self.alpha = alpha
        self.h = (h_a+h_f)*self.m
        self.h_a = h_a
        self.h_f = h_f
        self.pitch = self.m * PI
        self.rb = self.rp * np.cos(self.alpha*DEGREES)

        # for inner teeth, the top / bottom extensions are reversed
        if inner_teeth:
            self.ra = self.rp + self.m * h_f
            self.rf = self.rp - self.m * h_a
        else:
            self.ra = self.rp + self.m * h_a
            self.rf = self.rp - self.m * h_f
        self.inner_teeth = inner_teeth



        # note: points are created by the 'generate_points' function, which is called by some of the supers upon init
        super().__init__(**kwargs)

    def generate_points(self):

        # find t-range for the involute that lies inside the rf-ra range
        undercut = False
        res = fsolve(lambda u: involute_height_func(u,self.rb)-(self.ra-self.rb) , self.alpha * DEGREES)
        tmax = res[0]
        if(self.rf>self.rb):
            res = fsolve(lambda u: involute_height_func(u,self.rb)-(self.rf-self.rb) , self.alpha * DEGREES)
            tmin = res[0]
        else:
            tmin=0

        # involute starts at 0 angle at rb, but it should be at 0 on rp, so need an offset angle
        angle_base = fsolve(lambda u: involute_height_func(u, self.rb) - (self.rp - self.rb), self.alpha * DEGREES)
        angle_ofs =  angle_base[0]

        ucut_amount = (self.rf / np.cos(self.alpha * DEGREES) - self.rb)
        v_loc = (self.rb + ucut_amount) * RIGHT
        v_loc_2 = rotate_vector(v_loc, -self.alpha * DEGREES)
        ofs_vector = -self.rp * RIGHT + v_loc_2
        rad_ucut = ofs_vector[0]
        tan_ucut = ofs_vector[1]

        def undercut_func(t):
            return involute_func(t, self.rp, rad_offs=rad_ucut, tan_offs=tan_ucut)

        if self.z < 2 / (np.sin(self.alpha * DEGREES) ** 2):
            undercut = True
            def diff_cost_func(t):
                invo_val = rotate_vector(involute_func(-t[1], self.rb),- self.alpha*DEGREES)
                ucut_val = undercut_func(t[0])
                diff = ucut_val - invo_val
                return np.linalg.norm(diff)
                # return diff
            tres_ucut = fmin_powell(diff_cost_func, np.array([0, 0.1]))
            tmin = tres_ucut[1]
            tmax_ucut = tres_ucut[0]


            [tmin_ucut] = fsolve(lambda t: np.linalg.norm(undercut_func(t))-self.rf,0)
            t_step_ucut = (tmax_ucut-tmin_ucut) / 5
            undercut_curve = ParametricFunction(undercut_func,t_range=[tmin_ucut,tmax_ucut,t_step_ucut])
        elif self.rf < self.rb :
            [tmin_ucut] = fsolve(lambda t: np.linalg.norm(undercut_func(t)) - self.rf, 0)
            [tmax_ucut] = fsolve(lambda t: np.linalg.norm(undercut_func(t)) - self.rb, 0)
            t_step_ucut = (tmax_ucut-tmin_ucut) / 5
            undercut_curve = ParametricFunction(undercut_func, t_range=[tmin_ucut, tmax_ucut, t_step_ucut])
            undercut = True

        trange_step = (tmax - tmin) / 5
        involute_curve = ParametricFunction(lambda u: rotate_vector(involute_func(u, self.rb),- self.alpha*DEGREES), t_range=[-tmax,-tmin,trange_step])

        if undercut:
            undercut_curve.reverse_direction()
            mid_point = (undercut_curve.points[1,:] + involute_curve.points[-2,:])/2
            undercut_curve.points[0, :] = mid_point
            involute_curve.points[-1, :] = mid_point
            involute_curve.append_points(undercut_curve.points)

        # rotate back to construction position
        involute_curve.rotate(angle=-self.pitch/self.rp/4 + angle_ofs,
                              about_point=ORIGIN)

        involute_curve2 = involute_curve.copy().flip(axis=RIGHT, about_point=ORIGIN)

        arc_bot = ArcBetweenPoints(radius=self.rf,
                                   start=involute_curve.points[-1],
                                   end=involute_curve2.points[-1])
        arc_top = ArcBetweenPoints(radius=self.rf,
                                   start=rotate_vector(involute_curve2.points[0],-self.pitch/self.rp),
                                   end=involute_curve.points[0])

        involute_curve2.reverse_direction()

        # pre-smoothing joining parts
        mid_point = (involute_curve.points[1, :] + arc_top.points[-2, :]) / 2
        involute_curve.points[0, :] = mid_point
        arc_top.points[-1, :] = mid_point

        mid_point = (arc_bot.points[1, :] + involute_curve.points[-2, :]) / 2
        arc_bot.points[0, :] = mid_point
        involute_curve.points[-1, :] = mid_point

        mid_point = (involute_curve2.points[1, :] + arc_bot.points[-2, :]) / 2
        involute_curve2.points[0, :] = mid_point
        arc_bot.points[-1, :] = mid_point

        mid_point = (rotate_vector(arc_top.points[1, :], self.pitch / self.rp)+ involute_curve2.points[-2, :]) / 2
        arc_top.points[0, :] = rotate_vector(mid_point,-self.pitch / self.rp)
        involute_curve2.points[-1, :] = mid_point

        tooth_curve_points = np.concatenate((
            arc_top.points,
            involute_curve.points,
            arc_bot.points,
            involute_curve2.points
        ))
        self.points = involute_curve.points

        self.points=tooth_curve_points
        for i in range(self.z-1):
            self.rotate(-self.pitch / self.rp, about_point=ORIGIN)
            self.points = np.concatenate((self.points,tooth_curve_points),0)

        self.make_smooth()
        self.rotate(-self.pitch/self.rp/4)