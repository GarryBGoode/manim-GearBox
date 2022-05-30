import numpy as np
from manim import *
from scipy.optimize import fsolve

__all__ = [
    "involute_func",
    "invo_height_func",
    "Gear"
]

def involute_func(t, r, a=0):
    '''
    Returns the x-y-z values of the involute function.
    t: input angle
    r: base circle radius
    a: offset angle
    '''
    def involute_val(val):
        x = r * (np.cos(val) + (val-a) * np.sin(val-a))
        y = r * (np.sin(val) - (val-a) * np.cos(val-a))
        z = 0
        return np.array((x,y,z))
    if hasattr(t, '__iter__'):
        return np.array([involute_val(u) for u in t])
    else:
        return involute_val(t)

def invo_height_func(k, r):
    '''
    Returns the radial height of the involute compared to the base circle.
    '''
    return np.linalg.norm(involute_func(k, r)) - r



class Gear(VMobject):
    def __init__(self, modulus, z_teeth, alpha=20, h_a=1,h_f=1.25, inner_teeth = False, **kwargs):

        self.z = z_teeth
        self.m = modulus
        self.rp = modulus*self.z/2
        self.alpha = alpha
        self.h = (h_a+h_f)*self.m
        self.pitch = self.m * PI
        self.rb = self.rp * np.cos(self.alpha*DEGREES)
        self.ra = self.rp + self.m * h_a
        self.rf = self.rp - self.m * h_f

        if inner_teeth:
            self.ra = self.rp + self.m * h_f
            self.rf = self.rp - self.m * h_a
        self.inner_teeth = inner_teeth


        super().__init__(**kwargs)

    def generate_points(self):
        tmax = fsolve(lambda u: invo_height_func(u,self.rb)-(self.ra-self.rb) , self.alpha * DEGREES)
        if(self.rf>self.rb):
            tmin = fsolve(lambda u: invo_height_func(u,self.rb)-(self.rf-self.rb) , self.alpha * DEGREES)
        else:
            tmin=[0]
        trange_step = (tmax[0]-tmin[0])/10

        involute_curve = ParametricFunction(lambda u: involute_func(u,self.rb), t_range=[tmin[0],tmax[0],trange_step])
        angle_ofs = fsolve(lambda u: invo_height_func(u,self.rb)-(self.rp-self.rb) , self.alpha * DEGREES)-self.alpha*DEGREES
        involute_curve_2 = involute_curve.copy().flip(RIGHT,about_point=ORIGIN).rotate(self.pitch/self.rp/2+angle_ofs[0],about_point=ORIGIN)
        involute_curve.rotate(-angle_ofs[0], about_point=ORIGIN)
        self.rp_tooth_ofs = angle_ofs[0]
        fi = self.pitch / self.rp
        if self.inner_teeth:
            tanv1 = involute_curve.points[-2, :] - involute_curve.points[-1, :]
            tanv1 = tanv1 / np.linalg.norm(tanv1)
            tanv2 = involute_curve_2.points[-2, :] - involute_curve_2.points[-1, :]
            tanv2 = tanv2 / np.linalg.norm(tanv2)
            angle_arc = np.arctan2(np.linalg.norm(np.cross(tanv1, tanv2)), np.dot(tanv1, tanv2))
            arc_top = ArcBetweenPoints(involute_curve.points[-1, :], involute_curve_2.points[-1, :], angle=PI-angle_arc)
            arc_bot = ArcBetweenPoints(involute_curve_2.points[0, :],
                                       involute_curve.copy().rotate(fi, about_point=ORIGIN).points[0, :],
                                       radius=self.rf)
            #smoothing sharp corner at the bottom of the tooth
            ref_point_1 = (involute_curve_2.points[1,:]+arc_bot.points[1,:])/2
            ref_point_2 = (involute_curve.copy().rotate(fi, about_point=ORIGIN).points[1,:]+arc_bot.points[-2,:])/2
            arc_bot.points[0,:] = ref_point_1
            arc_bot.points[-1,:] = ref_point_2
            involute_curve_2.points[0, :] = ref_point_1
            involute_curve.points[0,:] = arc_bot.copy().rotate(-fi, about_point=ORIGIN).points[-1,:]

            tooth_curve_points = np.concatenate(
                [involute_curve.points, arc_top.points, involute_curve_2.reverse_points().points, arc_bot.points], 0)
        else:
            tanv1 = involute_curve.points[1, :] - involute_curve.points[0, :]
            tanv1 = tanv1 / np.linalg.norm(tanv1)
            tanv2 = involute_curve_2.points[1, :] - involute_curve_2.points[0, :]
            tanv2 = tanv2 / np.linalg.norm(tanv2)
            angle_arc = np.arctan2(np.linalg.norm(np.cross(tanv1, tanv2)), np.dot(tanv1, tanv2))
            arc_top = ArcBetweenPoints(involute_curve.points[-1, :],
                                       involute_curve_2.points[-1, :],
                                       radius=self.ra)
            arc_bot = ArcBetweenPoints(involute_curve_2.points[0, :],
                                       involute_curve.copy().rotate(fi, about_point=ORIGIN).points[0, :],
                                       angle=-PI+angle_arc)
            # smoothing sharp corner at the bottom of the tooth
            ref_point_1 = (involute_curve.points[-2, :] + arc_top.points[1, :]) / 2
            ref_point_2 = (involute_curve_2.points[-2, :] + arc_top.points[-2,:]) / 2
            arc_top.points[0, :] = ref_point_1
            arc_top.points[-1, :] = ref_point_2
            involute_curve.points[-1, :] = ref_point_1
            involute_curve_2.points[-1, :] = ref_point_2
            tooth_curve_points = np.concatenate(
                [involute_curve.points, arc_top.points, involute_curve_2.reverse_points().points, arc_bot.points], 0)

        for i in range(self.z):
            self.points = np.concatenate((self.points,tooth_curve_points),0)
            self.rotate(-self.pitch/self.rp,about_point=ORIGIN)
        #
        self.make_smooth()
