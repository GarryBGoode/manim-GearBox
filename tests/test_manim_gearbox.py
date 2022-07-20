from manim import *
from manim_gearbox import *

class gear_example(Scene):
    def construct(self):
        # small gear
        gear1=Gear(15, module=1, stroke_opacity=0, fill_color=WHITE,fill_opacity=1)
        # larger gear
        gear2=Gear(27, module=1, stroke_opacity=0, fill_color=RED, fill_opacity=1)
        # shifting gear1 away from center
        gear1.shift(gear1.rp * 1 *LEFT)
        # position gear2 next to gear1 so that they mesh together
        gear2.mesh_to(gear1,offset=0.2*gear1.m)

        # self.play(Create(gear1), Create(gear2))
        self.add(gear1,gear2)
        self.play(Rotate(gear1, gear1.pitch_angle, rate_func=linear),
                  Rotate(gear2, - gear2.pitch_angle, rate_func=linear),
                  run_time=4)

class gear_example_inner(Scene):
    def construct(self):
        # smaller gear
        gear1 = Gear(12, module=1, profile_shift=0.3, stroke_opacity=0, fill_color=WHITE,fill_opacity=1)
        # larger gear with inner teeth
        gear2 = Gear(36, module=1, inner_teeth=True, profile_shift=0.1, stroke_opacity=0, fill_color=RED, fill_opacity=1)
        gear1.shift(gear1.rp * UP)
        gear2.shift(gear2.rp * UP)
        gear2.mesh_to(gear1,offset=0.15,bias=-1)

        self.add(gear1)
        self.add(gear2)
        self.play(Rotate(gear1, gear1.pitch_angle, rate_func=linear),
                  Rotate(gear2, gear2.pitch_angle, rate_func=linear),
                  run_time=10)

class test_Gear(MovingCameraScene):
    def construct(self):
        gear1 = Gear(20, module=0.2, stroke_opacity=0.3, stroke_width=0.25, fill_color=WHITE, fill_opacity=0.3,
                     h_a=1,h_f=1.2, cutout_teeth_num=1, profile_shift=0)
        circ1 = Circle(radius=gear1.rf, stroke_opacity=0.3,stroke_width=1)
        circ2 = Circle(radius=gear1.ra, stroke_opacity=0.3,stroke_width=1)
        circ3 = Circle(radius=gear1.rp, stroke_opacity=0.3, stroke_color=GREEN, stroke_width=1)
        circ4 = Circle(radius=gear1.rb, stroke_opacity=0.3, stroke_color=BLUE, stroke_width=1)
        grp1 = VGroup(gear1, circ1, circ2, circ3, circ4)
        grid1 = NumberPlane()
        grp1.shift(gear1.rp * LEFT)
        # self.camera.frame.move_to(gear1)
        # self.camera.frame.set(width=2)
        self.add(circ1,circ2,circ3,circ4,grid1,gear1)

class test_Gear_mesh(MovingCameraScene):
    def construct(self):
        m=1
        x=0
        gear1 = Gear(10, module=m, stroke_opacity=0, stroke_width=0.5, fill_color=WHITE, fill_opacity=0.5,
                     h_a=1,h_f=1.2,profile_shift=x, nppc=20)
        gear2 = Gear(30, module=m, stroke_opacity=0, stroke_width=0.5, fill_color=WHITE, fill_opacity=0.5,
                     h_a=1, h_f=1.2, profile_shift=0.2,nppc=20,
                     inner_teeth=False)
        rack1 = Rack(21,module=m,h_a=1,h_f=1,
                     fill_color=WHITE, fill_opacity=0.5,
                     stroke_opacity=0, stroke_width=0.5,)
        rack1.rotate_about_origin(PI)
        circ1 = Circle(radius=gear1.rf, stroke_opacity=0.3,stroke_width=1)
        circ2 = Circle(radius=gear1.ra, stroke_opacity=0.3,stroke_width=1)
        circ3 = Circle(radius=gear1.rp, stroke_opacity=0.3, stroke_color=GREEN, stroke_width=1)
        circ4 = Circle(radius=gear1.rb, stroke_opacity=0.3, stroke_color=BLUE, stroke_width=1)
        grp1 = VGroup(gear1,circ1,circ2,circ3,circ4)
        grid1 = NumberPlane()
        grp1.shift(gear1.rp*LEFT*1.2)
        gear2.shift(gear2.rp*LEFT)
        rack1.shift(x*m*RIGHT)

        gear1.rotate(gear1.pitch_angle*0.2)

        ofs_tr = ValueTracker(0.0)
        bias_tr = ValueTracker(1)
        gear1.mesh_to(gear2,offset=0.3*m,bias=1)


        gear1.add_updater(lambda mob: mob.mesh_to(gear2,offset=ofs_tr.get_value(),bias=bias_tr.get_value()))
        # self.camera.frame.move_to(gear1)
        # self.camera.frame.set(width=2)
        self.add(circ1,circ2,circ3,circ4,grid1,gear1,gear2)
        self.play(ofs_tr.animate.set_value(0.8))
        self.wait()
        self.play(bias_tr.animate.set_value(0))
        self.play(bias_tr.animate.set_value(-1))
        self.wait()
        self.play(ofs_tr.animate.set_value(0))
        self.wait()


class test_Rack(Scene):
    def construct(self):
        rack1 = Rack(10,module=1, h_a=1.2, stroke_opacity=0,fill_opacity=1,fill_color=RED)
        gear1 = Gear(10, module=1,stroke_opacity=0,fill_opacity=1,fill_color=WHITE,nppc=10)
        gear1.shift(UP*gear1.rp)
        gear1.rotate(PI/2)
        rack1.rotate(PI/2)
        self.add(gear1,rack1)
        self.play(Rotate(gear1, gear1.pitch_angle * 2),
                  rack1.animate.shift(RIGHT*rack1.pitch * 2),
                  rate_func=linear, run_time=10)


class test_Gear_inner(MovingCameraScene):
    def construct(self):
        gear1 = Gear(24, module=0.2, inner_teeth=True,
                     cutout_teeth_num=0, profile_shift=0.0,
                     stroke_opacity=1, stroke_width=1, fill_color=WHITE, fill_opacity=0.3, h_a=1,h_f=1)

        circ1 = Circle(radius=gear1.rf, stroke_opacity=0.3, stroke_width=1)
        circ2 = Circle(radius=gear1.ra, stroke_opacity=0.3, stroke_color=YELLOW, stroke_width=1)
        circ3 = Circle(radius=gear1.rp, stroke_opacity=0.3, stroke_color=GREEN, stroke_width=1)
        circ4 = Circle(radius=gear1.rb, stroke_opacity=0.3, stroke_color=BLUE, stroke_width=1)
        grid1 = NumberPlane()


        # self.camera.frame.move_to(RIGHT*gear1.ra)
        # self.camera.frame.set(width=2)
        self.add(circ1,circ2,circ3,circ4,grid1,gear1)


# brute force param combination tester
class test_gear_param(Scene):
    def construct(self):
        ha_params = np.array([3,  5,  7,  9, 10]) /10
        hf_params = ha_params * 1.5

        min_tooth = 8
        max_tooth = 40
        tooth_arr = np.arange(min_tooth,max_tooth,3)

        gear_grp = VGroup()

        for tooth in tooth_arr:
            for h_a in ha_params:
                for h_f in hf_params:
                    gear =  Gear(tooth, h_a=h_a,h_f=h_f,stroke_width=0.5)
                    print("ha=",h_a,"hf=",h_f,"tooth=",tooth)
                    # gear2 = Gear(tooth, h_a=h_a,h_f=h_f, inner_teeth=True)
                    gear_grp.add(gear)

        self.add(gear_grp)

# this part is used for debugging
# with tempconfig({"quality": "medium_quality", "disable_caching": True}):
#     scene = gear_example()
#     scene.render()
