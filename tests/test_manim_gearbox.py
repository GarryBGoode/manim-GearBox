from manim import *
from manim_gearbox import *

class gear_example(Scene):
    def construct(self):
        # 10 tooth gear
        gear1=Gear(14, module=0.5,stroke_opacity=0, fill_color=WHITE,fill_opacity=1,h_a=1, h_f=1)
        # 40 tooth gear
        gear2=Gear(40, module=0.5, stroke_opacity=0, fill_color=RED, fill_opacity=1,h_a=1, h_f=1)
        # shifting gears away from center
        gear2.shift(gear2.rp * RIGHT)
        gear1.shift(-gear1.rp * RIGHT)

        self.add(gear1)
        self.add(gear2)
        self.play(Rotate(gear1, 2*PI / gear1.z, rate_func=linear),
                  Rotate(gear2, - 2* PI / gear2.z, rate_func=linear), run_time=4)

class test_Gear(MovingCameraScene):
    def construct(self):
        gear1 = Gear(22, module=0.4, stroke_opacity=1, stroke_width=0.5, fill_color=WHITE, fill_opacity=0.3, h_a=1,h_f=1.2)
        circ1 = Circle(radius=gear1.rf, stroke_opacity=0.3,stroke_width=1)
        circ2 = Circle(radius=gear1.ra, stroke_opacity=0.3,stroke_width=1)
        circ3 = Circle(radius=gear1.rp, stroke_opacity=0.3, stroke_color=GREEN, stroke_width=1)
        circ4 = Circle(radius=gear1.rb, stroke_opacity=0.3, stroke_color=BLUE, stroke_width=1)
        grid1 = NumberPlane()
        # self.camera.frame.move_to(gear1)
        # self.camera.frame.set(width=2)
        self.add(circ1,circ2,circ3,circ4,grid1,gear1)

class test_Gear_inner(MovingCameraScene):
    def construct(self):
        gear1 = Gear(10, module=0.4, inner_teeth=True,
                      stroke_opacity=1, stroke_width=1, fill_color=WHITE, fill_opacity=0.3, h_a=1,h_f=1)
        circ_gear = Circle(gear1.ra+0.2)
        gear2 = Difference(circ_gear,gear1,stroke_opacity=1, stroke_width=0.3, fill_color=WHITE, fill_opacity=0.3,)
        circ1 = Circle(radius=gear1.rf, stroke_opacity=0.3, stroke_width=1)
        circ2 = Circle(radius=gear1.ra, stroke_opacity=0.3, stroke_color=YELLOW, stroke_width=1)
        circ3 = Circle(radius=gear1.rp, stroke_opacity=0.3, stroke_color=GREEN, stroke_width=1)
        circ4 = Circle(radius=gear1.rb, stroke_opacity=0.3, stroke_color=BLUE, stroke_width=1)
        grid1 = NumberPlane()


        # self.camera.frame.move_to(RIGHT*gear1.ra)
        # self.camera.frame.set(width=2)
        self.add(circ1,circ2,circ3,circ4,grid1,gear2)

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
                    gear =  Gear(tooth, h_a=h_a,h_f=h_f)
                    print("ha=",h_a,"hf=",h_f,"tooth=",tooth)
                    # gear2 = Gear(tooth, h_a=h_a,h_f=h_f, inner_teeth=True)
                    gear_grp.add(gear)

        self.add(gear_grp)
#
#
# with tempconfig({"quality": "medium_quality", "disable_caching": True}):
#     scene = test_gear_param()
#     scene.render()
