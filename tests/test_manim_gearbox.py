from manim_gearbox import *
from manim import *

class gear_example(Scene):
    def construct(self):
        # 20 tooth gear
        gear1=Gear(20, stroke_opacity=0, fill_color=WHITE,fill_opacity=1,h_a=0.9)
        # 40 tooth gear
        gear2=Gear(40, stroke_opacity=0, fill_color=RED, fill_opacity=1,h_a=0.9)
        # shifting gears away from center
        gear2.shift(gear2.rp * RIGHT)
        gear1.shift(-gear1.rp * RIGHT)

        self.add(gear1)
        self.add(gear2)
        self.play(Rotate(gear1, PI / 10, rate_func=linear),
                  Rotate(gear2, - PI / 10 / 2, rate_func=linear), run_time=4)

class test_Gear(Scene):
    def construct(self):
        gear1 = Gear(20, stroke_opacity=0, fill_color=WHITE, fill_opacity=1, h_a=0.9)
        self.add(gear1)




