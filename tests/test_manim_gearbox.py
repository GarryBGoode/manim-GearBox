from manim_gearbox import *
from manim import *

def test_version():
    assert __version__ == '0.1.0'


class gear_test(Scene):
    def construct(self):
        gear1=Gear(0.2,20,stroke_color=WHITE,stroke_width=1.5, stroke_opacity=0, fill_color=WHITE,fill_opacity=0.5)
        gear2=Gear(0.2,40,
                   stroke_color=RED,stroke_width=1.5,stroke_opacity=0,
                   inner_teeth=True,fill_color=RED,fill_opacity=1,h_a=0.7)
        # gear1.rotate(-gear1.rp_tooth_ofs)
        # gear2.rotate(-gear2.rp_tooth_ofs)
        circ = Circle(gear2.ra+0.3)

        gear3 = Difference(circ,gear2)
        gear2.shift((gear2.rp + gear1.rp) * RIGHT)
        gear3.shift((-gear2.rp+gear1.rp*0)*RIGHT)
        gear3.set_fill(color=RED,opacity=1)
        gear3.set_stroke(opacity=0)
        gear1.shift(-gear1.rp*RIGHT)

        self.add(gear1)
        # self.add(gear2)
        self.add(gear3)

        # self.play(Rotate(gear1, PI / 10, rate_func=linear), Rotate(gear2, -PI / 10 / 3, rate_func=linear), run_time=4)
        self.play(Rotate(gear1, PI / 10, rate_func=linear), Rotate(gear3, PI / 10 / 2, rate_func=linear), run_time=4)