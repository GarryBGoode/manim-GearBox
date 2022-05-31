from manim_gearbox import *
from manim import *

def test_version():
    assert __version__ == '0.1.0'


class gear_test(Scene):
    def construct(self):
        # 20 tooth gear with module = 0.2
        gear1=Gear(0.2, 20, stroke_opacity=0, fill_color=WHITE,fill_opacity=1)
        # 40 tooth gear with module = 0.2
        gear2=Gear(0.2, 40, stroke_opacity=0, fill_color=RED, fill_opacity=1)

        # shifting gears away from center
        gear2.shift(gear2.rp * RIGHT)
        gear1.shift(-gear1.rp * RIGHT)

        self.add(gear1)
        self.add(gear2)
        self.play(Rotate(gear1, PI / 10, rate_func=linear), Rotate(gear2, - PI / 10 / 2, rate_func=linear), run_time=4)