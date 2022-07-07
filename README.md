# manim-Gearbox
This is a plugin for Manim that enables you to draw realistic looking gears and mechanisms.
So far only involute gears are supported, with inside and outside gears.

Planned further development:
- Rack and pinion
- Cycloid gears, cycloid rack
- Sliced gears
- Animation helpers

# Installation
`manim-gearbox` is a package on pypi, and can be directly installed using pip:
```
pip install manim-gearbox
```
# Usage
Make sure include these two imports at the top of the .py file
```py
from manim import *
from manim_gearbox import *
```

#Examples

**2 basic gears**
```py
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
		
```
![involute_gear_example](/media/involute_gear_example.gif)

**inner gear**
```py
class gear_example_inner(Scene):
    def construct(self):
        # smaller gear
        gear1 = Gear(15, module=1, stroke_opacity=0, fill_color=WHITE,fill_opacity=1)
        # larger gear with inner teeth
        gear2 = Gear(36, module=1, inner_teeth=True, stroke_opacity=0, fill_color=RED, fill_opacity=1)
        gear1.shift(gear1.rp * UP)
        gear2.mesh_to(gear1)

        self.add(gear1)
        self.add(gear2)
        self.play(Rotate(gear1, gear1.pitch_angle, rate_func=linear),
                  Rotate(gear2, gear2.pitch_angle, rate_func=linear),
                  run_time=10)
		
```
![inner_gear_example](/media/inner_gear_example.gif)
