# manim-Gearbox
This is a plugin for Manim that enables you to draw realistic looking gears and mechanisms.
Mostly based on these tec-science pages:
[https://www.tec-science.com/mechanical-power-transmission/involute-gear/geometry-of-involute-gears/](https://www.tec-science.com/mechanical-power-transmission/involute-gear/geometry-of-involute-gears/)

Currently supported Involute gear features:
- Basic spur gears
- Inside ring-gears
- Basic rack
- Undercutting (gears with fewer than 17 teeth)
- Profile shifted gears
- Meshing calculation with distance variation

Planned further development:
- Cycloid gears, cycloid rack


## Installation
`manim-gearbox` is a package on pypi, and can be directly installed using pip:
```
pip install manim-gearbox
```
Note: `manim-gearbox` uses, and depends on SciPy and Manim.

## Usage
Make sure include these two imports at the top of the .py file
```py
from manim import *
from manim_gearbox import *
```
Add Gear objects, use mesh_to() method to position 2 gears into meshing.
I tend to use 'fill_opacity=1' and 'stroke_opacity=0' options because the stroke increases the gear size by a couple pixels, and gives the feeling of interference.

# Examples

## 2 basic gears

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

## inner gear

```py
class gear_example_inner(Scene):
    def construct(self):
        # smaller gear
        gear1 = Gear(12, module=1, profile_shift=0.3, stroke_opacity=0, fill_color=WHITE,fill_opacity=1)
        # larger gear with inner teeth
        gear2 = Gear(36, module=1, inner_teeth=True, profile_shift=0.1, stroke_opacity=0, fill_color=RED, fill_opacity=1)
        gear1.shift(gear1.rp * UP)
        gear2.shift(gear2.rp * UP)
		# mesh with 0.15*module larger distance than default
		# positive_bias param is used to define left or right tooth flank shall engage if there is offset and play
        gear2.mesh_to(gear1,offset=0.15,positive_bias=False)

        self.add(gear1)
        self.add(gear2)
        self.play(Rotate(gear1, gear1.pitch_angle, rate_func=linear),
                  Rotate(gear2, gear2.pitch_angle, rate_func=linear),
                  run_time=10)

```
![inner_gear_example](/media/inner_gear_example.gif)
