# manim-Gearbox
This is a plugin for Manim that enables you to draw realistic looking gears and mechanisms.
So far only involute gears are supported, with inside and outside gears.

Planned further development:
- Rack and pinion
- Cycloid gears, cycloid rack
- Sliced gears
- Animation helpers

#Installation
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

**Example**
```py
	class gear_example(Scene):
		def construct(self):
			# small gear
			gear1=Gear(15, stroke_opacity=0, fill_color=WHITE,fill_opacity=1)
			# larger gear
			gear2=Gear(25,  stroke_opacity=0, fill_color=RED, fill_opacity=1)
			# shifting gear away from center
			gear1.shift(-gear1.rp * 1.5 * RIGHT)
			gear2.mesh_to(gear1)
	
			self.add(gear1, gear2)
			self.play(Rotate(gear1, gear1.pitch_angle, rate_func=linear,about_point=gear1.get_center()),
					  Rotate(gear2, - gear2.pitch_angle, rate_func=linear,about_point=gear2.get_center()),
					  run_time=4)
		
```
![involute_gear_example](/media/involute_gear_example.gif)