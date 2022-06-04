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
		
```
![involute_gear_example](/media/involute_gear_example.gif)