# Copyright (c) 2019 Simón Flores (https://github.com/cheece)

# Permission is hereby granted, free of charge,
# to any person obtaining a copy of this software
# and associated documentation files (the "Software"),
# to deal in the Software without restriction,
# including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so,
# subject to the following conditions:The above copyright to_quaternion
# notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY
# OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
# EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
# AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

# This addon adds jiggle physics to bones.
# How to use:
#   Enable "Jiggle Scene" in the scene properties.
#   Enable "Jiggle Bone" on the bones.
#   In the Jiggle Scene settings, click "Initialize Bones".

# Based on the Position Based Dynamics paper by Müller et al. http://matthias-mueller-fischer.ch/publications/posBasedDyn.pdf

bl_info = {
	"name": "Jiggle Armature",
	"author": "Simón Flores",
	"version": (2, 2, 1),
	"blender": (2, 80, 0),
	"description": "Jiggle bone animation tool",
	"warning": "",
	"wiki_url": "",
	"category": "Animation",
}

import bpy
from . import operators
from . import physics
from . import ui

files = (
	operators,
	physics,
	ui,
)

def register():
	for f in files:
		f.register()

def unregister():
	for f in files:
		f.unregister()