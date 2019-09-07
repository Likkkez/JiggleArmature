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

# TODO

### REFACTOR ###
# Rename classes, functions, variables to be consistent, verbose and descriptive.
# Organize functions and shit properly.
# 	Split into separate files: __init__, operators, jiggle, utils, ui.
#	__init__
#	operators
# 		put bake() inside JARM_OT_bake
#		put initialize_bones() inside JARM_OT_init_bones
#	utils
#		Just slap in all the floating utility functions so I don't have to constantly scroll past them.
#		I feel like some of those functions might already be in the Matrix or Quaternion classes from mathutils.
#	ui
#	jiggle (can find a better name)
#		JiggleBone/Armature/Scene classes and stuff related to it.
# 			Pack all jiggle bone properties into a CollectionProperty (need a new JiggleProperties class) and append that to bpy.types.Bone.
#				This would only be possible if all those properties could be custom properties, but that's not the case with matrices.
#			Let JiggleBones simulate themselves via step(). This would probably be quite tricky. But imo it would be a lot nicer.
#				Then the scene would have a jiggle.step() which calles every armatures' jiggle.step() which calls each of its bones' jiggle.step().

### FEATURES ###
# The whole idea of having to enable the jiggle in the scene, then in the armature, then in the bone, seems crazy to me. We should be able to simply enable it in the bone, and have some settings for it in the scene.
# There should be an option to simply use the scene's framerate as the simulation framerate. In this case, the simulation would be expected to run "faster" - That is to say, we shouldn't compensate with the assumption that 1 second is always 24 frames. 1 second is always fps frames.
# Let jiggle bones inherit scale of their parent bone, just like they inherit scale of the armature object.
# Would be nice to clean up the physics code so we could actually make changes to it and know wtf is actually happening.
# Currently, the simulation will progress forward in time, no matter where you move on the timeline.
#	We could cache animation either into internal variable or an Action.

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
from mathutils import *
from math import sqrt
from bpy.types import Menu, Panel, Operator
from bpy_extras.io_utils import ImportHelper, ExportHelper
from bpy.props import *
from collections import defaultdict
from time import time

class JiggleArmature(bpy.types.PropertyGroup):
	enabled: BoolProperty(default=True, 
		name="Enable physics bones on this armature")
	fps: FloatProperty(name = "simulation fps", default = 24)
	time_acc: FloatProperty(default = 0.0)

class JiggleScene(bpy.types.PropertyGroup):
	iterations: IntProperty(
		name = "Iterations",
		description="Higher values result in slower but higher quality simulation",
		min=1,
		default = 4)
	last_frame: IntProperty()

class JARM_PT_armature(bpy.types.Panel):
	bl_idname = "ARMATURE_PT_jiggle"
	bl_label = "Jiggle Armature"
	bl_space_type = 'PROPERTIES'
	bl_region_type = 'WINDOW'
	bl_context = "data"
	bl_options = {'DEFAULT_CLOSED'}

	def draw_header(self, context):
		layout = self.layout
		arm = context.object
		layout.prop(context.object.data.jiggle, "enabled", text="",)

	@classmethod
	def poll(cls, context):
		return (context.object is not None and context.object.type == 'ARMATURE')

	def draw(self, context):
		layout = self.layout
		col = layout.column()
		col.prop(context.object.data.jiggle, "fps")

class JARM_PT_scene(bpy.types.Panel):
	bl_idname = "SCENE_PT_jiggle"
	bl_label = "Jiggle Scene"
	bl_space_type = 'PROPERTIES'
	bl_region_type = 'WINDOW'
	bl_context = "scene"
	bl_options = {'DEFAULT_CLOSED'}

	def draw(self, context):
		layout = self.layout
		col = layout.column()
		col.prop(context.scene.jiggle,"iterations")
		col.operator("jiggle.initialize", text="Initialize Bones")
		col.operator("jiggle.bake", text="Bake Selected").a = False
		col.operator("jiggle.bake", text="Bake All").a = True

def setq(om, m):
	for i in range(4):
		om[i]= m[i]

class JARM_OT_set_rest(bpy.types.Operator):
	"""Set jiggle rest pose"""
	bl_idname = "jiggle.set_rest"
	bl_label = "Set Rest"

	def execute(self, context):
		scene = context.scene
		for o in scene.objects:
			if(o.select_get() and o.type == 'ARMATURE' ):
				arm_matrix = o.matrix_world
				scale = maxis(arm_matrix,0).length
				iow = arm_matrix.inverted()
				i=0
				for b in o.pose.bones:
					if(b.bone.select):
						M = b.parent.matrix.inverted() @ b.matrix #arm_matrix*Sbp.wmat* Sb.rmat #im
						db = b.bone
						setq(db.jiggle_rest, M.to_quaternion().normalized())
						db.jiggle_use_custom_rest = True

		return {'FINISHED'}

class JARM_PT_bone(bpy.types.Panel):
	bl_idname = "BONE_PT_jiggle_bone"
	bl_label = "Jiggle Bone"
	bl_space_type = 'PROPERTIES'
	bl_region_type = 'WINDOW'
	bl_context = "bone"
	bl_options = {'DEFAULT_CLOSED'}

	@classmethod
	def poll(cls, context):
		return ( context.bone is not None and context.object is not None and context.object.type == 'ARMATURE')

	def draw_header(self, context):
		layout = self.layout
		bon = context.bone
		layout.prop(bon , "jiggle_enabled", text="")

	def draw(self, context):
		layout = self.layout
		armature = context.object.data

		bon = context.bone
		col = layout.column()
		layout.enabled = armature.jiggle.enabled

		if(not armature.jiggle.enabled):
			col.label(text= "JiggleArmature is disabled for the armature, see the armature properties")

		if(bon.jiggle_enabled):
			col.prop(context.bone,"jiggle_Ks")
			col.prop(bon,"jiggle_Kd")
			col.prop(bon,"jiggle_Kld")
			col.prop(bon,"jiggle_mass")
			col.prop(bon,"gravity_multiplier")
			col.prop_search(bon,"jiggle_control_object",bpy.data,"objects")
			if(bon.jiggle_control_object in bpy.data.objects):
				o = bpy.data.objects[bon.jiggle_control_object]
				if(o.type == 'ARMATURE'):
					 col.prop_search(bon,"jiggle_control_bone",o.data,"bones")
				col.prop(bon,"jiggle_control")

			col.prop(bon,"jiggle_use_custom_rest")
			if(bon.jiggle_use_custom_rest):
				col.prop(bon,"jiggle_rest")
			col.operator("jiggle.set_rest")
			if(bon.parent==None):
				col.label(text= "warning: jibblebones without parent will fall",icon='COLOR_RED')

class JiggleBone:
	def __init__(self, pbone, matrix, parent):
		self.M = matrix.copy()
		self.length = pbone.bone.length*maxis(matrix,0).length
		self.b = pbone
		self.parent = parent	# Type: JiggleBone
		self.rest = None		# Type: Matrix
		self.rest_w = None		# Type: Matrix
		self.w = 0				# Type: Float - Inverse of mass?
		self.Kc = 0				# Type: Float - ?
		self.cQ = None			# Type: Quaternion
		self.X = None			# Type: Float
		self.P = None			# Type: Float
		self.R = None			# Type: Float
		self.Q = None			# Type: Float
		self.iI = Matrix.Identity(3) #first naive approach
		self.iIw = self.iI

	def computeI(self):
		self.iI = Matrix.Identity(3)*(self.w/(self.l*self.l)*5.0/2.0)

	def updateIW(self):
		rot = self.Q.to_matrix()
		self.iIw = rot@self.iI @ rot.transposed()
	
	def locSpring(self):

		Q0 = self.parent.Q
		Q1 = self.Q
		w0 = self.parent.w
		w1 = self.w

		v0 = self.rest_p

		P0 = self.parent.P
		P1 = self.P
		lf = self.l*0.5

		self.updateIW()
		self.parent.updateIW()

		connector0 = self.parent.P+self.parent.Q @ v0
		connector1 = self.P+self.Q @ Vector((0,-lf,0))

		K1 = self.computeMatrixK(connector0, w0, P0, self.parent.iIw)
		K2 = self.computeMatrixK(connector1, w1, P1, self.iIw)

		Kinv = (K1 + K2).inverted()

		pt = Kinv @ (connector1 - connector0)
		if (w0 != 0.0):
			r0 = connector0 - P0
			self.parent.P += w0*pt
			ot = (self.parent.iIw @ (r0.cross(pt)))

			otQ = Quaternion()
			otQ.x =ot[0]
			otQ.y = ot[1]
			otQ.z = ot[2]
			otQ.w = 0
			self.parent.Q = qadd(self.parent.Q, otQ @ self.parent.Q*0.5).normalized()

		if (w1 != 0.0):
			r1 = connector1 - P1
			self.P += -w1*pt
			ot = (self.iIw @ (r1.cross(-pt)))

			otQ = Quaternion()
			otQ.x =ot[0]
			otQ.y = ot[1]
			otQ.z = ot[2]
			otQ.w = 0
			self.Q = qadd(self.Q, otQ @ self.Q*0.5).normalized()
	
	def quatSpring(self, r=None, k=None):
		Q0 = self.parent.Q
		Q1 = self.Q
		w0 = self.parent.w
		w1 = self.w
		if(r==None):
			r = self.rest.to_quaternion()
		if(k==None):
			k = self.k

		ra = Q0.inverted() @ Q1
		if ra.dot(r) < 0:
			r = -r

		c, dQ0x,dQ0y,dQ0z,dQ0w,dQ1x,dQ1y,dQ1z,dQ1w = quatSpringGradient2(Q0,Q1,r)

		div = dQ0x*dQ0x*w0 + \
			dQ0y*dQ0y*w0 + \
			dQ0z*dQ0z*w0 + \
			dQ0w*dQ0w*w0 + \
			dQ1x*dQ1x*w1 + \
			dQ1y*dQ1y*w1 + \
			dQ1z*dQ1z*w1 + \
			dQ1w*dQ1w*w1

		if(div> 1e-8):
			s = -c/div
			if(w0>0.0):

				Q0.x+=dQ0x*s*w0*k
				Q0.y+=dQ0y*s*w0*k
				Q0.z+=dQ0z*s*w0*k
				Q0.w+=dQ0w*s*w0*k
				self.parent.Q = Q0.normalized()

			Q1.x+=dQ1x*s*w1*k
			Q1.y+=dQ1y*s*w1*k
			Q1.z+=dQ1z*s*w1*k
			Q1.w+=dQ1w*s*w1*k
			self.Q = Q1.normalized()

	@staticmethod
	def computeMatrixK(connector, invMass, x, inertiaInverseW):
		#adapted from https://github.com/InteractiveComputerGraphics/PositionBasedDynamics/blob/master/PositionBasedDynamics/PositionBasedRigidBodyDynamics.cpp
		K = Matrix().to_3x3()

		if (invMass == 0.0): 
			K.zero()
			return K

		v = connector - x
		a = v[0]
		b = v[1]
		c = v[2]
		j11 = inertiaInverseW[0][0]
		j12 = inertiaInverseW[1][0]
		j13 = inertiaInverseW[2][0]
		j22 = inertiaInverseW[1][1]
		j23 = inertiaInverseW[2][1]
		j33 = inertiaInverseW[2][2]

		K[0][0] = c*c*j22 - b*c*(j23 + j23) + b*b*j33 + invMass
		K[1][0] = -(c*c*j12) + a*c*j23 + b*c*j13 - a*b*j33
		K[2][0] = b*c*j12 - a*c*j22 - b*b*j13 + a*b*j23
		K[0][1] = K[1][0]
		K[1][1] = c*c*j11 - a*c*(j13 + j13) + a*a*j33 + invMass
		K[2][1] = -(b*c*j11) + a*c*j12 + a*b*j13 - a*a*j23
		K[0][2] = K[2][0]
		K[1][2] = K[2][1]
		K[2][2] = b*b*j11 - a*b*(j12 + j12) + a*a*j22 + invMass
			
		return K

def get_jiggle_children(arm_matrix, pbone, lst, parent, children_of_bone):
	""" Recursive function to build a flat list of JiggleBone objects. """
	jiggle_bone = JiggleBone(pbone, arm_matrix @ pbone.matrix, parent)
	lst.append(jiggle_bone)
	for c in children_of_bone[pbone]:
		get_jiggle_children(arm_matrix, c, lst, jiggle_bone, children_of_bone)

def maxis(M, i):
	return Vector((
		M[0][i], 
		M[1][i], 
		M[2][i])
	)
def saxis(M, i, v):
	M[0][i] = v[0]
	M[1][i] = v[1]
	M[2][i] = v[2]
def qadd(a, b):
	return Quaternion((
		a[0] + b[0], 
		a[1] + b[1], 
		a[2] + b[2], 
		a[3] + b[3])
	)

#NOTE: the following gradient computation implementation was automatically generated, if possible, it should be change for a clearer implementation
def quatSpringGradient2(Q0, Q1, r):
	"""Returns the gradient of C = |Q0*r - Q1|^2 wrt Q0 and Q1"""
	Q0x = Q0.x
	Q0y = Q0.y
	Q0z = Q0.z
	Q0w = Q0.w
	Q1x = Q1.x
	Q1y = Q1.y
	Q1z = Q1.z
	Q1w = Q1.w
	rx  =  r.x
	ry  =  r.y
	rz  =  r.z
	rw  =  r.w

	tmp0 = sqrt(((((((((-(Q0x*Q1w)-(Q0y*Q1z))+(Q0w*Q1x))+(Q0z*Q1y))-rx)*((((-(Q0x*Q1w)-(Q0y*Q1z))+(Q0w*Q1x))+(Q0z*Q1y))-rx))+(((((-(Q0x*Q1y)-(Q0z*Q1w))+(Q0w*Q1z))+(Q0y*Q1x))-rz)*((((-(Q0x*Q1y)-(Q0z*Q1w))+(Q0w*Q1z))+(Q0y*Q1x))-rz)))+(((((-(Q0y*Q1w)-(Q0z*Q1x))+(Q0w*Q1y))+(Q0x*Q1z))-ry)*((((-(Q0y*Q1w)-(Q0z*Q1x))+(Q0w*Q1y))+(Q0x*Q1z))-ry)))+((((((Q0w*Q1w)+(Q0x*Q1x))+(Q0y*Q1y))+(Q0z*Q1z))-rw)*(((((Q0w*Q1w)+(Q0x*Q1x))+(Q0y*Q1y))+(Q0z*Q1z))-rw))))
	tmp1 = 1.0/tmp0*Q0w*Q0y
	tmp2 = 1.0/tmp0*Q0w*Q1x
	tmp3 = 1.0/tmp0*Q0w*Q0x
	tmp4 = 1.0/tmp0*Q0x*Q1w
	tmp5 = 1.0/tmp0*Q0w*Q1w
	tmp6 = 1.0/tmp0*Q0y*Q1w
	tmp7 = 1.0/tmp0*Q0w*Q0z
	tmp8 = 1.0/tmp0*Q0x*Q1x
	tmp9 = 1.0/tmp0*Q0y*Q1x
	tmp10 = 1.0/tmp0*Q0x*Q0y
	tmp11 = 1.0/tmp0*Q0x*Q0z
	tmp12 = 1.0/tmp0*Q0z*Q1w
	tmp13 = 1.0/tmp0*Q0z*Q1x
	tmp14 = 1.0/tmp0*Q0y*Q0z
	tmp15 = 1.0/tmp0*Q0w*Q0w
	tmp16 = Q1w*Q1w
	tmp17 = 1.0/tmp0*Q0x*Q0x
	tmp18 = Q1x*Q1x
	tmp19 = 1.0/tmp0*Q0y*Q0y
	tmp20 = 1.0/tmp0
	tmp21 = Q1y*Q1y
	tmp22 = tmp20*Q0z*Q0z
	tmp23 = Q1z*Q1z
	tmp24 = tmp20*Q0x
	tmp25 = tmp20*Q0y
	tmp26 = tmp4*Q1x
	tmp27 = tmp24*Q1y*Q1z
	tmp28 = tmp3*Q1y
	tmp29 = tmp20*Q0z
	tmp30 = tmp3*Q1z
	tmp31 = tmp3*Q1w
	tmp32 = tmp5*Q1y
	tmp33 = tmp5*Q1z
	tmp34 = tmp1*Q1z
	tmp35 = tmp5*Q1x
	tmp36 = tmp1*Q1x
	tmp37 = tmp1*Q1w
	tmp38 = tmp6*Q1y
	tmp39 = tmp7*Q1y
	tmp40 = tmp2*Q1z
	tmp41 = tmp7*Q1x
	tmp42 = tmp9*Q1z
	tmp43 = tmp2*Q1y
	tmp44 = tmp3*Q1x
	tmp45 = tmp7*Q1w
	tmp46 = tmp20*Q0w*Q1y*Q1z
	tmp47 = tmp10*Q1x
	tmp48 = tmp4*Q1z
	tmp49 = tmp10*Q1y
	tmp50 = tmp10*Q1w
	tmp51 = tmp6*Q1z
	tmp52 = tmp4*Q1y
	tmp53 = tmp1*Q1y
	tmp54 = tmp12*Q1z
	tmp55 = -Q0x*Q1w-Q0y*Q1z+Q0w*Q1x+Q0z*Q1y-rx
	tmp56 = tmp20*Q1w
	tmp57 = tmp11*Q1z
	tmp58 = tmp9*Q1y
	tmp59 = tmp7*Q1z
	tmp60 = tmp11*Q1x
	tmp61 = tmp8*Q1y
	tmp62 = tmp13*Q1y
	tmp63 = tmp11*Q1w
	tmp64 = tmp8*Q1z
	tmp65 = -Q0x*Q1y-Q0z*Q1w+Q0w*Q1z+Q0y*Q1x-rz
	tmp66 = tmp14*Q1y
	tmp67 = tmp14*Q1w
	tmp68 = tmp12*Q1y
	tmp69 = tmp14*Q1z
	tmp70 = tmp20*Q1x
	tmp71 = -Q0y*Q1w-Q0z*Q1x+Q0w*Q1y+Q0x*Q1z-ry
	tmp72 = tmp6*Q1x
	tmp73 = tmp10*Q1z
	tmp74 = tmp12*Q1x
	tmp75 = tmp20*Q1y
	tmp76 = Q0w*Q1w+Q0x*Q1x+Q0y*Q1y+Q0z*Q1z-rw
	tmp77 = tmp29*Q1y*Q1z
	tmp78 = tmp25*Q1y*Q1z
	tmp79 = tmp13*Q1z
	tmp80 = tmp11*Q1y
	tmp81 = tmp20*Q0w
	tmp82 = tmp20*Q1z
	tmp83 = tmp14*Q1x
	c = tmp0
	dQ0x = tmp35+tmp46+tmp51+tmp58+tmp68+tmp79+tmp24*tmp16+tmp24*tmp18+tmp24*tmp21+tmp24*tmp23+tmp56*rx+tmp75*rz-tmp35-tmp46-tmp51-tmp58-tmp68-tmp79-tmp70*rw-tmp82*ry
	dQ0y = tmp32+tmp40+tmp48+tmp61+tmp74+tmp77+tmp25*tmp16+tmp25*tmp18+tmp25*tmp21+tmp25*tmp23+tmp56*ry+tmp82*rx-tmp32-tmp40-tmp48-tmp61-tmp74-tmp77-tmp70*rz-tmp75*rw
	dQ0z = tmp33+tmp43+tmp52+tmp64+tmp72+tmp78+tmp29*tmp16+tmp29*tmp18+tmp29*tmp21+tmp29*tmp23+tmp56*rz+tmp70*ry-tmp33-tmp43-tmp52-tmp64-tmp72-tmp78-tmp75*rx-tmp82*rw
	dQ0w = tmp26+tmp27+tmp38+tmp42+tmp54+tmp62+tmp81*tmp16+tmp81*tmp18+tmp81*tmp21+tmp81*tmp23-tmp26-tmp27-tmp38-tmp42-tmp54-tmp62-tmp56*rw-tmp70*rx-tmp75*ry-tmp82*rz
	dQ1x = tmp31+tmp34+tmp39+tmp49+tmp57+tmp67+tmp15*Q1x+tmp17*Q1x+tmp19*Q1x+tmp22*Q1x+tmp29*ry-tmp31-tmp34-tmp39-tmp49-tmp57-tmp67-tmp81*rx-tmp24*rw-tmp25*rz
	dQ1y = tmp30+tmp37+tmp41+tmp47+tmp63+tmp69+tmp15*Q1y+tmp17*Q1y+tmp19*Q1y+tmp22*Q1y+tmp24*rz-tmp30-tmp37-tmp41-tmp47-tmp63-tmp69-tmp81*ry-tmp25*rw-tmp29*rx
	dQ1z = tmp28+tmp36+tmp45+tmp50+tmp60+tmp66+tmp15*Q1z+tmp17*Q1z+tmp19*Q1z+tmp22*Q1z+tmp25*rx-tmp28-tmp36-tmp45-tmp50-tmp60-tmp66-tmp81*rz-tmp24*ry-tmp29*rw
	dQ1w = tmp44+tmp53+tmp59+tmp73+tmp80+tmp83+tmp15*Q1w+tmp17*Q1w+tmp19*Q1w+tmp22*Q1w+tmp24*rx+tmp25*ry+tmp29*rz-tmp44-tmp53-tmp59-tmp73-tmp80-tmp83-tmp81*rw

	return c, dQ0x,dQ0y,dQ0z,dQ0w,dQ1x,dQ1y,dQ1z,dQ1w

def step(scene=bpy.context.scene):
	dt = 1.0/(scene.render.fps)

	jiggle_armatures = [o for o in scene.objects if o.type=='ARMATURE' and o.data.jiggle.enabled]
	for arm_obj in jiggle_armatures:
		arm_data = arm_obj.data
		arm_matrix = arm_obj.matrix_world.copy()
		scale =maxis(arm_matrix,0).length
		
		children_of_bone = defaultdict(list)
		for bone in arm_obj.pose.bones:
			if(bone.parent):
				children_of_bone[bone.parent].append(bone)

		arm_data.jiggle.time_acc += dt * arm_data.jiggle.fps
		while arm_data.jiggle.time_acc > 1:
			arm_data.jiggle.time_acc -= 1

			jiggle_bones = []	# List that will store our JiggleBone objects. (We create a JiggleBone for every bone, not just jigglebones... which makes no sense, so TODO)
			for b in arm_obj.pose.bones:
				if(b.parent==None):
					get_jiggle_children(arm_matrix, b, jiggle_bones, None, children_of_bone)

			bl2 = []	# TODO How does this differ from jiggle_bones?
			for jb in jiggle_bones:
				b = jb.b
				jb.rest_w = b.bone.matrix_local.copy()
				saxis(jb.rest_w, 3, maxis(jb.rest_w, 3) * scale)
				saxis(jb.rest_w, 3, maxis(jb.rest_w, 3)+maxis(jb.rest_w, 1) * b.bone.length * 0.5 * scale)

			for jb in jiggle_bones:
				b = jb.b

				jb.restW = b.bone.matrix_local.copy() * scale
				saxis(jb.restW, 3, maxis(jb.restW, 3) * scale)

				M = jb.M
				if(b.bone.jiggle_enabled):
					db = b.bone
					jb.X = jb.P = db.jiggle_P
					jb.R = jb.Q = db.jiggle_R
					jb.rest = jb.rest_w
					if(b.parent):
						jb.rest = jb.parent.rest_w.inverted() @ jb.rest_w

					jb.rest_base = b.bone.matrix_local
					if(b.parent):
						jb.rest_base = b.parent.bone.matrix_local.inverted() @ jb.rest_base

					jb.rest_p = jb.parent.rest_w.inverted() @ (maxis(jb.rest_w, 3) - maxis(jb.rest_w, 1) * b.bone.length * 0.5 * scale) # mpos(jb.rest)

					jb.l = b.bone.length * scale
					jb.w = 1.0/db.jiggle_mass
					jb.k = 1 - pow(1 - db.jiggle_Ks, 1/scene.jiggle.iterations)
					db.jiggle_V *= 1.0 - db.jiggle_Kld
					db.jiggle_V += scene.gravity * db.gravity_multiplier * dt
					db.jiggle_W *= 1.0 - db.jiggle_Kd
					qv = Quaternion()
					qv.x = db.jiggle_W[0]
					qv.y = db.jiggle_W[1]
					qv.z = db.jiggle_W[2]
					qv.w = 0

					jb.Q = qadd(jb.Q, qv @ jb.Q * dt * 0.5).normalized() #removing * 0.5 seems to make single bone jiggle better but breaks arrays of bones

					jb.P = jb.X + db.jiggle_V * dt
					jb.computeI()

					#control object/bone constraint
					if(db.jiggle_control_object in bpy.data.objects):

						target_object = bpy.data.objects[db.jiggle_control_object]

						target_matrix = target_object.matrix_local

						if(target_object.type =='ARMATURE' and db.jiggle_control_bone in target_object.pose.bones):
							control_pb = target_object.pose.bones[db.jiggle_control_bone]
							target_matrix = control_pb.matrix
							if(control_pb.parent):
								target_matrix = control_pb.parent.matrix.inverted() @ target_matrix

						jb.cQ = target_matrix.to_quaternion().normalized()
						jb.Kc = 1- pow(1-db.jiggle_control, 1.0/scene.jiggle.iterations)

					bl2.append(jb)
				else:
					jb.w = 0
					jb.X = jb.P = M.translation + maxis(M, 1) * b.bone.length * 0.5
					jb.R = jb.Q = M.to_quaternion().normalized()

					M = arm_matrix @ b.matrix

					db = b.bone
					setq(db.jiggle_R, M.to_quaternion().normalized())
					db.jiggle_V = Vector((0,0,0))
					db.jiggle_P = M.translation + maxis(M, 1) * b.bone.length * 0.5
					db.jiggle_W = Vector((0,0,0))

			for i in range(scene.jiggle.iterations):
				#parent constraint
				for jb in bl2:
					if(jb.b.parent==None):
						continue
					jb.locSpring()
				#spring constraint
				for jb in bl2:
					if(jb.b.parent==None):
						continue
					jb.quatSpring(db.jiggle_rest if db.jiggle_use_custom_rest else jb.rest.to_quaternion().normalized())
					if(jb.cQ):
						jb.quatSpring(jb.cQ, jb.Kc)

			for jb in bl2:
				db = jb.b.bone

				jb.Q = jb.Q.normalized()
				m = jb.Q.to_matrix()
				for i in range(3):
					for j in range(3):
						jb.M[i][j] = m[i][j]*scale
				jb.M[3][3]=1

				db.jiggle_V = (jb.P - jb.X)/dt
				db.jiggle_P = jb.P.copy()
				qv = jb.Q @ db.jiggle_R.conjugated()
				db.jiggle_W = Vector((qv.x,qv.y,qv.z))*(2/dt)
				db.jiggle_R = jb.Q

				cp = db.jiggle_P - maxis(jb.M,1) * db.length * 0.5	# TODO: What's this?

				jb.M[0][3] = cp[0]
				jb.M[1][3] = cp[1]
				jb.M[2][3] = cp[2]

			for jb in bl2:
				pb = jb.b
				pM = arm_matrix
				if(pb.parent):
					pM = jb.parent.M
				mb = (pM @ jb.rest_base).inverted() @ jb.M

				pb.matrix_basis = mb

	scene.jiggle.last_frame+= 1

def bake(bake_all):
	print("Baking " + ("all" if(bake_all) else "selected") + "...")

	scene = bpy.context.scene
	scene.frame_set(scene.frame_start)
	
	jiggle_armatures = [o for o in scene.objects if o.type=='ARMATURE' and (o.select_get() or bake_all)]
	
	for i in range(scene.frame_start, scene.frame_end+1):
		print("Frame " + str(i))
		scene.frame_set(i)
		for arm in jiggle_armatures:
			for b in arm.pose.bones:
				b.bone.select = (b.bone.select or bake_all) and b.bone.jiggle_enabled
				if(b.bone.select):
					M = arm.matrix_world @ b.matrix
					db = b.bone
					setq(db.jiggle_R, M.to_quaternion().normalized())
					db.jiggle_V = Vector((0,0,0))
					db.jiggle_P = M.translation + maxis(M, 1) * db.length * 0.5
					db.jiggle_W = Vector((0,0,0))
			bpy.ops.anim.keyframe_insert_menu(type='LocRotScale')

class JARM_OT_bake(bpy.types.Operator):
	"""Bake jiggle bone motion to keyframes"""
	a: BoolProperty()
	bl_idname = "jiggle.bake"
	bl_label = "Bake Jiggle Physics"

	def execute(self, context):
		bake(self.a)
		return {'FINISHED'}

classes = (
	JARM_PT_armature,
	JiggleScene,
	JARM_PT_scene,
	JiggleArmature,
	JARM_OT_bake,
	JARM_OT_set_rest,
	JARM_PT_bone
)

def register():
	from bpy.utils import register_class
	for cls in classes:
		register_class(cls)

	bpy.app.handlers.frame_change_post.append(step)

	bpy.types.Scene.jiggle = PointerProperty(type = JiggleScene)

	bpy.types.Armature.jiggle = PointerProperty(type = JiggleArmature,options={'ANIMATABLE'})

	bpy.types.Bone.jiggle_enabled = BoolProperty(
		name = "Enable Jiggle",
		description = "Enable Jiggle Physics",
		default = False)
	bpy.types.Bone.jiggle_Kld = FloatProperty(
		name = "Linear Damping",
		min = 0.0,
		max = 1.0,
		default = 0.01)
	bpy.types.Bone.jiggle_Kd = FloatProperty(
		name = "Angular Damping",
		min = 0.0,
		max = 1.0,
		default = 0.01)
	bpy.types.Bone.jiggle_Ks = FloatProperty(
		name = "Stiffness",
		min = 0.0,
		max = 1.0,
		default = 0.8)
	bpy.types.Bone.gravity_multiplier = FloatProperty(
		name = "Gravity Multiplier",
		min = 0.0,
		default = 1.0)
	bpy.types.Bone.jiggle_mass = FloatProperty(
		name = "Mass",
		min = 0.0001,
		default = 1.0)
	bpy.types.Bone.jiggle_R = FloatVectorProperty(
		name = "Rotation",
		size = 4,
		subtype = 'QUATERNION')

	bpy.types.Bone.jiggle_W = FloatVectorProperty(size=3,subtype='XYZ') #angular velocity
	bpy.types.Bone.jiggle_P = FloatVectorProperty(size=3,subtype='XYZ')
	bpy.types.Bone.jiggle_V = FloatVectorProperty(size=3,subtype='XYZ') #linear velocity, ok?

	bpy.types.Bone.jiggle_use_custom_rest =BoolProperty(
		name = "Use Custom Rest Pose",
		default = False)
	bpy.types.Bone.jiggle_rest = FloatVectorProperty(
		name = "Rotation",
		size = 4,
		subtype = 'QUATERNION')
	bpy.types.Bone.jiggle_control = FloatProperty(
		name = "Control",
		min = 0.0,
		max = 1.0,
		default = 1)
	bpy.types.Bone.jiggle_control_object = StringProperty(name = "Control Object")
	bpy.types.Bone.jiggle_control_bone = StringProperty(name = "Control Bone")

def unregister():
	from bpy.utils import unregister_class
	for cls in reversed(classes):
		unregister_class(cls)
	bpy.app.handlers.frame_change_post.remove(step)

if __name__ == '__main__':
	register()
