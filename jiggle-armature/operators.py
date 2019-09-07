import bpy
from mathutils import *
from bpy.props import *
from .utils import *

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

class JARM_OT_bake(bpy.types.Operator):
	"""Bake jiggle simulation to keyframes"""

	bl_idname = "jiggle.bake"
	bl_label = "Bake Jiggle Physics"

	bake_all: BoolProperty()

	def execute(self, context):
		bake_all = self.bake_all
		print("Baking " + ("all" if(bake_all) else "selected") + "...")

		scene = context.scene
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

		return {'FINISHED'}

classes = (
	JARM_OT_bake,
	JARM_OT_set_rest,
)

def register():
	from bpy.utils import register_class
	for cls in classes:
		register_class(cls)

def unregister():
	from bpy.utils import unregister_class
	for cls in classes:
		unregister_class(cls)