import bpy
from bpy.props import *

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
		col.prop(context.object.data.jiggle,"iterations")
		col.operator("jiggle.bake", text="Bake Selected").bake_all = False
		col.operator("jiggle.bake", text="Bake All").bake_all = True

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

classes = (
	JARM_PT_armature,
	JARM_PT_bone,
)

def register():
	from bpy.utils import register_class
	for cls in classes:
		register_class(cls)

def unregister():
	from bpy.utils import unregister_class
	for cls in classes:
		unregister_class(cls)