from mathutils import *

def setq(om, m):
	for i in range(4):
		om[i]= m[i]

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