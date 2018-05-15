#!BPY
""" 
Name: 'Minimal Surface (by SJSimmons 2008)' 
Blender: 263 
Group: 'Mesh' 
Tooltip: 'computes minimal surface' 
""" 
numtimestosubdivide = 3
subdivmeth = 1  # triangulation method:  1  subdivides each face into 6 faces. 
                #                        2  subdivides each face into 4 faces. 
omega = 1      # reduce this (divide by 5, say) if surface blows up.
zerotol = 0.05 # reduce this for more accuracy in the final approximation.
showsurfarea = 0  # If 1, shows total surface area of each approximation.
setmesh = 1

# NOTE: This doesn't run (by a long shot) in newer versions fo Blender  
#
# Todo:
#    - remove removedoubles - done
#    - remove errant edges - done
#    - fix and finish description
#    - futher subdivide  - done
#    - futher subdivide up to some threshhold 
#    - try to automatically set threshold 
#    - write your own initial triangulation routine
#    - iterate until only negligible change - done (uses zerotol on max of derivative
#    - handle boundaries that are not connected (and maybe not even simply connected)
#
# Some References [1] Numerical Solution of Partial Differential Equations  (2nd Ed.)
#			          by K.W.Morton and D.F.Mayers.
#                 [2] Numerical Solution of Plateau's Problem byh a Finite Element 
#                     Method by Masahiro Hinata, Masaaki Shimasaki and Takeshi Kiyono.
#		          [3] Quadratic Functions, Optimization, and Quadratic Forms - Robert 
#                     M. Fruend
#
# Requires a python distribution along with the additional numpy python package.
#
# Description of method:
#
#     Take a selected closed edge path (which currently needs to be isomorphic to S1
#     in such a way that its projection onto the xy-plane is also isomophic to S1). 
#     Triangulate said projection (currently using a built-in function).  The
#     projection is called pme.  Give pme the original boundary data. 
#
#     Then cleanly subdivide, where cleanly means adding only the necessary vertices.  
#     Also, after each subdivision the faces and edges from the previous refinement 
#     level are removed. 
#
#     Note:  The number of subdivisions is controled by setting numtimestosubdivide 
#            below.          
#     

#from Blender import Scene, Mesh, Window, sys 
#from Blender import * 
#import Blender
#import BPyMessages 
import bpy
import mathutils
from numpy import *
#from copy import copy

Vector= mathutils.Vector

def vertedgeskel(mesh,v):  # number edges containing a given vertex v in mesh mesh
	count = 0
	for e in mesh.edges:
		if v == e.v1 or v == e.v2: 	
			count+=1
	return count	

def vertfaceskel(mesh,v):  # number faces containing a given vertex v in mesh mesh
	count = 0
	for f in mesh.faces:
		if v in f.verts: 	
			count+=1
	return count	
 
def my_mesh_util(me): 
	# This function runs out of editmode with a mesh error cases are already checked for 

	if setmesh != 1:
		#Find the most center (projected) vertex of the starting loop.  Not sure if this is always 'center most'.
		mind=10000
		minindex=0
		indexcount=0
		for v in me.vertices:
			d = 0.0
			for w in me.vertices:
				#Compute sum of distances of projection of v to projection of w
				d += ((v.co.x-w.co.x)**2+(v.co.y-w.co.y)**2)**0.5
			if d < mind:
				mind=d
				mindindex=indexcount
			indexcount+=1

		print('Finding projection in domain.')

		nextx=me.vertices[minindex].co.x
		nexty=me.vertices[minindex].co.y
		nextz=me.vertices[minindex].co.z

		print('Finding polyline in domain.')

		polyline=[Vector((nextx,nexty,0))]

		edgecount=1
		print('Original number of edges:',len(me.edges))

		while edgecount < len(me.edges):
			for e in me.edges:
				if e.select and edgecount < len(me.edges):
					print(e.v1.co)
					if e.vertices[1].co == Vector((nextx,nexty,nextz)) or e.vertices[2].co == Vector((nextx,nexty,nextz)): 
						if edgecount == 1:
							if e.v1.co == Vector(nextx,nexty,nextz):	
								polyline= polyline+[Vector(e.v2.co.x,e.v2.co.y,0)]
								nextx=e.v2.co.x
								nexty=e.v2.co.y
								nextz=e.v2.co.z
							else:
								polyline= polyline+[Vector(e.v1.co.x,e.v1.co.y,0)]
								nextx=e.v1.co.x
								nexty=e.v1.co.y
								nextz=e.v1.co.z
							edgecount +=1
							laste=e
						elif laste != e:
							if e.v1.co == Vector(nextx,nexty,nextz):	
								polyline= polyline+[Vector(e.v2.co.x,e.v2.co.y,0)]
								nextx=e.v2.co.x
								nexty=e.v2.co.y
								nextz=e.v2.co.z
							else:
								polyline= polyline+[Vector(e.v1.co.x,e.v1.co.y,0)]
								nextx=e.v1.co.x
								nexty=e.v1.co.y
								nextz=e.v1.co.z
							edgecount += 1
							laste=e


		print('Polyfilling (builtin)')
		fill= Blender.Geometry.PolyFill([polyline])  #triangulate using built in funcion.
		print('  done')

		pme=Blender.Mesh.New()   # pme is the projection projection 
		pme.verts.extend(polyline)
		pme.faces.extend(fill)

		scn = Blender.Scene.GetCurrent()
		ob = scn.objects.new(pme)
		Blender.Redraw()
	else:
		pme=Blender.Mesh.New()   # pme is the projection projection 

		#Could just do pme=me?
		pme=me
		#pme.verts.extend([v.co for v in me.verts])
		#pme.faces.extend([[f.v[0].index, f.v[1].index, f.v[2].index] for f in me.faces])

		scn = Blender.Scene.GetCurrent()
		ob = scn.objects.new(pme)
		Blender.Redraw()


	print('Before subdividing: ', len(pme.verts), 'verts,', len(pme.edges), 'edges and ', len(pme.faces), 'faces.')

	def subdiv(pme):
		print('Subdividing.')
		numfaces=len(pme.faces)
		counter=1
		if subdivmeth == 1:
			for f in pme.faces:
				if counter <= numfaces:
					v1=tuple((f.v[0].co + f.v[1].co)/2)
					v2=tuple((f.v[1].co + f.v[2].co)/2)
					v3=tuple((f.v[2].co + f.v[0].co)/2)
					f0=f.v[0].index
					f1=f.v[1].index
					f2=f.v[2].index
					fc=f.cent
					vertsset=[tuple(pme.verts[i].co) for i in range(len(pme.verts))]
					if v1 not in vertsset and v2 not in vertsset and v3 not in vertsset:
						pme.verts.extend([v1,v2,v3,fc])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,f0,vind-4])
						pme.faces.extend([vind-1,vind-4,f1])
						pme.faces.extend([vind-1,vind-3,f1])
						pme.faces.extend([vind-1,f2,vind-3])
						pme.faces.extend([vind-1,f2,vind-2])
						pme.faces.extend([vind-1,vind-2,f0])
					elif v1 in vertsset and v2 not in vertsset and v3 not in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						pme.verts.extend([v2,v3,fc])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,f0,v1])
						pme.faces.extend([vind-1,v1,f1])
						pme.faces.extend([vind-1,vind-3,f1])
						pme.faces.extend([vind-1,f2,vind-3])
						pme.faces.extend([vind-1,f2,vind-2])
						pme.faces.extend([vind-1,vind-2,f0])
					elif v1 not in vertsset and v2 in vertsset and v3 not in vertsset:
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						pme.verts.extend([v1,v3,fc])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,f0,vind-3])
						pme.faces.extend([vind-1,vind-3,f1])
						pme.faces.extend([vind-1,v2,f1])
						pme.faces.extend([vind-1,f2,v2])
						pme.faces.extend([vind-1,f2,vind-2])
						pme.faces.extend([vind-1,vind-2,f0])
					elif v1 not in vertsset and v2 not in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.extend([v1,v2,fc])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,f0,vind-3])
						pme.faces.extend([vind-1,vind-3,f1])
						pme.faces.extend([vind-1,vind-2,f1])
						pme.faces.extend([vind-1,f2,vind-2])
						pme.faces.extend([vind-1,f2,v3])
						pme.faces.extend([vind-1,v3,f0])
					elif v1 in vertsset and v2 in vertsset and v3 not in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						pme.verts.extend([v3,fc])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,f0,v1])
						pme.faces.extend([vind-1,v1,f1])
						pme.faces.extend([vind-1,v2,f1])
						pme.faces.extend([vind-1,f2,v2])
						pme.faces.extend([vind-1,f2,vind-2])
						pme.faces.extend([vind-1,vind-2,f0])
					elif v1 in vertsset and v2 not in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.extend([v2,fc])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,f0,v1])
						pme.faces.extend([vind-1,v1,f1])
						pme.faces.extend([vind-1,vind-2,f1])
						pme.faces.extend([vind-1,f2,vind-2])
						pme.faces.extend([vind-1,f2,v3])
						pme.faces.extend([vind-1,v3,f0])
					elif v1 not in vertsset and v2 in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.extend([v1,fc])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,f0,vind-2])
						pme.faces.extend([vind-1,vind-2,f1])
						pme.faces.extend([vind-1,v2,f1])
						pme.faces.extend([vind-1,f2,v2])
						pme.faces.extend([vind-1,f2,v3])
						pme.faces.extend([vind-1,v3,f0])
					elif v1 in vertsset and v2 in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.extend([fc])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,f0,v1])
						pme.faces.extend([vind-1,v1,f1])
						pme.faces.extend([vind-1,v2,f1])
						pme.faces.extend([vind-1,f2,v2])
						pme.faces.extend([vind-1,f2,v3])
						pme.faces.extend([vind-1,v3,f0])
					counter += 1
					Blender.Redraw()
		elif subdivmeth == 2:
			for f in pme.faces:
				if counter <= numfaces:
					v1=tuple((f.v[0].co + f.v[1].co)/2)
					v2=tuple((f.v[1].co + f.v[2].co)/2)
					v3=tuple((f.v[2].co + f.v[0].co)/2)
					f0=f.v[0].index
					f1=f.v[1].index
					f2=f.v[2].index
					vertsset=[tuple(pme.verts[i].co) for i in range(len(pme.verts))]
					if v1 not in vertsset and v2 not in vertsset and v3 not in vertsset:
						pme.verts.extend([v1,v2,v3])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,vind-2,vind-3])
						pme.faces.extend([f0,vind-1,vind-3])
						pme.faces.extend([f1,vind-2,vind-3])
						pme.faces.extend([f2,vind-1,vind-2])
					elif v1 in vertsset and v2 not in vertsset and v3 not in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						pme.verts.extend([v2,v3])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,vind-2,v1])
						pme.faces.extend([f0,v1,vind-1])
						pme.faces.extend([f1,v1,vind-2])
						pme.faces.extend([f2,vind-2,vind-1])
					elif v1 not in vertsset and v2 in vertsset and v3 not in vertsset:
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						pme.verts.extend([v1,v3])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,v2,vind-2])
						pme.faces.extend([f0,vind-1,vind-2])
						pme.faces.extend([f1,vind-2,v2])
						pme.faces.extend([f2,vind-1,v2])
					elif v1 not in vertsset and v2 not in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.extend([v1,v2])
						vind=len(pme.verts)
						pme.faces.extend([vind-1,vind-2,v3])
						pme.faces.extend([f0,vind-2,v3])
						pme.faces.extend([f1,vind-2,vind-1])
						pme.faces.extend([f2,v3,vind-1])
					elif v1 in vertsset and v2 in vertsset and v3 not in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						pme.verts.extend([v3])
						vind=len(pme.verts)
						pme.faces.extend([v1,v2,vind-1])
						pme.faces.extend([f0,v1,vind-1])
						pme.faces.extend([f1,v1,v2])
						pme.faces.extend([f2,vind-1,v2])
					elif v1 in vertsset and v2 not in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.extend([v2])
						vind=len(pme.verts)
						pme.faces.extend([v1,v3,vind-1])
						pme.faces.extend([f0,v1,v3])
						pme.faces.extend([f1,v1,vind-1])
						pme.faces.extend([f2,v3,vind-1])
					elif v1 not in vertsset and v2 in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.extend([v1])
						vind=len(pme.verts)
						pme.faces.extend([v2,v3,vind-1])
						pme.faces.extend([f0,vind-1,v3])
						pme.faces.extend([f1,v2,vind-1])
						pme.faces.extend([f2,v3,v2])
					elif v1 in vertsset and v2 in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.faces.extend([v2,v3,v1])
						pme.faces.extend([f0,v1,v3])
						pme.faces.extend([f1,v2,v1])
						pme.faces.extend([f2,v3,v2])
					else:
						print('error')
					counter += 1
					Blender.Redraw()
		else:
			print('subdiv method', subdivmeth, ' undefined')


		print(' After subdividing:',len(pme.verts),'verts',len(pme.edges),'edges',len(pme.faces),'faces.') 
		pme.faces.delete(1,[i for i in range(numfaces)])  #Removes the original faces (and edges)
		print(' Removing',numfaces,'(previous) faces. Now',len(pme.verts),'verts',len(pme.edges),'edges',len(pme.faces),'faces.') 
	# end of subdiv

	if setmesh != 1:
		for v in pme.verts:
			for i in range(len(me.verts)):     # Give pme.verts the original boundary data
				if v.co.x == me.verts[i].co.x and v.co.y == me.verts[i].co.y:
					v.co.z = me.verts[i].co.z

	for i in range(numtimestosubdivide):
		subdiv(pme)
		print('Redrawing.')
		Blender.Redraw()  #draw the projection 

	Blender.Redraw()  #draw the projection 

	# Start integrating

	print('Integrating to find starting approximation.')

	print('Finding boundary verts ...',)
	bdverts=[] 
	for v in pme.verts:
		if vertedgeskel(pme,v) > vertfaceskel(pme,v):
			bdverts=bdverts+[tuple(v.co)]  #tuples only
	print(' done.')

	bdvals=[]  # get boundary data for later.
	for v in bdverts:
		bdvals=bdvals+[v[2]]
	xbdvals=array(bdvals)

	dbdverts=[] # projection of boundary data vertices 
	for v in bdverts:
		dbdverts=dbdverts+[[v[0],v[1],0]]

	print('Finding interior verts.')

	dintverts=[] # (projection of) interior vertices
	for v in pme.verts:
		if vertedgeskel(pme,v) == vertfaceskel(pme,v):
			dintverts=dintverts+[[v.co.x,v.co.y,0]]
	print(' done.')

	print(' There are',len(bdverts),'boundary vertices and',len(dintverts),'interior verts.')

	dverts=dintverts+dbdverts  # vector version of pme.verts except with interior verts before bd verts.

	A = zeros((len(dverts),len(dverts)),dtype=float32)

	print('generating A ...',)            # for A, see ref [1], section 6.7,  6.114


	#Best way to generate A: 								

	for f in pme.faces:
		x0=f.v[0].co.x
		y0=f.v[0].co.y
		x1=f.v[1].co.x
		y1=f.v[1].co.y
		x2=f.v[2].co.x
		y2=f.v[2].co.y
		farea=2.0*abs(x0*(y2-y1)+x1*(y0-y2)+x2*(y1-y0))  # This is 4 times area of face.
		A[dverts.index([x0,y0,0]),dverts.index([x0,y0,0])]+=((y1-y2)**2+(x1-x2)**2)/farea
		A[dverts.index([x1,y1,0]),dverts.index([x1,y1,0])]+=((y0-y2)**2+(x2-x0)**2)/farea
		A[dverts.index([x2,y2,0]),dverts.index([x2,y2,0])]+=((y0-y1)**2+(x1-x0)**2)/farea
		A[dverts.index([x0,y0,0]),dverts.index([x1,y1,0])]+=((y1-y2)*(y2-y0)+(x2-x1)*(x0-x2))/farea
		A[dverts.index([x1,y1,0]),dverts.index([x2,y2,0])]+=((y2-y0)*(y0-y1)+(x0-x2)*(x1-x0))/farea
		A[dverts.index([x2,y2,0]),dverts.index([x0,y0,0])]+=((y1-y2)*(y0-y1)+(x2-x1)*(x1-x0))/farea
		A[dverts.index([x1,y1,0]),dverts.index([x0,y0,0])]=A[dverts.index([x0,y0,0]),dverts.index([x1,y1,0])]
		A[dverts.index([x2,y2,0]),dverts.index([x1,y1,0])]=A[dverts.index([x1,y1,0]),dverts.index([x2,y2,0])]
		A[dverts.index([x0,y0,0]),dverts.index([x2,y2,0])]=A[dverts.index([x2,y2,0]),dverts.index([x0,y0,0])]
	print(' done.')


	Aul=A[0:len(dintverts),0:len(dintverts)]             #A upper left
	Aur=A[0:len(dintverts),len(dintverts):len(dverts)]   #A upper right
	All=A[len(dintverts):len(dverts),0:len(dintverts)]   #A lower left

	set_printoptions(linewidth=150,precision=5)
	print('size of Aul', Aul.shape)
	#print 'size of Aur', Aur.shape
	#print 'size of All', All.shape
	#print 'size of xbdvals', xbdvals.shape
	#print 'xbdvals',xbdvals
#	print Aul
#	print 'diagonal ',
#	for i in range(len(dintverts)):
#		print Aul[i,i],
#	print ''

	tempsum = transpose(Aur)+All
	D = mat(tempsum.copy())
	k = mat(xbdvals.copy())
	b= k*D*0.5   #important: This 1/2 has to be here. 

	print('inverting A ... ',)
#	print 'det A =', linalg.det(A)
#	print 'det Aul =', linalg.det(Aul)
	xintvals=linalg.solve(Aul,-transpose(b))
	print(' done.')


	finalverts=[]  #final verts for the first approx to minimal surface
	for i in range(0,len(dintverts)):
		dintverts[i][2]=xintvals[i]
	for v in dintverts:
		finalverts=finalverts+[v]
	for v in bdverts:
		finalverts=finalverts+[v]

	for v in pme.verts:
		for vv in finalverts: 
			if v.co.x == vv[0] and v.co.y == vv[1]:
				v.co.z = vv[2]

	Blender.Redraw()  #redraw the first approximation
#	sys.sleep(1000)

	def newton():

		print('Minimizing surface area.') 

		# Let S(v_1,...,v_n) be the surface area of the current mesh
		# pme.  Compute partial dS of S w/r to v_i where v_i changes in
		# only the z-direction:

		def sumsq(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2): #This is what's under the square root in 2.8
			return (x0*(y1-y2)-x1*(y0-y2)+x2*(y0-y1))**2+((u0+z0)*(y1-y2)-(u1+z1)*(y0-y2)+(u2+z2)*(y0-y1))**2+(x0*((u1+z1)-(u2+z2))-x1*((u0+z0)-(u2+z2))+x2*((u0+z0)-(u1+z1)))**2

		def chain(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2): #chain rule part
			return (u0+z0)*((y1-y2)**2+(x1-x2)**2)+(u1+z1)*((y2-y0)*(y1-y2)+(x0-x2)*(x2-x1))+(u2+z2)*((y0-y1)*(y1-y2)+(x1-x0)*(x2-x1))

		def secder(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2):
			blahsumsq=sumsq(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2)
			blahchain=chain(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2)
			return -0.25*(1/blahsumsq)**1.5*blahchain**2+0.5*(1/blahsumsq)**0.5*((y1-y2)**2+(x2-x1)**2)
		
		def dS(z,q):
			sumt = 0.0
			firstder = 0.0
			secondder = 0.0
			v=intverts[q]  # differentiate w/r to v
#			print 'differentiating w/r to index ',q,' which is vertex',v.co.x,v.co.y,v.co.z,'  z=',z
			temp=v
			for f in pme.faces:	
				if v in f:	
					#always set (x0,y0,u0) to the vertex you are taking partial w/r to.
					if f.verts[0] == temp: 
						x0=f.verts[0].co.x
						y0=f.verts[0].co.y
						u0=f.verts[0].co.z
						z0=z[q]
						x1=f.verts[1].co.x
						y1=f.verts[1].co.y
						u1=f.verts[1].co.z
						z1=0
						cc=0
						for w in intverts:
							if w == f.verts[1]:
								z1+=z[cc]		
							cc+=1
						x2=f.verts[2].co.x
						y2=f.verts[2].co.y
						u2=f.verts[2].co.z
						z2=0
						cc=0
						for w in intverts:
							if w == f.verts[2]:
								z2+=z[cc]		
							cc+=1
					elif f.verts[1] == temp:
						x0=f.verts[1].co.x
						y0=f.verts[1].co.y
						u0=f.verts[1].co.z
						z0=z[q]
						x1=f.verts[0].co.x
						y1=f.verts[0].co.y
						u1=f.verts[0].co.z
						z1=0
						cc=0
						for w in intverts:
							if w == f.verts[0]:
								z1+=z[cc]		
							cc+=1
						x2=f.verts[2].co.x
						y2=f.verts[2].co.y
						u2=f.verts[2].co.z
						z2=0
						cc=0
						for w in intverts:
							if w == f.verts[2]:
								z2+=z[cc]		
							cc+=1
					else:
						x0=f.verts[2].co.x
						y0=f.verts[2].co.y
						u0=f.verts[2].co.z
						z0=z[q]
						x1=f.verts[0].co.x
						y1=f.verts[0].co.y
						u1=f.verts[0].co.z
						z1=0
						cc=0
						for w in intverts:
							if w == f.verts[0]:
								z1+=z[cc]		
							cc+=1
						x2=f.verts[1].co.x
						y2=f.verts[1].co.y
						u2=f.verts[1].co.z
						z2=0
						cc=0
						for w in intverts:
							if w == f.verts[1]:
								z2+=z[cc]		
							cc+=1
	#				sumsqt=sumsq(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2)
	#				chaint=chain(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2)
	#				secdert=secder(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2)
					firstder += 0.5*(1/sumsq(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2)**(0.5))*chain(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2)
					secondder += secder(x0,y0,u0,z0,x1,y1,u1,z1,x2,y2,u2,z2)
			sumt= firstder/secondder 
#			print 'derivative',firstder,', sumt=',sumt
			return [sumt,firstder]	
		# end of ds

		z=[]
		intlen=len(intverts)
		print('length intverts', intlen)

		for i in range(intlen):  #Initialize z to be the correct length 
			z=z+[0.0] 

		maxder=1.0
		numiterations=0
		znext=z
		print('Iterating')
		while maxder > zerotol:
			numiterations+=1
			print(numiterations,maxder,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b',)
			z=znext
			maxder=0.0
			for q in range(intlen):
				[dds,fder]=dS(znext[0:q]+z[q:intlen],q)
				if abs(fder) > maxder:
					maxder = abs(fder)
				znext[q]=z[q]-omega*dds # generate next z vector
			for v in pme.verts:
				if v in intverts:  # set pme.verts to the new approximation
					cc=0
					for w in intverts:
						if w == v:
							v.co.z=origu[cc]+znext[cc]
						cc += 1
			if showsurfarea == 1:     # print surface area for this approximation
				surfarea=0.0
				for f in pme.faces:
					surfarea+=f.area
				print('surface area',numiterations,surfarea,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b',)
			Blender.Redraw()           # redraw with this approximation
			for v in pme.verts:
				if v in intverts:  # reset pme to first approximation but don't redraw ...
					cc=0
					for w in intverts:
						if w == v:
							v.co.z=origu[cc]
						cc += 1
			pme.update()
		for v in pme.verts:
			if v in intverts:  # set pme.verts to the new approximation
				cc=0
				for w in intverts:
					if w == v:
						v.co.z=origu[cc]+znext[cc]
					cc += 1
		if showsurfarea == 1:     # print surface area for this approximation
			surfarea=0.0
			for f in pme.faces:
				surfarea+=f.area
			print('final surface area',numiterations,surfarea,)
		print(numiterations,maxder)
		Blender.Redraw()           # redraw with this approximation
	#end of newton

	bdverts2=[] # (for newton) boundary data vertices. (not in terms of vectors only)
	for v in pme.verts:
		if vertedgeskel(pme,v) > vertfaceskel(pme,v):
			bdverts=bdverts+[v]
	origu=[]

	intverts=[] # interior data vertices    
	for v in pme.verts:
		if vertedgeskel(pme,v) == vertfaceskel(pme,v):
			intverts=intverts+[v]

	for v in intverts:
		origu=origu+[v.co.z]
#	for q in range(len(intverts)):
#		origu=origu+[intverts[q].co.z]

	print('Newton',)
	newton()


def main(): 

	print('-- jminsurf --')
	
	if bpy.app.version[0] < 2 or bpy.app.version[1] < 62:
		raise Exception("Only for Blenderr 2.62 and above.")
		
	#for ob in bpy.data.objects:
	#	print(ob.name)
		
	# print all scene names in a list
	#print(bpy.data.scenes.keys())

	
	# remove mesh Cube
	#if "Cube" in bpy.data.meshes:
	#	mesh = bpy.data.meshes["Cube"]
	#	print("removing mesh", mesh)
	#	bpy.data.meshes.remove(mesh)
	
	# bpy.ops.object.mode_set(mode='OBJECT')  # Can't access corrdinate data in edit mode
	
	# me = bpy.data.meshes['Cube']
	
	## Gets the current scene, there can be many scenes in 1 blend file. 
	#sce = bpy.data.scenes.active 
	 
	## Get the active object, there can only ever be 1 
	## and the active object is always the editmode object. 
	#ob_act = sce.objects.active 
	 
	#if not ob_act or ob_act.type != 'Mesh': 
	#	BPyMessages.Error_NoMeshActive() 
	#	return  
	 
	## Saves the editmode state and go's out of  
	## editmode if its enabled, we cant make 
	## changes to the mesh data while in editmode. 
	#is_editmode = Window.EditMode() 
	#if is_editmode: Window.EditMode(0) 

	#t = sys.time() 

	#if setmesh != 1:
	#	Window.WaitCursor(1) 
	#	me = ob_act.getData(mesh=1) # old NMesh api is default 

	if setmesh == 1:
		print('Using hard-coded mesh.')	
		polyline1= [Vector((-3.0, 3.0, 0.0)), Vector((-3.0, -3.0, 0.0)), Vector((3.0, -3.0, 0.0)), Vector((3.0, 3.0, 0.0))] 
		polyline2= [Vector((-2.0, 2.0, 0.0)), Vector((-2.0, -2.0, 0.0)), Vector((2.0, -2.0, 0.0)), Vector((2.0, 2.0, 0.0))] 
		fill= Blender.Geometry.PolyFill([polyline1, polyline2])

		me= Blender.Mesh.New()
		me.verts.extend(polyline1)
		me.verts.extend(polyline2)
		me.faces.extend(fill)

		me= bpy.data.meshes.new(name="blah")
		polyline1= [Vector((-8.0, 8.0, -3.0)), Vector((-8.0, -8.0, 0.0)), Vector((8.0, -8.0, -3.0)), Vector((8.0, 8.0, 0.0))] 
		polyline2= [Vector((-3.0, 0.0, 1.0)), Vector((0.0, -3.0, 3.0)), Vector((3.0, 0.0, 5.0)), Vector((0.0, 3.5, 4.0))] 
		me.verts.extend(polyline1+polyline2)
		me.faces.extend([[0,1,4],[1,2,5],[4,5,1],[5,2,6],[3,6,7],[7,0,4],[3,7,0],[3,2,6]]) 
		scn = Blender.Scene.GetCurrent() 
		ob = scn.objects.new(me)
		Blender.Redraw()
	
	
	my_mesh_util(me) 
	 
	 
	# Restore editmode if it was enabled 
	#if is_editmode: Window.EditMode(1) 
	 
	# Timing the script is a good way to be aware on any speed hits when scripting 
	#print('Script finished in %.2f seconds' % (sys.time()-t)) 
	#Window.WaitCursor(0) 
	 
	 
# This lets you can import the script without running it 
if __name__ == '__main__': 
	main() 
