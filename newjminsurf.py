# jminsurf   
#
# generate minimal surface
#
# This was written to be run in blender
# Copyright Scott Simmons 2012 

import bpy, bmesh, math, mathutils, numpy

def MakeInitFreiOtto1(me):  # Frei Otto 1
	#make verts
	srad=2
	lrad=6
	height=2
	nptscirc=8   # use even number here
	angl=2*math.pi/nptscirc
	bignptscirc=2*nptscirc
	bigangl=2*math.pi/bignptscirc
	verts=[]
	for i in range(nptscirc): # top small circle
		verts.append(mathutils.Vector((lrad/2+srad*math.cos(i*angl),srad*math.sin(i*angl),height)))
	for i in range(nptscirc): # bottom small circle
		verts.append(mathutils.Vector((-lrad/2+srad*math.cos(i*angl),srad*math.sin(i*angl),-height)))
	for i in range(bignptscirc): # large circle
		verts.append(mathutils.Vector((lrad*math.cos(i*bigangl),lrad*math.sin(i*bigangl),0)))
	edges=[]	#make edges
	for i in range(nptscirc-1): # top small circle
		edges.append([i,i+1])
	edges.append([nptscirc-1,0])# complete the top circle
	for i in range(nptscirc-1): # bottom small circle
		i+=nptscirc
		edges.append([i,i+1])
	edges.append([2*nptscirc-1,nptscirc])# complete the bottom circle
	for i in range(bignptscirc-1): # large circle
		i+=2*nptscirc
		edges.append([i,i+1])
	edges.append([bignptscirc+2*nptscirc-1,2*nptscirc])# complete the big circle
	
	me.from_pydata(verts,edges,[])

	#Get a BMesh representation
	bm=bmesh.new()  #Create empty BMesh
	bm.from_mesh(me) #fill it in from a Mesh

	if hasattr(bm.verts, "ensure_lookup_table"):
		bm.verts.ensure_lookup_table()
	for i in range(int(nptscirc/2)):  # Add faces
		bm.faces.new((bm.verts[i],bm.verts[i+1],bm.verts[2*nptscirc+i]))
		bm.faces.new((bm.verts[i+1],bm.verts[2*nptscirc+i],bm.verts[2*nptscirc+1+i]))
		bm.faces.new((bm.verts[nptscirc+i],bm.verts[nptscirc+1+i],bm.verts[2*nptscirc+int(nptscirc/2)+i]))
		bm.faces.new((bm.verts[nptscirc+1+i],bm.verts[2*nptscirc+int(nptscirc/2)+i],bm.verts[2*nptscirc+int(nptscirc/2)+1+i]))
		if i < int(nptscirc/2)-1:
			bm.faces.new((bm.verts[int(3*nptscirc/2)+i],bm.verts[int(3*nptscirc/2)+1+i],bm.verts[3*nptscirc+i]))
		else:
			bm.faces.new((bm.verts[2*nptscirc-1],bm.verts[nptscirc],bm.verts[int(7*nptscirc/2)-1]))
		if i < int(nptscirc/2)-1:	
			bm.faces.new((bm.verts[int(3*nptscirc/2)+1+i],bm.verts[3*nptscirc+i],bm.verts[3*nptscirc+1+i]))
		else:
			bm.faces.new((bm.verts[nptscirc],bm.verts[int(7*nptscirc/2)-1],bm.verts[int(7*nptscirc/2)]))
		if i < int(nptscirc/2)-1:	
			bm.faces.new((bm.verts[int(nptscirc/2)+i],bm.verts[int(nptscirc/2)+1+i],bm.verts[int(7*nptscirc/2)+i]))
		else:
			bm.faces.new((bm.verts[nptscirc-1],bm.verts[4*nptscirc-1],bm.verts[0]))
		if i < int(nptscirc/2)-1:
			bm.faces.new((bm.verts[int(nptscirc/2)+1+i],bm.verts[int(7*nptscirc/2)+i],bm.verts[int(7*nptscirc/2)+1+i]))
		else:
			bm.faces.new((bm.verts[2*nptscirc],bm.verts[4*nptscirc-1],bm.verts[0]))
	#add the last few faces	
	bm.faces.new((bm.verts[int(nptscirc/2)],bm.verts[nptscirc],bm.verts[int(5*nptscirc/2)]))
	bm.faces.new((bm.verts[int(nptscirc/2)],bm.verts[nptscirc],bm.verts[int(7*nptscirc/2)]))

	bm.to_mesh(me) 
	me.update()
	return me

def subdiv(pme,subdivmeth):
		originalfaces=[]
		for f in pme.faces:
			originalfaces.append(f)
		originaledges=[]
		for e in pme.edges:
			originaledges.append(e)
		print('Subdividing. Intially', len(pme.verts),'verts,',len(pme.edges),'edges, and',len(pme.faces), 'faces.')
		counter=1
		if subdivmeth == 1:
			for f in pme.faces:
				if counter <= len(originalfaces):
					v1=tuple((f.verts[0].co + f.verts[1].co)/2)
					v2=tuple((f.verts[1].co + f.verts[2].co)/2)
					v3=tuple((f.verts[2].co + f.verts[0].co)/2)
					f0=f.verts[0].index
					f1=f.verts[1].index
					f2=f.verts[2].index
					#fc=f.cent
					fc=f.calc_center_median()
					if hasattr(pme.verts, "ensure_lookup_table"):
						pme.verts.ensure_lookup_table()
					vertsset=[tuple(pme.verts[i].co) for i in range(len(pme.verts))]
					if hasattr(pme.verts, "ensure_lookup_table"):
						pme.verts.ensure_lookup_table()
					if v1 not in vertsset and v2 not in vertsset and v3 not in vertsset:
						pme.verts.new(v1)
						pme.verts.new(v2)
						pme.verts.new(v3)
						pme.verts.new(fc)
						pme.verts.index_update()
						vind=len(pme.verts)
						pme.verts.ensure_lookup_table()
						pme.faces.new((pme.verts[vind-1],pme.verts[f0],pme.verts[vind-4]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-4],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-3],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[vind-3]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[vind-2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-2],pme.verts[f0]))
					elif v1 in vertsset and v2 not in vertsset and v3 not in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index					
						pme.verts.new(v2)
						pme.verts.new(v3)
						pme.verts.new(fc)
						pme.verts.index_update()
						vind=len(pme.verts)
						pme.verts.ensure_lookup_table()
						pme.faces.new((pme.verts[vind-1],pme.verts[f0],pme.verts[v1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v1],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-3],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[vind-3]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[vind-2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-2],pme.verts[f0]))
					elif v1 not in vertsset and v2 in vertsset and v3 not in vertsset:
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						pme.verts.new(v1)
						pme.verts.new(v3)
						pme.verts.new(fc)
						pme.verts.index_update()
						vind=len(pme.verts)
						pme.verts.ensure_lookup_table()
						pme.faces.new((pme.verts[vind-1],pme.verts[f0],pme.verts[vind-3]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-3],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v2],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[v2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[vind-2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-2],pme.verts[f0]))
					elif v1 not in vertsset and v2 not in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.new(v1)
						pme.verts.new(v2)
						pme.verts.new(fc)
						pme.verts.index_update()
						vind=len(pme.verts)
						pme.verts.ensure_lookup_table()
						pme.faces.new((pme.verts[vind-1],pme.verts[f0],pme.verts[vind-3]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-3],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-2],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[vind-2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[v3]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v3],pme.verts[f0]))
					elif v1 in vertsset and v2 in vertsset and v3 not in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						pme.verts.new(v3)
						pme.verts.new(fc)
						pme.verts.index_update()
						vind=len(pme.verts)
						pme.verts.ensure_lookup_table()
						pme.faces.new((pme.verts[vind-1],pme.verts[f0],pme.verts[v1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v1],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v2],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[v2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[vind-2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-2],pme.verts[f0]))
					elif v1 in vertsset and v2 not in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v1 == tuple(w.co):
								v1=w.index
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.new(v2)
						pme.verts.new(fc)
						pme.verts.index_update()
						vind=len(pme.verts)
						pme.verts.ensure_lookup_table()
						pme.faces.new((pme.verts[vind-1],pme.verts[f0],pme.verts[v1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v1],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-2],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[vind-2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[v3]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v3],pme.verts[f0]))
					elif v1 not in vertsset and v2 in vertsset and v3 in vertsset:
						for w in pme.verts:
							if v2 == tuple(w.co):
								v2=w.index
						for w in pme.verts:
							if v3 == tuple(w.co):
								v3=w.index
						pme.verts.new(v1)
						pme.verts.new(fc)
						pme.verts.index_update()
						vind=len(pme.verts)
						pme.verts.ensure_lookup_table()
						pme.faces.new((pme.verts[vind-1],pme.verts[f0],pme.verts[vind-2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[vind-2],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v2],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[v2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[v3]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v3],pme.verts[f0]))
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
						pme.verts.new(fc)
						pme.verts.index_update()
						vind=len(pme.verts)
						pme.verts.ensure_lookup_table()
						pme.faces.new((pme.verts[vind-1],pme.verts[f0],pme.verts[v1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v1],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v2],pme.verts[f1]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[v2]))
						pme.faces.new((pme.verts[vind-1],pme.verts[f2],pme.verts[v3]))
						pme.faces.new((pme.verts[vind-1],pme.verts[v3],pme.verts[f0]))
					counter += 1
					#Blender.Redraw()
		elif subdivmeth == 2:                #THIS NEEDS TO BE FIXED, ALONG THE LINES OF subdivmeth=1
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
		
		print(' After subdividing: ',len(pme.verts),'verts,',len(pme.edges),'edges,',len(pme.faces),'faces.')
		for f in originalfaces:
			pme.faces.remove(f)
		print(' Romoved original faces. Now',len(pme.verts),'verts,',len(pme.edges),'edges,',len(pme.faces),'faces.')

		for e in originaledges:
			pme.edges.remove(e)
		print(' Removed original edges. Now',len(pme.verts),'verts,',len(pme.edges),'edges,',len(pme.faces),'faces.')
		return(pme)
	# end of subdiv

def minim(bbm):	
	numverts=len(bbm.verts)
	S=numpy.zeros((numverts,numverts))
						  # Build S from Pinkall and Polthier section 3
	for e in bbm.edges:   # Best way to generate S.
		print("	edge", e) 
		print("	number of faces linked to edge:",len(e.link_faces))
		for f in e.link_faces:  # Pick a face.  Notice that this automatically handles both boundary and interior edges.				
			for l in f.loops:   # Pick a loop of that face.
				if l.vert not in e.verts:
					print("	angle",l.calc_angle(),"adding ",-.125*(1/numpy.tan(l.calc_angle())))
					S[e.verts[0].index,e.verts[1].index]= S[e.verts[0].index,e.verts[1].index]-.125*(1/numpy.tan(l.calc_angle()))
		S=.5*(S+numpy.transpose(S))
	for i in range(numverts):
		for e in bbm.verts[i].link_edges:
			S[bbm.verts[i].index,bbm.verts[i].index]=S[bbm.verts[i].index,bbm.verts[i].index]-S[e.verts[0].index,e.verts[1].index]
	print(" Size of S",S.shape)	
	print(S)
					
	#Find boundary and interior verts
	bdvertsindices=[]
	intvertsindices=[]
	for v in bbm.verts:
		if len(v.link_edges) > len(v.link_faces):
			bdvertsindices.append(v.index)
		else:
			intvertsindices.append(v.index)
	
	#Build the matrix S_int for the system that needs	to be inverted
	numintverts=len(intvertsindices)
	S_int=numpy.zeros((numintverts,numintverts))
	for i in range(numintverts):
		for j in range(numintverts):
			S_int[i,j] = S[intvertsindices[i],intvertsindices[j]]
	print(" Number of interior vertices",numintverts,"Size of S_int",S_int.shape)
	print(S_int)
	
	#Build the matrix S_bd to generate the right-side vector of the nonhomogeneous matrix equation to solve.
	numbdverts=len(bdvertsindices)
	S_bd=numpy.zeros((numintverts,numbdverts))
	for i in range(numintverts):
		for j in range(numbdverts):
			S_bd[i,j] = S[intvertsindices[i],bdvertsindices[j]]
	print(" Number of boundary vertices",numbdverts,"Size of S_bd",S_bd.shape)
	print(S_bd)
	
	#Compute the right-side vector
	bx=numpy.zeros((numintverts,1))
	by=numpy.zeros((numintverts,1))
	bz=numpy.zeros((numintverts,1))
	
	for i in range(numintverts):
		for j in range(numbdverts):
			bx[i]=bx[i]+S_bd[i,j]*bbm.verts[bdvertsindices[j]].co.x
			
	newx=numpy.linalg.solve(S_int,-bx)
	print("newx",sep=" ")
	for i in range(numintverts):
		print(newx[i],sep=" ")
		
	for i in range(numintverts):
		bbm.verts[intvertsindices[i]].co.x=newx[i]
	
	for i in range(numintverts):
		for j in range(numbdverts):
			by[i]=by[i]+S_bd[i,j]*bbm.verts[bdvertsindices[j]].co.y
			
	newy=numpy.linalg.solve(S_int,-by)
	print("newy",sep=" ")
	for i in range(numintverts):
		print(newy[i],sep=" ")

		
	for i in range(numintverts):
		bbm.verts[intvertsindices[i]].co.y=newy[i]
		
	for i in range(numintverts):
		for j in range(numbdverts):
			bz[i]=bz[i]+S_bd[i,j]*bbm.verts[bdvertsindices[j]].co.z
			
	newz=numpy.linalg.solve(S_int,-bz)
	print("newz",sep=" ")
	for i in range(numintverts):
		print(newz[i],sep=" ")
		
	for i in range(numintverts):
		bbm.verts[intvertsindices[i]].co.z=newz[i]
		
	return(bbm)	


def main(): 
	print('-- jminsurf --')	
	if bpy.app.version[0] < 2 or bpy.app.version[1] < 62:
		raise Exception("Only for Blender 2.62 and above.")
	if "Cube" in bpy.data.meshes: #This might be a way to delete the cube.  Might still show in the outliner's datablock
		cub=bpy.data.objects['Cube']
		cu=cub.data
		sce=bpy.data.scenes[0]
		sce.objects.unlink(cub)
		bpy.data.objects.remove(cub)
		bpy.data.meshes.remove(cu)

	me=bpy.data.meshes.new('myMesh')  #Create mesh
	ob=bpy.data.objects.new('myObject',me)  #Create object
	#ob.location=origin
	ob.show_name=True
	bpy.context.scene.objects.link(ob) #link to scene

	me=MakeInitFreiOtto1(me)   # Pass this a mesh - but the function creates and uses a bmesh, then returns a mesh version of that.
	
	#Get a BMesh representation
	bbm=bmesh.new()  #Create empty BMesh
	bbm.from_mesh(me) #fill it in from a Mesh
	
	simple=0
	verysimple=0
	if simple==1:
		#make something simple
		vertex1=bbm.verts.new((-6.0,0.0,6.0))
		vertex2=bbm.verts.new((0.0,6.0,0.0))
		vertex3=bbm.verts.new((6.0,0.0,6.0))
		vertex4=bbm.verts.new((0.0,-6.0,0.0))
		bbm.edges.new( (vertex1, vertex2) )
		bbm.edges.new( (vertex2, vertex3) )
		bbm.edges.new( (vertex3, vertex4) )
		bbm.edges.new( (vertex4, vertex1) )
		bbm.edges.new( (vertex1, vertex3) )
		bbm.faces.new((vertex1,vertex2,vertex3))
		bbm.faces.new((vertex3,vertex4,vertex1))
		bbm.to_mesh(me)
	if verysimple==1:
		#make something simple
		vertex1=bbm.verts.new((0.0,numpy.sqrt(3)/4,0.0))
		vertex2=bbm.verts.new((.5,-numpy.sqrt(3)/4,0.0))
		vertex3=bbm.verts.new((-.5,-numpy.sqrt(3)/4,0.0))
		#vertex4=bbm.verts.new((0.0,0.0,1.0))
		vertex4=bbm.verts.new((0.0,-1/(4*numpy.sqrt(3)),numpy.sqrt(3)/2))  # y coord here is -1.44337...  so the centriod.  I.e. regular tetrahedron
		bbm.edges.new( (vertex1, vertex2) )
		bbm.edges.new( (vertex2, vertex3) )
		bbm.edges.new( (vertex3, vertex1) )
		bbm.edges.new( (vertex1, vertex4) )
		bbm.edges.new( (vertex2, vertex4) )
		bbm.edges.new( (vertex3, vertex4) )
		bbm.faces.new((vertex1,vertex2,vertex4))
		bbm.faces.new((vertex2,vertex3,vertex4))
		bbm.faces.new((vertex3,vertex1,vertex4))
		bbm.to_mesh(me)

	#print(len(bbm.verts),'verts,',len(bbm.edges),'edges, and',len(bbm.faces), 'faces.')
	numtimestosubdivide = 2
	subdivmeth = 1  # triangulation method:  1  subdivides each face into 6 faces. 
					#                        2  subdivides each face into 4 faces. 
					
	for i in range(numtimestosubdivide):
		print("Subdividing",i+1," time")
		bbm=subdiv(bbm,subdivmeth)   # Pass this a bmesh. Returns a bmesh.
		bbm.to_mesh(me)
	
	for i in range(2):
		print("Minimizing",i+1," time:")	
		bbm=minim(bbm)
		bbm.to_mesh(me)
	
if __name__ == '__main__':      # This lets you import the script without running it 
	main() 
	
