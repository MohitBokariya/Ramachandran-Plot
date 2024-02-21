#To Find phi psi and omega from pdb file 

#MOHIT BOKARIYA
#18/may/2019

import sys 
import os
from math import*
import numpy as np
import matplotlib.pyplot as plt

pdb=sys.argv[1]  #for line command argument

dirc=os.getcwd()   #directory of program

E=[]
#x,y,z coordinates for Nitrogen list
xn=[]
yn=[]
zn=[]

#x,y,z coordinates for Alpha carbon list
xca=[]
yca=[]
zca=[]

#x,y,z coordinates for Carbonyl Carbon list
xc=[]
yc=[]
zc=[]

#x,y,z coordinates for Oxygen list
xo=[]
yo=[]
zo=[]

#phi, psi, omega angles list
P=[]
PS=[]
W=[]

#coordinates list
N2=[] # Nitrogen 
CA2=[] #Alpha carbon 
C2=[]  # Carbonyl Carbon
O2=[]	#Oxygen

#Residue number list 

RON=[]
ROAC=[]
ROC=[]
ROO=[]

with open(os.path.expanduser(r"{:s}/{:s}".format(dirc,pdb)),'r') as f:

	data=f.readlines()
	for line in data:
		list1=line.split()
		ID=list1[0]
		if ID == 'ATOM':
			ID1=list1[2]  # for Nitrogen
			if ID1=='N':
					
			
					
					position=list1[6:9]

					p1=[position[0],position[1],position[2]]

					a1=float(p1[0])
					b1=float(p1[1])
					c1=float(p1[2])
					#x,y,z coordinates of Nitrogen					
					xn.append(a1) 
					yn.append(b1)
					zn.append(c1)
					RON.append(float(list1[1]))
					
			#for Alpha carbon
			
			ID2=list1[2]
			if ID2=='CA':
				
					position=list1[6:9]

					p2=[position[0],position[1],position[2]]

					a2=float(p2[0])
					b2=float(p2[1])
					c2=float(p2[2])
					#x,y,z coordinates						
					xca.append(a2)
					yca.append(b2)
					zca.append(c2)
					ROAC.append(float(list1[1]))

#for carbonyl carbon 

			ID3=list1[2]
			if ID3=='C':
					
					position=list1[6:9]

					p3=[position[0],position[1],position[2]]

					a3=float(p3[0])
					b3=float(p3[1])
					c3=float(p3[2])
					#x,y,z coordinates
					xc.append(a3)
					yc.append(b3)
					zc.append(c3)
					ROC.append(float(list1[1]))

#for  Oxygen

			ID4=list1[2]
			if ID4=='O':
					
					position=list1[6:9]

					p4=[position[0],position[1],position[2]]

					a4=float(p4[0])
					b4=float(p4[1])
					c4=float(p4[2])
					#x,y,z coordinates					
					xo.append(a4)
					yo.append(b4)
					zo.append(c4)
					ROO.append(float(list1[1]))					
					
for i in range(len(xn)):

	N1=(xn[i],yn[i],zn[i])
	N2.append(N1)


	CA1=(xca[i],yca[i],zca[i])
	CA2.append(CA1)

	C1=(xc[i],yc[i],zc[i])

	C2.append(C1)

	O1=(xo[i],yo[i],zo[i])
	O2.append(O1)

N=np.array(N2)   # All nitrogen  coordinates in an array
CA=np.array(CA2)
C=np.array(C2)
O=np.array(O2)
		


#phi(C1,N2,CA2,C2)

for i in range(len(xn)-1):

	v01=C[i]-N[i+1]
	v12=N[i+1]-CA[i+1]
	v32=C[i+1]-CA[i+1]

	v1=np.cross(v12,v01)
	v2=np.cross(v12,v32)
	m0=(np.sum(v1**2))**(1/2)
	m3=(np.sum(v2**2))**(1/2)
	p=np.dot(v1,v2)/(m0*m3)	
	phi=acos(p)
	phi=phi*(180/pi)	
	#print(phi)
	P.append(phi)


#Psi(N2,CA2,C2,N3)
for i in range(len(xn)-2):
	
	ve01=N[i+1]-CA[i+1]
	ve12=CA[i+1]-C[i+1]
	ve32=N[i+2]-C[i+1]

	v3=np.cross(ve12,ve01)
	v4=np.cross(ve12,ve32)
	m3=(np.sum(v3**2))**(1/2)
	m4=(np.sum(v4**2))**(1/2)
	ps=np.dot(v3,v4)/(m3*m4)	
	psi=acos(ps)
	psi=psi*(180/pi)	
	#print(psi)
	PS.append(psi)



#omega(CA2,C2,N3,CA3)

	vec01=CA[i+1]-C[i+1]
	vec12=C[i+1]-N[i+2]
	vec32=CA[i+2]-N[i+2]

	v5=np.cross(vec12,vec01)
	v6=np.cross(vec12,vec32)
	m5=(np.sum(v5**2))**(1/2)
	m6=(np.sum(v6**2))**(1/2)
	o=np.dot(v5,v6)/(m5*m6)	
	omega=acos(o)
	omega=omega*(180/pi)	
	#print(omega)
	W.append(omega)
	
	print("Angle phi,psi is :",'{:.4f}'.format(P[i]),';','{:.4f}'.format(PS[i]),';','{:.4f}'.format(W[i]),'in Degree')

plt.title("OMEGA vs Residue number")
plt.xlabel('Residue numbers')
plt.ylabel("OMEGA")
		
#plt.plot(P,PS,'r*')
#plt.show()
#for e in range(len(P)):
	#E.append(e)
#plt.plot(E,W,'b*')
#plt.show()

