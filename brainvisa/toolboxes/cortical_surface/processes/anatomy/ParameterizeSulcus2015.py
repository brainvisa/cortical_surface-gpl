#  This software and supporting documentation are distributed by
#      Institut Federatif de Recherche 49
#      CEA/NeuroSpin, Batiment 145,
#      91191 Gif-sur-Yvette cedex
#      France
#
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the 
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license version 2 and that you accept its terms.

from brainvisa.processes import *
import shfjGlobals
import math
from numpy import *
from scipy import sparse
import scipy.sparse.linalg as alg
from scipy.sparse.linalg import spsolve
from soma import aims

name = 'Sulcus Parameterization 2015'

userLevel = 0

signature = Signature(
    'graph', ReadDiskItem( 'Cortical folds graph', 'Graph' ),
    'mri', ReadDiskItem( 'T1 MRI Bias Corrected', 'Aims readable volume formats' ),
    'label_attributes', Choice( 'label', 'name' ),
    'label_values', String(),
    'orientation', Choice( 'Top->Bottom', 'Front->Back' ),
    'sulcus_mesh', WriteDiskItem( 'Sulcus mesh', 'Aims mesh formats' ),
    'texture_param1', WriteDiskItem( 'Sulcus y coordinate texture', 'Texture' ),
    'coordinates_grid', WriteDiskItem( 'Sulcus coordinate grid mesh', 'Aims mesh formats' ),
    'depth_profile', WriteDiskItem( 'Sulcus depth profile', 'Text file' ),
    'dilation', Float(),
)

def initialization( self ):
     self.linkParameters( 'mri', 'graph' )
     self.linkParameters( 'sulcus_mesh', 'graph' )
     self.linkParameters( 'texture_param1', 'sulcus_mesh' )
     self.linkParameters( 'coordinates_grid', 'sulcus_mesh' )
     self.linkParameters( 'depth_profile', 'sulcus_mesh' )
     self.label_attributes = 'name'
     self.dilation = 1.0


##################################################################
# compute_mesh_weight - compute a weight matrix
#
#   W = compute_mesh_weight(vertex,face,type,options);
#
#   W is sparse weight matrix and W(i,j)=0 is vertex i and vertex j are not
#   connected in the mesh.
#
#   type is either 
#       'combinatorial': W(i,j)=1 is vertex i is conntected to vertex j.
#       'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
#           i and j.
#       'conformal': W(i,j) = cot(alpha_ij)+cot(beta_ij) where alpha_ij and
#           beta_ij are the adjacent angle to edge (i,j)
#
#   If options.normalize=1, the the rows of W are normalize to sum to 1.
#
####################################################################

def computeMeshWeights( mesh ):
     
     print '    Computing mesh weights'
     vert=array(mesh.vertex())
     poly=array(mesh.polygon())
     
     Nv=vert.shape[0]
     Np=poly.shape[0]
     
     W=sparse.lil_matrix((Nv,Nv))
     #W=zeros((Nv,Nv))
     # this old numpy array representation cannot handle big meshes in memory
     
     for i in range(3):
          i1=mod(i,3)
          i2=mod(i+1,3)
          i3=mod(i+2,3)
          print '    ', 3-i1
          pp=vert[poly[:,i2], :] - vert[poly[:,i1], :]
          qq=vert[poly[:,i3], :] - vert[poly[:,i1], :]
          np=apply_along_axis(linalg.norm, 1, pp)
          nq=apply_along_axis(linalg.norm, 1, qq)
          pp=pp / vstack( (np,vstack((np,np))) ).transpose()
          qq=qq / vstack( (nq,vstack((nq,nq))) ).transpose()
          ang = arccos(sum(pp*qq,1))
          
          for j in range(Np):
               ind1=poly[j,i1]
               ind2=poly[j,i2]
               ind3=poly[j,i3]
               W[ind2, ind3]=W[ind2, ind3]+1/tan(ang[j])
               W[ind3, ind2]=W[ind3, ind2]+1/tan(ang[j])
               
     print '    OK'
               
     return W
               
####################################################################
# 
# compute laplacian of a mesh
#
####################################################################
               
def computeMeshLaplacian( mesh ):
     print '    Computing Laplacian'
                    
     weights=computeMeshWeights( mesh )
     N=weights.shape[0]
     s=weights.sum(axis=1)
     dia=sparse.lil_matrix((N,N))
     dia.setdiag(s)
     L = dia - weights
                    
     print '    OK'
                    
     return L
     
####################################################################
# 
# compute comformal mapping of the mesh to a sphere
#
####################################################################


def sphereConformalMapping( mesh, lap, ps, pn, radius ):

     print '    Spherical mapping'
     L=sparse.lil_matrix(lap)
     #print 'Laplacian : ', L
     
     L[ps,:]=0
     L[ps,ps]=1
     L[pn,:]=0
     L[pn,pn]=1
     
     L=L.tocsr()

     vert=array(mesh.vertex())
     Nv=array(mesh.vertex()).shape[0]
     mesh.updateNormals()
     Nor=array(mesh.normal())
     Nor[ps,0]=0
     Nor[ps,1]=0
     Nor[ps,2]=-1
     Nor[pn,0]=0
     Nor[pn,1]=0
     Nor[pn,2]=1
     Rx=sparse.lil_matrix(Nor[:,0]).tocsr()
     Ry=sparse.lil_matrix(Nor[:,1]).tocsr()
     Rz=sparse.lil_matrix(Nor[:,2]).tocsr()

     print '    Solving linear system'

     x=spsolve(L,Rx)
     y=spsolve(L,Ry)
     z=spsolve(L,Rz)     
     
     print '    OK'
     vv=aims.vector_POINT3DF()
     for i in range(Nv):
          no=sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i])
          vv.append([x[i]/no, y[i]/no, z[i]/no])

     sphere=aims.AimsTimeSurface_3_VOID()
     sphere.vertex().assign(vv)
     sphere.polygon().assign(mesh.polygon())
     sphere.updateNormals()
     return(sphere)

####################################################################
# 
# compute iso-parameter lines on the mesh
#
####################################################################  

def meshIsoLine(mesh, tex, val):
     #print 'Looking for isoLine'
     points=array(mesh.vertex())
     values=array(tex[0])
     #print('isoLine: points:', points.shape)
     #print('isoLine: values:', values.shape)
     sign=zeros(values.size)
     sign[where(values < val)[0]]=10
     sign[where(values >= val)[0]]=20
     line=aims.AimsTimeSurface_2_VOID()
     isoV=aims.vector_POINT3DF()
     isoP=aims.vector_AimsVector_U32_2()
     
     triangles=array(mesh.polygon())
     for tr in triangles:
          count=sign[tr].sum()
          if (count==40):
               if (sign[tr[0]]==20):
                    v1=interpolateVertices(tr[0], tr[1], points, values, val)
                    v2=interpolateVertices(tr[0], tr[2], points, values, val)
               elif (sign[tr[1]]==20):
                    v1=interpolateVertices(tr[1], tr[0], points, values, val)
                    v2=interpolateVertices(tr[1], tr[2], points, values, val)
               elif (sign[tr[2]]==20):
                    v1=interpolateVertices(tr[2], tr[0], points, values, val)
                    v2=interpolateVertices(tr[2], tr[1], points, values, val)
               isoV, isoP = addSegment(v1, v2, isoV, isoP)       
          elif (count==50):
               if (sign[tr[0]]==10):
                    v1=interpolateVertices(tr[0], tr[1], points, values, val)
                    v2=interpolateVertices(tr[0], tr[2], points, values, val)
               elif (sign[tr[1]]==10):
                    v1=interpolateVertices(tr[1], tr[0], points, values, val)
                    v2=interpolateVertices(tr[1], tr[2], points, values, val)
               elif (sign[tr[2]]==10):
                    v1=interpolateVertices(tr[2], tr[0], points, values, val)
                    v2=interpolateVertices(tr[2], tr[1], points, values, val)
               isoV, isoP=addSegment(v1, v2, isoV, isoP)

     line.vertex().assign(isoV)
     line.polygon().assign(isoP)
     return line

# THIS ONE DOES NOT INTERPOLATE ON EDGES: 
# it sends back a list of vertex indices closest to the iso-line

def meshAlmostIsoLine(mesh, tex, val):
     #print 'Looking for isoLine'
     points=array(mesh.vertex())
     values=array(tex[0])
     #print('isoLine: points:', points.shape)
     #print('isoLine: values:', values.shape)
     sign=zeros(values.size)
     sign[where(values < val)[0]]=10
     sign[where(values >= val)[0]]=20
     
     si=set()
     
     triangles=array(mesh.polygon())
     for tr in triangles:
          count=sign[tr].sum()
          if (count==40):
               if (sign[tr[0]]==20):
                    i1=bestVertex(tr[0], tr[1], values, val)
                    i2=bestVertex(tr[0], tr[2], values, val)
               elif (sign[tr[1]]==20):
                    i1=bestVertex(tr[1], tr[0], values, val)
                    i2=bestVertex(tr[1], tr[2], values, val)
               elif (sign[tr[2]]==20):
                    i1=bestVertex(tr[2], tr[0], values, val)
                    i2=bestVertex(tr[2], tr[1], values, val)
               si.add(i1)
               si.add(i2)
          elif (count==50):
               if (sign[tr[0]]==10):
                    i1=bestVertex(tr[0], tr[1], values, val)
                    i2=bestVertex(tr[0], tr[2], values, val)
               elif (sign[tr[1]]==10):
                    i1=bestVertex(tr[1], tr[0], values, val)
                    i2=bestVertex(tr[1], tr[2], values, val)
               elif (sign[tr[2]]==10):
                    i1=bestVertex(tr[2], tr[0], values, val)
                    i2=bestVertex(tr[2], tr[1], values, val)
               si.add(i1)
               si.add(i2)
     return array(list(si))
               
def interpolateVertices(v1, v2, points, texture, valeur):
     t=abs(valeur - texture[v1])/abs(texture[v1] - texture[v2])
     new=(1-t)*points[v1] + t*points[v2]
     return new
     
def bestVertex(v1, v2, texture, valeur):
     t=abs(valeur - texture[v1])/abs(texture[v1] - texture[v2])
     if t<0.5:
          return v1
     else:
          return v2
     
def addSegment(v1, v2, vert, seg):
     nv=vert.size()
     a1=0
     a2=0
     
     for i in range(nv):
          v=vert[i]
          d1=v-v1
          d2=v-v2
          if (sqrt(d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]))<0.0001:
               i1=i
               a1=1
          if (sqrt(d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]))<0.0001:
               i2=i
               a2=1
 
     if (a1==0):
          vert.append(v1)
          i1=vert.size()-1
     if (a2==0):
          vert.append(v2)
          i2=vert.size()-1
     seg.append([i1, i2])
     return vert, seg 
     

     
####################################################################
# 
# find the vertex of a mesh closest to a point p
#
####################################################################

def closerVert(p, mesh):
     diff=p-mesh
     return(dot(diff,diff.transpose()).diagonal().argmin())

######################################################################

def execution( self, context ):
     sulcusIm=context.temporary( 'GIS image' )
     closedIm=context.temporary( 'GIS image' )
     hullIm=context.temporary(  'GIS image' )
#     dilHull=context.temporary(  'GIS image' )
     bottomIm=context.temporary(  'GIS image' )
#     dilBottom=context.temporary(  'GIS image' )
     simplesurf=context.temporary( 'GIS image' )
     dilatedIm=context.temporary(  'GIS image' )
     isoL=context.temporary( 'MESH mesh')
     meshNonDec=context.temporary( 'MESH mesh')
     
     distToPlan=context.temporary( 'Texture' )
     tempParam=context.temporary( 'Texture' )
     transform=''

     if (self.orientation=='Top->Bottom'):
          orient=0
     else :
          orient=1

     context.write('Extracting sulcus and buckets from graph')

     context.runProcess( 'Create Sulcus Label Volume',
                         graph = self.graph,
                         mri = self.mri,
                         sulci =  sulcusIm.fullPath(),
                         simple_surface = simplesurf.fullPath(),
                         compress= 'No',
                         bucket= 'Sulci',
                         label_attributes = self.label_attributes,
                         label_values = self.label_values,
                         node_edge_types='All' )
     context.runProcess( 'Create Sulcus Label Volume',
                         graph = self.graph,
                         mri = self.mri,
                         sulci =  sulcusIm.fullPath(),
                         simple_surface = simplesurf.fullPath(),
                         bottom = bottomIm.fullPath(),
                         compress= 'No',
                         bucket= 'Bottoms',
                         label_attributes = self.label_attributes,
                         label_values = self.label_values,
                         node_edge_types='All' )
     context.runProcess( 'Create Sulcus Label Volume',
                         graph = self.graph,
                         mri = self.mri,
                         sulci =  sulcusIm.fullPath(),
                         simple_surface = simplesurf.fullPath(),
                         hull_junction = hullIm.fullPath(),
                         compress= 'No',
                         bucket= 'Junctions with brain hull',
                         label_attributes = self.label_attributes,
                         label_values = self.label_values,
                         node_edge_types='All' )

     dilating = [ 'AimsDilation',
                 '-i', sulcusIm.fullPath(),
                 '-o', dilatedIm.fullPath(),
                 '-e', self.dilation ]

     apply( context.system, dilating )

     closing = [ 'AimsClosing',
                 '-i', dilatedIm.fullPath(),
                 '-o', closedIm.fullPath(),
                 '-r', 2 ]

     apply( context.system, closing )

     context.write('Remeshing sulcus')

     meshing = [ 'AimsMesh',
                 '-i', closedIm.fullPath(),
                 '-o', self.sulcus_mesh.fullPath(),
                 '-l', '32767',
                 '--smooth',
                 '--smoothIt', '20' ]
     apply( context.system, meshing )

     test=self.sulcus_mesh.fullName()
     sulcusMname=test + '_32767_0.mesh'

     shelltools.mv(sulcusMname, self.sulcus_mesh.fullPath())
     shelltools.mv(sulcusMname + '.minf', self.sulcus_mesh.fullPath() + '.minf')

     #############################################################
     #
     #     This is a new parameterization along the y axis only
     #     using the fiedler vector
     #
     #############################################################
     context.write('Starting parameterization')
     
     re=aims.Reader()  

     mesh=re.read(self.sulcus_mesh.fullPath())
     vert=array(mesh.vertex())
     N= mesh.vertex().size()
     parameter=zeros(N)
     
     print 'Computing Laplacian'
     L=computeMeshLaplacian(mesh)
     Lap=0.5*(L+L.transpose())
     print 'Computing fiedler vector'
     w,v=alg.eigsh(Lap, 4, which='SM')
 
     texOut=aims.TimeTexture_FLOAT(1,N)
     isoParam=aims.TimeTexture_FLOAT(1,N)
     fiedler=v[:,1]
     
     # gestion de l'orientation
     
     imin=fiedler.argmin()
     imax=fiedler.argmax()
     
     if (orient==1):
          cmin=vert[imin][1]
          cmax=vert[imax][1]
     elif (orient==0):
          cmin=vert[imin][2]
          cmax=vert[imax][2]
     
     if (cmin<cmax):
          top=imin
          bot=imax
     else:
          top=imax
          bot=imin
          
     ftop=fiedler[top]
     fbot=fiedler[bot]
          
     # rescaling de la texture antre 0 et 100
     
     a=100.0/float(fbot - ftop)
     b=-100.0*ftop/float(fbot - ftop)
     
     for i in range(N):
          texOut[0][i]=(fiedler[i]*a) + b
 
     #print 'Writing texture', self.texture_param1.fullName()
     ws=aims.Writer()
 
     ws.write( texOut, tempParam.fullPath() )
     
     context.write('Recomputing isometric parameterization')
     
     ########################################################
     # Computing medial axis et reparametrisation isometrique
     # (par morceaux) a partir de sa longueur
     ########################################################
     axis=aims.AimsSurfaceTriangle()
     gen=aims.SurfaceGenerator()
     vx=0
     vy=0
     vz=0

     points=zeros((101,3))
     dist=zeros(101)
     newCoord=zeros(101)
     a=zeros(101)
     b=zeros(101)
     points[0]=vert[top]
     points[100]=vert[bot]

     for i in range(1,100):
          cont= meshIsoLine(mesh, texOut, i)
          vertC=array(cont.vertex())
          bary=vertC.mean(axis=0)
          points[i]=bary
          vec=points[i] - points[i-1]
          dist[i]=dist[i-1]+sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] )
     vec=points[100] - points[99]
     dist[100]=dist[99]+sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] )
     
     #print 'Dist:', dist
     newCoord[0]=0.0
     newCoord[100]=100.0
     for i in range(1, 100):
          newCoord[i]=(dist[i]/dist[100]) *100.0
     #print 'newCoord:', newCoord

     for i in range(0, 100):
          a[i]=(newCoord[i+1]-newCoord[i])
          b[i]=(newCoord[i]*(i+1)-newCoord[i+1]*i)
     for i in range(N):
          ind=math.floor(texOut[0][i])
          isoParam[0][i]=(texOut[0][i]*a[ind]) + b[ind]
          parameter[i]=(texOut[0][i]*a[ind]) + b[ind]
     ws.write( isoParam, self.texture_param1.fullPath() )
     
     context.write('Isometric reparameterization done')
     
     ########################################################
     context.write('Computing coordinate grid')

     i=0;
     iso = [ 'AimsMeshIsoLine',
              '-i', self.sulcus_mesh.fullPath(),
              '-t', self.texture_param1.fullPath(),
              '-o', isoL.fullPath(),
              '-v', i ]
     conc = [ 'AimsZCat',
              '-i', isoL.fullPath(),
              '-o', self.coordinates_grid.fullPath() ]
     apply(context.system, iso)
     apply(context.system, conc)
     
     i=i+5
     while (i<=100):
          iso = [ 'AimsMeshIsoLine',
                   '-i', self.sulcus_mesh.fullPath(),
                   '-t', self.texture_param1.fullPath(),
                   '-o', isoL.fullPath(),
                   '-v', i ]
          conc = [ 'AimsZCat',
                   '-i', isoL.fullPath(), self.coordinates_grid.fullPath(),
                   '-o', self.coordinates_grid.fullPath() ]
          apply(context.system, iso)
          apply(context.system, conc)
          i=i+5
          
     read=aims.Reader()

     # here we compute the signed distance to the inertia plane
     # inertia plane is computed with a PCA of the vertices.
     # extremities (defined with 'offset') are removed.
     
     offset=10
     vert2=vert[where((parameter >= offset) & (parameter <= (100-offset)))].copy()
     bary=mean(vert2, axis=0)
     vert2=vert2-bary
     vert=vert-bary
     tvert2=vert2.transpose()
     coord=dot(tvert2,vert2)
     
 
     val,vect=linalg.eig(coord)
     i=argmin(val)
     k=argmax(val)
     for t in range(3):
          if (t!=i) and (t!=k):
               j=t
        
     u1=vect[:,i]
     u2=vect[:,j]
     u3=vect[:,k]

     texturex=aims.TimeTexture_FLOAT() 
     nn= vert.shape[0]
     for i in range(nn):
          texturex.push_back(dot(vert[i], -u1))

     dist=array(texturex[0])    
     context.write('Computing morphological curves (depth and profile)')
     
     
     # Computing curves for depth and profile
     
     morpho=zeros((101, 3))
     vert=array(mesh.vertex())
     for i in range(101):
          line=meshIsoLine(mesh, isoParam, i)
          vline=meshAlmostIsoLine(mesh, isoParam, i)
          depth=0
          prof=0
          prof2=0
          vl=array(line.vertex())
          pl=array(line.polygon())
          for p in pl:
               i1=p[0]
               i2=p[1]
               depth+=sqrt(dot((vl[i2]-vl[i1]), (vl[i2]-vl[i1])))
          for v in vline:
               prof+=dist[v]
          depth=depth/2.0
          if (depth==0):
               prof=0
               prof2=0
          else:
               prof2=prof/vl.size
               prof=prof/(2.0*depth)    
          morpho[i]=array([i, depth, prof])
     savetxt(self.depth_profile.fullPath(), morpho, delimiter='\t')

     context.write('Finished')

