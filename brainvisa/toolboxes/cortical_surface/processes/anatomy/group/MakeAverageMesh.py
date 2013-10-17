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
from soma import aims
import numpy as np

    
name = 'Make Average Mesh From a Set of Subjects'

userLevel = 0

signature = Signature(
                      
    'remeshed_meshes',ListOf( ReadDiskItem( 'Remeshed mesh', shfjGlobals.aimsMeshFormats) ),
    'side', Choice('left', 'right'),
    'CACP_tranfos',ListOf( ReadDiskItem( 'Transform Raw T1 MRI to Talairach-AC/PC-Anatomist','Transformation matrix' ) ), 
    'average_mesh',WriteDiskItem( 'Mesh', shfjGlobals.aimsMeshFormats),
)

def initialization( self ):
    def linkSide( proc, dummy ):
        if proc.remeshed_meshes is not None \
                and len( proc.remeshed_meshes ) != 0:
            return proc.remeshed_meshes[0].get( 'side' )
    self.linkParameters( 'side', 'remeshed_meshes', linkSide )
    self.linkParameters( 'CACP_tranfos','remeshed_meshes')
    
def execution( self, context ):

    re = aims.Reader()
    nb_mesh = len(self.remeshed_meshes)
    tmp_mesh = context.temporary(  'Mesh mesh' )

    for ind_mesh,r_mesh in enumerate(self.remeshed_meshes):
        context.write('working on mesh nb: ',ind_mesh+1)
        context.system('AimsMeshTransform', '-i',r_mesh.fullPath(),'-o',tmp_mesh.fullPath(),'-t',self.CACP_tranfos[ind_mesh].fullPath())
        mesh = re.read(tmp_mesh.fullPath())
        context.write('nb of vertices in the mesh '+str(len(mesh.vertex())))
        if ind_mesh==0:
            avg_verts = np.array(mesh.vertex())
        else:
            if len(mesh.vertex()) == avg_verts.shape[0]:
                avg_verts += np.array(mesh.vertex())
            else:
                raise Exception('mesh # '+str(ind_mesh)+' do not have the good number of vertices')
    avg_verts = avg_verts / nb_mesh
    avg_mesh = aims.AimsTimeSurface_3()
    avg_mesh.vertex().assign([ aims.Point3df(x) for x in avg_verts ])
    if self.side == 'right':
        poly = np.array(mesh.polygon())
        poly_tmp = poly.copy()
        poly_tmp[:,0] = poly[:,1]
        poly_tmp[:,1] = poly[:,0]
        pp = aims.vector_AimsVector_U32_3()
        for i in poly_tmp:
            pp.append(i)
        avg_mesh.polygon().assign(pp)
    else:
        avg_mesh.polygon().assign(mesh.polygon())
    # ensure there is no Nan in the vertex of outmesh, should never be the case because it was checked in remeshFromSphere
    verts_Nan = np.where(np.isnan(avg_verts))[0]
    if len(verts_Nan)>0:
        neigh = aims.SurfaceManip.surfaceNeighbours(avg_mesh)
        for vert_ind in verts_Nan:
            neig_vert = avg_verts[neigh[vert_ind].list(), :]
            avg_verts[vert_ind, :] = np.mean(neig_vert, axis = 0)
            avg_mesh.vertex().assign([ aims.Point3df(x) for x in avg_verts ])
           
    avg_mesh.updateNormals()
    
    ws = aims.Writer()
    ws.write(avg_mesh, self.average_mesh.fullPath())
