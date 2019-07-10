
from brainvisa.processes import *

name = 'Surface sulci extraction from DPF'
userLevel = 0

signature = Signature(
    'white_mesh', ReadDiskItem('Hemisphere white mesh',
                               'aims mesh formats'),
    'dpf', ReadDiskItem('DPF texture', 'aims texture formats'),
    'sulci_texture', WriteDiskItem('Sulci white texture',
                                   'aims texture formats'),
    'dpf_threshold', Float(),
    'min_size', Integer(),
    'min_branch_size', Integer(),
)


def link_texture(self, proc, dummy):
    # Sulci white texture in in folds/<graph_version> in the hierarchy,
    # thus needs a graph_version. For an unknown reason we could not
    # place it somewhere else not depending on graph_version (a bug in
    # findValue()) so we needed this trick.
    if self.white_mesh is not None:
        print('link_texture:', self.signature['sulci_texture'].findValue(
            self.white_mesh, requiredAttributes={'graph_version': '3.1', 'labelled': 'No'}))
        return self.signature['sulci_texture'].findValue(
            self.white_mesh, requiredAttributes={'graph_version': '3.1', 'labelled': 'No'})

def initialization(self):
    self.linkParameters('dpf', 'white_mesh')
    self.linkParameters('sulci_texture', 'white_mesh', self.link_texture)
    self.dpf_threshold = -0.5
    self.min_size = 10
    self.min_branch_size = 10


def execution(self, context):
    from soma.aimsalgo import mesh_skeleton
    from soma import aims
    import numpy as np

    mesh = aims.read(self.white_mesh.fullPath())
    tex = aims.read(self.dpf.fullPath())
    texture = aims.TimeTexture((np.asarray(
        tex[0]) >= self.dpf_threshold).astype(np.int32))

    stex = mesh_skeleton.mesh_skeleton(mesh, texture, dist_tex=tex,
                                       min_cc_size=self.min_size,
                                       min_branch_size=self.min_branch_size)
    aims.write(stex, self.sulci_texture.fullPath())


