
from __future__ import print_function

from brainvisa.processes import *

name = 'Surface sulci extraction from DPF and labels projection'
userLevel = 0

signature = Signature(
    'sulci_graph', ReadDiskItem('Labelled cortical folds graph',
                                'Graph and data'),
    'white_mesh', ReadDiskItem('Hemisphere white mesh',
                               'aims mesh formats'),
    'dpf', WriteDiskItem('DPF texture', 'aims texture formats'),
    'sulci_texture', WriteDiskItem('Sulci white texture',
                                   'aims texture formats'),
    'dpf_threshold', Float(),
    'min_size', Integer(),
    'min_branch_size', Integer(),
)



def initialization(self):

    def link_mask_textures(voronoi, sulci):
        return [voronoi, sulci]

    enode = SerialExecutionNode(self.name, parameterized=self)
    enode.addChild('dpf',
                   ProcessExecutionNode('DepthPotentialFunction', optional=1))
    enode.addChild('surface_sulci',
                   ProcessExecutionNode('surface_sulci_dpf', optional=1))
    enode.addChild('sulci_voronoi',
                   ProcessExecutionNode('sulciGraphVoronoi', optional=1))
    enode.addChild('masking',
                   ProcessExecutionNode('texture_operation', optional=1))

    enode.surface_sulci.removeLink('dpf', 'white_mesh')
    #enode.surface_sulci.removeLink('sulci_texture', 'white_mesh')
    #enode.sulci_voronoi.removeLink('white_mesh', 'graph')

    enode.addDoubleLink('sulci_graph', 'sulci_voronoi.graph')
    enode.addDoubleLink('white_mesh', 'dpf.input_mesh')
    enode.addDoubleLink('dpf', 'dpf.DPF_texture')
    enode.addDoubleLink('white_mesh', 'surface_sulci.white_mesh')
    enode.addDoubleLink('dpf', 'surface_sulci.dpf')
    #enode.addDoubleLink('sulci_texture', 'surface_sulci.sulci_texture')
    enode.addDoubleLink('sulci_graph', 'sulci_voronoi.graph')
    enode.addDoubleLink('white_mesh', 'sulci_voronoi.white_mesh')
    enode.addLink('masking.textures',
                  ['sulci_voronoi.sulci_voronoi',
                   'surface_sulci.sulci_texture'], link_mask_textures)
    enode.addDoubleLink('sulci_texture', 'masking.output_texture')
    self.linkParameters('sulci_texture', 'sulci_graph')
    enode.addDoubleLink('dpf_threshold', 'surface_sulci.dpf_threshold')
    enode.addDoubleLink('min_size', 'surface_sulci.min_size')
    enode.addDoubleLink('min_branch_size', 'surface_sulci.min_branch_size')

    enode.masking.formula = 'NT0[NT1==0] = 0'
    enode.masking.output_is_first_input = True
    enode.operation = 'Formula'
    self.dpf_threshold = enode.surface_sulci.dpf_threshold
    self.min_size= enode.surface_sulci.min_size
    self.min_branch_size = enode.surface_sulci.min_branch_size

    self.setExecutionNode(enode)

