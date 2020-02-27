# -*- coding: utf-8 -*-
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

from __future__ import absolute_import
from brainvisa.processes import *
import brainvisa.tools.aimsGlobals as shfjGlobals
from brainvisa.data.neuroHierarchy import databases
from brainvisa import registration

name = 'Gyral Parcellation'

userLevel = 2

signature = Signature(
    'Side', Choice("Both","Left","Right"),
    'Projection',Choice("Yes","No"),
    'Parcellation',Choice("White Matter Surface","Cortical Ribbon"),
    'Lgraph', ReadDiskItem( 'Labelled Cortical folds graph', 'Graph',
                            requiredAttributes={ 'side': 'left' } ),
    'Rgraph', ReadDiskItem( 'Labelled Cortical folds graph', 'Graph',
                            requiredAttributes={ 'side': 'right' } ),
    'mri_corrected', ReadDiskItem( 'T1 MRI Bias Corrected', 'aims readable volume formats' ),
    'sulcus_identification',Choice('label','name'),
    'translation',ReadDiskItem('Label Translation','Label Translation' ),
    'gyri_model',ReadDiskItem('Gyri Model','Gyri Model' ),
    'sulci_label_to_sulci_name',WriteDiskItem( 'Sulci To White Texture Translation', 'Text File'),
    'left_gyri_label_to_gyri_name',WriteDiskItem( 'Gyri To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'left' }),
    'right_gyri_label_to_gyri_name',WriteDiskItem( 'Gyri To White Texture Translation', 'Text File',requiredAttributes={ 'side': 'right' }),
    'left_white_mesh',ReadDiskItem( 'Left Hemisphere White Mesh' ,
                                    shfjGlobals.aimsMeshFormats),
    'right_white_mesh',ReadDiskItem( 'Right Hemisphere White Mesh' ,
                                     shfjGlobals.aimsMeshFormats),
    'left_hemi_mesh', ReadDiskItem( 'Left Hemisphere Mesh', 'aims mesh formats' ),
    'right_hemi_mesh', ReadDiskItem( 'Right Hemisphere Mesh', 'aims mesh formats' ),
    'left_grey_white', ReadDiskItem( 'Left Grey White Mask', 'aims readable volume formats' ),
    'right_grey_white', ReadDiskItem( 'Right Grey White Mask', 'aims readable volume formats' ),
    'left_white_sulci',WriteDiskItem( 'Sulci White Texture' ,'aims texture formats',
                                      requiredAttributes={ 'side': 'left' } ),
    'right_white_sulci',WriteDiskItem( 'Sulci White Texture',
                                       'aims texture formats',
                                       requiredAttributes= \
                                       { 'side': 'right' } ),
    'left_white_gyri',WriteDiskItem( 'Gyri White Texture',
                                     'aims texture formats',
                                     requiredAttributes={ 'side': 'left' } ),
    'right_white_gyri',WriteDiskItem( 'Gyri White Texture',
                                      'aims texture formats',
                                      requiredAttributes={ 'side': 'right' } ),
    'left_gyri_graph',WriteDiskItem( 'Gyri Graph' ,'Graph',
                                           requiredAttributes= \
                                           { 'side': 'left' } ),
    'right_gyri_graph',WriteDiskItem( 'Gyri Graph' ,'Graph',
                                            requiredAttributes= \
                                            { 'side': 'right' } ),
    'left_white_gyri_volume',WriteDiskItem( 'Gyri White Volume',
                                            'aims writable volume formats',
                                            requiredAttributes= \
                                            { 'side': 'left' } ),
    'right_white_gyri_volume',WriteDiskItem( 'Gyri White Volume',
                                             'aims writable volume formats',
                                             requiredAttributes= \
                                             { 'side': 'right' } ),
)

def initialization( self ):
    def linkLabelAtt( self, dummy ):
        if self.Lgraph is not None:
            m = self.Lgraph.get( 'manually_labelled' )
            if m and m == 'Yes':
                return 'name'
        return 'label'
    self.linkParameters( 'left_white_mesh', 'Lgraph' )
    self.linkParameters( 'left_hemi_mesh', 'Lgraph' )
    self.linkParameters( 'Rgraph', 'Lgraph' )
    self.linkParameters( 'mri_corrected', 'Lgraph' )
    self.linkParameters( 'right_white_mesh', 'Rgraph' )
    self.linkParameters( 'right_hemi_mesh', 'Rgraph' )
    self.linkParameters( 'left_white_sulci', 'Lgraph' )
    self.linkParameters( 'right_white_sulci', 'Rgraph' )
    self.linkParameters( 'left_white_gyri', 'Lgraph' )
    self.linkParameters( 'right_white_gyri', 'Rgraph' )
    self.linkParameters( 'right_gyri_graph', 'Rgraph' )
    self.linkParameters( 'left_gyri_graph', 'Lgraph' )
    self.linkParameters( 'left_white_gyri_volume', 'Lgraph' )
    self.linkParameters( 'right_white_gyri_volume', 'Rgraph' )
    self.linkParameters( 'left_grey_white', 'Lgraph' )
    self.linkParameters( 'right_grey_white', 'Rgraph' )
    self.linkParameters( 'sulci_label_to_sulci_name', 'Lgraph' )
    self.linkParameters( 'left_gyri_label_to_gyri_name', 'Lgraph' )
    self.linkParameters( 'right_gyri_label_to_gyri_name', 'Rgraph' )
    self.sulcus_identification = 'label'
    self.linkParameters( 'sulcus_identification', 'Lgraph', linkLabelAtt )
    self.setOptional('right_grey_white','left_grey_white',
                    'left_white_gyri_volume','right_white_gyri_volume',
                    'Rgraph', 'Lgraph',
                    'left_white_mesh','right_white_mesh',
                    'left_white_gyri','right_white_gyri',
                    'left_hemi_mesh','right_hemi_mesh',
                    'left_gyri_graph','right_gyri_graph',
                    'left_white_sulci','right_white_sulci' )
    self.translation = ReadDiskItem('Label Translation','Label Translation' ).findValue( { 'filename_variable' : 'gyri' } )
    try:
      self.gyri_model = databases.getDiskItemFromUuid( '172c4168-a9d3-dc41-464c-1226ad07c19c' )
    except: pass
    self.Projection = 'Yes'

def execution( self, context ): 

     call_list = ['siParcellation',
                  '-m', self.gyri_model.fullPath() ]

     if self.Parcellation == 'Cortical Ribbon':
          call_list += ['--3D']

     if self.Side in ('Left','Both'):
          if ( self.Projection == 'Yes' ):
               context.runProcess( 'surface_sulci_projection_nearest',
                                   Side = 'Left',
                                   Lgraph = self.Lgraph,
                                   Rgraph = self.Rgraph,
                                   left_white_mesh = self.left_white_mesh,
                                   right_white_mesh = self.right_white_mesh,
                                   left_white_sulci = self.left_white_sulci,
                                   right_white_sulci = self.right_white_sulci,
                                   translation = self.translation,
                                   left_sulci_label_to_sulci_name = self.sulci_label_to_sulci_name,
                                   right_sulci_label_to_sulci_name = self.sulci_label_to_sulci_name,
                                   gyri_model = self.gyri_model,
                                   mri_corrected = self.mri_corrected,
                                   sulcus_identification = self.sulcus_identification )
          io = ['-i', self.left_white_mesh.fullPath() ,
                '-s' ,self.left_white_sulci.fullPath() ,
                '-o' ,self.left_white_gyri.fullPath() ,
                '-g',self.left_gyri_graph.fullPath(),
                '-p',self.left_white_gyri_volume.fullPath(),
                '-b', self.left_hemi_mesh.fullPath(),
                '--sulcitraduction', self.sulci_label_to_sulci_name.fullPath() ,
                '--gyritraduction', self.left_gyri_label_to_gyri_name.fullPath() ]

          if self.Parcellation == 'Cortical Ribbon':
              if self.left_grey_white is None \
                or self.left_white_gyri_volume is None:
                raise ValueError( _t_( 'In Cortical Ribbon mode, ' \
                'left_grey_white and left_white_gyri_volume parameters ' \
                'are mandatory' ) )
              io += ['-V',self.left_grey_white.fullPath()]

          context.system(*call_list+io)
          tm = registration.getTransformationManager()
          tm.copyReferential( self.Lgraph, self.left_gyri_graph )

     if self.Side in ('Right','Both'):
          if ( self.Projection == 'Yes' ):
               context.runProcess( 'surface_sulci_projection_nearest',
                                   Side = 'Right',
                                   Lgraph = self.Lgraph,
                                   Rgraph = self.Rgraph,
                                   left_white_mesh = self.left_white_mesh,
                                   right_white_mesh = self.right_white_mesh,
                                   left_white_sulci = self.left_white_sulci,
                                   right_white_sulci = self.right_white_sulci,
                                   translation = self.translation,
                                   left_sulci_label_to_sulci_name = self.sulci_label_to_sulci_name,
                                   right_sulci_label_to_sulci_name = self.sulci_label_to_sulci_name,
                                   gyri_model = self.gyri_model,
                                   mri_corrected = self.mri_corrected,
                                   sulcus_identification = self.sulcus_identification)
          io = ['-i', self.right_white_mesh.fullPath() ,
                "-s" ,self.right_white_sulci.fullPath() ,
                '-o' ,self.right_white_gyri.fullPath() ,
                '-g',self.right_gyri_graph.fullPath(),
                '-p',self.right_white_gyri_volume.fullPath(),
                '-b', self.right_hemi_mesh.fullPath(),
                '--sulcitraduction', self.sulci_label_to_sulci_name.fullPath() ,
                '--gyritraduction', self.right_gyri_label_to_gyri_name.fullPath() ]

          if self.Parcellation == 'Cortical Ribbon':
              if self.left_grey_white is None \
                or self.left_white_gyri_volume is None:
                raise ValueError( _t_( 'In Cortical Ribbon mode, ' \
                'right_grey_white and right_white_gyri_volume parameters ' \
                'are mandatory' ) )
              io += ['-V',self.right_grey_white.fullPath()]

          context.system(*call_list+io)
          tm = registration.getTransformationManager()
          tm.copyReferential( self.Rgraph, self.right_gyri_graph )
