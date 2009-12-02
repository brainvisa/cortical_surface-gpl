# Copyright CEA and IFR 49 (2000-2005)
#
#  This software and supporting documentation were developed by
#      CEA/DSV/SHFJ and IFR 49
#      4 place du General Leclerc
#      91401 Orsay cedex
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


from neuroProcesses import *
import shfjGlobals
import csv
import numpy as np
import os

name = 'Plot Curvature Connex Components'

userLevel = 3

signature = Signature(
    'white', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh'),
    'scale_space', ReadDiskItem( 'Scale Space White Curvature Texture', 'Texture'),
    'csv_file', WriteDiskItem('Text File', 'Text File'),
    'gyri', ReadDiskItem( 'Hemisphere gyri parcellation texture', 'Texture')
    )

def count_label_max(texture):
  label_max = 0
  for i in xrange(len(texture)):
    if (texture[i] > label_max):
      label_max = texture[i]
  return label_max

def initialization( self ):
  self.linkParameters( 'gyri', 'white' )
  self.linkParameters( 'scale_space', 'white' )

def execution( self, context ):
    sujet = self.white.get( 'subject')
    context.write(sujet)

    # on lit la texture d'espace echelle et on cree un fichier temporaire
    scale_space = aims.read(str(self.scale_space))
    curvature_tmp_texfile = context.temporary("Texture")
    context.write( "curvature_tmp_texfile : " + str (curvature_tmp_texfile) )
    nb_scales = scale_space.size()
    context.write( str(nb_scales) + " scales" )

    #on calcule la texture d'intersection de gyri
    intersection_texfile = context.temporary("Texture")
    context.write("intersection_texfile : " + str(intersection_texfile))
    command = [ 'AimsTextureLabelIntersection', '-m', self.white, '-t', self.gyri, '-o', str(intersection_texfile)]
    apply(context.system, command)
    intersection_texture = aims.read(str(intersection_texfile))


    data = []
    data.append(["sujet", "scale_number", "threshold", "cc", "nb_cc_inters", "nb_cc_noninters"])
    # pour chaque niveau d'echelle
    for scale_number in xrange(nb_scales) :
      context.write( "scale :" + str(scale_number) )
      threshold = 0.0
      nb_cc = 1

      all_thresholds_cc_tex = []
      iteration = 0
      
      while (nb_cc != 0):
        
        # on recupere la texture du niveau d'echelle courant
        scale_texture = aims.TimeTexture_FLOAT(1,len(scale_space[scale_number]))
        scale_texture[0] = scale_space[scale_number]
        wtex = aims.Writer()
        wtex.write(scale_texture, str(curvature_tmp_texfile) )


        # on calcule une carte seuillee
        command = [ 'AimsTextureThreshold', '-i', str(curvature_tmp_texfile), '-o', str(curvature_tmp_texfile), '-m', "le", '-t', str(threshold) ]
        threshold -= 0.05
        apply(context.system, command)

        # on convertit en FLOAT texture
        command = [ 'AimsFileConvert', '-i', str(curvature_tmp_texfile), '-o', str(curvature_tmp_texfile), '-t', 'FLOAT' ]
        apply(context.system, command)

        # on calcule une premiere fois les composantes connexes
        command = [ 'AimsMeshConnectedComponent', '-t', str(curvature_tmp_texfile), '-i', self.white, '-o', str(curvature_tmp_texfile), '-m', '1', '-T', '0']
        apply(context.system, command)

        # on compte combien de composantes connexes
        compconned_texture_to_be_filtered = aims.read(str(curvature_tmp_texfile))
        nb_cc = count_label_max(compconned_texture_to_be_filtered[0])

        # mais d'abord on filtre par la texture d'intersection
        compconn = {}
        
        for i in xrange(len(compconned_texture_to_be_filtered[0])):
          label = compconned_texture_to_be_filtered[0][i]
          if (label != -1):
            if (compconn.has_key(label) is False):
              compconn[label] = []
            compconn[label].append(i)

        nb_inters_labels = 0
        nb_noninters_labels = 0
        for label in compconn.keys():
          is_one_inters = False
          for node in compconn[label]:
            if (intersection_texture[0][node] != 0):
              is_one_inters = True
          if (is_one_inters):
            nb_inters_labels = nb_inters_labels + 1
          else :
            nb_noninters_labels = nb_noninters_labels + 1
            # on efface la compconn de la texture
            for node in compconn[label]:
              compconned_texture_to_be_filtered[0][node] = -1

        all_thresholds_cc_tex.append(compconned_texture_to_be_filtered)

        context.write(str(nb_cc) + " " + str(nb_inters_labels) + " " + str(nb_noninters_labels))

        #context.write( "scale :" + str(scale_number) + " thresh : "+ str(threshold) +
           #" label_max : " + str (label_max) )
        

        item = []
        item.append(sujet)
        item.append(scale_number)
        item.append(-threshold)
        item.append(nb_cc)
        item.append(nb_inters_labels)
        item.append(nb_noninters_labels)
        data.append(item)
        iteration = iteration + 1

      # ecriture de la texture all thresholds cc
      temptexture = aims.TimeTexture_FLOAT(len(all_thresholds_cc_tex), 
         len(scale_space[scale_number]))
      for i in xrange(len(all_thresholds_cc_tex)):
        temptexture[i] = all_thresholds_cc_tex[i][0]
      wtex.write(temptexture, "/volatile/operto/temp/temp_tex/%s_%i.tex" % (sujet,scale_number))

    
    csv_file = self.csv_file
    writer = csv.writer( open(str(csv_file),'w') )
    writer.writerows(data)
    #data = np.recfromcsv( str(csv_file) )
    #print "test"
    #for scale, max_list in reader:
      #for i in max_list:
        #print ("echelle:" + str(scale) + "max:" + str( i))

      

     
     #call_list = [ 'surfLabelsTex2Graph',
                #'-m', self.white,
                #'-t', self.texture,                   
                #'-g', self.graph,
                #'--lat', self.latitude,
                #'--lon', self.longitude,
                #'-s', s,
                #'--repM', self.repMesh]
     #context.write ( "Executing process.." )
     #apply(context.system, call_list)
     #context.write('Finished')

