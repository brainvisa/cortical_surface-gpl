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
    'csv_file', WriteDiskItem('Text File', 'Text File')
    )

def initialization( self ):
  pass

def execution( self, context ):
    sujet = self.white.get( 'subject')
    context.write(sujet)

    # on lit la texture d'espace echelle et on cree un fichier temporaire
    scale_space = aims.read(str(self.scale_space))
    curvature_tmp_texfile = context.temporary("Texture")
    context.write( "temporary : " + str (curvature_tmp_texfile) )
    nb_scales = scale_space.size()
    context.write( str(nb_scales) + " scales" )


    data = []
    data.append(["scale_number", "threshold", "cc"])
    # pour chaque niveau d'echelle
    for scale_number in xrange(nb_scales) :
      context.write( "scale :" + str(scale_number) )
      threshold = 0.0
      label_max = 1

      all_thresholds_cc_tex = []
      iteration = 0
      
      while (label_max != 0):
        
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
        
        # on calcule les composantes connexes
        command = [ 'AimsMeshConnectedComponent', '-t', str(curvature_tmp_texfile), '-i', self.white, '-o', str(curvature_tmp_texfile), '-m', '1', '-T', '0']
        apply(context.system, command)

        ##on calcule un voronoi pour tester
        #command = [ 'AimsTextureVoronoi', '-m', self.white, '-t', str(curvature_tmp_texfile), '-o', "/volatile/operto/voronoi.tex"]
        #apply(context.system, command)

        
        # on compte combien de composantes connexes
        compconn_texture = aims.read(str(curvature_tmp_texfile))

        all_thresholds_cc_tex.append(compconn_texture)
        
        label_max = 0
        for i in xrange(len(compconn_texture[0])):
          if (compconn_texture[0][i] > label_max):
            label_max = compconn_texture[0][i]

        #context.write( "scale :" + str(scale_number) + " thresh : "+ str(threshold) +
           #" label_max : " + str (label_max) )
        

        item = []
        item.append(scale_number)
        item.append(-threshold)
        item.append(label_max)
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

