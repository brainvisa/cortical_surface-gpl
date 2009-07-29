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


name = 'Surface-based statistical parametric map'

userLevel = 2

signature = Signature(
  'projection_texture', ReadDiskItem( 'Texture', 'Texture'),
  'protocol_text', ReadDiskItem( 'Text File', 'Text File' ),
  'beta', WriteDiskItem('Texture', 'Texture'),
  'spmt_texture', WriteDiskItem( 'Texture', 'Texture'),
  'contraste', String()
)

         
def processtag(reader):
   testa = str(reader.Value())
   j=0
   k=1
   for test in string.split(testa,'\n'):
      if test == "" or test == "None":
         continue
      for i in test:
         if i not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']:
            j=j+1
            continue
         if i in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']:
            k=k+1
            continue
         break
      testi = int((test[j:k+j]))
      #print reader.Name()
      return testi

def streamFile(filename):
    import libxml2
    try:
        reader = libxml2.newTextReaderFilename(filename)
    except Exception, e:
        print e
        print "unable to open %s" % (filename)
        return
    TR=24
    t_time = []
    t_type = []
    hrf = []
    ret = reader.Read()
    while ret == 1:
        #processNode(reader)
        if reader.Name() == "time" and int(reader.NodeType()) == 1:
           ret = reader.Read()
           t_time.append(int(processtag(reader)))
        elif reader.Name() == "type" and int(reader.NodeType()) == 1:
           ret = reader.Read()
           t_type.append(int(processtag(reader)))
        elif reader.Name() == "tr" and int(reader.NodeType()) == 1:
          ret = reader.Read()
          TR = int(processtag(reader))
        else :
           ret = reader.Read()
    
    if ret != 0:
        print "%s : failed to parse" % (filename)
        
    return (TR,t_type,t_time)

        
def LogGamma(x): 
    import math
    coeff0 =  7.61800917300e+1
    coeff1 = -8.65053203300e+1
    coeff2 =  2.40140982200e+1
    coeff3 = -1.23173951600e+0
    coeff4 =  1.20858003000e-3
    coeff5 = -5.36382000000e-6
    stp    =  2.50662827465e+0
    half   =  5.00000000000e-1
    fourpf =  4.50000000000e+0
    one    =  1.00000000000e+0
    two    =  2.00000000000e+0 
    three  =  3.00000000000e+0
    four   =  4.00000000000e+0 
    five   =  5.00000000000e+0
    r = coeff0 / ( x ) + coeff1 / ( x + one   ) + coeff2 / ( x + two  ) + coeff3 / ( x + three ) + coeff4 / ( x + four ) + coeff5 / ( x + five  )
    s = x + fourpf
    t = ( x - half ) * math.log( s ) - s
    return t + math.log( stp * ( r + one ) )

def gammapdf(in1, g):
  import numpy as N
  import math
  terme1 = N.log(N.array(in1))* (g-1)
  terme2 =  N.array(in1) + LogGamma(g)
  out = terme1-terme2
  out = N.exp(out)
  return out
    
def get_hrf(TR, longueur):
   import numpy as N
   
   tp1 = 6.0
   tp2 = 12.5
   fwhm1 = 1.0
   fwhm2 = 1.0
   alp = .16
   
   dxA = N.arange(1,longueur+1,1)
   dxA = dxA*TR
   dxB = N.arange(1,longueur+1,1)
   dxB = dxB*TR
   #print dxA
   
   A = gammapdf(dxA, tp1)
   B = gammapdf(dxB, tp2)
   
   maxA = max(A)
   maxB = max(B)

   hf = A/maxA - alp*B/maxB
   return hf
  
def creer_prereg(condition, maxtime, hrfduration, t_type, t_time):
  import numpy as N
  prereg = N.zeros(int((maxtime + hrfduration)/100), float)
  #print len(prereg)
  for j in range(0,len(t_type)):
     if int(t_type[j]) == int(condition)+1:
        prereg[int(t_time[j]/100)] = 1
  return prereg

def initialization( self ):
  #self.projection_texture = "/home/olivier/time_Lwhite.tex"
  #self.protocol_text = "/home/olivier/localizer.txt"
  #self.spmt_texture = "/home/olivier/t_tex.tex"
  self.setOptional('beta')

  self.contraste = "0 0 1 1 -1 -1 1 -1 -1 1"

def execution( self, context ):
  import fff.glm as G
  import numpy as N
  from soma import aims
  import sys,os,string 

  #########################
  
  reader2 = aims.Reader()
  texfilename = self.projection_texture
  texture = reader2.read(str(texfilename))
  nb_nodes = int(texture[0].nItem())
  nb_scans = int(texture.size())
 
  tab = N.arange(float(nb_nodes*nb_scans))
  k=0
  baseline = N.zeros(nb_scans)
  
  for i in range(0,nb_nodes):
     for j in range(0,nb_scans):
        tab[k] = texture[j][i]
        k=k+1 
        baseline[j] += texture[j][i]
     
  baseline /= nb_nodes
  
  tab = tab.reshape(nb_nodes,nb_scans)
  
  (TR,t_type,t_time) = streamFile(str(self.protocol_text))
  #print t_type
  #print t_time
  #print TR
  nb_cond = int(N.max(t_type))
  maxtime = int(N.max(t_time))
  lentype = len(t_type)
  #print lentype
  #print str(len(t_time))
  
  hrf = get_hrf(0.1, 250)
  #print 'hrf'
  #print hrf
  #print ' '
  reg = N.zeros((nb_cond+1)*nb_scans, float)
  reg = reg.reshape(nb_scans, (nb_cond+1))
  
  for i in range(0,nb_cond):
     prereg = creer_prereg(i,maxtime, 250, t_type, t_time)
     hrf_aux = N.convolve(prereg,hrf)
     #print 'prereg' + str(len(prereg))
     #for j in range(0, len(prereg)):
       #print prereg[j]
     #print 'hrf_aux' + str(len(hrf_aux))
     #for j in range(0, len(hrf_aux)):
       #print hrf_aux[j]
     #print int(len(hrf_aux))
     #print 'cond' + str(i)
     for j in range(0,nb_scans):
         ##print str(j) + ' ' +str(int(j*TR))
         #print str(float(hrf_aux[int(j*TR)]))
         reg[j,i] = float(hrf_aux[int(j*TR)])
   
  reg[:,nb_cond] = baseline
  #for i in range(0,nb_cond):
    #print 'reg' + str(i)
    #print reg[:,i]
  #B, nVB, s2, dof = G.KF_fit(tab, reg, axis=1)
  #model = G.glm(tab, reg, axis=1)
  B, nVB, s2, a, dof = G.RKF_fit(tab, reg, axis=1)

  #model.fit(method='kalman.ar1') ## default is 'ols'
  #print model.beta[:,1].size
  #print model.beta[1].size
  c = N.zeros(int(nb_cond+1))
  j=0
  for i in string.split(str(self.contraste)):
    c[j] = int(i)
    j=j+1


  cB, VcB, t = G.t_contrast(c, B, nVB, s2, axis=1)

  #tcon = model.contrast(c) ## recognizes a t contrast
  #t, p, z = tcon.test(zscore=True)
  #cB, VcB, t = G.t_contrast(c, B, nVB, s2, axis=1 )
  
  writer = aims.Writer()
  textur = aims.TimeTexture_FLOAT()
  textur[0] = aims.Texture_FLOAT(int(t.size))
  
  for i in range(0,t.size):
     textur[0][i] = float(t[i])
  writer.write(textur, str(self.spmt_texture))

  #if self.beta is not '' :
    #betatex = aims.TimeTexture_FLOAT()
    #for i in range(0,int(B.size)):
      #betatex[i] = aims.Texture_FLOAT(int(B[:,1].size))
      #for j in range(0, B[:,1].size):
        #betatex[i][j] = float(B[j][i])
    #writer2 = aims.Writer()
    #writer2.write(betatex, str(self.beta))
  
    
