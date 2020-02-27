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
from numpy import *
from soma import aims
import six


name = 'Average Sulcus Morphological Curves'

userLevel = 0

signature = Signature(
    'depth_profiles', ListOf(ReadDiskItem('Sulcus depth profile', 'Text file' )),
    'average_profiles', WriteDiskItem('Text file','Text file'),
)

def initialization( self ):
    pass

def execution( self, context ):
    ns=len(self.depth_profiles)
    ldepth=zeros((101,ns))
    lprof=zeros((101, ns))
    context.write("ldepth shape", ldepth.shape)
    context.write("lprof shape", lprof.shape)
    context.write("Found ", ns, " profiles")
    for i in six.moves.xrange(ns) :
        context.write('i= ', i)
        context.write("Reading ", self.depth_profiles[i].fullPath())
        data=loadtxt(self.depth_profiles[i].fullPath())
        depth=data[:,1]
        profile=data[:,2]
        context.write("profile shape:", profile.shape)
        context.write("depth shape:", depth.shape)
        x=data[:,0]
        context.write("ldepth shape", ldepth.shape)
        context.write("lprof shape", lprof.shape)
        ldepth[:, i]=depth
        lprof[:, i]=profile
    mdepth=ldepth.mean(axis=1)
    mprof=lprof.mean(axis=1)
    sdepth=ldepth.std(axis=1)
    sprof=lprof.std(axis=1)
    mx=x

    context.write('Computing stats')

    mdata=zeros((101,5))
    mdata[:,0]=mx
    mdata[:,1]=mdepth
    mdata[:,2]=sdepth
    mdata[:,3]=mprof
    mdata[:,4]=sprof

    savetxt(self.average_profiles.fullPath(), mdata, delimiter='\t' )
    context.write('Finished')

