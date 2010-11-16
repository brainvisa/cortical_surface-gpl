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
from neuroProcesses import *
import shfjGlobals     
import numpy as np
from soma import aims    
import sys, os, imp#, kalman
#from nipy.neurospin.glm.glm import glm
#from nipy.neurospin.glm import kalman

name = '1 - Create Surface-Based Statistical Parametric Maps'
userLevel = 2

signature = Signature(  #'meshes', ListOf(ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh' )),
    'boldtextures', ListOf(ReadDiskItem('Functional Time Texture', 'Texture')),
    'protocolfile', ReadDiskItem( 'Text File', 'Text File' ),
    'contrast', String(),
    #'contrastname', String(),
    'betamaps', ListOf(WriteDiskItem('Surface-Based Beta Map', 'Texture')),
    'spmtmaps', ListOf(WriteDiskItem( 'Surface-Based SPMt Map', 'Texture')),
  )

DEF_TINY = 1e-50
DEF_DOFMAX = 1e10

models = {'spherical':['ols']}
#models =  {'spherical':['kalman']}

def ols ( Y, X, axis = 0 ) :
    """ beta, nvbeta, s2, dof = ols(Y, X, axis=0)
    Essentially, compute pinv(X)*Y """
    ndims = len(Y.shape)
    pX = np.linalg.pinv(X)
    beta = np.rollaxis(np.inner(pX, np.rollaxis(Y, axis, ndims)), 0, axis+1)
    nvbeta = np.inner(pX, pX)
    res = Y - np.rollaxis(np.inner(X, np.rollaxis(beta, axis, ndims)), 0, axis+1)
    n = res.shape[axis]
    s2 = (res**2).sum(axis) / float(n-X.shape[1])
    dof = float(X.shape[0] - X.shape[1])
    return beta, nvbeta, s2, dof

class glm:
    def __init__ ( self, Y = None, X = None, formula = None, axis = 0,
              model = 'spherical', method = None, niter = 2 ) :

        # Check dimensions
        if Y == None:
            return
        else:
            self.fit(Y, X, formula, axis, model, method, niter)

    def fit ( self, Y, X, formula = None, axis = 0, model = 'spherical', method = None, niter = 2 ) :

        if Y.shape[axis] != X.shape[0]:
            raise ValueError, 'Response and predictors are inconsistent'

        # Find model type
        self._axis = axis
        if isinstance(formula, str):
            model = 'mfx'
        if models.has_key(model):
            self.model = model
            if method == None:
                self.method = models[model][0]
            elif models[model].count(method):
                self.method = method
            else:
                raise ValueError, 'Unknown method'
        else:
            raise ValueError, 'Unknown model'

        # Initialize fields
        constants = []
        a = 0

        # Switch on models / methods
        if self.model == 'spherical':
            constants = ['nvbeta', 'a']
            if self.method == 'ols':
                out = ols(Y, X, axis=axis)
            #elif self.method == 'kalman':
                #out = kalman.ols(Y, X, axis=axis)
        #elif self.model == 'ar1':
            #constants = ['a']
            #out = kalman.ar1(Y, X, axis=axis, niter=niter)
            #a = out[4]
            #out = out[0:4]


        # Finalize
        self.beta, self.nvbeta, self.s2, self.dof = out
        self.s2 = self.s2.squeeze()
        self.a = a
        self._constants = constants

    """
    Save fit into a .npz file
    """
    def save(self, file):
        np.savez(file,
            beta=self.beta,
            nvbeta=self.nvbeta,
            s2=self.s2,
            dof=self.dof,
            a=self.a,
            model=self.model,
            method=self.method,
            axis=self._axis,
            constants=self._constants)


    """
    c must be a numpy.ndarray (or anything that numpy.asarray
    can cast to a ndarray).
    For a F contrast, c must be q x p where q is the number of contrast vectors and
    p is the total number of regressors.
    """
    def contrast(self, c, type='t', tiny=DEF_TINY, dofmax=DEF_DOFMAX):
        import numpy as np
        c = np.asarray(c)
        #dim = len(c.shape)
        if c.ndim == 1:
            dim=1
        else:
            dim = c.shape[0]
        axis = self._axis
        ndims = len(self.beta.shape)

        # Compute the contrast estimate: c*B
        B = np.rollaxis(self.beta, axis, ndims)
        con = np.inner(c, B) # q, X

        # Compute the variance of the contrast estimate: s2 * (c' * nvbeta * c)
        # Two cases are considered: either the input effect variance
        # is position-dependent (output by RKF_fit), or it is a global
        # one (output by KF_fit)
        s2 = self.s2.squeeze()
        nvbeta = self.nvbeta
        if not 'nvbeta' in self._constants:
            nvbeta = np.rollaxis(nvbeta, axis, ndims+1)
            nvbeta = np.rollaxis(nvbeta, axis, ndims+1) # X, p, p
        if dim == 1:
            vcon = np.inner(c, np.inner(c, nvbeta))
            vcon = vcon.squeeze()*s2
        else:
            vcon = np.dot(c, np.inner(nvbeta, c)) # q, X, q or q, q
            if not 'nvbeta' in self._constants:
                vcon = np.rollaxis(vcon, ndims, 1)*s2 # q, q, X
            else:
                aux = vcon.shape # q, q
                vcon = np.resize(vcon, s2.shape+aux) # X, q, q
                vcon = vcon.transpose().reshape(aux+(s2.size,))*s2.reshape((s2.size,)) # q, q, Xflat
                vcon = vcon.reshape(aux+s2.shape) # q, q, X
        # Create contrast instance
        c = contrast(dim, type, tiny, dofmax)
        c.effect = con
        c.variance = vcon
        c.dof = self.dof
        return c

class contrast:

    def __init__(self, dim, type='t', tiny=DEF_TINY, dofmax=DEF_DOFMAX):
        """
        tiny is a numerical constant for computations.
        """
        self.dim = dim
        self.effect = None
        self.variance = None
        self.dof = None
        if dim > 1:
            if type is 't':
                type = 'F'
        self.type = type
        self._stat = None
        self._pvalue = None
        self._baseline = 0
        self._tiny = tiny
        self._dofmax = dofmax

    def summary(self):
        """
        Return a dictionary containing the estimated contrast effect,
        the associated ReML-based estimation variance, and the estimated
        degrees	of freedom (variance of the variance).
        """
        return {'effect':self.effect, 'variance':self.variance, 'dof':self.dof}

    def stat(self, baseline=0.0):
        import numpy as np
        """
        Return the decision statistic associated with the test of the
        null hypothesis: (H0) 'contrast equals baseline'
        """
        self._baseline = baseline

        # Case: one-dimensional contrast ==> t or t**2
        if self.dim == 1:
            # avoids division by zero
            t = (self.effect-baseline) / np.sqrt(np.maximum(self.variance, self._tiny))
            if self.type == 'F':
                t = t**2
        # Unknwon stat
        else:
            raise ValueError, 'Unknown statistic type'

        self._stat = t
        return t


    def pvalue(self, baseline=0.0):
        import scipy.stats as sps
        import numpy as np

        """
        Return a parametric approximation of the p-value associated
        with the null hypothesis: (H0) 'contrast equals baseline'
        """
        if self._stat == None or not self._baseline == baseline:
            self._stat = self.stat(baseline)
        # Valid conjunction as in Nichols et al, Neuroimage 25, 2005.
        if self.type in ['t', 'tmin']:
            p = sps.t.sf(self._stat, np.minimum(self.dof, self._dofmax))
        elif self.type == 'F':
            p = sps.f.sf(self._stat, self.dim, np.minimum(self.dof, self._dofmax))
        else:
            raise ValueError, 'Unknown statistic type'
        self._pvalue = p
        return p

    def __add__(self, other):
        if self.dim != other.dim:
            return None
        con = contrast(self.dim)
        con.type = self.type
        con.effect = self.effect + other.effect
        con.variance = self.variance + other.variance
        con.dof = self.dof + other.dof
        return con

    def __rmul__(self, other):
        k = float(other)
        con = contrast(self.dim)
        con.type = self.type
        con.effect = k * self.effect
        con.variance = k**2 * self.variance
        con.dof = self.dof
        return con

    __mul__ = __rmul__

    def __div__(self, other):
        return self.__rmul__(1/float(other))

def getContrastName(self, data):
    if (self.contrast_name is not None and len(self.meshes) == len(self.BOLD_textures) and len(self.meshes)>0):
        result = []
        for mesh in self.meshes:
            attributes = mesh.hierarchyAttributes()
            attributes[ 'contrast' ] = str(self.contrast_name)
            print attributes[ 'contrast' ]
            res = self.signature[ 'spmt_maps' ].findValue( attributes )
            print res[0]
            result.append( res[0] )
        return result
    return None

def getBetaName(self, data):
    if (self.contrast_name is not None and len(self.meshes) == len(self.BOLD_textures) and len(self.meshes)>0):
        result = []
        for mesh in self.meshes:
            attributes = mesh.hierarchyAttributes()
            attributes[ 'contrast' ] = str(self.contrast_name)
            print attributes[ 'contrast' ]
            res = self.signature[ 'beta' ].findValue( attributes )
            print res[0]
            result.append( res[0] )
        return result
    return None

def initialization( self ):
    pass
    #self.setOptional ( 'beta_maps' )
    #self.setOptional ( 'meshes' )
    #self.linkParameters ( 'boldtextures', 'meshes' )
    #self.addLink ( 'spmt_maps', 'meshes', self.getContrastName )
    #self.addLink ( 'spmt_maps', 'contrast_name', self.getContrastName )
    #self.addLink ( 'beta_maps', 'meshes', self.getBetaName )
    #self.addLink ( 'beta_maps', 'contrast_name', self.getBetaName )

def LogGamma ( x ) : 
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

def GammaPdf ( in1, g ) :
    terme1 = np.log ( np.array(in1) ) * (g-1)
    terme2 =  np.array(in1) + LogGamma(g)
    out = terme1 - terme2
    out = np.exp ( out )
    return out

def HrfFunction ( sampling_rate ) :
    ''' get_hrf returns an array with estimated samples of a canonic Haemodynamic
        Response Function (HRF)
        sampling_rate is the time interval between two samples
        number_of_samples is the size of the output array, hence allowing to tune the
            time of sampling'''
    tp1 = 6.0
    tp2 = 12.5
    fwhm1 = 1.0
    fwhm2 = 1.0
    alp = .16
    duration_hrf = 25.0
    number_of_samples = duration_hrf / sampling_rate

    dxA = np.linspace ( 1.0, duration_hrf + 1, number_of_samples, endpoint=False)
    dxB = np.linspace ( 1.0, duration_hrf + 1, number_of_samples, endpoint=False)
    print dxA, dxB
    A = GammaPdf ( dxA, tp1 )
    B = GammaPdf ( dxB, tp2 )

    maxA = max ( A )
    maxB = max ( B )
    
    hf = A / maxA - alp * B / maxB
    return hf

def PreRegressor ( condition, types, times, conversion_rate ):
    ''' builds and returns a pre-version of regressor (made of boxcars or diracs) according to :
        - condition : a given condition index
        - types : the sequence of stimulations, by types
        - times : the sequence of stimulations, by times
        - conversion_rate gives the conversion factor to apply to onsets/durations
            given in the protocol file. (for instance, if onsets are in ms and the
            preregressor should be sampled in s, conversion_rate equals 0.001 )'''
    assert (len(times)==len(types))
    assert (len(times) > 0)
    
    try:
        if ( len(times[0]) == 2 ) :
            mode = 'EPOCH'
    except TypeError :
        mode = 'EVENT'
    print 'mode:', mode
    ''' A specific mode is triggered according to the syntax used in the protocol file
        times = [(onset1,duration1), ...] triggers the EPOCH omde
        times = [onset1, onset2, ...] triggers the EVENT mode'''
    
 
    if ( mode == 'EVENT' ) :
        prereg = np.zeros ( int ( times[-1] * conversion_rate ) + 1 , float )
        for j in xrange( len(types) ):
            if int(types[j]) == int(condition) :
                prereg [ int(times[j] * conversion_rate ) ] = 1
                
    elif ( mode == 'EPOCH' ):
        onsets = [ each[0] for each in times ]
        durations = [ each[1] for each in times ]
        size_prereg = int ( ( onsets[-1] + durations[-1] ) * conversion_rate ) + 1
        prereg = np.zeros ( size_prereg, float )
        for j in xrange ( len(types) ):
            if int(types[j]) == int(condition) :
                prereg [ int(onsets[j] * conversion_rate ) ] = 1
                for k in xrange ( int( durations[j] * conversion_rate ) ) :
                    prereg [ int(onsets[j] * conversion_rate + k ) ] = 1
    return prereg
        

def execution ( self, context ) :

    execfile ( self.protocolfile.fullPath(), locals(), globals() )
    nptimes = np.array(times)
    nptypes = np.array(types)
    
    context.write ( 'TR (must be in ms):', TR )
    context.write ( 'times (must be in ms):', nptimes )
    context.write ( 'types:', nptypes )

    for index in xrange (len(self.boldtextures) ) :
        print 'texture number ', index
        #meshpath = self.meshes[index].fullPath()
        boldtexpath = self.boldtextures[index].fullPath()
        betapath = self.betamaps[index].fullPath()
        spmtpath = self.spmtmaps[index].fullPath()

        texture = aims.read( str( boldtexpath ) )
        nb_nodes = int ( texture[0].nItem() )
        nb_scans = int ( texture.size() )

        tab = np.arange ( float ( nb_nodes * nb_scans ) )
        k = 0
        baseline = np.zeros ( nb_scans )

        for i in range ( 0, nb_nodes ) :
            for j in range ( 0, nb_scans ) :
                tab[k] = texture[j][i]
                k = k + 1
                baseline[j] += texture[j][i]

        baseline /= nb_nodes

        tab = tab.reshape ( nb_nodes, nb_scans )
        nb_cond = nptypes.max()
        lentype = len(types)
        print nb_cond
        print lentype
        sampling_rate = 0.1 #in (10x)seconds, depends on how the onsets are defined
        # multiplying onsets by this factor should still give integers. If not,
        # then the sampling_rate is too high.
        conversion_rate = 0.001 / sampling_rate

        hrf = HrfFunction ( sampling_rate )
        ''' hrf contains the canonic HRF function sampled at sampling_rate (given in s)'''

        reg = np.zeros ( ( nb_cond + 1 ) * nb_scans, float )
        reg = reg.reshape ( nb_scans, ( nb_cond + 1 ) )

        for condition_index in xrange ( nb_cond ):
            prereg = PreRegressor ( condition_index + 1, types, times, conversion_rate )
            '''prereg contains a pre-version of the regressor made of boxcars or diracs'''
            hrf_aux = np.convolve ( hrf, prereg )
            '''hrf_aux contains a convoluted version of prereg'''
            aux_x = np.linspace ( 0.0, len(hrf_aux), len(hrf_aux) )
            '''aux_x is the sampling space of hrf_aux (normally [0, 1, 2, 3, ...,'''
            reg_x = np.linspace ( 0.0, nb_scans * TR * conversion_rate, nb_scans, endpoint=False )
            reg_aux = np.interp ( reg_x, aux_x, hrf_aux ).tolist()
            reg[:, condition_index] = reg_aux

        reg[:,nb_cond] = baseline

        data = tab
        data = data.transpose()

        m = glm ( data, reg, axis=0 )
        v = m.s2
        b = m.beta
        tex = aims.TimeTexture_FLOAT()
        for j in xrange(b.shape[0]):
            t = tex[j]
            t.assign(list(b[j,:]))
        aims.write ( tex, betapath )

        print self.contrast
        c = [int(i) for i in string.split(str(self.contrast))]
        print len(c), nb_cond
        c.extend ( np.zeros ( max( nb_cond + 1 - len(c), 0 ) ) )
        assert(len(c) == nb_cond + 1)
        print len(c), nb_cond
        tcon = m.contrast(c)

        tex = aims.TimeTexture_FLOAT()
        tmap = tcon.stat()
        t = tex[0]
        t.assign( tmap )

        aims.write ( tex, spmtpath )
        context.write("Finished")
