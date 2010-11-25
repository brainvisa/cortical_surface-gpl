#!/usr/bin/env python

import os, string, shutil, sys 
from soma import aims
from brainvisa.cortical_surface.multiprocessing import MultiProcExecute, os_system



def IterateAnalysis ( groupgraphpath, tmpdir, analyze, nb_iter, params, number_of_proc = 2 ) :
    ''' Specific function allowing to run various structural analyses on various
    processors '''
    [ddweight, intrapsweight, simweight, lsweight, ddx1, ddx2, simx1, simx2, ddh, globalweight ] = [float(params[i]) for i in range(0,10)]

    if ( not os.path.exists( tmpdir ) ):
        os.makedirs( tmpdir )

    energies = []
    jobs = []
    tmppathes = [ str(os.path.join ( tmpdir, 'tmp_analysisgraph_%i.arg'%i ) ) for i in xrange(nb_iter) ]
    
    for i in xrange ( nb_iter ) :
        ana = [each for each in analyze]
        tmppath = tmppathes[i]
        ana += [ '-i', str(groupgraphpath) ]
        ana += [ '-o', str(tmppath) ]
        ana += [ '--run', 'True' ]

        analyze_command = string.join( [ "'"+str(each)+"'" for each in ana ], ' ' )
        print analyze_command
        jobs.append ( (analyze_command, i, 'test') )

    print "JOBS:", len(jobs)
    results = MultiProcExecute ( os_system, jobs, number_of_proc )

    return tmppathes

