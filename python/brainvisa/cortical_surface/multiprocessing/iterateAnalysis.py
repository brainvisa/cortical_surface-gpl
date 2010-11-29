#!/usr/bin/env python

import os, string, shutil, sys 
from soma import aims
from brainvisa.cortical_surface.multiprocessing import MultiProcExecute, os_system



def IterateAnalysis ( groupgraphpath, tmppathes, analyze, params, number_of_proc = 2 ) :
    ''' Specific function allowing to run various structural analyses on various
    processors '''
    nb_iter = len(tmppathes)
    [ddweight, intrapsweight, simweight, lsweight, ddx1, ddx2, simx1, simx2, ddh, globalweight ] = [float(params[i]) for i in range(0,10)]

    if ( not os.path.exists( tmpdir ) ):
        os.makedirs( tmpdir )

    energies = []
    jobs = []
    
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


