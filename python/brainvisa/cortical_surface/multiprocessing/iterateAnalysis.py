#!/usr/bin/env python

from __future__ import print_function

import os, string, shutil, sys
from soma import aims
from brainvisa.cortical_surface.multiprocessing.mproc import MultiProcExecute, os_system


def IterateAnalysis ( groupgraphpath, tmppathes, analyze, number_of_proc = 2 ) :
    ''' Specific function allowing to run various structural analyses on various
    processors '''
    nb_iter = len(tmppathes)
    #[ddweight, intrapsweight, simweight, lsweight, ddx1, ddx2, simx1, simx2, ddh ] = [float(params[i]) for i in range(0,9)]

    #if ( not os.path.exists( tmpdir ) ):
        #os.makedirs( tmpdir )

    energies = []
    jobs = []
    
    for i in xrange ( nb_iter ) :
        ana = [each for each in analyze]
        tmppath = tmppathes[i]
        ana += [ '-i', str(groupgraphpath) ]
        ana += [ '-o', str(tmppath) ]
        ana += [ '--run', 'True' ]

        analyze_command = string.join( [ "'"+str(each)+"'" for each in ana ], ' ' )
        print(analyze_command)
        jobs.append ( (analyze_command, i, 'test') )

    print("JOBS:", len(jobs))
    results = MultiProcExecute ( os_system, jobs, number_of_proc )

def DistributeIndividualGraphs ( indiv_commands, number_of_proc = 2 ) :
    ''' Specific function allowing to distribute the computation of individual
    graphs on various processors '''
    jobs = []
    for i, indiv in enumerate( indiv_commands ) :
        indiv_command = string.join( [ "'"+str(each)+"'" for each in indiv ], ' ' )
        print(indiv_command)
        jobs.append ( (indiv_command, i, 'test') )

    print("JOBS:", len(jobs))
    results = MultiProcExecute ( os_system, jobs, number_of_proc )



