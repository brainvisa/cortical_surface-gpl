#!/usr/bin/env python

import os, sys
from soma import aims
import brainvisa.cortical_surface.structural.blobManip as oo
from brainvisa.cortical_surface.shell.inputParameters import *
from brainvisa.cortical_surface.multiprocessing import MultiProcExecute, os_system


def getfMRIPathes ( db, contrast, subjects = None ) :
    """Builds and returns the pathes according to the data base path and a
        contrast"""
    pathes = {}
    if subjects is None :
        subjects = oo.getSubjects ( db )
    for sujet in subjects :
        pathes[sujet] = {}
        pathes[sujet]['fMRI'] = os.path.join ( db, sujet, 'fMRI/default_acquisition/loc1/spm_analysis/spm_analysis_Norm_S/spmT_00%i.img' % contrast )
        pathes[sujet]['mask'] = os.path.join ( db, sujet, 'fMRI/default_acquisition/loc1/spm_analysis/spm_analysis_Norm_S/mask.img' )
        pathes[sujet]['primal'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/%s_spmT_00%i_primal.arg' % (sujet, contrast) )
        pathes[sujet]['scale'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/%s_spmT_00%i_ss.ima' % (sujet, contrast) )

    return pathes

messages_defaults = [ ('db path', '/home/go224932/data/structural/database_miccai09'),
                      ('contrast', '24'),
                      ('number of proc', '8') ]              

if __name__ == '__main__':

    params = InputParameters( sys.argv[1:], messages_defaults )
    if (raw_input('OK? y/n') != 'y'):
        sys.exit(0)

    db = str(params[0])
    contrast = int(params[1])
    number_of_proc = int(params[2])

    sujets = oo.getSubjects( db )

    pathes = getfMRIPathes ( db, contrast )
    jobs = []
    res = {}

    for sujet in sujets:
        res[sujet] = {}
        primal_command = 'volIma2Graph -i  '+ str(pathes[sujet]['fMRI']) +'  -m '+ str(pathes[sujet]['mask']) +' -g ' \
            + str(pathes[sujet]['primal']) +' --subject '+str(sujet)+' --ss '+ str(pathes[sujet]['scale']) +' --recover False '
        print(primal_command)
        jobs.append ( (primal_command, sujet, 'test') )

    print("JOBS:", len(jobs))
    results = MultiProcExecute ( os_system, jobs, number_of_proc )

    for each in results :
        (keyword, sujet, resultat) = each
        res[sujet] = resultat
