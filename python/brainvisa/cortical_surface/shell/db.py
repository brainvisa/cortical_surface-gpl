import os

def SubjectsSubset ( subjects, nb ) :
    '''Returns a random subset of size nb from a set of subjects names'''
    from random import randrange
    new_data = [d for d in subjects]
    assert ( nb <= len(subjects) )
    comp = []
    for n in xrange ( nb ) :
        pos = randrange ( len(new_data) )
        elem = new_data[pos]
        new_data[pos] = new_data[-1]
        del new_data[-1]
        comp.append( elem )
    return comp

def getSubjects ( db ) :
    """Returns the list of subjects (dirnames) contained in the data base directory
         (ignores any existing 'tmp')"""
    sujets = []
    for direc in os.listdir ( db ):
        if os.path.isdir( os.path.join(db, direc) ) and direc != 'tmp':
            sujets.append ( direc )
    return sujets

def getPathes ( db, hemis_side = 'L' , data_type = 'curv', subjects = None ) :
    """Builds and returns the pathes according to the data base path, hemisphere side
         and data type (curvature or depth)"""
    pathes = {}
    if subjects is None :
        subjects = getSubjects ( db )
    
    for sujet in subjects :
        pathes[sujet] = {}
        
        pathes[sujet]['mesh'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/mesh/%s_%swhite.mesh' % ( sujet, hemis_side ))
        pathes[sujet]['tex'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/%s_%swhite_%s.tex' % ( sujet, hemis_side, data_type ) )
        pathes[sujet]['primal'] = os.path.join( db, sujet, 'surface/%s_%swhite_primal_%s.arg' % ( sujet, hemis_side, data_type ) )
        
        pathes[sujet]['blobs'] = os.path.join( db, sujet, 'surface/%s_%swhite_%s_blobs.tex' % ( sujet, hemis_side, data_type ) )
        pathes[sujet]['scalespace'] = os.path.join ( db, sujet, 'surface/%s_%swhite_%s_ss.tex' % ( sujet, hemis_side, data_type ) )
        pathes[sujet]['lat'] = os.path.join ( db, sujet, 'surface/%s_%s_lat.tex' % ( sujet, hemis_side ) )
        pathes[sujet]['lon'] = os.path.join ( db, sujet, 'surface/%s_%s_lon.tex' % ( sujet, hemis_side ) )
        pathes[sujet]['gyri'] = os.path.join ( db, sujet, 'surface/%s_%s_gyri.tex' % ( sujet, hemis_side ) )
        
    return pathes

def getfMRIPathes ( db, contrast, subjects = None ) :
    """Builds and returns the pathes according to the data base path and a
        contrast"""
    pathes = {}
    if subjects is None :
        subjects = getSubjects ( db )
    for sujet in subjects :
        pathes[sujet] = {}
        pathes[sujet]['fMRI'] = os.path.join ( db, sujet, 'fMRI/default_acquisition/loc1/spm_analysis/spm_analysis_Norm_S/spmT_00%i.img' % contrast )
        pathes[sujet]['mask'] = os.path.join ( db, sujet, 'fMRI/default_acquisition/loc1/spm_analysis/spm_analysis_Norm_S/mask.img' )
        pathes[sujet]['primal'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/%s_spmT_00%i_primal.arg' % (sujet, contrast) )
        pathes[sujet]['scale'] = os.path.join ( db, sujet, 't1mri/t1/default_analysis/segmentation/%s_spmT_00%i_ss.ima' % (sujet, contrast) )

    return pathes