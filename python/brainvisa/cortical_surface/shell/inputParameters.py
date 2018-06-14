from __future__ import print_function

import os

def InputParameters( arguments, messages_defaults ):
    params = []
    nb_arguments = len( arguments )
    messages = []
    defaults = {}
    for each in messages_defaults:
        messages.append(each[0])
        defaults[each[0]] = each[1]
    assert(len(messages) == len(defaults))

    params = [ arguments[i] for i in xrange(nb_arguments) ]
    doAsk = True

    if ( nb_arguments < len (defaults) ):
        inputs = []
        for i in xrange( len(defaults) - nb_arguments) :
            if (doAsk):
                inp = raw_input(messages[i + nb_arguments] + str( ' ? (%s) '%defaults[messages[ i + nb_arguments ] ] ) )
            else :
                inp = ''
            if inp in ['OK', 'ok', 'Ok']:
                doAsk = False
                inp = ''
            inputs.append ( inp  )
        for i, p in enumerate(inputs):
            if p != '':
                params.append(p)
            else :
                params.append(defaults[ messages[i + nb_arguments] ])
    print(params)

    print('\n params')
    for i, p in enumerate(params):
        print(messages[i], p)
        defaults[messages[i]] = p

    return params
