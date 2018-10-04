#!/bin/env python

import os
import sys
import argparse

def main() :

    parser = argparse.ArgumentParser(description = 'From SusyNtuple\'s TriggerList, generate Superflow code for getting the variables')
    parser.add_argument('--trigger-list', required = True, help = 'Provide TriggerList.h')
    args = parser.parse_args()

    if not os.path.isfile(args.trigger_list) :
        print 'ERROR TriggerList file provided (=%s) not found' % args.trigger_list
        sys.exit()
    full_file = os.path.abspath(args.trigger_list)

    trigger_names = []
    with open(full_file, 'r') as trigger_file :
        for line in trigger_file :
            line = line.strip()
            if line.startswith('//') : continue
            if '"HLT_' not in line : continue
            remove = 'HLT_'
            trig_name = line[line.find(remove) + len(remove) : ]
            trig_name = trig_name.split('"')[0]
            trigger_names.append(trig_name)

    # declarations
    for trigger in trigger_names :
        print 'bool p_%s;' % trigger

    # testing
    print '*cutflow << [&](Superlink* sl, var_void*) {'
    for trigger in trigger_names :
        print '\tp_%s = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, \"HLT_%s\");' % (trigger, trigger)
    print '};'

    # fill vars
    for trigger in trigger_names :
        print '*cutflow << NewVar("pass %s"); {' % trigger
        print '\t*cutflow << HFTname("trig_%s");' % trigger
        print '\t*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {'
        print '\t\treturn p_%s;' % trigger
        print '\t};'
        print '\t*cutflow << SaveVar();'
        print '}'

if __name__ == '__main__' :
    main()
