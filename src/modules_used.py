#------------------------------------------------------------------------------------------------------------
"""modules_used.py - go through all f90 files in current folder and print includes and modules used."""
#
#------------------------------------------------------------------------------------------------------------

import re
import os

#------------------------------------------------------------------------------------------------------------
# global definitions/configuration

outfile      = 'modules.incl'
ignore_files = ['contact.f90', 'contact_addon.f90', 'caddon_license.f90', 'mkl_dfti.f90', 
                'test_caddon.f90', 'test_float.f90', 'test_mbench.f90', 'test_memory.f90',
                'test_table.f90', 'usetab_table.f90' ]
ignore_mods  = ['iso_c_binding', 'ifport', 'omp_lib' ]

# end global definitions

#------------------------------------------------------------------------------------------------------------

def parse_fortran_file( fname, f_out, idebug=3 ):
    """function [ ] = parse_fortran_file( fname, f_out, [idebug] )
                        open file, parse contents, print includes and modules used"""

    print('Parsing contents of "%s"' % fname)

    incl_used = []
    mods_used = []

    # parse contents of the file

    f = open(fname, 'r')
    iline = 0
    for line in f:
        iline = iline + 1
#       if (idebug>=3):
#           print('Line %4d:' % iline, line.strip('\n'))

        if (re.search('^ *include ',line) or re.search('^# *include',line)):
            print('Line %4d:' % iline, line.strip('\n'))
            m = re.search('[\'"](.*)[\'"]',line)
            inc = line[m.start(1):m.end(1)]
            # print('found include: i=[%d:%d], inc=%s' % (m.start(1), m.end(1), inc))

            # check if include should be ignored

            if (inc in ignore_mods):
                is_ignored = 1
                # print('include-file %s is in ignore-list' % inc)
            else:
                is_ignored = 0

            if (not is_ignored):
                incl_used.append( line[m.start(1):m.end(1)] )

        # search all lines starting with 'submodule'
        if (re.search('^ *submodule', line)):
            print('Line %4d:' % iline, line.strip('\n'))

            # extract the module name
            _, mod, sub = line.split()
            mod = mod[1:-1] # discard parentheses
            print('         module = %s, submodule = %s' % (mod, sub))
            mods_used.append( mod )

        # search all lines starting with 'use ' or 'use, ', avoid 'use_plast ='
        if (re.search('^ *use[ ,]', line)):

            # recognize intrinsic modules, esp. iso_c_binding

            m = re.search(', *intrinsic *::',line)
            if (m):
                is_intrinsic = 1
                line = line[0:m.start()] + line[m.end():-1]
                # print('found intrinsic: i=[%d:%d]' % (m.start(), m.end()))
                # print('stripped:', line.strip('\n'))
            else:
                is_intrinsic = 0

            # remove any ', only' clause

            m = re.search(', *only',line)
            if (m):
                line = line[0:m.start()]
                # print('found only: i=[%d:%d]' % (m.start(), m.end()))
                # print('stripped:', line.strip('\n'))

            # extract the module name

            _, mod = line.split()

            # check if module should be ignored

            if (mod in ignore_mods):
                is_ignored = 1
                # print('module %s is in ignore-list' % mod)
            else:
                is_ignored = 0

            # add module to list of modules used for this file

            if (not is_ignored):
                mods_used.append( mod )
                # print('Line %4d:' % iline, line.strip('\n'))

        # endif (use module)

    # end for (line in f)

    f.close()

    space = '                             '
    # print('File %s uses %d modules:' % (fname, len(mods_used)))
    # if (len(mods_used)>0):
    #     print(mods_used)

    if (len(incl_used)+len(mods_used)>0):
        base, ext = os.path.splitext( fname )
        l = len(base)
        f_out.write('$(P)%s.mod %s $(P)%s.$O: %s' % (base, space[0:16-l], base, space[0:16-l] ))
        for m in incl_used:
            l = len(m)
            f_out.write(' \\\n                                                    %s %s' % (m, space[0:24-l]))
        for m in mods_used:
            l = len(m)
            f_out.write(' \\\n                                                    $(P)%s.mod %s' % (m, space[0:16-l]))
        f_out.write('\n')

    return

# end function parse_fortran_file

#------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    f_out = open(outfile, 'w')

    for file in sorted(os.listdir(".")):
        if (file in ignore_files):
            print('Ignoring file %s (ignore list)' % file)
        elif (file.endswith(".f90")):
            parse_fortran_file(file, f_out)
        else:
            print('Ignoring file %s (no Fortran)' % file)

    f_out.close()

# end main program "modules_used

#------------------------------------------------------------------------------------------------------------
