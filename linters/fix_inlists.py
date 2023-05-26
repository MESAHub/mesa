#!/usr/bin/env python

import os
import sys
import re
import glob

MESA_DIR = "../"



def set_min_value(line,num='10'):
    if line.strip().startswith('!'):
        return line

    l = line.split('!')[0]
    interval = l.split('=')[1]
    if int(interval) < 10:
        line = line.replace(interval,f" {num} \n")

    return line

def set_flag(line,new_flag='.false.',old_flag='.true.'):
    if line.strip().startswith('!'):
        return line

    line = line.replace(old_flag,new_flag)

    return line

def fix_line(line):

    # Turn off solver
    if 'report_solver_progress' in line:
        line = set_flag(line)

    # Limit output
    if 'terminal_interval' in line:
        line = set_min_value(line)
        
    # Limit output
    if 'profile_interval' in line:
        line = set_min_value(line)
    return line


def fix_inlist(file):
    try:
        with open(file,'r') as f:
            lines = f.readlines()
    except IsADirectoryError:
        return

    lines = [fix_line(line) for line in lines]

    with open(file,'w') as f:
        f.writelines(lines)

def fix_inlists(*files):
    for file in files:
        print(file)
        fix_inlist(file)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        fix_inlists(sys.argv[1:])
    else:
        inlists = glob.glob(os.path.join(MESA_DIR,'star','test_suite','*','inlist*'))
        fix_inlists(*inlists)



