# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 11:08:25 2016

@author: cra
"""

#from __future__ import print_function

##--- import explicitly to help rope
#import numpy as np
#import fenics_helpers as fh
#import particle_helpers as ph
#import sys
#
#execfile('/home/cra/dolfin/useful_f/init_cra.py')
#
#
#import time
#
#
#Np = 10
#
#cnt = 0
#print str(cnt),
#while True:
#    cnt += 1
#    time.sleep(0.1)
#    sys.stdout.write("\rDoing thing %i" % cnt)
#    sys.stdout.flush()


import time
import curses

stdscr = curses.initscr()

#stdscr.addstr(0, 0, "Hello")
#stdscr.refresh()
#
#time.sleep(1)
#
#stdscr.addstr(0, 0, "World! (with curses)")
#stdscr.refresh()
cnt=1
while True:
    cnt += 1
    time.sleep(0.1)
    stdscr.addstr(0, 0, "p: %d"%(cnt))
    stdscr.refresh()