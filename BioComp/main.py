#!/usr/bin/python

import random
import sys 

# main script 

print("Hello World!")

nb_iterations = sys.argv[1]

global goal_profile
goal_profile = sys.argv[2]

global percent_zyup
percent_zyup = sys.argv[3]

Metropolis(nb_iterations)



