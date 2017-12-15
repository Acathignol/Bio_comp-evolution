#!/usr/bin/python

import random
import sys 

# main script 

print("Hello World!")

nb_iterations = 1 # sys.argv[1]

global goal_profile
goal_profile = [25,250,400,22,300,44,500,230,145,957] # sys.argv[2] 

global percent_zyup
percent_zyup = 1 # sys.argv[3] # 1 = 100%

Metropolis(nb_iterations)



