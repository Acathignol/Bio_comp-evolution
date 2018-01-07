#!/usr/bin/python

import numpy as np
import os
import time
import random
import sys 
import TestMain

# main script 

print("Hello World!")

nb_iterations = 1 # sys.argv[1]

global goal_profile
goal_profile = [25,250,400,22,300,44,500,230,145,957] # sys.argv[2] 

global percent_zyup
percent_zyup = 1 # sys.argv[3] # 1 = 100% = only indel

global Genome_fitness
Genome_fitness = 0
# TestMain.calc_fitness(Transcrit) #problem : il est en global dans l'autre, antoine ? help !

TestMain.Metropolis(nb_iterations)

print(Genome_pos)



