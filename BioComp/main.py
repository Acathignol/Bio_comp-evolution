#!/usr/bin/python

import numpy as np
import os
import time
import random
import sys 
import TestMain


# main script 

print("Hello World!")

nb_iterations = 1000 # sys.argv[1]
nb_it_fitness = 10 

goal_profile =  [200, 250, 300, 300, 250, 200, 250, 300, 200, 250] # sys.argv[2] 

#global percent_zyup
percent_zyup = 0.5 # sys.argv[3] # 1 = 100% = only indel
# 1, 0.75, 0.5 (et 0)

#~ global Genome_fitness
#~ Genome_fitness = 0
# TestMain.calc_fitness(Transcrit) #problem : il est en global dans l'autre, antoine ? help !

TestMain.Metropolis(nb_iterations, nb_it_fitness, goal_profile, percent_zyup)

