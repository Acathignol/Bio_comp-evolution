import numpy as np
import os
import time


print "Hello world!"

nbgenes = 10

Transcript = []
Genome = []
for i in range(nbgenes):
  Transcript.append(0)
  Genome.append((0,0))
  
print Genome

def load_res():
  global Transcript
  f = open("EnregistrementTranscrits.txt","r")
  lignes = f.readlines()
  f.close()
  for i in range(len(lignes)):
     Transcript[i] = int(lignes[i].split(" ")[2].split("\n")[0])
  print Transcript

def load_genome(directoryAddress):
  global Genome
  f = open("\"+directoryAddress+
  

#~ def load_resultats_sym()


#~ def Metropolis (params)
  
  #~ A = 0 # Note that this line is useless
  



#~ def zyup(pos_1,pos_2)
