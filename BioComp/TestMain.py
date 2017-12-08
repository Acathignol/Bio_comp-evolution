import numpy as np
import os
import time


print ("Hello world!")

SaveTranscript = "EnregistrementTranscrits.txt"

nbgenes = 10

Transcript = []
Genome_pos = []
Genome_sign = []
Range = (0,0)
for i in range(nbgenes):
  Transcript.append(0)
  Genome_pos.append((0,0))
  Genome_sign.append("+")
  
def load_res(SaveTranscript):
  global Transcript
  f = open(SaveTranscript,"r")
  lignes = f.readlines()
  f.close()
  for i in range(len(lignes)):
     Transcript[i] = int(lignes[i].split(" ")[2].split("\n")[0])
  

def load_genome(directoryAddress):
  global Genome_pos,Genome_sign,Range
  f = open(directoryAddress+'tousgenesidentiques.gff','r')
  lignes = f.readlines()
  f.close()
  for i in range(5,len(lignes)):
    l = lignes[i].split("\t")
    Genome_pos[i-5] = (int(l[3]),int(l[4]))
    Genome_sign[i-5] = l[6]
  ligneRange = lignes[4].split("\t")
  Range = (ligneRange[3],ligneRange[4])


load_res(SaveTranscript)
load_genome("tousgenesidentiques/")

print ("#############################################")
print (Transcript)
print (Genome_pos)
print (Genome_sign)
print (Range)
print ("#############################################")

def changeGenomeDir():
	dostuff = 0
	

#~ def Metropolis (params)
  
  #~ A = 0 # Note that this line is useless
  



#~ def zyup(pos_1,pos_2)

