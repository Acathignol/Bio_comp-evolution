import numpy as np
import os
import time
import random

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



def Metropolis (params){
  
  # il faut une grosse taille simulation pour estimer la fitness => peut dépendre architecture du génome
  # simulation avec que des invertions => génome constant au début 
  # poid relatif d'insertion délétion plutot que 
  #~ A = 0 # Note that this line is useless
  fitness = 0.5
  
  zyup()
  
}

def calc_fitness(Transcript) :
  goal_profile = [25,250,400,22,300,44,500,230,145,957]
  fitness = 0
  for i in xrange(len(Transcript)) : 
    fitness = fitness + (goal_profile[i] - Transcript[i])*(goal_profile[i] - Transcript[i])
  return fitness
  

def pic_interg_pos(nb_of_pos){
  pos_chosen = []
  # creating intergenic positions list
  interg_pos = []
  for pos in range(Range[0],Range[1]+1) : 
    is_interg = "TRUE"
    for g in Genome : 
      if pos in range(g[0],g[1]+1) :
        is_interg = "FALSE"
        break
    if is_interg=="TRUE" : 
      interg_pos.append(pos)
	
  # picking one position or two from it
  for i in range(nb_of_pos) :
    pos_chosen.append(random.choice(interg_pos))
  return pos_chosen
}

def zyup(Genome_pos_new, Genome_sign_new) : 
  (pos_1,pos_2) = pic_interg_pos(2)
  for i,a in enumerate(Genome_pos_new) :
    if (Genome_pos_new[i][0] >= pos_1 and Genome_pos_new[i][0] <= pos_2 and Genome_pos_new[i][1] >= pos_1 and Genome_pos_new[i][1] <= pos_2 ) :
      if (Genome_sign_new[i] == "+") :
        Genome_pos_new[i] = (pos_2 - (min(Genome_pos_new[i]) - pos_1), pos_1 + pos_2 - max(Genome_pos_new[i]))
      	Genome_sign_new[i] = "-"
      else : 
        Genome_pos_new[i] = ( pos_1 + pos_2 - max(Genome_pos_new[i]),pos_2 - (min(Genome_pos_new[i]) - pos_1))
	Genome_sign_new[i] = "+"


def zyop(Genome_pos_new) : 
  pos = pic_interg_pos(1)
  sign = random.choice(("+", "-"))
  for i,a in enumerate(Genome_pos_new) :
    if (min(Genome_pos_new[i]) > pos) :
      if (sign == "+") :
        Genome_pos_new[i] = (Genome_pos_new[i][0]+1, Genome_pos_new[i][1]+1)
      else : 
        Genome_pos_new[i] = (Genome_pos_new[i][0]-1, Genome_pos_new[i][1]-1)
  
