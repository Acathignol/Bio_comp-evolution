import numpy as np
import os
import time
import random
import math

print ("Hello world!")

SaveTranscript = "EnregistrementTranscrits.txt"
Dir_start = "tousgenesidentiques/"
Dir_curr = "current_profile/"
nbgenes = 10

Transcript = []
Genome_pos = []
Genome_sign = []
Genome_names = []
Prot_pos = []
Range = (0,0)

for i in range(nbgenes):
  Transcript.append(0)
  Genome_pos.append((0,0))
  Genome_sign.append("+")
  Genome_names.append("test")
  Prot_pos.append(0)

Transcript_new = Transcript
Genome_pos_new = Genome_pos
Genome_sign_new = Genome_sign
Genome_names_new = Genome_names
Prot_pos_new = Prot_pos
Range_new = Range
  
def load_res():
  global Transcript_new,SaveTranscript
  f = open(SaveTranscript,"r")
  lignes = f.readlines()
  f.close()
  for i in range(len(lignes)):
     Transcript_new[i] = int(lignes[i].split(" ")[2].split("\n")[0])
 
def load_genome(directoryAddress):
  global Genome_pos,Genome_sign,Prot_pos,Range,Genome_names,SaveTranscript
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  f = open(directoryAddress+'tousgenesidentiques.gff','r')
  lignes = f.readlines()
  f.close()
  for i in range(5,len(lignes)):
    l = lignes[i].split("\t")
    Genome_pos[i-5] = (int(l[3]),int(l[4]))
    Genome_sign[i-5] = l[6]
    Genome_names[i-5] = l[8].split("=")[2]
  ligneRange = lignes[4].split("\t")
  Range = (ligneRange[3],ligneRange[4])
  f = open(directoryAddress+'prot.dat','r')
  lignes = f.readlines()
  f.close()
  for i in range(1,len(lignes)):
	  Prot_pos[i-1] = int(lignes[i].split("\t")[1])
  Transcript_new = Transcript
  Genome_pos_new = Genome_pos
  Genome_sign_new = Genome_sign
  Genome_names_new = Genome_names
  Prot_pos_new = Prot_pos
  Range_new = Range
  ## Run once the code with starting params
  os.system('python3 start_simulation.py params.ini out.ouput')
  f = open(SaveTranscript,"r")
  lignes = f.readlines()
  f.close()
  for i in range(len(lignes)):
     Transcript[i] = int(lignes[i].split(" ")[2].split("\n")[0])
	
	  
def changeGenomeDir(Dir_curr):
  reorder()
  writeProt_dat(Dir_curr)
  write_gff(Dir_curr)
  writeTTS_dat(Dir_curr)
  writeTSS_dat(Dir_curr)


def reorder():
  dostuff=0

def writeProt_dat(directoryAddress):
  global Prot_pos_new
  newtab = []
  newtab.append("prot_name\tprot_pos\n")
  for i in range(len(Genome_pos_new)):
  	newtab.append("hns"+"\t"+str(Prot_pos_new[i])+"\n")
  f = open(directoryAddress+'prot.dat','w')
  f.writelines(newtab)
  f.close()

def write_gff(directoryAddress):
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  f = open(directoryAddress+"current_profile.gff","r")
  lignes = f.readlines()
  f.close()
  newtab = lignes[0:3]
  newtab.append(lignes[3].split(" ")[0]+" "+lignes[3].split(" ")[1]+" "+Range_new[0]+" "+Range_new[1]+"\n")
  newtab.append(lignes[4].split("\t")[0]+"\t"+lignes[4].split("\t")[1]+"\t"+lignes[4].split("\t")[2]+"\t"+Range_new[0]+"\t"+Range_new[1]+"\t.\t+\t.\tID=id0;Name=tousgenesidentiques\n")
  for i in range(len(Genome_pos_new)):
    s = "tousgenesidentiques\tRefSeq\tgene\t"
  newtab.append(s+str(Genome_pos_new[i][0])+"\t"+str(Genome_pos_new[i][1])+"\t.\t"+Genome_sign_new[i]+"\t.\tID=g1;Name="+Genome_names_new[i])
  f = open(directoryAddress+"current_profile.gff","w")
  f.writelines(newtab)
  f.close()
	
def writeTSS_dat(directoryAddress):
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  newtab = []
  newtab.append("TUindex\tTUorient\tTSS_pos\tTSS_strength\n")
  for i in range(len(Genome_pos_new)):
    newtab.append(str(i)+"\t"+Genome_sign_new[i]+"\t"+str(min(Genome_pos_new[i]))+"\t.2\n") # Note can't change TSS_strengh with this command yet !
  f = open(directoryAddress+"TSS.dat","w")
  f.writelines(newtab)
  f.close()
	
def writeTTS_dat(directoryAddress):
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  newtab = []
  newtab.append("TUindex\tTUorient\tTTS_pos\tTTS_proba_off\n")
  for i in range(len(Genome_pos_new)):
    newtab.append(str(i)+"\t"+Genome_sign_new[i]+"\t"+str(max(Genome_pos_new[i]))+"\t1.\n") # Note can't change TTS_proba_off with this command yet !
  f = open(directoryAddress+"TTS.dat","w")
  f.writelines(newtab)
  f.close()
	


def Metropolis (nb_iterations) :
  
  # il faut une grosse taille simulation pour estimer la fitness => peut dependre architecture du genome
  # simulation avec que des invertions => genome constant au debut 
  # poid relatif d'insertion deletion plutot que 
  #~ A = 0 # Note that this line is useless
  fitness = 0.5
  
  for i in xrange(nb_iterations) :
    
    Genome_pos_new = list(Genome_pos)
    Genome_sign_new = list(Genome_sign)
    
    fitness_prev = Genome_fitness
    fitness_new = []
    #ecrire les fichiers 
    zyup(Genome_pos_new, Genome_sign_new)
    changeGenomeDir(Dir_curr)
    
    for j in xrange(10) :
      os.system('python3 start_simulation.py params_new.ini out.ouput')
      
      #recuperer transcrits : Transcript_new
      load_res(SaveTranscript)
      fitness_new.append(calc_fitness(Transcript_new))
    
    fitness_new = np.mean(fitness_new)
    
    if fitness_new > Genome_fitness :
      Genome_fitness = fitness_new
      Genome_pos = list(Genome_pos_new)
      Genome_sign = list(Genome_sign_new)
      Genome_names = list(Genome_names_new)
      Prot_pos = list(Prot_pos_new)
      Transcript = list(Transcript_new)
      Range = list(Range_new)
      
    else : 
     
      p_fitness = (math.exp(fitness_new)-1)/(math.exp(1)-1)  # proba aceptation, depend de la fitness
      if random.random() > p_fitness : 
        Genome_fitness = fitness_new
        Genome_pos = list(Genome_pos_new)
        Genome_sign = list(Genome_sign_new)
        Genome_names = list(Genome_names_new)
        Prot_pos = list(Prot_pos_new)
        Transcript = list(Transcript_new)
        Range = list(Range_new)
    #enregistrer genome final
    changeGenomeDir(Dir_curr)
  



def calc_fitness(Transcript) :
  goal_profile = [25,250,400,22,300,44,500,230,145,957]
  fitness = 0
  for i in xrange(len(Transcript)) : 
    fitness = fitness + (goal_profile[i] - Transcript[i])*(goal_profile[i] - Transcript[i])
  return fitness
  

def pic_interg_pos(nb_of_pos):
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
  
load_res()
load_genome(Dir_start)
changeGenomeDir(Dir_curr)

print ("#############################################")
print (Transcript)
print (Genome_pos)
print (Genome_sign)
print (Genome_names)
print (Prot_pos)
print (Range)
print ("#############################################")
