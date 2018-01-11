import numpy as np
import os
import time
import random
import math
from tkinter import *

print ("Hello world!")

SaveTranscript = "EnregistrementTranscrits.txt"
Dir_start = "tousgenesidentiques/"
Dir_curr = "current_profile/"
nbgenes = 10

Transcript = []
Genome_pos = []
Genome_sign = []
Genome_names = []
Genome_fitness = 00000000
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
  change = 0
  for i in range(len(lignes)):
     if not Transcript_new[i] ==int(lignes[i].split(" ")[2].split("\n")[0]):
       change = 1
     Transcript_new[i] = int(lignes[i].split(" ")[2].split("\n")[0])
  return change
 
def load_genome(directoryAddress):
  global Genome_pos,Genome_sign,Prot_pos,Range,Genome_names,SaveTranscript
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  f = open(directoryAddress+'tousgenesidentiques.gff','r')
  lignes = f.readlines()
  f.close()
  for i in range(5,len(lignes)):
    l = lignes[i].split("\t")
    if len(l)>1:
      Genome_pos[i-5] = (int(l[3]),int(l[4]))
      Genome_sign[i-5] = l[6]
      Genome_names[i-5] = l[8].split("=")[2]
  ligneRange = lignes[4].split("\t")
  Range = (int(ligneRange[3]),int(ligneRange[4]))
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
  dostuff = 0
  #~ global Genome_pos_new,Genome_sign_new
  #~ for i in range(len(Genome_pos_new)):
    #~ if(Genome_pos_new[i][0]>Genome_pos_new[i][1]):
      #~ Genome_sign_new[i] = "-"
    #~ else:
      #~ Genome_sign_new[i] = "+"

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
  newtab.append(lignes[3].split(" ")[0]+" "+lignes[3].split(" ")[1]+" "+str(Range_new[0])+" "+str(Range_new[1])+"\n")
  newtab.append(lignes[4].split("\t")[0]+"\t"+lignes[4].split("\t")[1]+"\t"+lignes[4].split("\t")[2]+"\t"+str(Range_new[0])+"\t"+str(Range_new[1])+"\t.\t+\t.\tID=id0;Name=tousgenesidentiques\n")
  #~ print (Genome_pos_new)
  for i in range(len(Genome_pos_new)):
    s = "tousgenesidentiques\tRefSeq\tgene\t"
    newtab.append(s+str(Genome_pos_new[i][0])+"\t"+str(Genome_pos_new[i][1])+"\t.\t"+Genome_sign_new[i]+"\t.\tID=g1;Name="+Genome_names_new[i])
  newtab.append("\n")
  f = open(directoryAddress+"current_profile.gff","w")
  f.writelines(newtab)
  f.close()
	
def writeTSS_dat(directoryAddress):
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  newtab = []
  newtab.append("TUindex\tTUorient\tTSS_pos\tTSS_strength\n")
  for i in range(len(Genome_pos_new)):
    newtab.append(str(i)+"\t"+Genome_sign_new[i]+"\t"+str(Genome_pos_new[i][0])+"\t.2\n") # Note can't change TSS_strengh with this command yet !
  f = open(directoryAddress+"TSS.dat","w")
  f.writelines(newtab)
  f.close()
	
def writeTTS_dat(directoryAddress):
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  newtab = []
  newtab.append("TUindex\tTUorient\tTTS_pos\tTTS_proba_off\n")
  for i in range(len(Genome_pos_new)):
    newtab.append(str(i)+"\t"+Genome_sign_new[i]+"\t"+str(Genome_pos_new[i][1])+"\t1.\n") # Note can't change TTS_proba_off with this command yet !
  f = open(directoryAddress+"TTS.dat","w")
  f.writelines(newtab)
  f.close()
	
def colorPalette():
	
	col=[]
	col.append(['slategray','#6f8090',10])
	col.append(['cornsilk3','#cdc7b0',10])
	#~ col.append(['burlywood2','#edc591',10])
	col.append(['lightpink','#ffadb8',10])
	#~ col.append(['palegreen','#97fa97',10])
	col.append(['paleturquoise','#baffff',10])
	#~ col.append(['thistle2','#cac3ca',10])
	col.append(['slateblue3','#6958cd',10])
	col.append(['steelblue3','#4f93cd',10])
	col.append(['tomato','#ff6246',10])
	col.append(['orangered1','#ff4500',10])
	col.append(['orange','#ffa400',10])
	col.append(['firebrick3','#cd2525',10])
	#~ col.append(['brown2','#ed3a3a',10])
	#~ col.append(['yellow1','#ffff00',10])
	col.append(['gold2','#edc800',10])
	col.append(['deeppink2','#ed1288',10])
	col.append(['mediumpurple','#926fdb',10])
	col.append(['chartreuse2','#7eff00',10])
	col.append(['green2','#00cd00',10])
	#~ col.append(['olivedrab3','#9acd31',10])
	#~ col.append(['limegreen','#31cd31',10])
	col.append(['royalblue2','#436ded',10])
	col.append(['seagreen3','#43cd80',10])
	col.append(['skyblue3','#6ca6cd',10])
	col.append(['darkmagenta','#8a008a',10])
	col.append(['hotpink1','#ff6db4',10])
	
	return col

def Metropolis (nb_iterations, nb_it_fitness, goal_profile,percent_zyup) :
  global Genome_pos,Genome_sign,Prot_pos,Range,Genome_names,SaveTranscript
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  global Genome_fitness,Transcript_new,Transcript
  # il faut une grosse taille simulation pour estimer la fitness => peut dependre architecture du genome
  # simulation avec que des invertions => genome constant au debut 
  # poid relatif d'insertion deletion plutot que 
  #~ A = 0 # Note that this line is useless
  fitness = 0.5
  newtab = []
  i = 0 
  window = Tk()
  window.title('Metropolis')
  canvas = Canvas(window,bg='#353535',height=300,width=1200)
  canvas.pack(side=TOP)
  col = colorPalette()
  
  
  while i < nb_iterations :
    Genome_pos_new = list(Genome_pos)
    Genome_sign_new = list(Genome_sign)
    
    fitness_prev = Genome_fitness
    fitness_new = []
    change = 0
    #ecrire les fichiers 
    
    if random.random() <= percent_zyup :
      z = 0
      while z == 0:
        z = zyup()
    else : zyop()
		  
    changeGenomeDir(Dir_curr)
    
    for j in range(nb_it_fitness) :
      os.system('python3 start_simulation.py parnew.ini out.ouput')
      
      #recuperer transcrits : Transcript_new
      change = change + load_res()
      fitness_new.append(calc_fitness(Transcript_new, goal_profile))
      print (Transcript_new)
    
    fitness_new = np.mean(fitness_new)
    print(fitness_new)
    s = " "
    if change > 0:
      i+=1
      canvas.delete("all")
      canvas.create_rectangle(100,50,1100,250,fill="#000000")
      canvas.create_text(30,25,text=str(i),fill="#FFFFFF",anchor=NW)
      canvas.create_text(30,75,text=str(Genome_fitness),fill="#FFFFFF",anchor=NW)
      canvas.create_text(30,175,text=str(fitness_new),fill="#FFFFFF",anchor=NW)
      for k in range(len(Genome_pos)):
        canvas.create_rectangle(100+Genome_pos[k][0]/30,60,100+Genome_pos[k][1]/30,80,fill=col[k][1])
      for k in range(len(Prot_pos)):
        canvas.create_rectangle(100+Prot_pos[k]/30,90,110+Prot_pos[k]/30,110,fill=col[k][1])
      for k in range(len(Genome_pos)):
        canvas.create_rectangle(100+Genome_pos_new[k][0]/30,160,100+Genome_pos_new[k][1]/30,180,fill=col[k][1])
      for k in range(len(Prot_pos)):
        canvas.create_rectangle(100+Prot_pos_new[k]/30,190,110+Prot_pos_new[k]/30,210,fill=col[k][1])
      window.update()
      #~ print(" ")
      #~ print("####################################")
      #~ print("     Iterations number "+str(i))
      #~ print("####################################") 
      if fitness_new > Genome_fitness :
        Genome_fitness = fitness_new
        Genome_pos = list(Genome_pos_new)
        Genome_sign = list(Genome_sign_new)
        Genome_names = list(Genome_names_new)
        Prot_pos = list(Prot_pos_new)
        Transcript = list(Transcript_new)
        Range = list(Range_new)
      
      
      else :
        p_fitness = (math.exp(math.exp(math.exp(fitness_new/Genome_fitness)))-1)/(math.exp(math.exp(math.exp(1)))-1)  # proba aceptation, depend de la fitness   #### changes here
        #~ p_fitness = -0.1
        if random.random() < p_fitness : 
          Genome_fitness = fitness_new
          Genome_pos = list(Genome_pos_new)
          Genome_sign = list(Genome_sign_new)
          Genome_names = list(Genome_names_new)
          Prot_pos = list(Prot_pos_new)
          Transcript = list(Transcript_new)
          Range = list(Range_new)
        else :		
          Genome_pos_new = list(Genome_pos)
          Genome_sign_new = list(Genome_sign)
          Genome_names_new = list(Genome_names)
          Prot_pos_new = list(Prot_pos)
          Transcript_new = list(Transcript)
          Range_new = list(Range)
      #~ else : 
        #~ #p_fitness = (math.exp(fitness_new)-1)/(math.exp(1)-1)  # proba aceptation, depend de la fitness
        #~ p_fitness = 2
        #~ if random.random() > p_fitness : 
          #~ Genome_fitness = fitness_new
          #~ Genome_pos = list(Genome_pos_new)
          #~ Genome_sign = list(Genome_sign_new)
          #~ Genome_names = list(Genome_names_new)
          #~ Prot_pos = list(Prot_pos_new)
          #~ Transcript = list(Transcript_new)
          #~ Range = list(Range_new)
      #enregistrer genome final
      changeGenomeDir(Dir_curr)
      s += str(i)+"|"+str(Genome_fitness)+"|"+str(fitness_new)# +"|"+str(change)
      #~ print (Transcript_new)
      s += "|  "
      for k in range(len(Transcript_new)):
        s+= "|"+str(Transcript_new[k])
      s += "|  "
      for k in range(len(Prot_pos)):
        s+= "|"+str(Prot_pos[k])
      s += "|  "
      for k in range(len(Genome_pos)):
        s+="|"+str(Genome_pos[k][0])+";"+str(Genome_pos[k][1])
      newtab.append(s+'\n')
      f = open("SimulationResults.txt","w")
      f.writelines(newtab)
      f.close()



def calc_fitness(Transcript, goal_profile) :
  fitness = 1   #### changes here
  for i in range(len(Transcript)) : 
    fitness = fitness + (goal_profile[i] - Transcript[i])*(goal_profile[i] - Transcript[i])
  return 1/fitness   #### changes here
  

def pic_interg_pos():
  global Genome_pos
  pos_chosen = False
  pos = 0
  while pos_chosen == False:
    isok = True
    pos = np.random.randint(Range[0],Range[1]+1)
    for k in range(len(Genome_pos)):
      if pos > min(Genome_pos[k]):
        if pos < max(Genome_pos[k]):
          isok = False
    if isok == True:
      pos_chosen = True
    print (pos,isok)
  return pos
	  
  #~ global Genome_pos
  #~ pos_chosen = []
  #~ # creating intergenic positions list
  #~ interg_pos = []
  #~ for pos in range(Range[0],Range[1]+1) : 
    #~ is_interg = "TRUE"
    #~ for g in Genome_pos : 
      #~ if pos in range(g[0],g[1]+1) :
        #~ is_interg = "FALSE"
        #~ break
    #~ if is_interg=="TRUE" : 
      #~ interg_pos.append(pos)
	
  #~ # picking one position or two from it
  #~ for i in range(nb_of_pos) :
    #~ pos_chosen.append(random.choice(interg_pos))
  #~ return pos_chosen


def zyup() : 
  global Genome_pos,Genome_sign,Prot_pos,Range,Genome_names,SaveTranscript
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  global Genome_fitness,Transcript_new
  counter = 0
  (pos_1,pos_2) = (pic_interg_pos(),pic_interg_pos())
  for i in range(len(Genome_pos)) :
    if Genome_pos[i][0]>min(pos_1,pos_2):
      if Genome_pos[i][0]<max(pos_1,pos_2):
        if Genome_pos[i][1]>min(pos_1,pos_2):
          if Genome_pos[i][1]<max(pos_1,pos_2):
            counter +=1
            Genome_pos_new[i] = (pos_1+pos_2-Genome_pos[i][0],pos_1+pos_2-Genome_pos[i][1])
            if (Genome_sign[i] == "+") :
              Genome_sign_new[i] = "-"
            else : 
              Genome_sign_new[i] = "+"
  #~ for i,a in enumerate(Genome_pos) :
    #~ if (Genome_pos[i][0] >= pos_1 and Genome_pos[i][0] <= pos_2 and Genome_pos[i][1] >= pos_1 and Genome_pos[i][1] <= pos_2 ) :
      #~ if (Genome_sign[i] == "+") :
        #~ Genome_pos_new[i] = (pos_2 - (min(Genome_pos[i]) - pos_1), pos_1 + pos_2 - max(Genome_pos[i]))
        #~ Genome_sign_new[i] = "-"
      #~ else : 
        #~ Genome_pos_new[i] = ( pos_1 + pos_2 - max(Genome_pos[i]),pos_2 - (min(Genome_pos[i]) - pos_1))
        #~ Genome_sign_new[i] = "+"
  for i in range(len(Prot_pos)):
    if Prot_pos[i] > min(pos_1,pos_2):
      if Prot_pos[i] < max(pos_1,pos_2):
        counter +=1
        Prot_pos_new[i] =  pos_1 + pos_2 - Prot_pos[i]
  return counter

def zyop() :
  global Genome_pos,Genome_sign,Prot_pos,Range,Genome_names,SaveTranscript
  global Genome_pos_new,Genome_sign_new,Prot_pos_new,Range_new,Genome_names_new
  global Genome_fitness,Transcript_new 
  pos = pic_interg_pos()
  sign = random.choice(("+", "-"))
  for i,a in enumerate(Genome_pos) :
    if (min(Genome_pos[i]) > pos) :
      if (sign == "+") :
        Range_new[1] += 1 
        Genome_pos_new[i] = (Genome_pos[i][0]+1, Genome_pos[i][1]+1)
      else : 
        Range_new[1] -= 1 
        Genome_pos_new[i] = (Genome_pos[i][0]-1, Genome_pos[i][1]-1)
  
load_res()
load_genome(Dir_start)
changeGenomeDir(Dir_curr)



print ("#############################################")
print (Transcript,"Transcript")
print (Genome_pos,"Genome_pos")
print (Genome_sign,"Genome_sign")
print (Genome_names,"Genome_names")
print (Prot_pos,"Prot_pos")
print (Range,"Range")
print ("#############################################")
