# Bio_comp-evolution
Aim : To develop a simulation code to test the adaptation of transcriptional regulation in the absence of mutations in gene sequences and promotors

- We have a bacterial DNA with a set of genes and intergenic regions. At the end of a simulation with the code of old 5BIM, we obtain the level of genes expression (how much a gene is transcript). If we simulate just like that, we find an even expression for all the genes.
- We want to see a differiential expression that we will choose before. (so individuals must evolve to it).

FIRST PART OF THE CODE :
- zyup : function(intergenic_pos1, intergenic_pos2)
- zyop : function(intergenic_pos, sign) # + = insertion and - = deletion

SECOND PART OF THE CODE :
- metropolis algorithm : We have one individual per generation.
At each time we :
   1) take the individual
   2) measure its fitness (difference from the expected genes expression)
   3) do a random change to the individual (zyup or zyop) with a probability
   4) measure the new fitness of the individual
   5) if the new-individual's fitness is better : keep the new ???
      else : keep the new with a probability given. (prob decreases with the fitness).
