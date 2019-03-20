""" Adapted from code by Adam Weiner: https://github.com/adamcweiner/flu-vax/blob/master/fluv/sequence_processing.py
"""

from Bio import SeqIO
import os
import numpy as np
import random

def trim(seqFile):
    cwd = os.getcwd()
    homeDir = cwd[1:5]
    if (homeDir == 'home'):
        print('using path from server to load sequences')
        pathToFile = os.path.join("/home","zoekim","Desktop",str(seqFile)) #aretha server
    else:
        print('using path from mac machine to load sequences')
        pathToFile = os.path.join("/Users","zoekim","Desktop",str(seqFile))  #mac machine
    
    allSeqs = []
    allLabels = []
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
            allSeqs.append(seq_record.seq)
            allLabels.append(seq_record.id)
    
    seqMat = np.array(allSeqs)
    label = np.array(allLabels)
    
    sequence = seqMat[:, 0:317]
    
    # filtering out residues not included in PAM250 pymsa distance matrix (http://www.matrixscience.com/blog/non-standard-amino-acid-residues.html)
    for i in range(0, sequence.shape[0]):
        for j in range(0,sequence.shape[1]):
            if (sequence[i,j] == 'J'):
                sequence[i,j] = random.choice(['I', 'L'])
    
    return (label, sequence)
