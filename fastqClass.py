# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import pandas as pd
import gzip

# Quality value dictionary
qualityDict = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, '\'': 6, '(': 7, ')': 8, '*': 9, '+': 10,
               ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
               '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
               '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40,
               'J': 41, 'K': 42}

# Class to read fastq data files and place into a dataframe
class DataFile:
    
    
    def __init__(self, fastq):
       
        self.fastq = fastq
        
        # Takes the full file name cnad  shortens it to a readable name
        self.sample = '_'.join(fastq.split("_")[:2])
        
        # Splice R2 in place of R1 in fastq file name.
        fastq2 = fastq.split("_R1_")[0] + "_R2_" + fastq.split("_R1_")[1]
        
        # Read fastq with gzip library.
        with gzip.open(fastq, 'rt') as fastqFile:
           
            fqlines = fastqFile.readlines()
           
        # If there is a reverse file, read it and add it to fqlines
        if fastq2 in os.listdir():
        
            with gzip.open(fastq2, 'rt') as fastq2File:
                
                fq2lines = fastq2File.readlines()

            fileLinesList = []
        
            for i in range(0, len(fqlines), 4):
                fileLinesList.append(fqlines[i])
                fileLinesList.append(fqlines[i+1])
                fileLinesList.append(fqlines[i+2])
                fileLinesList.append(fqlines[i+3])
        
                fileLinesList.append(fq2lines[i])
                fileLinesList.append(fq2lines[i+1])
                fileLinesList.append(fq2lines[i+2])
                fileLinesList.append(fq2lines[i+3])
                
        # Create lists to hold FASTQ sequences and quality strings from fqlines
        fqseqList = []
        fqqualList = []
        
        # Starting at index 1, add every fourth line (nucleotide sequence) to fqseqList
        for i in range(1, len(fileLinesList), 4):
            fqseqList.append(fileLinesList[i])
        
        # Starting at index 3, add every fourth line (quality sequence) to fqqualList
        for j in range(3, len(fileLinesList), 4):
            fqqualList.append(fileLinesList[j])
        
        # Create a data frame "fqfiledf" from fqseqList and fqqualList
        self.fqfiledf = pd.DataFrame({'Seq': fqseqList, 'Qual': fqqualList})
        
        
    def HQ(self):
         
        # Create an empty column to track whether the SNP base meets quality threshold.
        self.fqfiledf['HQ'] = ''
    
    
    def qualAverage(self):
        
        
        quallist = self.fqfiledf['Qual']
        qualscores = []
        
        # Create a list for the quality strings and append the corresponding quality scores
        for item in quallist:
            qualstringlist = []
            
            for symbol in item.strip():
                qualstringlist.append(qualityDict[symbol])
            
            # Create a variable for the overall average quality of each line
            avequal = sum(qualstringlist)/float(len(qualstringlist))    
            qualscores.append(avequal)
        
        # Adds the quality scores in a coloumn to the dataframe
        self.fqfiledf['Avg Qual'] = qualscores
        
"""
    def seqStats(self):    
        
        # Store number of total sequence reads in file
        numSeqsInt = len(fqseqList)
"""     
        
        
# Class reads fasta files and creates dataframe
class FastAReader:
    
    def __init__(self, fasta):
        
        self.fasta = fasta
        
        # Read fastq with gzip library.
        with open(fasta, 'r') as fastaFile:
            
            falines = fastaFile.readlines()

        # Create lists to hold FASTQ sequences and quality strings from fqlines
        fanameList = []
        faseqList = []
        
        # Starting at index 1, add every other line to faseqList
        for i in range(0, len(falines), 2):
            fanameList.append(falines[i])

        for j in range(1, len(falines), 2):
            faseqList.append(falines[j])

        self.fafiledf = pd.DataFrame({'Name': fanameList, 'Seq': faseqList})
       
        
    def positionSNP(self):
        
        
        positionList = []
        
        # Initialize a counter called "i" with a value of 0
        i = 0
        
        faseqList = self.fafiledf['Seq']
        fanameList = self.fafiledf['Name']
               
        # Loop through the length of namesList while the counter is within the index of namesList
        while i < len(fanameList)-1:
            # If adjacent rows in namesList have the same beginning (and are thus different
            # versions/SNPs of the same sequence), set sequence variable seq1 and seq2 to
            # be the corresponding sequences. These should differ by 1 base (the SNP).
            if fanameList[i].split("_")[0] == fanameList[i+1].split("_")[0]:
                seq1 = faseqList[i]
                seq2 = faseqList[i+1]
                
                # Find the position of the SNP and add it to positionList at the index
                # of EACH sequence. If the SNP is at position 5, set BOTH indices of positionList
                # corresponding to the sequences to 5.
                self.position = [i for i in range(len(seq1)) if seq1[i] != seq2[i]][0]
                positionList[i] = self.position
                positionList[i+1] = self.position
                
            # Increment the counter to move forward through namesList
            i += 1
            
        # Create a data frame called "searchdf" from the name, sequence, and position lists, with an
        # empty column for counts that will be filled in later.
        print(len(positionList))
        if len(positionList) > 0:
            self.fafiledf["Position"] = positionList
        
    
    
    
    
    
    
    
    
    
        
        
        