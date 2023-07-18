######FLFB(Feature-based LLPS-Fibrils-Background classification model)########
import localcider
from localcider.sequenceParameters import SequenceParameters
from localcider import plots
import localcider
from localcider.sequenceParameters import SequenceParameters
from localcider import plots
import pandas as pd

import subprocess
from subprocess import PIPE
from pandas.core.frame import DataFrame
import argparse
import numpy as np
import joblib
from pathlib import Path
import os
import iupred3_lib
import sklearn
from tqdm import tqdm

def get_SCD_SHD_v(seqs):   ###.csv format###
    SHDs = []
    
    for i,seq in zip(tqdm(range(len(seqs)), desc='Caculating the SHD'),seqs):
        seq1 = []
        N = len(seq)
        
        for aa in seq:      
            if aa == 'A':     
                seq1.append(0.730)
            elif aa == 'R':
                seq1.append(0)
            elif aa == 'N':
                seq1.append(0.432)
            elif aa == 'D':
                seq1.append(0.378)
            elif aa == 'C':
                seq1.append(0.595)
            elif aa == 'Q':
                seq1.append(0.514)
            elif aa == 'E':
                seq1.append(0.459)
            elif aa == 'G':
                seq1.append(0.649)
            elif aa == 'H':
                seq1.append(0.514)
            elif aa == 'I':
                seq1.append(0.973)
            elif aa == 'L':
                seq1.append(0.973)
            elif aa == 'K':
                seq1.append(0.514)
            elif aa == 'M':
                seq1.append(0.838)
            elif aa == 'F' or aa == 'P':
                seq1.append(1.000)
            elif aa == 'S':
                seq1.append(0.595)
            elif aa == 'T':
                seq1.append(0.676)
            elif aa == 'W':
                seq1.append(0.946)
            elif aa == 'Y':
                seq1.append(0.865)
            elif aa == 'V':                            
                seq1.append(0.892)
                
        m = 1  
        n = 0   
        SHD = 0    
        while m < N:  
            while n < m:   
                S1 = (seq1[m] + seq1[n]) * (m - n) ** -1  # calculation of SHD
                n += 1 
                SHD = SHD + S1 
            m += 1 
            n = 0
        SHD = SHD/N
    
        SHDs.append(SHD)
          
    return SHDs

def get_localcider(seqs):
    N=[]
    W=[]
    L=[]
    C=[]
    AAs =[N,W,L,C]
    AAS =['N','W','L','C']
    for i,seq in zip(tqdm(range(len(seqs)), desc='Caculating N,W,L,C'),seqs):
        SeqOb = SequenceParameters(sequence=seq)
        AADIC = SeqOb.get_amino_acid_fractions()
        
        for i in range(len(AAs)):
            AAs[i].append(AADIC[AAS[i]])
    return N,W,L,C

def get_IDR(seq):

    iupred2_result = iupred3_lib.iupred(seq,'long', smoothing='medium')

    scores = []
    for pos, residue in enumerate(seq):
        # print('{}\t{}\t{:.4f}'.format(pos + 1, residue, iupred2_result[0][pos]), end="")
        scores.append(iupred2_result[0][pos])

    score_region = 0
    for i in range(len(scores)):
        if scores[i] >= 0.5:
            score_region = score_region + 1
    score_region = score_region / len(scores)

    return score_region

def get_IDRs(seqs):
    print('Calculating the IDR fraction')
    IDRs= []
    for seq in seqs :
        IDR = get_IDR(seq)
        IDRs.append(IDR)
    return IDRs

def fasta2seqs(filepath):
    """Function to read a fasta file into a {tag:seq} dictionary."""
    f = open(filepath, 'r')
    seqs = []
    names = []
    regions = []
    for line in f.readlines():
        if line[0]  == '>':
            seqs.append(''.join(regions))
            regions = []
            names.append(line.strip())
        elif line[0] != '>':
            regions.append(line.strip())
    seqs.append(''.join(regions))
    return names, seqs[1:]

def get_result(features):
    current_work_dir = os.path.dirname(__file__)
    model60 = joblib.load(os.path.join(current_work_dir,'model','6model0.pkl'))
    model61 = joblib.load(os.path.join(current_work_dir,'model','6model1.pkl'))
    model62 = joblib.load(os.path.join(current_work_dir,'model','6model2.pkl'))
    model63 = joblib.load(os.path.join(current_work_dir,'model','6model3.pkl'))
    model64 = joblib.load(os.path.join(current_work_dir,'model','6model4.pkl'))
    
    model_6 = [model60,model61,model62,model63,model64]
    
    num_seqs , num_features = features.shape
    results = np.zeros(3*num_seqs).reshape(num_seqs,3)
    
    for i in range(5):
        clf = model_6[i]
        result = clf.predict_proba(features)
        results = results + np.array(result)

    class_predicted = []
    results_ave = results/5
    for i in range(num_seqs):
        probablity = results[i]
        score = np.argmax(probablity)
        if score == 0:
            class_predicted.append('Background')
        if score == 1:
            class_predicted.append('LLPS')
        if score == 2:
            class_predicted.append('Fibrils')
            
    return class_predicted,results_ave

def df2file(df,file_name):
    path = ''.join([file_name,'.csv'])
    path = Path(path)
    df.T.to_csv(path)

def run_seq(args):
    protein_names,seqs = fasta2seqs(args.seqs_fasta)
    SHD = get_SCD_SHD_v(seqs)
    IDR = get_IDRs(seqs)
    N,W,L,C = get_localcider(seqs)

    features = np.array([C,L,N,SHD,W,IDR]).T
    predicted_class,results_proba = get_result(features)
    df_result1 = DataFrame(data = [protein_names,seqs,predicted_class],index = ['protein_name','seqs','Class'] )
    df_result2 = DataFrame(data = results_proba,columns = ['P_background','P_llps','P_fibrils']).T
    df_result = pd.concat([df_result1,df_result2],axis = 0)
    df2file(df_result,args.output_filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="run LLPS-Fibrils-Background model on fasta files")
    parser.add_argument('--input file name','-i', dest = 'seqs_fasta',help = ('input fasta file name'), type = str )
    parser.add_argument('--output filename','-o', dest = 'output_filename', help = ('output csv file name'), type = str)
    args = parser.parse_args()
    run_seq(args)
    print('OK')