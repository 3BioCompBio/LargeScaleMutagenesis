from config import lst_pth,aas_dir,nas_dir,ali_dir,clustalw2

#Extract & Translate AA seq
inc_Codon=0
x_Codon=0
def na_to_aa(fasta_file):
    global inc_Codon
    global x_Codon
    code={'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'TGT':'C','TGC':'C',
            'GAT':'D','GAC':'D',
            'GAA':'E','GAG':'E',
            'TTT':'F','TTC':'F',
            'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
            'CAT':'H','CAC':'H',
            'ATT':'I','ATC':'I','ATA':'I',
            'AAA':'K','AAG':'K',
            'TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
            'ATG':'M',
            'AAT':'N','AAC':'N',
            'TAG':'O',
            'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'CAA':'Q','CAG':'Q',
            'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
            'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S',
            'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
            'TGA':'U',
            'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
            'TGG':'W',
            'TAT':'Y','TAC':'Y'}
    #Extract sequence
    ffile=open(fasta_file,'r')
    na_fasta=ffile.read()
    ffile.close()
    na_seq=''
    aa_seq=''
    for lines in na_fasta.split('\n'):
        if '>' in lines or ';' in lines:
            header=lines+'\n'
        else:
            na_seq=na_seq+lines
    #Translate sequence
    na_seq=[na_seq[i:i+3] for i in range(0, len(na_seq), 3)]
    for codons in na_seq:
        if len(codons)!=3:
            print(fasta_file+' : incomplete codon')
            inc_Codon+=1
            continue
        if codons not in code:
            print(fasta_file+' : '+codons+' codon found')
            x_Codon+=1
            aa_seq=aa_seq+'X'
        else:
            aa_seq=aa_seq+code[codons]
    #Format and return        
    aa_seq='\n'.join([aa_seq[i:i+80] for i in range(0, len(aa_seq), 80)])
    return (str(header+aa_seq+'\n'))

import os
os.chdir(nas_dir)

#CLUSTALW2 alignment
from os import system
def get_alignment(fastafile):
    alignfile=fastafile.strip('_tr+pdb.txt')+'_align.fa'
    os.system(clustalw2+' -QUIET -OUTPUT=FASTA -OUTORDER=INPUT -TYPE=PROTEIN -GAPOPEN=10 -GAPEXT=1 -INFILE='+fastafile+' -OUTFILE='+alignfile+' &> /dev/null')

from os.path import exists
import glob
def align_list(pdblist):
    nopdbfasta=0
    for x in pdblist:
        pdbid=x[0:4]
        chain=x[5]

        ###GET PDB FASTA
        if not exists(aas_dir+pdbid+chain+'.txt'):
            print(pdbid+chain+' : aa_seq file not found')
            nopdbfasta+=1
        else:
            fpdb=open(aas_dir+pdbid+chain+'.txt','r')
            AAseq=fpdb.read()

        ###MATCH EMBL FASTA
        embflist=glob.glob(pdbid+chain+'*.txt')
        for files in embflist:
            trNAseq=na_to_aa(files)
            trFILE_path=ali_dir+files.strip('.txt')+'_tr+pdb.txt'
            trFILE=open(trFILE_path,'w')
            trFILE.write(trNAseq+AAseq)
            trFILE.close()
            if not exists(trFILE_path):
                print(files.strip('.txt')+'_tr+pdb.txt'+' : not created')
            get_alignment(trFILE_path)
            
        embflist=None
    print(str(inc_Codon)+' incomplete codons, '+str(x_Codon)+' unknown codons, '+str(nopdbfasta)+' missing pdb fasta')

listfile=open(lst_pth, 'r')
pdblist=listfile.readlines()
align_list(pdblist)
listfile.close()