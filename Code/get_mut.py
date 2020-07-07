from config import map_pth,nas_dir,ali_dir,mut_dir

def get_codons(alipath,embpath):
    #Loop to save the emb aaseq and the pdb aaseq of the alignment
    alifile=open(alipath,'r').read()
    header1=False #trNAseq
    header2=False #AAseq
    aaseq_emb=''
    aaseq_pdb=''
    for lines in alifile.split('\n'):
        if header1 and '>' in lines or ';' in lines:
            header2=True
            continue
        if '>' in lines or ';' in lines:
            header1=True
            continue
        if header1 and not header2:
            aaseq_emb=aaseq_emb+lines
        if header1 and header2:
            aaseq_pdb=aaseq_pdb+lines
    #Loop to save the emb naseq
    embfile=open(embpath,'r').read()
    naseq=''
    for lines in embfile.split('\n'):
        if '>' in lines or ';' in lines:
            continue
        else:
            naseq=naseq+lines
    #Attribute codons
    naseq=[naseq[i:i+3] for i in range(0, len(naseq), 3)]
    codons=[]
    j=0
    for i in range(len(aaseq_pdb)):
        if aaseq_pdb[i]==aaseq_emb[i]:
            codons.append(naseq[j])
            j+=1
        elif aaseq_pdb[i]=='-':
            j+=1
            continue
        elif aaseq_emb[i]=='-':
            codons.append('-') #j=j
            continue
        elif aaseq_pdb[i]!=aaseq_emb[i] and aaseq_pdb[i]!='-' and aaseq_emb[i]!='-':
            codons.append('-')
            j+=1
            continue
    
    return(codons)

    #outfile=open(mut_dir+pdbid+chain+'_'+embid+'_cod.txt','w')
    #outfile.write(codons)
    #outfile.close()

def codon_mut(codon):
    codonspace=set(['A','C','G','T'])
    c13=[]
    c12=[]
    c11=[]
    c23=[]
    c22=[]
    c21=[]
    c33=[]
    for n in range(0,len(codon)):  #codon variable of this function must now be a list !
        for i in range(0,3):
            c13.append(codon[n][0]+codon[n][1]+list(codonspace-set(codon[n][2]))[i])
            c12.append(codon[n][0]+list(codonspace-set(codon[n][1]))[i]+codon[n][2])
            c11.append(list(codonspace-set(codon[n][0]))[i]+codon[n][1]+codon[n][2])
            for j in range(0,3):
                c23.append(codon[n][0]+list(codonspace-set(codon[n][1]))[i]+list(codonspace-set(codon[n][2]))[j])
                c22.append(list(codonspace-set(codon[n][0]))[i]+codon[n][1]+list(codonspace-set(codon[n][2]))[j])
                c21.append(list(codonspace-set(codon[n][0]))[i]+list(codonspace-set(codon[n][1]))[j]+codon[n][2])
                for k in range(0,3):
                    c33.append(list(codonspace-set(codon[n][0]))[i]+list(codonspace-set(codon[n][1]))[j]+list(codonspace-set(codon[n][2]))[k])
    return([c13,c12,c11,c23,c22,c21,c33])

trans={'GCT':'ALA','GCC':'ALA','GCA':'ALA','GCG':'ALA',
       'TGT':'CYS','TGC':'CYS',
       'GAT':'ASP','GAC':'ASP',
       'GAA':'GLU','GAG':'GLU',
       'TTT':'PHE','TTC':'PHE',
       'GGT':'GLY','GGC':'GLY','GGA':'GLY','GGG':'GLY',
       'CAT':'HIS','CAC':'HIS',
       'ATT':'ILE','ATC':'ILE','ATA':'ILE',
       'AAA':'LYS','AAG':'LYS',
       'TTA':'LEU','TTG':'LEU','CTT':'LEU','CTC':'LEU','CTA':'LEU','CTG':'LEU',
       'ATG':'MET',
       'AAT':'ASN','AAC':'ASN',
       'CCT':'PRO','CCC':'PRO','CCA':'PRO','CCG':'PRO',
       'CAA':'GLN','CAG':'GLN',
       'CGT':'ARG','CGC':'ARG','CGA':'ARG','CGG':'ARG','AGA':'ARG','AGG':'ARG',
       'TCT':'SER','TCC':'SER','TCA':'SER','TCG':'SER','AGT':'SER','AGC':'SER',
       'ACT':'THR','ACC':'THR','ACA':'THR','ACG':'THR',
       'GTT':'VAL','GTC':'VAL','GTA':'VAL','GTG':'VAL',
       'TGG':'TRP',
       'TAT':'TYR','TAC':'TYR',
       'TAG':'STP','TGA':'STP','TAA':'STP'}

alls={'ALA':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'CYS':['HIS','ARG','ASP','ILE','ASN','GLN','ALA','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'ASP':['HIS','ARG','ALA','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'GLU':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','ALA','TYR','TRP','VAL','PRO','THR'],
      'PHE':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','ALA','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'GLY':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','ALA','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'HIS':['ALA','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'ILE':['HIS','ARG','ASP','ALA','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'LYS':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','ALA','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'LEU':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','ALA','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'MET':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','ALA','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'ASN':['HIS','ARG','ASP','ILE','ALA','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'PRO':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','ALA','THR'],
      'GLN':['HIS','ARG','ASP','ILE','ASN','ALA','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'ARG':['HIS','ALA','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','THR'],
      'SER':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','ALA','GLU','TYR','TRP','VAL','PRO','THR'],
      'THR':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','VAL','PRO','ALA'],
      'VAL':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','TRP','ALA','PRO','THR'],
      'TRP':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','TYR','ALA','VAL','PRO','THR'],
      'TYR':['HIS','ARG','ASP','ILE','ASN','GLN','CYS','PHE','LYS','GLY','LEU','MET','SER','GLU','ALA','TRP','VAL','PRO','THR']}

snps={'ALA':['VAL','ASP','GLU','GLY','THR','PRO','SER'],
      'CYS':['TYR','SER','PHE','TRP','ARG','GLY'],
      'ASP':['VAL','ALA','GLY','GLU','ASN','HIS','TYR'],
      'GLU':['VAL','ALA','GLY','ASP','LYS','GLN'],
      'PHE':['LEU','ILE','VAL','SER','TYR','CYS'],
      'GLY':['ASP','GLU','ALA','VAL','SER','ARG','CYS','TRP'],
      'HIS':['ARG','PRO','LEU','GLN','TYR','ASN','ASP'],
      'ILE':['PHE','LEU','MET','VAL','THR','ASN','LYS','SER','ARG'],
      'LYS':['ARG','THR','MET','ILE','ASN','GLN','GLU'],
      'LEU':['PHE','ILE','MET','VAL','SER','TRP','PRO','HIS','GLN','ARG'],
      'MET':['LEU','ILE','VAL','THR','LYS','ARG'],
      'ASN':['LYS','ASP','HIS','TYR','SER','THR','ILE'],
      'PRO':['LEU','HIS','GLN','ARG','SER','THR','ALA'],
      'GLN':['LEU','PRO','ARG','HIS','LYS','GLU'],
      'ARG':['HIS','GLN','PRO','LEU','CYS','TRP','SER','GLY','LYS','THR','MET','ILE'],
      'SER':['PHE','LEU','TYR','CYS','TRP','PRO','THR','ALA','ARG','GLY','ASN','ILE'],
      'THR':['ILE','MET','ASN','LYS','SER','ARG','ALA','PRO'],
      'VAL':['ILE','MET','LEU','PHE','ALA','ASP','GLU','GLY'],
      'TRP':['CYS','ARG','GLY','SER','LEU'],
      'TYR':['CYS','SER','PHE','HIS','ASN','ASP']}
  
mnps={'ALA':['HIS','ARG','ILE','TYR','ASN','GLN','CYS','TRP','PHE','LYS','LEU','MET'],
      'CYS':['HIS','ASP','ILE','GLU','ASN','GLN','VAL','LYS','LEU','MET','ALA','PRO','THR'],
      'ASP':['ARG','ILE','SER','GLN','CYS','TRP','PHE','LYS','LEU','MET','PRO','THR'],
      'GLU':['HIS','ARG','ILE','SER','TYR','ASN','CYS','TRP','PHE','LEU','MET','PRO','THR'],
      'PHE':['HIS','ARG','ASP','GLU','ASN','GLN','TRP','GLY','LYS','MET','ALA','PRO','THR'],
      'GLY':['HIS','ILE','TYR','ASN','GLN','PHE','LYS','LEU','MET','PRO','THR'],
      'HIS':['ILE','SER','GLU','CYS','TRP','VAL','PHE','LYS','GLY','MET','ALA','THR'],
      'ILE':['HIS','ASP','GLU','TYR','GLN','CYS','TRP','GLY','ALA','PRO'],
      'LYS':['HIS','ASP','SER','TYR','CYS','TRP','VAL','PHE','GLY','LEU','ALA','PRO'],
      'LEU':['ASP','GLU','TYR','ASN','CYS','GLY','LYS','ALA','THR'],
      'MET':['HIS','ASP','SER','GLU','TYR','ASN','GLN','CYS','TRP','PHE','GLY','ALA','PRO'],
      'ASN':['ARG','GLU','GLN','CYS','TRP','VAL','PHE','GLY','LEU','MET','ALA','PRO'],
      'PRO':['ASP','ILE','GLU','TYR','ASN','CYS','TRP','VAL','PHE','LYS','GLY','MET'],
      'GLN':['ASP','ILE','SER','TYR','ASN','CYS','TRP','VAL','PHE','GLY','MET','ALA','THR'],
      'ARG':['ASP','GLU','TYR','ASN','VAL','PHE','ALA'],
      'SER':['HIS','ASP','GLU','GLN','VAL','LYS','MET'],
      'THR':['HIS','ASP','GLU','TYR','GLN','CYS','TRP','VAL','PHE','GLY','LEU'],
      'VAL':['HIS','ARG','SER','TYR','ASN','GLN','CYS','TRP','LYS','PRO','THR'],
      'TRP':['HIS','ASP','ILE','GLU','TYR','ASN','GLN','VAL','PHE','LYS','MET','ALA','PRO','THR'],
      'TYR':['ARG','ILE','GLU','GLN','TRP','VAL','GLY','LYS','LEU','MET','ALA','PRO','THR']}

#mnps={}
#for key in snps.keys():
#    mnps[key]= list(set(alls[key])-set(snps[key]))
#print(mnps)

import pandas as pd
import os.path
def getmut(mylist):
    for x in mylist:
        pdbid=x[0:4]
        chain=x[5]
        embid=x.split('\t')[1]

        if os.path.isfile(mut_dir+pdbid+chain+'_'+embid+'_mut.txt'):
            continue
    
        alipath=ali_dir+pdbid+chain+'_'+embid+'_align.fa'
        embpath=nas_dir+pdbid+chain+'_'+embid+'.txt'
        try:
            codons=get_codons(alipath,embpath)
        except OSError as e:
            print(pdbid,chain,embid,e)
            continue
            
        output=pd.DataFrame(data={"codon":[],"res_type":[],"all":[],"snp":[],"mnp":[],
                                  "m13":[],"m12":[],"m11":[],"m23":[],"m22":[],"m21":[],"m33":[],
                                  "m13c":[],"m12c":[],"m11c":[],"m23c":[],"m22c":[],"m21c":[],"m33c":[]})
        for idx,codon in enumerate(codons):
            if codon=='-': #gap or mismatch
                output.loc[idx,:] = '-'
            elif codon not in trans: #IUPAC nucleotide code other than ATCG
                output.loc[idx,:] = 'UNK'
            elif trans[codon]=='STP': #stop codon
                output.loc[idx,:] = 'STP'
            else:
                clist=[key  for (key, value) in trans.items() if value==trans[codon]]
                clist.remove(codon)
                
                output.loc[idx,"codon"] = codon
                output.loc[idx,"res_type"] = trans[codon]
                
                output.loc[idx,"all"] = '_'.join(alls[trans[codon]])
                output.loc[idx,"snp"] = '_'.join(snps[trans[codon]])
                output.loc[idx,"mnp"] = '_'.join(mnps[trans[codon]])
                
                aamut=[trans[codon_mut([codon])[0][n]] for n in range(0,len(codon_mut([codon])[0]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m13"] = '_'.join(aamut)
                aamut=[trans[codon_mut([codon])[1][n]] for n in range(0,len(codon_mut([codon])[1]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m12"] = '_'.join(aamut)
                aamut=[trans[codon_mut([codon])[2][n]] for n in range(0,len(codon_mut([codon])[2]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m11"] = '_'.join(aamut)
                aamut=[trans[codon_mut([codon])[3][n]] for n in range(0,len(codon_mut([codon])[3]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m23"] = '_'.join(aamut)
                aamut=[trans[codon_mut([codon])[4][n]] for n in range(0,len(codon_mut([codon])[4]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m22"] = '_'.join(aamut)
                aamut=[trans[codon_mut([codon])[5][n]] for n in range(0,len(codon_mut([codon])[5]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m21"] = '_'.join(aamut)
                aamut=[trans[codon_mut([codon])[6][n]] for n in range(0,len(codon_mut([codon])[6]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m33"] = '_'.join(aamut)
                
                aamut=[trans[codon_mut(clist)[0][n]] for n in range(0,len(codon_mut(clist)[0]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m13c"] = '_'.join(aamut)
                aamut=[trans[codon_mut(clist)[1][n]] for n in range(0,len(codon_mut(clist)[1]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m12c"] = '_'.join(aamut)
                aamut=[trans[codon_mut(clist)[2][n]] for n in range(0,len(codon_mut(clist)[2]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m11c"] = '_'.join(aamut)
                aamut=[trans[codon_mut(clist)[3][n]] for n in range(0,len(codon_mut(clist)[3]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m23c"] = '_'.join(aamut)
                aamut=[trans[codon_mut(clist)[4][n]] for n in range(0,len(codon_mut(clist)[4]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m22c"] = '_'.join(aamut)
                aamut=[trans[codon_mut(clist)[5][n]] for n in range(0,len(codon_mut(clist)[5]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m21c"] = '_'.join(aamut)
                aamut=[trans[codon_mut(clist)[6][n]] for n in range(0,len(codon_mut(clist)[6]))]
                aamut=[e for e in aamut if e!=trans[codon]]
                output.loc[idx,"m33c"] = '_'.join(aamut)
        
        output=output[["codon","res_type","all","snp","mnp",
                       "m13","m12","m11","m23","m22","m21","m33",
                       "m13c","m12c","m11c","m23c","m22c","m21c","m33c"]]
        output.to_csv(mut_dir+pdbid+chain+'_'+embid+'_mut.txt',index=False,sep='\t')

listfile=open(map_pth, 'r')
pdblist=listfile.readlines()
getmut(pdblist)
listfile.close()  