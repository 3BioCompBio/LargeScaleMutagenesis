from config import ext,lst_pth,mus_dir,aas_dir

aa_dict={'ALA':'A','ARG':'R','ASN':'N',
         'ASP':'D','ASX':'B','CYS':'C',
         'GLU':'E','GLN':'Q','GLX':'Z',
         'GLY':'G','HIS':'H','ILE':'I',
         'LEU':'L','LYS':'K','MET':'M',
         'PHE':'F','PRO':'P','SER':'S',
         'THR':'T','TRP':'W','TYR':'Y',
         'VAL':'V','XAA':'X','UNK':'X'}

import pandas as pd
def fasta_mus(mylist):
    for x in mylist:
        pdbid=x[0:4]
        chain=x[5]
        filepath=mus_dir+pdbid+"."+ext
        
        #open file
        try:
            musfile=pd.read_csv(filepath,sep='\s+',comment='#',header=None)
        except Exception as e:
            print(pdbid,chain,ext,e)
            continue
            
        #focus on chain if applicable
        if chain!='-':
            musfile=musfile[musfile[0].astype(str)==chain]   
            if len(musfile)==0:
                print(pdbid,chain,ext,': chain not found')
                continue
                
        #record sequence
        seq=''
        for i in range(0,len(musfile)//19):
            aa3=str(musfile.iloc[i*19,2])
            aa1=aa_dict[aa3]
            seq+=aa1
        
        #output file
        output=open(aas_dir+pdbid+chain+'.txt','w')
        output.write(">"+pdbid+":"+chain+"|PDBID|CHAIN|SEQUENCE"+"\n")
        for j in range(0,len(seq)):
            if j!=0 and j%80==0:
                output.write(seq[j]+"\n")
            else:
                output.write(seq[j])
        output.close()

listfile=open(lst_pth, 'r')
pdblist=listfile.readlines()
fasta_mus(pdblist)