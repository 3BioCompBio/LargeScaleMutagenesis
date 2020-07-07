from config import lst_pth,map_pth,ali_dir

import os
os.chdir(ali_dir)

import glob
def grab_ali(pdbid,chain):
    fpathlist=glob.glob(pdbid+chain+'*.fa')
    return(fpathlist)

#make sure that if we have the same id% and same simil% we only select one of the two
def comp_seqid(files):
    tag=files[6:].split('_')[0]
    #Open alignment and extract sequences
    alignmentfile=open(files,'r')
    alignment=alignmentfile.read()
    header1=False
    header2=False
    seq1=''
    seq2=''
    for lines in alignment.split('\n'):
        if header1 and '>' in lines or ';' in lines:
            header2=True
            continue
        if '>' in lines or ';' in lines:
            header1=True
            continue
        if header1 and not header2:
            seq1=seq1+lines
        if header1 and header2:
            seq2=seq2+lines
    alignmentfile.close()
    #Compute sequences identity and similarity
    simlen=0
    simcount=0
    idlen=0
    idcount=0
    for aa in range(0,len(seq1)):
        simlen+=1
        if seq1[aa]==seq2[aa]:
            simcount+=1
        if seq1[aa]!='-' and seq2[aa]!='-':
            idlen+=1
            if seq1[aa]==seq2[aa]:
                idcount+=1   
    #Error test
    if simlen==0 or idlen==0 :
        identity=0
        similarity=0
        print(files+' : alignment missing')
    else:
        identity=idcount/idlen
        similarity=simcount/simlen
    return(tag,identity,similarity)

import pandas as pd
def ali_map(mylist):
    mydict={}
    for x in mylist:
        pdbid=x[0:4]
        chain=x[5]
        for ali in grab_ali(pdbid,chain):
            res=comp_seqid(ali)
            tag=res[0]
            ident=res[1]
            simil=res[2]
            if pdbid+'_'+chain in mydict:
                if ident>mydict[pdbid+'_'+chain][1]:
                    mydict.update({pdbid+'_'+chain:(tag,ident,simil)})
                if ident==mydict[pdbid+'_'+chain][1] and simil>mydict[pdbid+'_'+chain][2]:
                    mydict.update({pdbid+'_'+chain:(tag,ident,simil)})
            else:
                mydict.update({pdbid+'_'+chain:(tag,ident,simil)})
    return(pd.DataFrame.from_dict(mydict,orient='index'))

def IDthr(mydf,thr):
    mydf=mydf[mydf[1]>thr]
    return(mydf)

listfile=open(lst_pth, 'r')
pdblist=listfile.readlines()
listfile.close()
out=IDthr(ali_map(pdblist),0.95)
out.to_csv(map_pth,header=None,index=True,sep='\t')