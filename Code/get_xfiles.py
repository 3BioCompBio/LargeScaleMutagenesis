from config import ext,map_pth,mus_dir,mut_dir,mux_dir

import pandas as pd
import os.path
def xmut(mylist):
    for x in mylist:
        pdbid=x[0:4]
        chain=x[5]
        embid=x.split('\t')[1]

        if os.path.isfile(mux_dir+pdbid+chain+'_'+embid+'.'+ext+'x'):
            continue
            
        #open input files
        try:
            mutfile=pd.read_csv(mut_dir+pdbid+chain+'_'+embid+'_mut.txt',sep='\t',na_filter=False)
            musfile=pd.read_csv(mus_dir+pdbid+'.'+ext,sep='\s+',comment='#',header=None)
        except Exception as e:
            print(pdbid,chain,ext,e)
            continue
        #focus on chain
        try:
            musfile=musfile[musfile[0].astype(str)==chain]
        except Exception as e:
            print(pdbid,chain,ext,': chain not found')
            continue
            
        #check if the length are corresponding (for psy)
        if len(musfile)/19!=len(mutfile):
            print(pdbid,chain,ext,': incomplete file')
            continue
            
        #output
        output=musfile[[1,2,3,4,5,6]].copy()
        output.columns=['res_num','res_typ','res_mut','sec_str','sol_acc',ext]
        output=output.reset_index(drop=True)    
        for i in range(0,len(mutfile)):
            if mutfile.loc[i,'codon']=='-' or mutfile.loc[i,'codon']=='STP'or mutfile.loc[i,'codon']=='UNK':
                output.loc[i*19:i*19+18,'cod']=mutfile.loc[i,'codon']
                output.loc[i*19:i*19+18,'all']=float('nan')
                output.loc[i*19:i*19+18,'snp']=float('nan')
                output.loc[i*19:i*19+18,'mnp']=float('nan')
                output.loc[i*19:i*19+18,'m13']=float('nan')
                output.loc[i*19:i*19+18,'m12']=float('nan')
                output.loc[i*19:i*19+18,'m11']=float('nan')
                output.loc[i*19:i*19+18,'m23']=float('nan')
                output.loc[i*19:i*19+18,'m22']=float('nan')
                output.loc[i*19:i*19+18,'m21']=float('nan')
                output.loc[i*19:i*19+18,'m33']=float('nan')
                output.loc[i*19:i*19+18,'m13c']=float('nan')
                output.loc[i*19:i*19+18,'m12c']=float('nan')
                output.loc[i*19:i*19+18,'m11c']=float('nan')
                output.loc[i*19:i*19+18,'m23c']=float('nan')
                output.loc[i*19:i*19+18,'m22c']=float('nan')
                output.loc[i*19:i*19+18,'m21c']=float('nan')
                output.loc[i*19:i*19+18,'m33c']=float('nan')
            else:
                mutlist=['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
                wt=output.loc[i*19,'res_typ']
                mutlist.remove(wt)
                output.loc[i*19:i*19+18,'cod']=mutfile.loc[i,'codon']
                output.loc[i*19:i*19+18,'all']=[mutfile.loc[i,'all'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'snp']=[mutfile.loc[i,'snp'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'mnp']=[mutfile.loc[i,'mnp'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m13']=[mutfile.loc[i,'m13'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m12']=[mutfile.loc[i,'m12'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m11']=[mutfile.loc[i,'m11'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m23']=[mutfile.loc[i,'m23'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m22']=[mutfile.loc[i,'m22'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m21']=[mutfile.loc[i,'m21'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m33']=[mutfile.loc[i,'m33'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m13c']=[mutfile.loc[i,'m13c'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m12c']=[mutfile.loc[i,'m12c'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m11c']=[mutfile.loc[i,'m11c'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m23c']=[mutfile.loc[i,'m23c'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m22c']=[mutfile.loc[i,'m22c'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m21c']=[mutfile.loc[i,'m21c'].split('_').count(x) for x in mutlist]
                output.loc[i*19:i*19+18,'m33c']=[mutfile.loc[i,'m33c'].split('_').count(x) for x in mutlist]
        output.to_csv(mux_dir+pdbid+chain+'_'+embid+'.'+ext+'x',sep='\t',index=False)

listfile=open(map_pth, 'r')
pdblist=listfile.readlines()
xmut(pdblist)
listfile.close()