from config import lst_pth,nas_dir

import os
os.chdir(nas_dir)

import requests
def get_fasta_embl(pdblist):
    nouniid=0
    noembid=0
    nofasta=0
    notfound=0
    suppressed=0
    cofail=0
    
    for x in pdblist:
        pdbid=x[0:4]
        chain=x[5]

        ###GET UNIPROT ID
        if chain!='-':
            url='http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=pdb&id='+pdbid+'&chain='+chain
        else:
            url='http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=pdb&id='+pdbid
        try:
            r=requests.get(url)
        except requests.exceptions.RequestException:
            print(pdbid+chain+' Connection Failed')
            cofail+=1
            continue
        n=r.text.find('AC:')
        uniid=r.text[n+4:n+10]
        ###TEST UNIPROT ID
        if n==-1:
            print(pdbid+chain+' : uniprot id not found')
            nouniid=nouniid+1
            continue

        ###GET EMBL ID
        url='http://www.uniprot.org/uniprot/'+uniid+'.txt'
        try:
            r=requests.get(url)
        except requests.exceptions.RequestException:
            print(pdbid+chain+'_'+uniid+' Connection Failed')
            cofail+=1
            continue
        f=r.text
        ebiid=[]
        for item in f.split('\n'): 
            if 'DR   EMBL' in item and item.strip().split(';')[2][1:] != '-':                        # AND / OR ?
                ebiid.append(item.strip().split(';')[2][1:])
        ###TEST EMBL ID	
        if not ebiid:
            print(pdbid+chain+'_'+uniid+' : embl ids not found')
            noembid=noembid+1
            continue

        ###DOWNLOAD FASTA FILE
        fcount=0
        for i in ebiid:
            url='https://www.ebi.ac.uk/ena/data/view/'+i+'&display=fasta'
            try:
                r=requests.get(url)
            except requests.exceptions.RequestException:
                print(pdbid+chain+'_'+uniid+'_'+i+' Connection Failed')
                cofail+=1
                continue
            fasta=r.text
            fpath=pdbid+chain+'_'+i+'.txt'
            for item in fasta.split('\n'):
                if 'is not found' in item:                                                   #
                    print(pdbid+chain+'_'+uniid+'_'+i+' : FASTA file not found')
                    fcount=fcount+1
                    notfound=notfound+1
                    continue
                if 'suppressed' in item:                                                      #elif
                    print(pdbid+chain+'_'+uniid+'_'+i+' : FASTA file suppressed')
                    fcount=fcount+1
                    suppressed=suppressed+1
                    continue
                if '>' in item or ';' in item:                                                #elif
                    if i.split('.')[0] in item:
                        seq_file=open(fpath,'w')
                        seq_file.write(fasta)
                        seq_file.close()
        ###TEST FASTA FILE
        if fcount == len(ebiid):
            nofasta=nofasta+1

        ###CLEAR VARIABLES
        url=None
        r=None
        n=None
        f=None
        fasta=None
    tot=nouniid+noembid+nofasta
    print(str(tot)+' pdb with no results: '+str(nouniid)+' missing uniprot, '+str(noembid)+' with no embl ids, '+str(nofasta)+' with no fasta files ('+str(notfound)+' FASTA not found, '+str(suppressed)+' FASTA suppressed)')
    print('Moreover, there were '+str(cofail)+' failed connection')

listfile=open(lst_pth, 'r')
pdblist=listfile.readlines()
get_fasta_embl(pdblist)