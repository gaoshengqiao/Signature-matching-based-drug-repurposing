import pandas as pd
import numpy as np
import os
import scipy

def Jaccard(S_up,S_up2):
    INT =S_up.intersection(S_up2)
    UNI= S_up.union(S_up2)
    jaccard =len(INT)/len(UNI)
    return jaccard

data = 'LINCS database'
files =os.listdir(data)
query =pd.read_csv('./LPS input signatures.txt',sep=',')
gene = pd.read_csv('./metadata/Transcribed_Gene_Metadata.txt',sep='\t')
metadata =  pd.read_csv('./metadata/CD_signature_Metadata.csv',sep=',')
for num_sig in range(len(query)):
    final = pd.DataFrame()
    query_DE_up = pd.Index(set(query.loc[num_sig,'up_GENE'].split(",")))
    query_DE_up = pd.DataFrame(query_DE_up)
    query_DE_up2 = np.unique(query_DE_up.merge(gene,left_on=0,right_on='pr_gene_symbol')['pr_gene_id'])
    query_DE_up = pd.Index(query_DE_up2) 
    query_DE_down =  pd.Index(set(query.loc[num_sig,'down_GENE'].split(",")))
    query_DE_down = pd.DataFrame(query_DE_down)
    query_DE_down2 = np.unique(query_DE_down.merge(gene,left_on=0,right_on='pr_gene_symbol')['pr_gene_id'])
    query_DE_down = pd.Index(query_DE_down2)
    print(query.loc[num_sig,"drug"])
    n_up = len(query_DE_up)
    n_down = len(query_DE_down)
    for file in files:
        print(file)
        a = pd.read_csv('%s/'%data+file,dtype=str)
        out =np.zeros((len(a),1))
        out2=list() 
        for J in range(len(a)):

            DEG_up = pd.Index(set((a.loc[J,'up_GENE']).split(",")))
            DEG_down = pd.Index((set((a.loc[J,'down_GENE']).split(","))))
            SJ_score = Jaccard(query_DE_up,DEG_down)/2
            out2.append(SJ_score)                 
        out = pd.DataFrame(out,columns = ["Jaccard"])     
        out["Jaccard"]=pd.DataFrame(np.array(out2))
        out=pd.concat([a['sig_id'],out],axis=1)        
        final = pd.concat([final,out])
    final = final.merge(metadata,on='sig_id')
    final = final.sort_values(by="Jaccard",ascending=False)

    final.to_csv("output_%s.csv"%(query.loc[num_sig,'sig_id']),index=False)


