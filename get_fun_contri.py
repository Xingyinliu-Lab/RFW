import pandas as pd
import os
import sys
from multiprocessing import Pool
import numpy as np

fun_taxon_filename=sys.argv[1]
output_dir=sys.argv[2]
abu_rel_df_name=sys.argv[3]
processors=int(sys.argv[4])


abu_rel_df=pd.read_csv(abu_rel_df_name,header=0,sep=',')
fun_taxon=pd.read_csv(fun_taxon_filename,sep=',',header=0,converters={'nid':str})
all_fun_list=list(set(list(fun_taxon['e_id'])))

def cal_cl(cl):
    for c in cl:
        f=all_fun_list[c]
        if not os.path.exists(output_dir + '/taxon_contri/'+f+'.tsv'):
            nid_list=list(set(list(fun_taxon.loc[fun_taxon['e_id']==f,'nid'])))
            nid_list=[str(x) for x in nid_list]
            tmp_df=abu_rel_df[['sample']+nid_list]
            tmp_df.to_csv(output_dir + '/taxon_contri/'+f+'.tsv',index=False,sep='\t')

if __name__ == '__main__':
    processor = processors
    cl = np.array_split(np.asarray(np.asarray(range(len(all_fun_list)))), processor, axis=0)
    reslist = []
    p = Pool(processor)
    for i in range(processor):
        reslist.append(p.apply_async(cal_cl, args=(cl[i], )))
    p.close()
    p.join()