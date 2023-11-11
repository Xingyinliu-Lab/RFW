import pandas as pd
import sqlite3
import os
from multiprocessing import Pool
import numpy as np


import sys
output_dir=sys.argv[1]
db_dir=sys.argv[2]
processors=int(sys.argv[3])
fun_db_file=sys.argv[4]
anno_type=sys.argv[5]

mapping_taxa_df=pd.read_csv(output_dir+'/mapping_taxa_df.csv',sep=',',header=0)
taxonomyInfoDB_file= db_dir+'/taxonomyInfoDB.db'

mapping_taxa_dict = mapping_taxa_df.set_index('nid')['refseq_id'].to_dict()
def get_rep_species(refid):
    sql_str = "SELECT refseq_id,gtdb_genome_representative FROM taxonomy_info where refseq_id='{refid}';"
    str_value = {'refid': refid}
    sql_str = sql_str.format(**str_value)
    #print(sql_str)
    with sqlite3.connect(taxonomyInfoDB_file) as conn:
        g_df = pd.read_sql(sql_str, conn)
    rep_gtdb= g_df.loc[0,'gtdb_genome_representative']
    sql_str = "SELECT refseq_id,gtdb_accession FROM taxonomy_info where gtdb_accession='{rep_gtdb}';"
    str_value = {'rep_gtdb': rep_gtdb}
    sql_str = sql_str.format(**str_value)
    # print(sql_str)
    with sqlite3.connect(taxonomyInfoDB_file) as conn:
        g_df = pd.read_sql(sql_str, conn)
    rep_refid = g_df.loc[0, 'refseq_id']
    return rep_refid


def addc(xstr):
    return "'" + xstr + "'"

def get_fun(nid):
    refseq_id=mapping_taxa_dict.get(nid)
    rep_refseq_id=get_rep_species(refseq_id)
    inchi_str="select e_id,copynumber from taxa_fun where refseq_id= '{target}';"
    str_value = {'target': str(rep_refseq_id)}
    inchi_str = inchi_str.format(**str_value)
    with sqlite3.connect(fun_db_file) as conn:
        gid_fun = pd.read_sql(inchi_str, conn)
    if len(gid_fun)>0:
        gid_fun.to_csv(output_dir + '/taxon_fun/' + str(nid) + '_'+anno_type+'.tsv', index=False, sep='\t')

def get_fun_cl(cl):
    for nid in cl:
        if not os.path.exists(output_dir + '/taxon_fun/' + str(nid) + '_'+anno_type+'.tsv'):
            get_fun(nid)

if __name__ == '__main__':
    processor = processors
    cl = np.array_split(np.asarray(np.asarray(list(mapping_taxa_df['nid']))), processor, axis=0)
    reslist = []
    p = Pool(processor)
    for i in range(processor):
        reslist.append(p.apply_async(get_fun_cl, args=(cl[i], )))
    p.close()
    p.join()