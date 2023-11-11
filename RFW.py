import pandas as pd
import taxoniq
import sqlite3
import os
import numpy as np
import sys
import time
import argparse


parser = argparse.ArgumentParser(description='RFW argparse')
parser.add_argument('--processors', type=int, default=4)
parser.add_argument('--col_info_name', type=str, default='')
parser.add_argument('--abu_df_name', type=str, default='')
parser.add_argument('--output_dir', type=str, default='')
parser.add_argument('--anno_type', type=str, default='EC')# EC KO COG
parser.add_argument('--rfw_type', type=str, default='FSCA')# FA FSCA
parser.add_argument('--db_place', type=str, default='')
parser.add_argument('--taxonomy_id_type', type=str, default='')#NCBI_ID / GTDB_ANNO

args = parser.parse_args()

processors=args.processors
col_info_name=args.col_info_name
abu_df_name=args.abu_df_name
output_dir=args.output_dir
anno_type=args.anno_type
db_place=args.db_place
rfw_type=args.rfw_type
taxonomy_id_type=args.taxonomy_id_type



if taxonomy_id_type!='NCBI_ID' and taxonomy_id_type!='GTDB_ANNO':
    print('taxonomy_id_type type error.')
    time.sleep(15)
    sys.exit()

if rfw_type!='FA' and rfw_type!='FSCA':
    print('RFW type error.')
    time.sleep(15)
    sys.exit()

if anno_type!='EC' and anno_type!='KO' and anno_type!='COG':
    print('Anno type error.')
    time.sleep(15)
    sys.exit()

if not os.path.exists(col_info_name):
    print('Col info file error.')
    time.sleep(15)
    sys.exit()

if not os.path.exists(abu_df_name):
    print('Abu file error.')
    time.sleep(15)
    sys.exit()


taxonomyInfoDB_file = db_place+'/taxonomyInfoDB.db'
if anno_type=='EC':
    fun_db_file=db_place+'/fundb-EC.db'
if anno_type=='KO':
    fun_db_file=db_place+'/fundb-KO.db'
if anno_type=='COG':
    fun_db_file=db_place+'/fundb-COG.db'

if not os.path.exists(db_place+'/taxonomyInfoDB.db'):
    print('DB fileplace error or taxonomyInfoDB.db was not found.')
    time.sleep(15)
    sys.exit()

if not os.path.exists(fun_db_file):
    print('DB fileplace error or function db file was not found.')
    time.sleep(15)
    sys.exit()

if not os.path.exists(output_dir + '/'):
    os.mkdir(output_dir + '/')

if not os.path.exists(output_dir + '/taxon_fun/'):
    os.mkdir(output_dir + '/taxon_fun/')

if not os.path.exists(output_dir + '/taxon_contri/'):
    os.mkdir(output_dir + '/taxon_contri/')


print('Map taxonomy ...')
col_info=pd.read_csv(col_info_name,sep=',',header=0)

if taxonomy_id_type=='NCBI_ID':
    taxonomy_id_label=list(col_info.loc[col_info['type']=='NCBI_ID','col'])[0]
if taxonomy_id_type=='GTDB_ANNO':
    taxonomy_id_label = list(col_info.loc[col_info['type'] == 'GTDB_ANNO', 'col'])[0]


sample_list=list(col_info.loc[col_info['type']=='sample','col'])

abu_df=pd.read_csv(abu_df_name,sep=',',header=0,converters={taxonomy_id_label:str})
abu_df=abu_df[abu_df[taxonomy_id_label]!='-1']
abu_df.reset_index(drop=True,inplace=True)

taxon_list=list(set(list(abu_df[taxonomy_id_label])))
abu_df = abu_df.groupby([taxonomy_id_label])[sample_list].sum().reset_index()

if not os.path.exists(output_dir+'/mapping_taxa_df.csv'):
    if taxonomy_id_type == 'NCBI_ID':
        info_df=pd.DataFrame()
        i=0
        for nid in taxon_list:
            info_df.loc[i, 'nid'] = nid
            try:
                t=taxoniq.Taxon(nid)
                tsname=t.scientific_name
                info_df.loc[i, 'strain_name'] = tsname
                trank=str(t.rank).replace('Rank.','')
                s=t.ranked_lineage
                for m in s:
                    mname=m.scientific_name
                    mrank=str(m.rank).replace('Rank.','')
                    mid=str(m.tax_id)
                    if mrank=='species':
                        info_df.loc[i,'specie_id']=str(mid)
                        info_df.loc[i,'specie_name']=mname
                    if mrank=='superkingdom':
                        info_df.loc[i,'k']=mname
            except:
                #prnid
                pass
            i=i+1
        #info_df  nid strain_name  specie_id specie_name k
        sql_str = "SELECT refseq_id,ncbi_name,ncbi_taxid,ncbi_species,ncbi_species_id FROM taxonomy_info;"
        with sqlite3.connect(taxonomyInfoDB_file) as conn:
            all_taxon_info = pd.read_sql(sql_str, conn)
        all_taxon_info.fillna(value=-1,inplace=True)
        all_taxon_info['ncbi_taxid']=all_taxon_info['ncbi_taxid'].apply(lambda x: str(int(x)))
        all_taxon_info['ncbi_species_id'] = all_taxon_info['ncbi_species_id'].apply(lambda x: str(int(x)))

        tax_df1=pd.merge(info_df[~info_df['nid'].isna()],all_taxon_info[~all_taxon_info['ncbi_taxid'].isna()],left_on='nid',right_on='ncbi_taxid',how='left')
        tax_df1=tax_df1[['nid','refseq_id']]

        tax_df2=pd.merge(info_df[~info_df['strain_name'].isna()],all_taxon_info[~all_taxon_info['ncbi_name'].isna()],left_on='strain_name',right_on='ncbi_name',how='left')
        tax_df2=tax_df2[['nid','refseq_id']]

        tax_df3=pd.merge(info_df[~info_df['specie_id'].isna()],all_taxon_info[~all_taxon_info['ncbi_species_id'].isna()],left_on='specie_id',right_on='ncbi_species_id',how='left')
        tax_df3=tax_df3[['nid','refseq_id']]

        tax_df4=pd.merge(info_df[~info_df['specie_name'].isna()],all_taxon_info[~all_taxon_info['ncbi_species'].isna()],left_on='specie_name',right_on='ncbi_species',how='left')
        tax_df4=tax_df4[['nid','refseq_id']]

        tax_df5=pd.merge(info_df[~info_df['nid'].isna()],all_taxon_info[~all_taxon_info['ncbi_species_id'].isna()],left_on='nid',right_on='ncbi_species_id',how='left')
        tax_df5=tax_df5[['nid','refseq_id']]

        tax_df6=pd.merge(info_df[~info_df['strain_name'].isna()],all_taxon_info[~all_taxon_info['ncbi_species'].isna()],left_on='strain_name',right_on='ncbi_species',how='left')
        tax_df6=tax_df6[['nid','refseq_id']]

        tax_df=pd.concat([tax_df1,tax_df2,tax_df3,tax_df4,tax_df5,tax_df6])
        tax_df.drop_duplicates(inplace=True)
        tax_df=tax_df[~tax_df['refseq_id'].isna()]
        tax_df = tax_df.groupby(['nid'])[['refseq_id']].first().reset_index()

        mapping_taxa_df=pd.merge(info_df,tax_df,on='nid',how='inner')
        mapping_taxa_df.reset_index(drop=True,inplace=True)

        no_mapping_taxa_df=info_df[~info_df['nid'].isin(tax_df['nid'])]
        no_mapping_taxa_df.reset_index(drop=True,inplace=True)

    if taxonomy_id_type == 'GTDB_ANNO':
        info_df = pd.DataFrame(taxon_list,columns=['nid'])
        sql_str = "SELECT refseq_id,gtdb_anno,ncbi_name,ncbi_taxid,ncbi_species,ncbi_species_id FROM taxonomy_info;"
        with sqlite3.connect(taxonomyInfoDB_file) as conn:
            all_taxon_info = pd.read_sql(sql_str, conn)
        info_df2=pd.merge(info_df,all_taxon_info,left_on='nid',right_on='gtdb_anno',how='inner')
        info_df2 = info_df2[~info_df2['refseq_id'].isna()]
        info_df2 = info_df2.groupby(['nid'])[['refseq_id','ncbi_name','ncbi_taxid','ncbi_species','ncbi_species_id']].first().reset_index()
        mapping_taxa_df=info_df2
        mapping_taxa_df.reset_index(drop=True, inplace=True)
        no_mapping_taxa_df=info_df[~info_df['nid'].isin(mapping_taxa_df['nid'])]
        no_mapping_taxa_df.reset_index(drop=True,inplace=True)

    mapping_taxa_df.to_csv(output_dir+'/mapping_taxa_df.csv',sep=',',index=False)
    no_mapping_taxa_df.to_csv(output_dir+'/no_mapping_taxa_df.csv',sep=',',index=False)

    mapping_taxa_df = pd.read_csv(output_dir + '/mapping_taxa_df.csv', sep=',', header=0,converters={'nid':str})
    no_mapping_taxa_df = pd.read_csv(output_dir + '/no_mapping_taxa_df.csv', sep=',', header=0)
else:
    mapping_taxa_df=pd.read_csv(output_dir+'/mapping_taxa_df.csv',sep=',',header=0,converters={'nid':str})
    no_mapping_taxa_df=pd.read_csv(output_dir+'/no_mapping_taxa_df.csv',sep=',',header=0)

if len(mapping_taxa_df)<10:
    print('Mapped taxa lessed than 10.')
    time.sleep(15)
    sys.exit()

print('Retrieve taxonomy function ...')
cmdstr='python get_taxon_fun.py '+output_dir+' '+db_place+' '+str(processors)+' '+fun_db_file+' '+anno_type
#print(cmdstr)
os.system(cmdstr)
print('Calculate function abundance...')
if not os.path.exists(output_dir + '/taxon_fun/'+anno_type+'fun_taxon.csv'):
    fun_taxon=[]
    for i in mapping_taxa_df.index:
        nid=mapping_taxa_df.loc[i,'nid']
        if os.path.exists(output_dir + '/taxon_fun/' + str(nid) + '_'+anno_type+'.tsv'):
            tmp_df=pd.read_csv(output_dir + '/taxon_fun/' + str(nid) + '_'+anno_type+'.tsv',header=0,sep='\t')
            tmp_df['nid']=nid
            fun_taxon.append(tmp_df)
    fun_taxon=pd.concat(fun_taxon)
    fun_taxon.to_csv(output_dir + '/taxon_fun/'+anno_type+'fun_taxon.csv',sep=',',index=False)
    fun_taxon = pd.read_csv(output_dir + '/taxon_fun/' + anno_type + 'fun_taxon.csv', sep=',', header=0,converters={'nid':str})
else:
    fun_taxon=pd.read_csv(output_dir + '/taxon_fun/'+anno_type+'fun_taxon.csv',sep=',',header=0,converters={'nid':str})

all_fun_list=list(set(list(fun_taxon['e_id'])))

def addc(xstr):
    return "'" + xstr + "'"

if not os.path.exists(output_dir+'/fun_desc.csv'):
    sql_str = "select e_id,e_name from fun_info where e_id in ({fid});"
    fid = list(map(addc, all_fun_list))
    fid = ','.join(fid)
    str_value = {'fid': fid}
    sql_str = sql_str.format(**str_value)
    with sqlite3.connect(taxonomyInfoDB_file) as conn:
        f_info_df = pd.read_sql(sql_str, conn)
    f_info_df.to_csv(output_dir+'/fun_desc.csv',index=False,sep=',')


if not os.path.exists(output_dir + '/rel_abu.csv'):
    abu_df=abu_df[abu_df[taxonomy_id_label].isin(list(mapping_taxa_df['nid']))]
    abu_df.reset_index(drop=True,inplace=True)
    taxon_list=list(mapping_taxa_df['nid'])
    abu_df.index=abu_df[taxonomy_id_label]
    abu_dft=abu_df.T
    abu_dft['sample']=abu_dft.index
    abu_dft.drop(index=taxonomy_id_label, inplace=True)
    abu_dft.reset_index(drop=True, inplace=True)
    sumabundance = abu_dft[taxon_list].sum(axis=1)
    abu_dft[taxon_list] = abu_dft[taxon_list].div(sumabundance, axis='rows')
    abu_rel_df=abu_dft
    abu_rel_df.to_csv(output_dir + '/rel_abu.csv', index=False, sep=',')
    abu_rel_df = pd.read_csv(output_dir + '/rel_abu.csv', header=0, sep=',')
else:
    abu_rel_df=pd.read_csv(output_dir + '/rel_abu.csv',header=0,sep=',')

print('Calculate taxon contribution to functions...')

cmdstr='python get_fun_contri.py '+output_dir + '/taxon_fun/'+anno_type+'fun_taxon.csv'+' '+output_dir+' '+output_dir + '/rel_abu.csv'+' '+str(processors)
#print(cmdstr)
os.system(cmdstr)


if not os.path.exists(output_dir + '/fun_abu_df_'+anno_type+'_'+rfw_type+'.csv'):
    fun_abu_df=pd.DataFrame(columns=all_fun_list)
    fun_abu_df['sample']=sample_list
    fun_abu_df.index=sample_list
    fun_abu_df.fillna(value=0, inplace=True)
    abu_rel_df.fillna(value=0,inplace=True)
    count=0
    if rfw_type=='FSCA':
        for nid in list(mapping_taxa_df['nid']):
            count=count+1
            if count%300==0:
                print('\tFunction abu:',count/len(list(mapping_taxa_df['nid'])))
            nid_fun=pd.read_csv(output_dir + '/taxon_fun/' + str(nid) + '_'+anno_type+'.tsv', header=0, sep='\t')
            # ec copynumber
            #无copynumber修正
            nid_abu=list(abu_rel_df[str(nid)])
            for e in list(nid_fun['e_id']):
                fun_abu_df['tmp']=nid_abu
                fun_abu_df[e]= fun_abu_df[e] + fun_abu_df['tmp']
    if rfw_type == 'FA':
        for nid in list(mapping_taxa_df['nid']):
            count = count + 1
            if count % 300 == 0:
                print('\tFunction abu:',count / len(list(mapping_taxa_df['nid'])))
            nid_fun = pd.read_csv(output_dir +'/taxon_fun/' + str(nid) + '_'+anno_type+'.tsv', header=0, sep='\t')
            cp_dict = nid_fun.set_index('e_id')['copynumber'].to_dict()
            # ec copynumber
            # copynumber修正
            nid_abu = np.array(abu_rel_df[str(nid)])
            for e in list(nid_fun['e_id']):
                fun_abu_df['tmp'] = nid_abu * cp_dict.get(e)
                fun_abu_df[e] = fun_abu_df[e] + fun_abu_df['tmp']

    fun_abu_df=fun_abu_df[['sample'] + all_fun_list]
    fun_abu_df.to_csv(output_dir + '/fun_abu_df_'+anno_type+'_'+rfw_type+'.csv', index=False, sep=',')
    fun_abu_df = pd.read_csv(output_dir + '/fun_abu_df_' + anno_type + '_' + rfw_type + '.csv', header=0, sep=',')
else:
    fun_abu_df=pd.read_csv(output_dir + '/fun_abu_df_'+anno_type+'_'+rfw_type+'.csv', header=0, sep=',')


if anno_type=='EC':
    print('Calculate pathway abundance...')
    if not os.path.exists(output_dir + '/pathway_coverage.csv'):
        sql_str = "SELECT pathway_id,pathway_name,ec,ec_count FROM pathway_info;"
        with sqlite3.connect(taxonomyInfoDB_file) as conn:
            pathway = pd.read_sql(sql_str, conn)

        pathway.dropna(inplace=True)
        pw_ec_list = list(pathway['ec'])
        pw_ec_list = list(set(list(','.join(pw_ec_list).split(','))))

        # 全覆盖
        def ecList_to_pathway_by_percent(ec_list):
            covered_pathway_df = pd.DataFrame()
            count = 0
            for p in pathway.index:
                pid = pathway.loc[p, 'pathway_id']
                pec_list = pathway.loc[p, 'ec']
                pname = pathway.loc[p, 'pathway_name']
                pc = pathway.loc[p, 'ec_count']
                if '.' in str(pec_list):
                    pec_list = pec_list.split(',')
                    u = list(set(pec_list) & set(ec_list))
                    if len(u) > 0:
                        covered_pathway_df.loc[count, 'covered_ec_count'] = len(u)
                        covered_pathway_df.loc[count, 'pathway_id'] = pid
                        covered_pathway_df.loc[count, 'pathway_name'] = pname
                        covered_pathway_df.loc[count, 'total_ec_count'] = pc
                        covered_pathway_df.loc[count, 'pathway_ec'] = ','.join(pec_list)
                        covered_pathway_df.loc[count, 'covered_ec'] = ','.join(u)
                        count = count + 1
            covered_pathway_df['covered_ratio'] = covered_pathway_df['covered_ec_count'] / covered_pathway_df[
                'total_ec_count']
            return covered_pathway_df
        covered_pathway_df=ecList_to_pathway_by_percent(all_fun_list)
        covered_pathway_df.to_csv(output_dir + '/pathway_coverage.csv', index=False, sep=',')
        covered_pathway_df = pd.read_csv(output_dir + '/pathway_coverage.csv', header=0, sep=',')
    else:
        covered_pathway_df=pd.read_csv(output_dir + '/pathway_coverage.csv', header=0, sep=',')

    if not os.path.exists(output_dir + '/pathway_abu.csv'):
        pathway_abu=[]
        for i in covered_pathway_df.index:
            if i % 300 == 0:
                print('\tPathway abu:',i / len(covered_pathway_df))

            pid=covered_pathway_df.loc[i, 'pathway_id']
            covered_ec=covered_pathway_df.loc[i, 'covered_ec']
            covered_ec = covered_ec.split(',')
            nid_list=list(set(list(fun_taxon.loc[fun_taxon['e_id'].isin(covered_ec),'nid'])))
            nid_list=[str(x) for x in nid_list]
            sumabundance = abu_rel_df[nid_list].sum(axis=1)
            if not os.path.exists(output_dir + '/taxon_contri/' + pid + '.tsv'):
                pathway_contri = abu_rel_df[['sample'] + nid_list]
                pathway_contri.to_csv(output_dir+'/taxon_contri/'+pid+'.tsv',sep='\t',index=False)
            tmp_pathway_df=pd.DataFrame({'sample':list(abu_rel_df['sample']),pid:sumabundance})
            tmp_pathway_df.index=tmp_pathway_df['sample']
            #tmp_pathway_df=tmp_pathway_df[pid]
            tmp_pathway_dft=tmp_pathway_df.T
            tmp_pathway_dft['pathway_id']=tmp_pathway_dft.index
            tmp_pathway_dft=tmp_pathway_dft[tmp_pathway_dft['pathway_id']!='sample']
            tmp_pathway_dft.reset_index(drop=True,inplace=True)
            pathway_abu.append(tmp_pathway_dft)

        pathway_abu=pd.concat(pathway_abu)

        pathway_abu.index = pathway_abu['pathway_id']
        pathway_abut = pathway_abu.T
        pathway_abut['sample'] = pathway_abut.index
        pathway_abut.drop(index='pathway_id', inplace=True)
        pathway_abut.reset_index(drop=True, inplace=True)

        pathway_abu=pathway_abut
        pathway_abu.to_csv(output_dir + '/pathway_abu.csv', index=False, sep=',')
        pathway_abu = pd.read_csv(output_dir + '/pathway_abu.csv', header=0, sep=',')
    else:
        pathway_abu=pd.read_csv(output_dir + '/pathway_abu.csv', header=0, sep=',')
print('Done')