import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
import os
import sys
import time
import argparse


parser = argparse.ArgumentParser(description='DFSCA argparse')
parser.add_argument('--col_info_name', type=str, default='')
parser.add_argument('--rel_abu_df_name', type=str, default='')
parser.add_argument('--fun_abu_df_name', type=str, default='')
parser.add_argument('--prefix', type=str, default='')
parser.add_argument('--output_dir', type=str, default='')
parser.add_argument('--control_label', type=str, default='')
parser.add_argument('--treat_label', type=str, default='')
parser.add_argument('--minimum_taxa_abundance_control', type=float, default=0)
parser.add_argument('--minimum_taxa_detection_num', type=int, default=5)

args = parser.parse_args()

col_info_name=args.col_info_name
rel_abu_df_name=args.rel_abu_df_name
fun_abu_df_name=args.fun_abu_df_name
prefix=args.prefix
output_dir=args.output_dir
minimum_taxa_abundance_control=args.minimum_taxa_abundance_control
minimum_taxa_detection_num=args.minimum_taxa_detection_num
control_label=args.control_label
treat_label=args.treat_label



if not os.path.exists(col_info_name):
    print('Col info file error.')
    time.sleep(15)
    sys.exit()

if not os.path.exists(fun_abu_df_name):
    print('Fun abu file error.')
    time.sleep(15)
    sys.exit()

if not os.path.exists(rel_abu_df_name):
    print('Abu file error.')
    time.sleep(15)
    sys.exit()

if not os.path.exists(output_dir + '/'):
    os.mkdir(output_dir + '/')

col_info=pd.read_csv(col_info_name,sep=',',header=0)
if 'NCBI_ID' in list(col_info['type']):
    taxonomy_id_label=list(col_info.loc[col_info['type']=='NCBI_ID','col'])[0]
else:
    taxonomy_id_label = list(col_info.loc[col_info['type'] == 'GTDB_ANNO', 'col'])[0]
sample_list=list(col_info.loc[col_info['type']=='sample','col'])

rel_abu_df=pd.read_csv(rel_abu_df_name,sep=',',header=0)
tax_list=list(rel_abu_df.columns)
tax_list.remove('sample')
rel_abu_df=pd.merge(rel_abu_df,col_info[['col','group']],how='left',left_on='sample',right_on='col')
rel_abu_df=rel_abu_df[tax_list+['group']]

fun_abu_df=pd.read_csv(fun_abu_df_name,sep=',',header=0)
fun_list=list(fun_abu_df.columns)
fun_list.remove('sample')
fun_abu_df=pd.merge(fun_abu_df,col_info[['col','group']],how='left',left_on='sample',right_on='col')
fun_abu_df=fun_abu_df[fun_list+['group']]


def preprocess_qmd(data_df, minimum_taxa_detection_num, control_label, treat_label, group_label, minimum_taxa_abundance_control):

    minimum_taxa_detection_num = int(minimum_taxa_detection_num)
    minimum_taxa_abundance_control = float(minimum_taxa_abundance_control)
    control = control_label
    treat = treat_label
    genusdata = data_df
    genusdata = genusdata.loc[(genusdata[group_label] == control) | (
            genusdata[group_label] == treat),]
    tmpgenuslist = list(genusdata.columns)
    genuslist = []
    for g in tmpgenuslist:
        if g != group_label:
            genuslist.append(g)
    sumabundance = genusdata[genuslist].sum(axis=1)
    genusdata[genuslist] = genusdata[genuslist].div(sumabundance, axis='rows')
    genuslist = pd.DataFrame(genuslist)
    genuslist.columns = ['taxaId']
    controllen = sum(genusdata[group_label] == control)
    treatlen = sum(genusdata[group_label] == treat)
    res_detection_ratio = pd.DataFrame(columns=[
        'taxaId', 'controlDetectionNum', 'treatDetectionNum', 'controlDetectionRate', 'treatDetectionRate',
        'controlAbundance'])
    count = 0
    for i in genuslist.index:
        gid = genuslist.loc[i, 'taxaId']
        # print(i, gid)
        try:
            controlgdata = genusdata.loc[genusdata[group_label]
                                         == control, gid]
            treatgdata = genusdata.loc[genusdata[group_label] == treat, gid]
            controlobs = sum(controlgdata > 0)
            treatobs = sum(treatgdata > 0)
            res_detection_ratio.loc[count, 'taxaId'] = gid
            res_detection_ratio.loc[count, 'controlDetectionNum'] = controlobs
            res_detection_ratio.loc[count, 'treatDetectionNum'] = treatobs
            res_detection_ratio.loc[count,
            'controlDetectionRate'] = controlobs / controllen
            res_detection_ratio.loc[count,
            'treatDetectionRate'] = treatobs / treatlen
            if controlobs > 1:
                res_detection_ratio.loc[count, 'controlAbundance'] = np.nanmedian(
                    controlgdata)
            else:
                res_detection_ratio.loc[count, 'controlAbundance'] = 0
            count = count + 1
        except:
            pass


    genuslist = res_detection_ratio.loc[(res_detection_ratio['controlDetectionNum'] >= minimum_taxa_detection_num)
                                        & (res_detection_ratio['treatDetectionNum'] >= minimum_taxa_detection_num)
                                        & (res_detection_ratio['controlAbundance'] >= minimum_taxa_abundance_control),]

    genuslist = genuslist.reset_index(drop=True)
    genusIntoModel = list(genuslist['taxaId'])

    genusdata = genusdata[genusIntoModel + [group_label]]

    genusIntoModel = list(genuslist['taxaId'])
    genusdata = genusdata[genusIntoModel + [group_label]]
    sumabundance = genusdata[genusIntoModel].sum(axis=1)
    genusdata[genusIntoModel] = genusdata[genusIntoModel].div(
        sumabundance, axis='rows')

    for i in genuslist.index:
        gid = genuslist.loc[i, 'taxaId']
        try:
            controlgdata = genusdata.loc[genusdata[group_label]
                                         == control, gid]
            treatgdata = genusdata.loc[genusdata[group_label] == treat, gid]
            controlgdata = controlgdata[controlgdata > 0]
            treatgdata = treatgdata[treatgdata > 0]
            mean_t= np.nanmean(np.log2(treatgdata))
            mean_c= np.nanmean(np.log2(controlgdata))
            mean_diff = mean_t-mean_c
            genuslist.loc[i, 'treat_logged_mean'] = mean_t
            genuslist.loc[i, 'control_logged_mean'] = mean_c
            genuslist.loc[i, 'logged_mean_diff'] = mean_diff
        except:
            pass
    genuslist['ID'] = genuslist.index
    genuslist = genuslist.reset_index(drop=True)
    return genuslist


def qmd_optimization(data_df):
    def cal_cost(x, mod, detectionV):
        tmpmod = np.tile(mod, (x.shape[0], 1))
        tmpdetection = np.tile(detectionV, (x.shape[0], 1))
        delta = x
        tmpmod = np.abs(tmpmod + delta)
        tmpmod = tmpmod * tmpdetection
        j1 = tmpmod.sum(axis=1)
        j = j1 / tmpmod.shape[1]
        return j
    genuslist = data_df
    genuslist['D'] = 1 / 2 * \
                     (genuslist['controlDetectionRate'] + genuslist['treatDetectionRate'])
    detectionV = np.array(genuslist['D'])
    mod = np.asarray(genuslist['logged_mean_diff'])
    arr2 = np.array(range(2000)).reshape(2000, 1)
    arr2 = arr2 / 100 - 10
    res = cal_cost(arr2, mod, detectionV)
    minres = np.min(res)
    pdelta = arr2[res == minres]
    absdelta = np.abs(pdelta)
    pdelta = pdelta[np.argmin(absdelta)][0]
    return pdelta


def qmd(data_df, minimum_taxa_detection_num, control_label, treat_label, group_label, minimum_taxa_abundance_control):
    pre_df = preprocess_qmd(data_df, minimum_taxa_detection_num, control_label, treat_label, group_label, minimum_taxa_abundance_control)
    pdelta = qmd_optimization(pre_df)
    return pdelta


def QMD_DTAE_BC(data_df, control_label, treat_label, group_label, pdelta,qtype):
    genusdata = data_df
    taxalist = []
    control = control_label
    treat = treat_label
    res = pd.DataFrame()
    count = 0
    for taxa in genusdata.columns:
        if taxa == group_label:
            continue
        taxalist.append(taxa)
        try:
            controlgdata = genusdata.loc[genusdata[group_label] == control, taxa]
            treatgdata = genusdata.loc[genusdata[group_label] == treat, taxa]
            controlgdata = controlgdata[controlgdata > 0]
            treatgdata = treatgdata[treatgdata > 0]

            res.loc[count, 'control_relabu'] = np.mean(controlgdata)
            res.loc[count, 'treat_relabu'] = np.mean(treatgdata)

            loggedControl = np.log2(controlgdata)
            loggedTreat = np.log2(treatgdata)
            loggedControl = loggedControl[~np.isnan(loggedControl)]
            loggedTreat = loggedTreat[~np.isnan(loggedTreat)]
            mean_diff = np.mean(loggedTreat) - np.mean(loggedControl)
            _, qmd_pvalue = mannwhitneyu(
                loggedControl, loggedTreat + pdelta, alternative='two-sided')
            res.loc[count, qtype] = taxa
            res.loc[count, 'rel_diff'] = mean_diff
            res.loc[count, 'qmd_diff'] = pdelta + mean_diff
            res.loc[count, 'qmd_pvalue'] = qmd_pvalue
            _, rel_pvalue = mannwhitneyu(
                loggedControl, loggedTreat, alternative='two-sided')
            res.loc[count, 'rel_pvalue'] = rel_pvalue
            count = count + 1
        except:
            pass
    _, qvaluelist = fdrcorrection(res['qmd_pvalue'], alpha=0.05, method='indep', is_sorted=False)
    res['qmd_qvalue'] = qvaluelist
    _, qvaluelist = fdrcorrection(res['rel_pvalue'], alpha=0.05, method='indep', is_sorted=False)
    res['rel_qvalue'] = qvaluelist

    return res

pdelta = qmd(rel_abu_df, minimum_taxa_detection_num, control_label, treat_label, 'group', minimum_taxa_abundance_control)
qmd_abu=QMD_DTAE_BC(rel_abu_df, control_label, treat_label, 'group', pdelta,'Taxon')
qmd_tae=QMD_DTAE_BC(fun_abu_df, control_label, treat_label, 'group', pdelta,'FSCA')

qmd_abu.to_csv(output_dir+'/'+prefix+'_QMD_taxon_Abu_diff.csv',sep=',',index=False)
qmd_tae.to_csv(output_dir+'/'+prefix+'_DFSCA_BC_diff.csv',sep=',',index=False)