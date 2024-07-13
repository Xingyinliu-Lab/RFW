# RFW
 Reference based functional profile inference on WMS


## Installation recommendations

```python
conda create -n RFW python=3.7 
conda activate RFW
pip install -r requirements.txt
```

## Instructions for use

### Reference database preparation

Reference database contains following four files.

> 1. fundb-COG.db
> 2. fundb-EC.db
> 3. fundb-KO.db
> 4. taxonomyInfoDB.db

Download reference database from [figshare](https://doi.org/10.6084/m9.figshare.24541876) and unzip the files into one filefolder. This will be used as db_place in RFW arguments.
 




### Arguments for RFW
> 1. col_info_name. Metainfo file.
> 2. abu_df_name. Abundance data file. Both relative abundance or raw counts file are ok.
> 3. output_dir. Place to store the result.
> 4. anno_type. Microbial function annonation type. Accept EC, KO, COG.
> 5. rfw_type. RFW abundance type. Accept FSCA, FA
> 6. db_place. Database filefolder
> 7. taxonomy_id_type. Mapping type on taxon. Accept NCBI_ID GTDB_ANNO
> 8. processors. Jobs in multiprocessing


### Arguments for DFSCA-BC
> 1. col_info_name. Metainfo file
> 2. rel_abu_df_name. Relative abundance data file
> 3. fun_abu_df_name. RFW FSCA file.
> 4. prefix. Prefix added to the result file.
> 5. output_dir. Place to store the result.
> 6. control_label. Label of control group. 
> 7. treat_label. Label of treat group. 
> 8. minimum_taxa_abundance_control. Minimum taxa median abundance in control group. The value should be in [0,1). If the median relative abundance of taxa is less than the set value, the taxa will be filtered out. 
> 9. minimum_taxa_detection_num. Minimum taxa detection samples in each group. The value should be an integer. If a taxa is only detected in a few samples less than the set value in either the control group or the treatment group, it will be filtered out. 

### Recommendation for 'NCBI_ID' or 'GTDB_ANNO'
For Kraken series pipeline derived taxonomic profiling, NCBI_ID was recommended. For metaphlan series tools, we recommend to use command sgb_to_gtdb_profile for transfering Metaphlan SGB to GTDB Taxonomy identifier, hence, the GTDB_ANNO type.

### Input file template

#### Note: Only csv(Comma delimited) supported. Please avoid common illegal characters. Don’t start or end you're the prefix with a space, period. Avoid using spaces or non-alphanumeric characters.
#### Demo input file could be found at 

#### col_info_name: Metainfo file. This file records columns information in the abundance file. Following is a demo for Metainfo file. The 'col' should match your abundance file columns. The 'type' states columns attributes. Choose 'NCBI_ID' or 'GTDB_ANNO' correctly according to your taxon identifier. 


|col|type|note|group|
|-|-|-|-|
|taxonomy_id|NCBI_ID|||		
|S1|sample||control|
|S2|sample||control|
|S3|sample||treat|

#### abu_df_name: Abundance data file. Demo raw counts file with NCBI_ID identifier

|taxonomy_id|S1|S2|S3|
|-|-|-|-|
|546|	10305020|	7104138|	7635799|
|1639133|	1045175|	725191	|776858|
|67827|	173754|	119556|	128371|
|57706|	502237|	346897|	368860|

#### abu_df_name: Abundance data file. Demo raw counts file with GTDB_ANNO identifier

|taxonomy_id|S1|S2|S3|
|-|-|-|-|
|d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__Blautia_A wexlerae|	10305020|	7104138|	7635799|
|d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus salivarius|	1045175|	725191	|776858|
|d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__Ruminococcus_E sp003526955|	173754|	119556|	128371|
|d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium prausnitzii_C|	502237|	346897|	368860|

#### rel_abu_df_name: Relative abundance data file. 

|sample|546|	1639133|	67827|
|-|-|-|-|
|S1|	0|	0	|0.00093498|
|S2|	3.53E-06|	0	|1.69E-05|
|S3|	0	|0.212940924|	0|

#### fun_abu_df_name: RFW FSCA file. 

|sample|	3.4.11.2|	5.4.2.6|	4.2.2.10|
|-|-|-|-|
|ERR688505|	0.1388|	0.1199	|0.0001|
|ERR688506|	0.0966	|0.1134	|0.0326|
|ERR688508|	0.0279	|0.0794	|0.0059|


### Output of RFW

> 1. mapping_taxa_df.csv: taxon mapping result. 
#### Columns definition for NCBI_ID type RFW:
```
 nid: NCBI taxa ID
 strain_name: NCBI taxa strain name 
 specie_id: NCBI speies id
 specie_name: NCBI speies name
 k: Kingdom
 refseq_id: NCBI refseq id
```
#### Columns definition for GTDB_ANNO type RFW:
```
 nid: GTDB taxonomic annotation
 refseq_id: NCBI refseq id
 ncbi_name: NCBI taxa strain name 
 ncbi_taxid: NCBI taxa ID
 ncbi_species_id: NCBI speies id
 ncbi_species: NCBI speies name
```
> 2. no_mapping_taxa_df.csv: unmapped taxa list
#### Columns definition for NCBI_ID type RFW:
```
 nid: NCBI taxa ID
```
#### Columns definition for GTDB_ANNO type RFW:
```
 nid: GTDB taxonomic annotation
```
> 3. fun_abu_df_${anno_type}_${rfw_type}.csv: RFW abundance

#### Output RFW abundance tabel example (EC):
|sample|	3.1.3.16|	2.3.2.23|	...|
|-|-|-|-|
|S1|	0.1003|	0.4831	|...|
|S2|	0.7547	|0.1571	|...|
|...|	...	|...	|...|

#### Output RFW abundance tabel example (COG):
|sample|	COG3715|	COG0853|	...|
|-|-|-|-|
|S1|	0.9749|	0.4305	|...|
|S2|	0.5462	|0.3032	|...|
|...|	...	|...	|...|

#### Output RFW abundance tabel example (KO):
|sample|	K06423|	K02987|	...|
|-|-|-|-|
|S1|	0.2149|	0.1254	|...|
|S2|	0.4754	|0.1341	|...|
|...|	...	|...	|...|

> 4. fun_desc.csv: microbial function description
#### Columns definition:
```
 e_id: microbial function id(EC, KO, COG)
 e_name: microbial function name
```
> 5. rel_abu.csv: microbial relative abundance data
#### Output taxon relative abundance tabel example (NCBI_ID):
|sample|546|	1639133|	67827|
|-|-|-|-|
|S1|	0|	0	|0.00093498|
|S2|	3.53E-06|	0	|1.69E-05|
|S3|	0	|0.212940924|	0|

#### Output taxon relative abundance tabel example (GTDB_ANNO):
|sample|d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobacterium;s__Methanobacterium sp000499765|	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter_A;s__Methanobrevibacter_A oralis|	...|
|-|-|-|-|
|S1|	0|	0	|0.00093498|
|S2|	3.53E-06|	0	|1.69E-05|
|S3|	0	|0.212940924|	0|

> 6. pathway_abu.csv: microbial pathway abundance data if anno_type=='EC'
#### Output pathway abundance tabel example:
|sample|PWY-5695|	PWY-5692|	...|
|-|-|-|-|
|S1|	0.0447|	0.19232	|...|
|S2|	0.078479|	0.14596	|...|
|...|	...	|...|	...|


> 7. pathway_coverage.csv: microbial pathway coverage data if anno_type=='EC'. This file records microbial encoded or "covered" enzyme counts info for each pathway.
#### Columns definition:
```
 pathway_id: Metacyc pathway id
 pathway_name: Metacyc pathway Name
 pathway_ec: ECs involved in the pathway, seperated by comma
 total_ec_count: total number of ECs involved in the pathway
 covered_ec: microbial derived ECs in the pathway, seperated by comma
 covered_ec_count:total number of microbial derived ECs in the pathway
 covered_ratio:	covered_ec_count/total_ec_count
```

> 8. taxon_fun/*.tsv: microbial encoded function for each taxon
#### Columns definition:
```
 e_id: microbial function id(EC, KO, COG)
 copynumber: microbial function copynumber in genome
```

> 9. taxon_contri/*.tsv: sub-community info for each microbial function 
#### Output taxa members and their abundance for each microbial function :
|sample|d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobacterium;s__Methanobacterium sp000499765|	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter_A;s__Methanobrevibacter_A oralis|	...|
|-|-|-|-|
|S1|	0.0244|	0.0315	|...|
|S2|	3.53E-06|	0	|...|
|S3|	...	|...|	...|


### Output of DFSCA-BC

> 1. ${prefix}_QMD_taxon_Abu_diff.csv: estimated taxon absoulte abundance changes between groups
#### Columns definition:
```
 Taxon: microbial taxon id(NCBI_ID,GTDB_ANNO)
 control_relabu: average microbe abundance at control group
 treat_relabu: average microbe abundance at case group
 rel_diff: logged(2-based) relative abundance foldchange between groups 
 rel_pvalue: p value for relative abundance change
 rel_qvalue: q value for relative abundance change
 qmd_diff: logged(2-based) absolute abundance foldchange between groups 
 qmd_pvalue: p value for absolute abundance change
 qmd_qvalue: q value for absolute abundance change
```
> 2. ${prefix}_DFSCA_BC_diff.csv: estimated microbial function absoulte abundance changes 
between groups
#### Columns definition:
```
 FSCA: microbial function id(EC, KO, COG)
 control_relabu: average microbial function abundance at control group
 treat_relabu: average microbial function abundance at case group
 rel_diff: logged(2-based) relative abundance foldchange between groups 
 rel_pvalue: p value for relative abundance change
 rel_qvalue: q value for relative abundance change
 qmd_diff: logged(2-based) absolute abundance foldchange between groups 
 qmd_pvalue: p value for absolute abundance change
 qmd_qvalue: q value for absolute abundance change
```
## Examples

```python
python RFW.py --processors 60 --col_info_name demo_data/col_info_gtdb.csv --abu_df_name demo_data/demo_gtdb_anno_abudance_table.csv --output_dir demo_data/gtdb_res_EC --anno_type EC --rfw_type FA --db_place ../db2 --taxonomy_id_type GTDB_ANNO
```

```python
python RFW.py --processors 60 --col_info_name demo_data/col_info_gtdb.csv --abu_df_name demo_data/demo_gtdb_anno_abudance_table.csv --output_dir demo_data/gtdb_res_KO --anno_type KO --rfw_type FA --db_place ../db2 --taxonomy_id_type GTDB_ANNO
```


```python
python RFW.py --processors 60 --col_info_name demo_data/col_info_gtdb.csv --abu_df_name demo_data/demo_gtdb_anno_abudance_table.csv --output_dir demo_data/gtdb_res_COG --anno_type COG --rfw_type FA --db_place ../db2 --taxonomy_id_type GTDB_ANNO
```



```python
python RFW.py --processors 60 --col_info_name demo_data/col_info_ncbi.csv --abu_df_name demo_data/demo_ncbi_id_abudance_table.csv --output_dir demo_data/ncbi_res_KO --anno_type KO --rfw_type FSCA --db_place ../db2 --taxonomy_id_type NCBI_ID
```


```python
python RFW.py --processors 60 --col_info_name demo_data/col_info_ncbi.csv --abu_df_name demo_data/demo_ncbi_id_abudance_table.csv --output_dir demo_data/ncbi_res_EC --anno_type EC --rfw_type FSCA --db_place ../db2 --taxonomy_id_type NCBI_ID
```

```python
python RFW.py --processors 60 --col_info_name demo_data/col_info_ncbi.csv --abu_df_name demo_data/demo_ncbi_id_abudance_table.csv --output_dir demo_data/ncbi_res_COG --anno_type COG --rfw_type FSCA --db_place ../db2 --taxonomy_id_type NCBI_ID
```

```python
python DFSCA_BC.py --col_info_name demo_data/col_info_ncbi.csv --rel_abu_df_name demo_data/ncbi_res_EC/rel_abu.csv --fun_abu_df_name demo_data/ncbi_res_EC/pathway_abu.csv --prefix pathway_k --output_dir demo_data/ncbi_res_EC/qmd --control_label A --treat_label B --minimum_taxa_abundance_control 0 --minimum_taxa_detection_num 2
```

```python
python DFSCA_BC.py --col_info_name demo_data/col_info_ncbi.csv --rel_abu_df_name demo_data/ncbi_res_EC/rel_abu.csv --fun_abu_df_name demo_data/ncbi_res_EC/fun_abu_df_EC_FSCA.csv --prefix kraken_y --output_dir demo_data/ncbi_res_EC/qmd --control_label A --treat_label B --minimum_taxa_abundance_control 0 --minimum_taxa_detection_num 2
```


```python
python DFSCA_BC.py --col_info_name demo_data/col_info_ncbi.csv --rel_abu_df_name demo_data/ncbi_res_COG/rel_abu.csv --fun_abu_df_name demo_data/ncbi_res_COG/fun_abu_df_COG_FSCA.csv --prefix kraken_y --output_dir demo_data/ncbi_res_COG/qmd --control_label A --treat_label B --minimum_taxa_abundance_control 0 --minimum_taxa_detection_num 2
```

```python
python DFSCA_BC.py --col_info_name demo_data/col_info_ncbi.csv --rel_abu_df_name demo_data/ncbi_res_KO/rel_abu.csv --fun_abu_df_name demo_data/ncbi_res_KO/fun_abu_df_KO_FSCA.csv --prefix kraken_y --output_dir demo_data/ncbi_res_KO/qmd --control_label A --treat_label B --minimum_taxa_abundance_control 0 --minimum_taxa_detection_num 2
```
## License

The RFW software is licensed under the MIT license.