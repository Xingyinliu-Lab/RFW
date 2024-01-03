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
|ERR688505|	0.138806003|	0.11994831	|0.000120429|
|ERR688506|	0.096647547	|0.113485571	|0.032600066|
|ERR688508|	0.027975259	|0.07947178	|0.005954884|


### Output of RFW

> 1. mapping_taxa_df.csv: taxon mapping result
> 2. no_mapping_taxa_df.csv: unmapped taxa list
> 3. fun_abu_df_${anno_type}_${rfw_type}.csv: RFW abundance
> 4. fun_desc.csv: microbial function description
> 5. rel_abu.csv: microbial relative abundance data
> 6. pathway_abu.csv: microbial pathway abundance data if anno_type=='EC'
> 7. pathway_coverage.csv: microbial pathway coverage data if anno_type=='EC'. This file records microbial encoded or "covered" enzyme counts info for each pathway.
> 8. taxon_fun/*.tsv: microbial encoded function for each taxon
> 9. taxon_contri/*.tsv: sub-community info for each microbial function 

### Output of DFSCA-BC

> 1. ${prefix}_QMD_taxon_Abu_diff.csv: estimated taxon absoulte abundance changes between groups
> 2. ${prefix}_DFSCA_BC_diff.csv: estimated microbial function absoulte abundance changes between groups
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
