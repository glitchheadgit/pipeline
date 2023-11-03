# Pipeline info

## Requirements

### GTDB-Tk

1. GTDB-Tk requires ~84G of external data that needs to be downloaded and unarchived:

```bash
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz # or mirror - https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
```
2. Link GTDB to the "gtdbtk_data" folder in pipeline:

```bash
ln -s /path/to/gtdbtk /path/to/pipeline/gtdbtk_data
```

### CheckM

1. CheckM DB requires ~1.4G of external data that needs to be downloaded and unarchived:

```bash
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xvzf checkm_data_2015_01_16.tar.gz
```
2. Link CheckM database to the "checkm_data" folder in pipeline:

```bash
ln -s path/to/checkm_data /path/to/pipeline/checkm_data
```

### Kaiju

1. Kaiju requires ~102G of external data that needs to be downloaded and unarchived:

```bash
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_progenomes_2023-05-25
tar xvzf kaiju_db_progenomes_2023-05-25
```
2. Link Kaiju db to the "kaiju_data" folder in pipeline:

```bash
ln -s path/to/kaiju_data /path/to/pipeline/kaiju_data
```

### Installing metaWRAP environment

**Install it in the pipeline directory!**

```bash
mamba create -p /path/to/pipeline/mw-env -c ursky metawrap-mg=1.3.2
```


## Preprocessing
```bash
nextflow run pipeline/preprocess.nf \
--reads 'pipeline/fq_check/*{1,2}.fq'
```

