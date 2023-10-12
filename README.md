# Pipeline info

## Preprocessing
```bash
nextflow run pipeline/preprocess.nf \
--reads 'pipeline/fq_check/*{1,2}.fq' \
-params-file pipeline/config/preprocessing_params.json
```
