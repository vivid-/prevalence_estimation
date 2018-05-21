# prevalence_estimation
Complete version of prevalence estimation

## Downlad pathogenicity annotation from mutation database
For EGL database, the python script `extract_EGL.py` would be used to fetch information from the website.
``` bash
python ./fetch_database/extract_EGL.sh --gene [GENE] --out [OUT]
## for example: 
python ./fetch_database/extract_EGL.sh --gene DYSF --out ../DYSF_EGL.txt
```
Another python script `convert_HGVS.py` would be furtherly used to extract genomic coordinates for annotated variants in the EGL database.
```bash
python ./fetch_database/convert_HGVS.py --inpu [INPU] --out [OUT]
## for example:
python ./fetch_database/convert_HGVS.py --inpu ../DYSF_EGL.txt --out ../DYSF_EGL_loc.txt
```

For ClinVar database, we 
### 
