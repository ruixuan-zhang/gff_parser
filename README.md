# GFF Parser

This script basically did two things 

1. parse the input gff <- refer to [gff3_parser](https://github.com/McClain-Thiel/gff3_parser/tree/main)
2. extract sequences surrounding start codon 

Usage 

```{python}
run_all(genome_path, gff_path, left_shift, right_shfit, out_folder, "seq_name")
```

This script will create 3 subfolders named `cds` , `n_term_long`, and `after_start_long` in the specified `out_folder` 

sequences `cds`, `[left: right]` and `[start_codon, right_shift]` will be saved into corresponding folder with the name specified in the last argument. 