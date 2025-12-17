# :sunflower: flowR: genomic environment visualization tool (R-based)

## :bulb: Bash Quick Examples

Example 1: Basic Plot

Basic plot using GFF file and YAML configuration and specifying name of gene.

```{bash, eval = FALSE}
Rscript flowR.R --gff_file data.gff --yaml config.yaml --gene gene_name --plot output_plot.png
```

```{bash, eval = FALSE}
Rscript flowR.R -g </path/of/file.gff> -y <config.yaml> -g <gene_name> -p <output_plot.png>
```

<!-- If you have any genome data, put this data in TSV or CSV format and use GFF to map the chromosome.

Example 2: Fill data to Plot

```{bash, eval = FALSE}
Rscript chromR.R --gff_file data.gff --fill_file data_to_fill_information.tsv/.csv --format tsv/csv
```

```{bash, eval = FALSE}
Rscript chromR.R -g </path/of/file.gff> -ff <path/of/data_to_fill> -f <tsv/csv>
``` -->