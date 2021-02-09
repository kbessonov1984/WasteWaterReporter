# Waste Water Reporter

This simple reporter tool will generate reports on Variants of Concern
from mixed waste water samples


# Requirements
* pysam
* pandas
* python > 3.6

# Install
As any other Python library.

# Run
Get any BED file (e.g. `21OTT5Barrie_Articv3_combined-MN908947.3.trim.bam`)

# Usage
```bash
usage: Generate VOCs Waste Water report
 [-h] -i INPUT -o OUTPUT_NAME

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        BAM file(s)
  -o OUTPUT_NAME, --output_name OUTPUT_NAME
                        Output summary file name
```