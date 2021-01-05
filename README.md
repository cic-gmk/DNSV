### DNSV
A tool for potential de novo variation detections from structural variations.
### Quick start

Download the files and set the appropriate work path.

The vcf/csv SVs files of a trio are needed. The details of data format are referred in docs.

python DNSV.py father.vcf mother.vcf son.vcf -o dnsv.csv --statistics True


### Parameters

```
Usage: python DNSV.py [options] father.vcf mother.vcf son.vcf -o dnsv.csv

positional arguments:
  father_SVs            Input father's SVs
  mother_SVs            Input mother's SVs
  son_SVs               Input son's SVs

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output Name (son_sv_DNSV.csv)
  --statistics STATISTICS
                        Whether the statistics of DNSVs is needed (False)

Comparison Threshold Arguments:
  -r REFDIST, --refdist REFDIST
                        Max reference location distance (200bp)
  -t TYPEIGNORE, --typeignore TYPEIGNORE
                        Variant types don't need to match to compare (False)
  -O OVERLAP_RATE, --overlap_rate OVERLAP_RATE
                        Reciprocal overlaps with the reference SVs (0.5)

Filtering Arguments:
  -p PRECISIONLIMIT, --precisionlimit PRECISIONLIMIT
                        Limit the precision description in INFO column
  -s SIZEMIN, --sizemin SIZEMIN
                        Minimum DNSV size (50bp)
  -l TYPELIMIT, --typelimit TYPELIMIT
                        Limit DNSV types in ['INS','DEL','DUP','INV']
```
