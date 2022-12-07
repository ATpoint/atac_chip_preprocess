# BLACKLISTS

We provide both the ENCODE NGS blacklist as well as the merge of the ENCODE blacklist and the ATAC-seq blacklist for mitochondrial
homologs in the genome.

## Mouse

```bash

mkdir mm10

# ENCODE:
wget https://github.com/Boyle-Lab/Blacklist/archive/v2.0.tar.gz && tar zxf v2.0.tar.gz
gzip -c -d Blacklist-2.0/lists/mm10-blacklist.v2.bed.gz | cut -f1-3  > ./mm10/mm10_encode_blacklist_v2.bed

# Mitochondrial
curl -o - -s https://raw.githubusercontent.com/buenrostrolab/mitoblacklist/master/peaks/mm10_peaks.narrowPeak \
| cut -f1-3 \
| tee ./mm10/mm10_mito_blacklist.bed \
| cat - ./mm10/mm10_encode_blacklist_v2.bed \
| cut -f1-3 \
| sort -k1,1 -k2,2n \
| bedtools merge -i - > ./mm10/mm10_combined_blacklist.bed

```

## Human

```bash

mkdir hg38

# ENCODE:
gzip -c -d Blacklist-2.0/lists/hg38-blacklist.v2.bed.gz | cut -f1-3 > ./hg38/hg38_encode_blacklist_v2.bed

# Mitochondrial
curl -o - -s https://raw.githubusercontent.com/buenrostrolab/mitoblacklist/master/peaks/hg38_peaks.narrowPeak \
| cut -f1-3 \
| tee ./hg38/hg38_mito_blacklist.bed \
| cat - ./hg38/hg38_encode_blacklist_v2.bed \
| cut -f1-3 \
| sort -k1,1 -k2,2n \
| bedtools merge -i - > ./hg38/hg38_combined_blacklist_v2.bed

rm -r Blacklist-2.0 v2.0.tar.gz

```
