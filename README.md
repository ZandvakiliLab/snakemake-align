# snakemake-align

This is a fork of nicely designed [snakemake-simple-mapping](https://github.com/MPUSP/snakemake-simple-mapping) repository. This fork is focused on genome alignment and adds the following features

1. Download genomes using NCBI datasets based on genome assembly accession or Genbank id 
2. Align to an arbitrary number of genomes
3. Deduplication 
4. Support for cram files
5. Support scaling coverage files from spike-in genomes [TO DO]

## Authors

- Arya Zandvakili MD PhD
  - Affiliation: University of Cincinnati College of Medicine
  - ORCID profile: https://orcid.org/my-orcid?orcid=0000-0001-8031-8067
  - github page: https://github.com/aryazand

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.