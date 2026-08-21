
## [jumble_count](https://github.com/ClinSeq/jumble)
Create read count files from bam files used for Jumble panel of normal creation.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__jumble__jumble_count#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__jumble__jumble_count#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__jumble_count#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__jumble_count#

---

## [jumble_reference](https://github.com/ClinSeq/jumble)
Creates a panel of normal for the CNV caller Jumble.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__jumble__jumble_reference#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__jumble__jumble_reference#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__jumble_reference#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__jumble_reference#

---

## [ichorcna_offtarget_read_counter](https://github.com/shahcompbio/hmmcopy_utils)
Create per-sample read count wig files from bam files, used for ichorCNA off-target panel of normals creation.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__ichorcna_offtarget__ichorcna_offtarget_read_counter#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__ichorcna_offtarget__ichorcna_offtarget_read_counter#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__ichorcna_offtarget_read_counter#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__ichorcna_offtarget_read_counter#

---

## [ichorcna_offtarget_wig_list](https://github.com/genomic-medicine-sweden/cnv_sv)
Make a wig file list used as input for the ichorCNA off-target panel of normals.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__ichorcna_offtarget__ichorcna_offtarget_wig_list#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__ichorcna_offtarget__ichorcna_offtarget_wig_list#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__ichorcna_offtarget_wig_list#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__ichorcna_offtarget_wig_list#

---

## [ichorcna_offtarget_panel_of_normals](https://github.com/broadinstitute/ichorCNA)
Creates a panel of normals for ichorCNA off-target tumor fraction estimation from a cohort of normal samples' read count wig files.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__ichorcna_offtarget__ichorcna_offtarget_panel_of_normals#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__ichorcna_offtarget__ichorcna_offtarget_panel_of_normals#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__ichorcna_offtarget_panel_of_normals#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__ichorcna_offtarget_panel_of_normals#

---

## [deepsomatic_deepsomatic_pon](https://github.com/google/deepsomatic)
Calls small variants from normal BAM files

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deepsomatic__deepsomatic_deepsomatic_pon#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deepsomatic__deepsomatic_deepsomatic_pon#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deepsomatic_deepsomatic_pon#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deepsomatic_deepsomatic_pon#
