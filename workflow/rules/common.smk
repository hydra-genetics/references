# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import os
import re
import sys
import typing

import pandas as pd
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics.utils.misc import get_input_aligned_bam

min_version("9.0.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "type"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


# Sentinel separating "no default given, so this entry is required" from a default of
# None, [] or "", each of which is a value a caller may legitimately want back.
_REQUIRED = object()


def get_config_value(*keys, default=_REQUIRED):
    """
    Fetch a value from the config, failing with a message that names the missing entry.

    Defaulting to "" is not usable here: an empty string reaches Snakemake either as
    a rule input, where it aborts with a MissingInputException that lists no file, or
    as a params value, where it silently produces a malformed shell command. Call this
    from an input/params function so the check stays lazy -- a workflow that never uses
    the rule does not have to configure it.

    Pass default=[] for an input file that the rule can run without. Snakemake reads an
    empty list as "no file", which is what "" was never able to express. Without a
    default the entry is required, and a missing or blank one raises.
    """
    value = config
    for i, key in enumerate(keys):
        if not isinstance(value, dict) or key not in value:
            if default is not _REQUIRED:
                return default
            missing = ":".join(keys[: i + 1])
            raise WorkflowError(f"references: missing config entry '{missing}', required by the rule being run")
        value = value[key]

    if not isinstance(value, str) or not value.strip():
        if default is not _REQUIRED:
            return default
        name = ":".join(keys)
        raise WorkflowError(f"references: config entry '{name}' must be a non-empty string, got {repr(value)}")

    return value


def design_bed_basename():
    """
    Basename of the design bed, used to name the PoN artefacts.

    Several rules build output, log and benchmark paths from this, and those must resolve
    at parse time, so this cannot raise the way get_config_value does -- a workflow that
    never builds a PoN must still be able to parse those rules. An unset design_bed
    therefore still yields "" here; config.schema.yaml constrains the value, and the
    matching input/params entries go through get_config_value and fail loudly.
    """
    return config.get("reference", {}).get("design_bed", "").split("/")[-1]


def get_bams(units: pd.DataFrame) -> typing.List[str]:
    """
    function used to extract all bam files found in units.tsv
    Args:
        units: DataFrame generate by importing a file following schema definition
               found in pre-alignment/workflow/schemas/units.schema.tsv
    Returns:
        List of strings with all bam file names and path
    """
    return get_units_column(units, "bam")


def get_bais(units: pd.DataFrame) -> typing.List[str]:
    """
    function used to extract all bam files found in units.tsv and add .bai to the filename
    Args:
        units: DataFrame generate by importing a file following schema definition
               found in pre-alignment/workflow/schemas/units.schema.tsv
    Returns:
        List of strings with all bam.bai file names and path
    """
    return [f"{bam_string}.bai" for bam_string in get_units_column(units, "bam")]


def get_counts(samples, units):
    count_list = [
        "references/jumble_count/%s_%s.bam.counts.RDS" % (sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    return count_list


def get_ichorcna_wigs(samples, units):
    wig_list = [
        "references/ichorcna_offtarget_read_counter/%s_%s.wig" % (sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    return wig_list


def get_coverage_files(samples, units):
    coverage_list = [
        "references/purecn_coverage/%s_%s_coverage_loess.txt.gz" % (sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    return coverage_list


def get_gvcfs(units: pd.DataFrame) -> typing.List[str]:
    """
    function used to extract all gvcf files found in units.tsv
    Args:
        units: DataFrame generate by importing a file following schema definition
               found in pre-alignment/workflow/schemas/units.schema.tsv
    Returns:
        List of strings with all gvcf file names and path
    """
    return get_units_column(units, "gvcf")


def get_vcfs(units: pd.DataFrame) -> typing.List[str]:
    """
    function used to extract all vcf files found in units.tsv
    Args:
        units: DataFrame generate by importing a file following schema definition
               found in pre-alignment/workflow/schemas/units.schema.tsv
    Returns:
        List of strings with all vcf file names and path
    """
    return get_units_column(units, "vcf")


def get_cnv_vcfs(units: pd.DataFrame) -> typing.List[str]:
    """
    function used to extract all cnv.vcf files found in units.tsv
    Args:
        units: DataFrame generate by importing a file following schema definition
               found in pre-alignment/workflow/schemas/units.schema.tsv
    Returns:
        List of strings with all cnv.vcf file names and path
    """
    return get_units_column(units, "cnv_vcf")


def get_units_column(units: pd.DataFrame, column: str) -> typing.List[str]:
    """
    extract a column from units.tsv
    Args:
        units: DataFrame generated by importing a file following schema definition
               found in workflow/schema/units.schema.yaml
        column: the name of the requested column
    Returns:
        List of strings representing the content of the requested column, or an
        empty list if the column does not exist, or if all values are missing/empty
    """
    if column not in units.columns:
        return []
    return list(set(units[column][units[column].notna()]))


wildcard_constraints:
    sample="|".join(re.escape(s) for s in samples.index),
    type="N|T|R",


# Output files commented out as they do not work in integration testing using small files
def compile_output_list(wildcards):
    return [
        "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",
        # "references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5",
        "references/exomedepth_reference/RefCount.Rdata",
        "references/msisensor_pro_baseline/Msisensor_pro_reference.list_baseline",
        "references/create_background_file/background_panel.tsv",
        "references/create_artifact_file/artifact_panel.tsv",
        "references/svdb_export/svdb_cnv.vcf",
        # "references/purecn_normal_db/output/normalDB.rds",
        # "references/purecn_normal_db/output/mapping_bias.rds",
        # "references/purecn_interval_file/targets_intervals.txt",
        # "references/jumble_reference/twist_DNA_solid.chr1.annotated.bed.reference.RDS",
        # "references/ichorcna_offtarget_panel_of_normals/ichorcna_offtarget_PoN_median.rds",
    ]
