# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("6.8.0")

### Set and validate config file


configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type"], drop=False)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


def get_units2(units: pandas.DataFrame, wildcards: snakemake.io.Wildcards, type: str = None) -> pandas.DataFrame:
    """
    function used to extract one or more units from units.tsv
    Args:
        units: DataFrame generate by importing a file following schema defintion
               found in pre-alignment/workflow/schemas/units.schema.tsv
        wildcards: wildcards object with at least the following wildcard names
               sample and type (optional, can also be passed as an argument).
        type: N,T or R
    Returns:
        all units from the DataFrame that can be filtereted out using sample name
        and unit type (N,T,R)
    Raises:
        raises an exception (KeyError) if no unit(s) can be extracted from the Dataframe
    """
    if type is None:
        files = units.loc[(wildcards.sample, wildcards.type)].dropna()
    else:
        file = units.loc[(wildcards.sample, types)].dropna()
    if files is not None:
        if isinstance(files, pandas.Series):
            files = pandas.DataFrame([[f[1] for f in files.iteritems()], ], columns=[f[0] for f in files.iteritems()]).set_index(units.index.names)
    return [file for file in files.itertuples()]


def get_units3(units: pandas.DataFrame, wildcards: snakemake.io.Wildcards, type: str = None) -> pandas.DataFrame:
    """
    function used to extract one or more units from units.tsv
    Args:
        units: DataFrame generate by importing a file following schema defintion
               found in pre-alignment/workflow/schemas/units.schema.tsv
        wildcards: wildcards object with at least the following wildcard names
               sample and type (optional, can also be passed as an argument).
        type: N,T or R
    Returns:
        all units from the DataFrame that can be filtereted out using sample name
        and unit type (N,T,R)
    Raises:
        raises an exception (KeyError) if no unit(s) can be extracted from the Dataframe
    """
    if type is None:
        files = units.loc[(wildcards.sample, wildcards.type)].dropna()
    else:
        files = units.loc[(wildcards.sample, type)].dropna()
    return [file for file in files.iteritems()]


def get_bams(units: pandas.DataFrame) -> typing.List[str]:
    """
    function used to extract all bam files found in units.tsv
    Args:
        units: DataFrame generate by importing a file following schema defintion
               found in pre-alignment/workflow/schemas/units.schema.tsv
    Returns:
        List of strings with all bam file names and path
    """
    return [unit.bam for unit in units.itertuples()]


print(type(units))
#print(get_units(units, snakemake.io.Wildcards(fromdict={'sample': 'HD832.HES45', "type": "T"})))
files = units.loc[("HD832.HES45", "T")].dropna()
print(type(files))
#print(get_units2(units, snakemake.io.Wildcards(fromdict={'sample': 'HD832.HES45', "type": "T"}))[0].bam)


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def compile_output_list(wildcards):
    return [
        "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",
        "references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5",
    ]
