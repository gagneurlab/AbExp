import os
import sys
import glob
import yaml
import json
import pathlib

from typing import Dict, Callable
from snakemake.io import glob_wildcards
import re
import numpy

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# from snakemk_util import recursive_format
def recursive_format(data, params, fail_on_unknown=False):
    if isinstance(data, str):
        return data.format(**params)
    elif isinstance(data, dict):
        return {k: recursive_format(v, params) for k, v in data.items()}
    elif isinstance(data, list):
        return [recursive_format(v, params) for v in data]
    else:
        if fail_on_unknown:
            raise ValueError("Handling of data type not implemented: %s" % type(data))
        else:
            return data


class SafeDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'


def glob_output(
    input_pattern: str,
    output_pattern: str,
    wildcard_constraint_regex: Dict[str, str]=None
) -> Callable:
    """
    Glob for wildcards in a given input_pattern and formats some output pattern with the detected wildcards.
    Generates an input function for Snakemake rules that globs for wildcards in a given input_pattern
    and formats some output pattern with the detected wildcards.
    Allows to restrict wildcards with a dictionary of wildcard constraints.
    
    :param input_pattern: string with some wildcards to glob for, e.g. 'my_input/{vcf_file}'
    :param output_pattern: string with some wildcards which should be formatted, e.g. 'my_output/{vcf_file}.gz'
    :param wildcard_constraint_regex: dictionary of wildcard constraints, e.g.: `{"vcf_file": ".+\.vcf$"}`
    :returns: callable taking wildcards as first argument; can be used as input function to snakemake rules
    """
    if wildcard_constraint_regex is None:
        wildcard_constraint_regex = {}
    def input_fn(
        wildcards,
        input_pattern=input_pattern,
        output_pattern=output_pattern,
        wildcard_constraint_regex=wildcard_constraint_regex,
    ):
        # print("begin globbing...")
        input_pattern = input_pattern.format_map(
            SafeDict(**wildcards)
        )
        output_pattern = output_pattern.format_map(
            SafeDict(**wildcards)
        )
        # print(f"input_pattern: {input_pattern}")
        # print(f"output_pattern: {output_pattern}")

        matched_wildcards = glob_wildcards(input_pattern)._asdict()
        # ensure that each wildcard matches the constraints:
        mask = None
        for k, constraint in wildcard_constraint_regex.items():
            if k not in matched_wildcards:
                continue

            pattern = re.compile(constraint)
            wildcard_list = matched_wildcards[k]
            if mask is None:
                mask = np.ones(shape=(len(wildcard_list), ), dtype="bool")

            for idx, w in enumerate(wildcard_list):
                if not pattern.match(w):
                    # print(f"'{w}' does not match {constraint}")
                    mask[idx] = False

        if mask is not None:
            filtered_wildcards = {}
            for k, wildcard_list in matched_wildcards.items():
                wildcard_list = np.asarray(wildcard_list, dtype=object)[mask]
                filtered_wildcards[k] = wildcard_list
        else:
            filtered_wildcards = matched_wildcards
        # print(filtered_wildcards)

        # compute output files
        output_files = expand(output_pattern, zip, **filtered_wildcards)

        return output_files
    
    return input_fn


checkpoint file_depends:
    output: 
        file=touch("{file}.depends"),
    input:
        file="{file}",
    shell:
        "echo '{wildcards.file}'"


def require(file):
    """
    Makes sure that a certain input file exists before the rest of the rule is being executed.
    Technically, returns a function that takes wildcards as arguments.
    The resulting dependency adds ".depends" to the file path.
    Usage example:
        ```
        def read_config():
            [...]
        
        rule test:
            input:
                config_depends = require("config.yaml") # will effectively require the file "config.yaml.depends"
                config = lambda wildcards: read_config("config.yaml")
        ```
    This script does not fail, as `config_depends` has to exist before `config` is being executed.
    """
    return lambda wildcards: checkpoints.file_depends.get(file=file).output.file


def read_yaml_input(file, wildcards):
    """
    Reads the "input" section from a config yaml
    """
    import yaml
    with open(file, "r") as fd:
        config = yaml.safe_load(fd)

    retval = dict()
    if "snakemake" in config:
        if "input" in config["snakemake"]:
            retval = config["snakemake"]["input"]
            retval = recursive_format(retval, wildcards)
            
            # make sure that the output is a dictionary
            if isinstance(retval, list):
                retval = dict(yaml_input=retval)
    
    return retval


def require_yaml_input(file):
    """
    Reads the "input" section from a config yaml and adds it as additional input rules.
    Also, makes sure that the config file exists beforehand.
    See also `require`
    """
    def retval(wildcards):
        file_fmt = file.format(**wildcards)
        # make sure that yaml file exists:
        config_depends = checkpoints.file_depends.get(file=file_fmt).output

        return read_yaml_input(file_fmt, wildcards)
    return retval



localrules: file_depends
