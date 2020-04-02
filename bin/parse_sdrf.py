#!/usr/bin/env python3
import re
import sys
import logging
import os

import click
import pandas as pd

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
  """This is the main tool that give access to all commands and options provided by the mslookup and dleamse algorithms"""


def openms_ify_mods(sdrf_mods):
  oms_mods = list()

  for m in sdrf_mods:
    if "AC=UNIMOD" not in m:
      raise Exception("only UNIMOD modifications supported.")

    name = re.search("NT=(.+?)(;|$)", m).group(1)
    name = name.capitalize()

    # workaround for missing PP in some sdrf TODO: fix in sdrf spec?
    if re.search("PP=(.+?)[;$]", m) is None:
      pp = "Anywhere"
    else:
      pp = re.search("PP=(.+?)(;|$)", m).group(
        1)  # one of [Anywhere, Protein N-term, Protein C-term, Any N-term, Any C-term

    if re.search("TA=(.+?)(;|$)", m) is None:  # TODO: missing in sdrf.
      print("Warning no TA= specified. Setting to N-term or C-term if possible: " + m)
      if "C-term" in pp:
        ta = "C-term"
      elif "N-term" in pp:
        ta = "N-term"
      else:
        print("Reassignment not possible. skipping")
        pass
    else:
      ta = re.search("TA=(.+?)(;|$)", m).group(1)  # target amino-acid
    aa = ta.split(",")  # multiply target site e.g., S,T,Y including potentially termini "C-term"

    if pp == "Protein N-term" or pp == "Protein C-term":
      for a in aa:
        if a == "C-term" or a == "N-term":  # no site specificity
          oms_mods.append(name + " (" + pp + ")")  # any Protein N/C-term
        else:
          oms_mods.append(name + " (" + pp + " " + a + ")")  # specific Protein N/C-term
    elif pp == "Any N-term" or pp == "Any C-term":
      pp = pp.replace("Any ", "")  # in OpenMS we just use N-term and C-term
      for a in aa:
        if a == "C-term" or aa == "N-term":  # no site specificity
          oms_mods.append(name + " (" + pp + ")")  # any N/C-term
        else:
          oms_mods.append(name + " (" + pp + " " + a + ")")  # specific N/C-term
    else:  # Anywhere in the peptide
      for a in aa:
        oms_mods.append(name + " (" + a + ")")  # specific site in peptide

  return ",".join(oms_mods)


def openms_convert(sdrf_file: str = None):
  sdrf = pd.read_table(sdrf_file)
  sdrf.columns = map(str.lower, sdrf.columns)  # convert column names to lower-case

  # map filename to tuple of [fixed, variable] mods
  mod_cols = [c for ind, c in enumerate(sdrf) if
              c.startswith('comment[modification parameters')]  # columns with modification parameters

  # get factor columns (except constant ones)
  factor_cols = [c for ind, c in enumerate(sdrf) if c.startswith('factor value[') and len(sdrf[c].unique()) > 1]

  file2mods = dict()
  file2pctol = dict()
  file2pctolunit = dict()
  file2fragtol = dict()
  file2fragtolunit = dict()
  file2diss = dict()
  file2enzyme = dict()
  file2fraction = dict()
  file2label = dict()
  file2source = dict()
  source_name_list = list()
  file2combined_factors = dict()
  file2technical_rep = dict()
  for index, row in sdrf.iterrows():
    ## extract mods
    all_mods = list(row[mod_cols])
    var_mods = [m for m in all_mods if 'MT=variable' in m or 'MT=Variable' in m]  # workaround for capitalization
    var_mods.sort()
    fixed_mods = [m for m in all_mods if 'MT=fixed' in m or 'MT=Fixed' in m]  # workaround for capitalization
    fixed_mods.sort()
    raw = row['comment[data file]']
    fixed_mods_string = ""
    if fixed_mods is not None:
      fixed_mods_string = openms_ify_mods(fixed_mods)

    variable_mods_string = ""
    if var_mods is not None:
      variable_mods_string = openms_ify_mods(var_mods)

    file2mods[raw] = (fixed_mods_string, variable_mods_string)

    source_name = row['source name']
    file2source[raw] = source_name
    if not source_name in source_name_list:
      source_name_list.append(source_name)

    if 'comment[precursor mass tolerance]' in row:
      pc_tol_str = row['comment[precursor mass tolerance]']
      if "ppm" in pc_tol_str or "Da" in pc_tol_str:
        pc_tmp = pc_tol_str.split(" ")
        file2pctol[raw] = pc_tmp[0]
        file2pctolunit[raw] = pc_tmp[1]
      else:
        print("Invalid precursor mass tolerance set. Assuming 10 ppm.")
        file2pctol[raw] = "10"
        file2pctolunit[raw] = "ppm"
    else:
      print("No precursor mass tolerance set. Assuming 10 ppm.")
      file2pctol[raw] = "10"
      file2pctolunit[raw] = "ppm"

    if 'comment[fragment mass tolerance]' in row:
      f_tol_str = row['comment[fragment mass tolerance]']
      f_tol_str.replace("PPM", "ppm")  # workaround
      if "ppm" in f_tol_str or "Da" in f_tol_str:
        f_tmp = f_tol_str.split(" ")
        file2fragtol[raw] = f_tmp[0]
        file2fragtolunit[raw] = f_tmp[1]
      else:
        print("Invalid fragment mass tolerance set. Assuming 10 ppm.")
        file2fragtol[raw] = "10"
        file2fragtolunit[raw] = "ppm"
    else:
      print("No fragment mass tolerance set. Assuming 10 ppm.")
      file2fragtol[raw] = "20"
      file2fragtolunit[raw] = "ppm"

    if 'comment[dissociation method]' in row:
      diss_method = row['comment[dissociation method]']
      file2diss[raw] = diss_method.toUpper()
    else:
      print("No dissociation method provided. Assuming HCD.")
      file2diss[raw] = 'HCD'

    if 'comment[technical replicate]' in row:
      file2technical_rep[raw] = str(row['comment[technical replicate]'])
    else:
      file2technical_rep[raw] = "1"

    enzyme = re.search("NT=(.+?)(;|$)", row['comment[cleavage agent details]']).group(1)
    enzyme = enzyme.capitalize()
    if "Trypsin/p" in enzyme:  # workaround
      enzyme = "Trypsin/P"
    file2enzyme[raw] = enzyme
    file2fraction[raw] = str(row['comment[fraction identifier]'])
    label = re.search("NT=(.+?)(;|$)", row['comment[label]']).group(1)
    file2label[raw] = label

    ## extract factors
    all_factors = list(row[factor_cols])
    combined_factors = "|".join(all_factors)
    if combined_factors == "":
      print("No factors specified. Adding dummy factor used as condition.")
      combined_factors = "none"

    file2combined_factors[raw] = combined_factors
    # print(combined_factors)

  ##################### only label-free supported right now

  # output of search settings
  f = open("openms.tsv", "w+")
  open_ms_search_settings_header = ["Filename", "FixedModifications", "VariableModifications", "Label",
                                    "PrecursorMassTolerance", "PrecursorMassToleranceUnit", "FragmentMassTolerance",
                                    "FragmentMassToleranceUnit", "DissociationMethod", "Enzyme"]
  f.write("\t".join(open_ms_search_settings_header) + "\n")
  for index, row in sdrf.iterrows():  # does only work for label-free not for multiplexed. TODO
    raw = row["comment[data file]"]
    f.write(raw + "\t" + file2mods[raw][0] + "\t" + file2mods[raw][1] + "\t" + file2label[raw] + "\t" + file2pctol[
      raw] + "\t" + file2pctolunit[raw] + "\t" + file2fragtol[raw] + "\t" + file2fragtolunit[raw] + "\t" + file2diss[
              raw] + "\t" + file2enzyme[raw] + "\n")
  f.close()

  # output of experimental design
  f = open("experimental_design.tsv", "w+")
  open_ms_experimental_design_header = ["Fraction_Group", "Fraction", "Spectra_Filepath", "Label", "Sample",
                                        "MSstats_Condition", "MSstats_BioReplicate"]
  f.write("\t".join(open_ms_experimental_design_header) + "\n")
  for index, row in sdrf.iterrows():  # does only work for label-free not for multiplexed. TODO
    raw = row["comment[data file]"]
    fraction_group = str(1 + source_name_list.index(row["source name"]))  # extract fraction group
    sample = fraction_group
    if 'none' in file2combined_factors[raw]:
      # no factor defined use sample as condition
      condition = sample
    else:
      condition = file2combined_factors[raw]
    replicate = file2technical_rep[raw]
    label = file2label[raw]
    if "label free sample" in label:
      label = "1"
    f.write(fraction_group + "\t" + file2fraction[
      raw] + "\t" + raw + "\t" + label + "\t" + sample + "\t" + condition + "\t" + replicate + "\n")
  f.close()


@click.command('convert-openms', short_help='convert sdrf to openms file output')
@click.option('--sdrf', '-s', help='SDRF file')
@click.pass_context
def openms_from_sdrf(ctx, sdrf: str):
  if sdrf is None:
    help()
  openms_convert(sdrf)


cli.add_command(openms_from_sdrf)

if __name__ == "__main__":
  cli()
