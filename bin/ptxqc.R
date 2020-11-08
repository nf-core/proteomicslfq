#!/usr/bin/env Rscript

options(encoding = 'UTF-8')

args = commandArgs(trailingOnly = TRUE)
require(PTXQC)
mztab_file = args[1]
fn = PTXQC:::getReportFilenames(dirname(mztab_file), TRUE, mztab_file)

yaml_obj = list()
if (length(args) >= 2) {
  if (file.exists(args[2])) {
  	yaml_obj = yaml::yaml.load_file(args[2])
  } 
}
createReport(txt_folder = NULL, mztab_file = mztab_file,  yaml_obj)

