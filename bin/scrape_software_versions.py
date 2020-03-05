#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

openms_version_regex = r"[0-9][.][0-9][.][0-9]"

regexes = {
    'nf-core/proteomicslfq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'ThermorawfileParser': ['v_thermorawfileparser.txt', r"(\S+)"],
}
results = OrderedDict()
results['nf-core/proteomicslfq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['ThermorawfileParser'] = '<span style="color:#999999;\">N/A</span>'
results['FileConverter'] = '<span style="color:#999999;\">N/A</span>'
results['DecoyDatabase'] = '<span style="color:#999999;\">N/A</span>'
results['MSGFPlusAdapter'] = '<span style="color:#999999;\">N/A</span>'
results['MSGFPlus'] = '<span style="color:#999999;\">N/A</span>'
results['CometAdapter'] = '<span style="color:#999999;\">N/A</span>'
results['Comet'] = '<span style="color:#999999;\">N/A</span>'
results['PeptideIndexer'] = '<span style="color:#999999;\">N/A</span>'
results['PSMFeatureExtractor'] = '<span style="color:#999999;\">N/A</span>'
results['PercolatorAdapter'] = '<span style="color:#999999;\">N/A</span>'
results['Percolator'] = '<span style="color:#999999;\">N/A</span>'
results['IDFilter'] = '<span style="color:#999999;\">N/A</span>'
results['IDScoreSwitcher'] = '<span style="color:#999999;\">N/A</span>'
results['FalseDiscoveryRate'] = '<span style="color:#999999;\">N/A</span>'
results['IDPosteriorErrorProbability'] = '<span style="color:#999999;\">N/A</span>'
results['IDFilter'] = '<span style="color:#999999;\">N/A</span>'
results['ProteomicsLFQ'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/proteomicslfq Software Versions'
section_href: 'https://github.com/nf-core/proteomicslfq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
