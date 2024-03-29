#!/usr/bin/env python3
#generate_versions.py
#*
#* --------------------------------------------------------------------------
#* Licensed under MIT (https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/blob/14a1c222e53f59391d96a2a2e1fd4995474c0d15/LICENSE)
#* --------------------------------------------------------------------------
#*

'''Make YAML of software versions.'''

from __future__ import print_function
from collections import OrderedDict
import re
import os
import logging
import glob
import argparse
import numpy as np

EPILOG = '''
For more details:
        %(prog)s --help
'''

# SETTINGS

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = False
logger.setLevel(logging.INFO)

SOFTWARE_REGEX = {
    'Python': ['version_python.txt', r"Python (\S+)"],
    'DERIVA': ['version_deriva.txt', r"(\S+)"],
    'BDBag': ['version_bdbag.txt', r"BDBag (\S+) \(Bagit \S+\)"],
    'Trim Galore!': ['version_trimgalore.txt', r"version (\S+)"],
    'HISAT2': ['version_hisat2.txt', r"version (\S+)"],
    'Samtools': ['version_samtools.txt', r"samtools (\S+)"],
    'picard (MarkDuplicates)': ['version_markdups.txt', r"Version:(\S+)"],
    'featureCounts': ['version_featurecounts.txt', r"featureCounts v(\S+)"],
    'deepTools': ['version_deeptools.txt', r"deeptools (\S+)"],
    'Seqtk': ['version_seqtk.txt', r"Version: (\S+)"],
    'R': ['version_r.txt', r"R version (\S+)"],
    'FastQC': ['version_fastqc.txt', r"FastQC v(\S+)"],
    'SeqWho': ['version_seqwho.txt', r"Version (\S+)"],
    'RSeQC': ['version_rseqc.txt', r"infer_experiment.py (\S+)"],
    'MultiQC': ['version_multiqc.txt', r"multiqc, version (\S+)"],
    'Pipeline Version': ['./nextflow.config', r"version = 'v(\S+)'"]
}


def get_args():
    '''Define arguments.'''

    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-o', '--output',
                        help="The out file name.",
                        required=True)

    parser.add_argument('-t', '--test',
                        help='Used for testing purposes',
                        default=False,
                        action='store_true')

    args = parser.parse_args()
    return args


def check_files(files, test):
    '''Check if version files are found.'''

    logger.info("Running file check.")

    software_files = np.array(list(SOFTWARE_REGEX.values()))[:,0]

    extra_files =  set(files) - set(software_files)

    if len(extra_files) > 0 and test:
            logger.error('Missing regex: %s', list(extra_files))
            raise Exception("Missing regex: %s" % list(extra_files))


def main():
    args = get_args()
    output = args.output
    test = args.test

    out_filename = output + '_mqc.yaml'

    results = OrderedDict()
    results['Python'] = '<span style="color:#999999;\">Not Run</span>'
    results['DERIVA'] = '<span style="color:#999999;\">Not Run</span>'
    results['BDBag'] = '<span style="color:#999999;\">Not Run</span>'
    results['Trim Galore!'] = '<span style="color:#999999;\">Not Run</span>'
    results['HISAT2'] = '<span style="color:#999999;\">Not Run</span>'
    results['Samtools'] = '<span style="color:#999999;\">Not Run</span>'
    results['picard (MarkDuplicates)'] = '<span style="color:#999999;\">Not Run</span>'
    results['featureCounts'] = '<span style="color:#999999;\">Not Run</span>'
    results['deepTools'] = '<span style="color:#999999;\">Not Run</span>'
    results['Seqtk'] = '<span style="color:#999999;\">Not Run</span>'
    results['R'] = '<span style="color:#999999;\">Not Run</span>'
    results['FastQC'] = '<span style="color:#999999;\">Not Run</span>'
    results['SeqWho'] = '<span style="color:#999999;\">Not Run</span>'
    results['RSeQC'] = '<span style="color:#999999;\">Not Run</span>'
    results['MultiQC'] = '<span style="color:#999999;\">Not Run</span>'
    results['Pipeline Version'] = '<span style="color:#999999;\">Not Run</span>'

    # list all files
    files = glob.glob('**/*.txt', recursive=True)

    # Check for version files:
    check_files(files, test)

    # Search each file using its regex
    for k, v in SOFTWARE_REGEX.items():
        if os.path.isfile(v[0]):
            with open(v[0]) as x:
                versions = x.read()
                match = re.search(v[1], versions)
                if match:
                    results[k] = "v{}".format(match.group(1))

    # Dump to YAML
    print(
        '''
        id: 'software_versions'
        section_name: 'Software Versions'
        section_href: 'https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/wikis/Pipeline/Tool-Versions'
        plot_type: 'html'
        description: 'are collected for pipeline version.'
        data: |
            <dl class="dl-horizontal">
        '''
    , file = open(out_filename, "w"))

    for k, v in results.items():
        print("            <dt>{}</dt><dd>{}</dd>".format(k, v), file = open(out_filename, "a"))
    print("            </dl>", file = open(out_filename, "a"))


if __name__ == '__main__':
    main()
