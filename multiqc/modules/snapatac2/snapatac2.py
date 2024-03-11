""" MultiQC module to parse output from snapatac2 """


import logging
import os
import re
import csv

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

def get_frac(x) -> str:
    x = round(float(x)*100,2)
    return str(x)+'%'

def get_int(x) -> str:
    return str(int(x))

def get_float(x) -> str:
    return str(float(x))


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="snapATAC2",
            anchor="snapatac2",
            href="https://github.com/kaizhang/SnapATAC2",
            info="Single-cell epigenomics analysis tools.",
        )

        # Find and load any STAR reports
        self.data = dict()
        for f in self.find_log_files("snapatac2", filehandles=True):
            parsed_data = self.parse_log(f["f"])
            if parsed_data is not None:
                s_name = f['s_name']
                if s_name in self.data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f)
                self.data[s_name] = parsed_data

        if len(self.data) > 0:
            log.info(f"Found {len(self.data)} reports")
            # Write parsed report data to a file
            self.write_data_file(self.data, "multiqc_snapatac2")

            # Basic Stats Table
            self.general_stats_table()


    def parse_log(self, f):
        parsed_data = {}
        reader = csv.reader(f)
        for row in reader:
            parsed_data[row[0]] = row[1]
        return parsed_data
    
    def general_stats_table(self):
        headers = {
            "Raw Cell Number": {
                'title': 'N Cells',
                'description': 'Raw Cell number ',
                'format': get_int,
            },
            'Valid Cell Number': {
                'title': 'Valid Cells',
                'description': 'Cell with TSS enrichment score>5 and not a doublet',
                'format': get_int,
                "cond_formatting_rules": {
                    "pass": [{"gt": 500}],
                    "warn": [{"lt": 500}, {"gt": 20000}],
                },
            },
            'Median TSS Enrichment Score': {
                'title': 'TSS score',
                'description': "Median TSS enrichment score of all cells",
                'format': get_float,
                "cond_formatting_rules": {
                    "pass": [{"gt": 5}],
                    "fail": [{"lt": 5}],
                },
            },
        }
        self.general_stats_addcols(self.data, headers=headers)

