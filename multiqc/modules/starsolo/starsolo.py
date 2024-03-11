""" MultiQC module to parse output from STARSolo """


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


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="STARSolo",
            anchor="starsolo",
            href="https://github.com/alexdobin/STAR",
            info="is an ultrafast universal RNA-seq aligner.",
            doi="10.1093/bioinformatics/bts635",
        )

        # Find and load any STAR reports
        self.starsolo_data = dict()
        for f in self.find_log_files("starsolo", filehandles=True):
            parsed_data = self.parse_starsolo_report(f["f"])
            if parsed_data is not None:
                solo_out_dir = os.path.dirname(f["root"])
                s_name = self.clean_s_name(os.path.basename(solo_out_dir), f, root=os.path.dirname(solo_out_dir))
                if s_name in self.starsolo_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(s_name=s_name, source=os.path.join(f["root"], f['fn']), section="summary")
                self.starsolo_data[s_name] = parsed_data

        if len(self.starsolo_data) > 0:
            log.info(f"Found {len(self.starsolo_data)} reports")
            # Write parsed report data to a file
            self.write_data_file(self.starsolo_data, "multiqc_starsolo")

            # Basic Stats Table
            self.general_stats_table()

            self.add_section(
                name='Barcode Rank plot',
                description='description here',
                helptext="""
                help here
                """
            )

    def parse_starsolo_report(self, f):
        parsed_data = {}
        reader = csv.reader(f)
        for row in reader:
            parsed_data[row[0]] = row[1]
        return parsed_data
    
    def general_stats_table(self):
        headers = {
            'Reads With Valid Barcodes': {
                'title': '% Valid Barcodes',
                'description': 'Fraction of reads with valid barcodes matching whitelist',
                'format': get_frac,
            },
            'Estimated Number of Cells': {
                'title': 'N Cells',
                'description': 'Estimated number of cells',
                'format': get_int,
            },
            'Fraction of Unique Reads in Cells': {
                'title': '% Reads in Cells',
                'description': 'Fraction of unique reads in cells',
                'format': get_frac,
            },
            'Median GeneFull_Ex50pAS per Cell': {
                'title': 'Median Genes',
                'description': 'Median genes per cell',
                'format': get_int,
            },
            'Mean Reads per Cell': {
                'title': 'Mean Reads',
                'description': 'Mean Reads per Cell',
                'format': get_int,
            },
            'Mean UMI per Cell': {
                'title': 'Mean UMI',
                'description': 'Mean UMI per Cell',
                'format': get_int,
            },
            'Sequencing Saturation': {
                'title': 'Saturation',
                'description': 'Sequencing Saturation',
                'format': get_frac,
            },
        }
        self.general_stats_addcols(self.starsolo_data, headers=headers)

