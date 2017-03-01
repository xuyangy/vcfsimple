#!/bin/env python

"""
Creates Excel/tsv reports from a vcf file. Valid genes are filtered using
a transcripts file from a genepanel.

It uses VcfReader for the data parsing/converint,
picking out the fields to include in the report.
"""

import argparse
import logging
import re
from collections import OrderedDict

from report.interpretation.vcfreader import VcfReader

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('ExomeReport')


class ExomeReport(object):

    FILTER_CRITERIAS = [
        ('inDB_alleleFreq', 0.05),
        ('ExAC_TOT', 0.01),
        ('1000g', 0.01)
    ]

    HGVS_INTRON_RE = re.compile('.*c\.[0-9]+?([\-\+])([0-9]+)')

    def __init__(self, vcf_path, transcripts_path):
        self.vcf_path = vcf_path
        self.transcripts_path = transcripts_path
        self.excel_writer = None
        self.tsv_fd = None

    def check_hgvs_intronic(self, hgvs):
        """
        Checks whether a variant is +6/-20 bp into intronic region.
        """
        match = re.match(ExomeReport.HGVS_INTRON_RE, hgvs)
        criterias = {
            '+': 6,
            '-': 20
        }
        if match:
            sign, num = match.groups()
            num = int(num)
            return sign in criterias and num > criterias[sign]
        return False

    def get_rows(self):

        reader = VcfReader(self.vcf_path, self.transcripts_path)

        empty_fields = [
            'IGV',
            'ACMG kriterier',
            'Klasse',
            'Resultatvurdering kommentar',
            'Kontroll kommentar',
            'Variantvurdering'
        ]

        rename_fields = {
            'inDB_noMutInd': 'inDB_noMutInd'
        }

        # Fields to include and their order
        fields = [
            'CHR',
            'POS',
            'REF',
            'ALT',
            'Gene',
            'Inheritance',
            'HGVSc',
            'HGVSp',
            'Consequence',
            'IGV',
            'ACMG kriterier',
            'Klasse',
            'Resultatvurdering kommentar',
            'Kontroll kommentar',
            'Variantvurdering',
            'Genotype',
            'ExAC_URL',
            'HGMD_URL',
            'AlamutColumn',
            'inDB_alleleFreq',
            'inDB_noMutInd',
            'inDB_indications',
            'sanger_verify',
            'VCF_FILTER',
            'VCF_QUAL',
            'VCF_AD_AlleleDepth',
            'VCF_DP_depth',
            'ExAC_TOT',
            'ExAC_NFE',
            'ExAC_FIN',
            'ExAC_AFR',
            'ExAC_AMR',
            'ExAC_EAS',
            'ExAC_SAS',
            'ExAC_OTH',
            '1000g',
            'inDB_genotypeFreq',
            'inDB_noTotal',
            'inDB_filter',
            'dbSNP',
            'repeatMasker',
        ]

        for entry in reader.get_data():
            inDB_noTotal = entry.get('inDB_noTotal', 'N/A')
            try:
                inDB_noTotal = int(inDB_noTotal)
                rename_fields['inDB_noMutInd'] = 'inDB_noMutInd ({})'.format(inDB_noTotal)
                break
            except ValueError:
                pass

        report_row = list()
        for entry in reader.get_data():
            for field in fields:
                # Rename field if in map
                fieldname = rename_fields.get(field, field)

                # If intronic, mark Klasse as 'u'
                if field == 'Klasse':
                    if self.check_hgvs_intronic(entry['HGVSc']):
                        report_row.append((fieldname, 'u'))
                    else:
                        report_row.append((fieldname, ''))
                elif field in empty_fields:
                    report_row.append((fieldname, ''))
                else:
                    report_row.append((fieldname, entry[field]))
            yield OrderedDict(report_row)

    @staticmethod
    def _filter_line(line):

        criteria_check = dict()
        for name, threshold in ExomeReport.FILTER_CRITERIAS:
            val = line[name]
            if val == 'N/A':
                criteria_check[name] = False
                continue
            try:
                val = float(val)
            except ValueError:
                log.debug("Value {} is not float for criteria {}".format(val, name))
                criteria_check[name] = False
                continue
            criteria_check[name] = val > threshold
        return any(criteria_check.values())

    def write(self, excel=None, tsv=None, igv=None, header=None):
        if not any(a for a in [excel, tsv, igv]):
            raise RuntimeError("You must at least provide one output type.")

        lines = sorted(self.get_rows(), key=lambda x: (x['Gene'],
                                                       x['Consequence'],
                                                       x['ExAC_TOT']))

        filtered_lines = [line for line in lines if not ExomeReport._filter_line(line)]

        if excel:
            # Lazy import to avoid unneccessary dependencies
            from report.interpretation.excelwriter import ExcelWriter, Worksheet, SectionedWorksheet
            from report.interpretation.custom_excel import VariantHelpWorksheet, VariantRowWriter, ChecklistWriter, VersionWorksheet

            self.excel_writer = ExcelWriter(excel)

            all_ws = Worksheet('All variants', freeze_panes=(1, 6))
            all_ws.set_writer(VariantRowWriter(lines))
            self.excel_writer.add_worksheet(all_ws)

            filtered_ws = SectionedWorksheet('Filtered variants', active=True, freeze_panes=(1, 6))
            filtered_ws.add_section(header, VariantRowWriter(filtered_lines, color_interpretation=True))
            filtered_ws.add_section(None, ChecklistWriter(offset=6))

            self.excel_writer.add_worksheet(filtered_ws)
            self.excel_writer.add_worksheet(VariantHelpWorksheet("Help"))
            self.excel_writer.add_worksheet(VersionWorksheet("Version"))

            self.excel_writer.save()

        if tsv:
            self.tsv_fd = open(tsv, 'w')
            self.tsv_fd.write('\t'.join(lines[0].keys()))
            self.tsv_fd.write('\n')
            for line in lines:
                if self.tsv_fd:
                    self.tsv_fd.write('\t'.join([str(f) for f in line.values()]))
                    self.tsv_fd.write('\n')
            self.tsv_fd.close()

        if igv:
            igv_data = list()
            for line in filtered_lines:
                igv_data.append(
                    [
                        line['CHR'],
                        line['POS']-1,
                        line['POS'],
                        line['Gene'],
                    ]
                )

            with open(igv, 'w') as igv_fd:
                for line in igv_data:
                    igv_fd.write('\t'.join(str(l) for l in line) + '\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Creates report tsv file from VEP annotation for use in HTS-interpretation.')
    parser.add_argument('-i', required=True, help='VCF file', dest='input')
    parser.add_argument('--transcripts', required=True, help='Transcript file from genepanel', dest='transcripts')
    parser.add_argument('--excel', required=False, help='Output to excel file name', dest='excel')
    parser.add_argument('--tsv', required=False, help='Output to tsv file name', dest='tsv')
    parser.add_argument('--igv', required=False, help='Output to igv file name', dest='igv')
    parser.add_argument('--header', required=False, help='Header to display in excel report (must be single word)', dest='header')

    args = parser.parse_args()

    ssr = ExomeReport(
        args.input,
        args.transcripts
    )

    ssr.write(excel=args.excel, tsv=args.tsv, igv=args.igv, header=args.header)
