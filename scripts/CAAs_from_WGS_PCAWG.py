#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import collections
import csv
import io
import math
import pathlib
import re

import pandas as pd


class Barcode:

    """PCAWG Barcodes
    """
    regex = re.compile(r'''(?P<project>\w+)-
                           (?P<tss>\d+)-
                           (?P<participant>\w+)-
                           (?P<sample>\w+)''',
                       re.VERBOSE)

    def __init__(self, barcode):
        self.barcode = barcode
        self.match = Barcode.regex.search(barcode)
        if self.match:
            self.groupdict = self.match.groupdict()
        else:
            raise ValueError(barcode)

    @property
    def project(self):
        return self.groupdict['project']

    @property
    def tss(self):
        return self.groupdict['tss']

    @property
    def participant(self):
        return self.groupdict['participant']

    @property
    def sample(self):
        """Return the sample code.

        Tumour types range from 01-09, normal types from 10-29 and control samples from 20-29.

        """
        return self.groupdict['sample']
 
class Centromere:

    @staticmethod
    def region(chrom):
        centromere_regions = {
            '1': (121236957, 123476957), '2': (91689898, 94689898),
            '3': (90587544, 93487544), '4': (49354874, 52354874),
            '5': (46441398, 49441398), '6': (58938125, 61938125),
            '7': (58058273, 61058273), '8': (43958052, 46958052),
            '9': (47107499, 50107499), '10': (39244941, 41624941),
            '11': (51450781, 54450781), '12': (34747961, 36142961),
            '13': (16000000, 17868000), '14': (15070000, 18070000),
            '15': (15260000, 18260000), '16': (35143302, 36943302),
            '17': (22187133, 22287133), '18': (15400898, 16764896),
            '19': (26923622, 29923622), '20': (26267569, 28033230),
            '21': (10260000, 13260000), '22': (11330000, 14330000)
        }
        return centromere_regions[chrom]


class Segment:
    def __init__(self, sample, chrom, start, end, num_probes, mean):
        # coordinates should be 0-based
        self.sample = sample
        self.chrom = chrom
        self.start = int(float(start))
        self.end = int(float(end))
        self.num_probes = int(float(num_probes))   # some firehose files have 1e+05 which will cause int to fail
        self.mean = float(mean)

    def __len__(self):
        #  return self.end - (self.start - 1)
        return self.end - self.start


class SegmentFile:

    @staticmethod
    def parse(file_name):
        with open(file_name, 'rt') as in_handle:
            reader = csv.DictReader(in_handle, delimiter=',')
            for row in reader:
                segment = Segment(row['Sample'], row['Chromosome'], row['Start'], row['End'],
                                  row['Num_Probes'], row['Segment_Mean'])
                if segment.chrom not in ['X', 'Y', '23', '24']:
                    yield segment


class Chromosomes:

    #  in_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    in_order = [str(i) for i in range(1, 23)]


def summarise_sample(segments):
    seg_ns = {}
    seg_lens = {}
    seg_means = {}
    arm_lengths = {}

    # Initialise all data structures.

    for chrom in Chromosomes.in_order:
        for arm in ['p', 'q']:
            arm_lengths[(chrom, arm)] = 0
            for direction in ['amp', 'del']:
                seg_ns[(chrom, arm, direction)] = 0
                seg_lens[(chrom, arm, direction)] = 0
                seg_means[(chrom, arm, direction)] = 0

    for segment in segments:
        if segment.chrom in ['X', 'Y']:
            continue

        centromere_start, centromere_end = Centromere.region(segment.chrom)

        if segment.end < centromere_start:
            arm = 'p'
            arm_lengths[(segment.chrom, arm)] += len(segment)
        elif segment.start > centromere_end:
            arm = 'q'
            arm_lengths[(segment.chrom, arm)] += len(segment)
        else:
            # Segment intersects centromere, skip
            continue

        if segment.mean > 0.2:
            direction = 'amp'
        elif segment.mean < -0.2:
            direction = 'del'
        else:
            # This segment hasn't reached the required threshold, skip.
            continue

        key = (segment.chrom, arm, direction)
        seg_ns[key] += 1
        seg_lens[key] += len(segment)
        seg_means[key] += segment.mean

    return (seg_ns, seg_lens, seg_means, arm_lengths)


def header():
    """Return the header row for the output."""
    yield 'Tumour'
    yield 'Sample'
    for chrom in Chromosomes.in_order:
        for end in ['p amp', 'p del', 'q amp', 'q del']:
            for middle in ['Num_Segments', 'Segments_Length', 'Frac_Length', 'Segments_Mean']:
                yield '{} {} {}'.format(chrom, middle, end)
   #  for column in ['OS_STATUS', 'OS_MONTHS', 'DFS_STATUS', 'DFS_MONTHS']:
   #     yield column


def format_output(tumour, sample_id, seg_ns, seg_lens, seg_means, arm_lengths):
    row = [tumour, sample_id]
    for chrom in Chromosomes.in_order:
        for arm, direction in [('p', 'amp'), ('p', 'del'), ('q', 'amp'), ('q', 'del')]:
            key = (chrom, arm, direction)

            if seg_lens[key] == 0 and arm_lengths[(chrom, arm)] == 0:
                frac = 0
            elif arm_lengths[(chrom, arm)] > 0:
                frac = seg_lens[key] / arm_lengths[(chrom, arm)]
            else:
                raise ValueError('have segment length without arm length, {}{} {}'.format(
                    chrom, arm, direction))

            row.extend([seg_ns[key],
                        seg_lens[key],
                        frac,
                        seg_means[key]])
    return row


def process_single_tumour(tumour, seg_files):


    segments = {}
    for seg_file in seg_files:
        for segment in SegmentFile.parse(seg_file):
            segments.setdefault(segment.sample, []).append(segment)

    out_handle = io.StringIO()
    writer = csv.writer(out_handle)
    writer.writerow(list(header()))

    for sample_id, sample_segments in segments.items():
        barcode = Barcode(sample_id)

        (seg_ns, seg_lens, seg_means, arm_lengths) = summarise_sample(sample_segments)
        writer.writerow(format_output(tumour, sample_id, seg_ns, seg_lens, seg_means, arm_lengths))

    out_handle.seek(0)
    dat = pd.read_csv(out_handle)
    return dat


def find_segment_files(tumour):
    result = []
    for path in pathlib.Path('cnv_data').glob('*.seg.txt'):
        if path.name.startswith(tumour):
            result.append(str(path))
    if result == []:
        raise ValueError("can not find segment file for {}".format(tumour))
    else:
        return result




def process_tumours(output_file):
    cohorts = ['msk_impact_2017_data_cna_hg19']

    with pd.ExcelWriter(output_file) as writer:
        for tumour in cohorts:
            seg_files = find_segment_files(tumour)
       #     survival_file = find_survival_file(tumour)
            dat = process_single_tumour(tumour, seg_files)
            dat.to_excel(writer, sheet_name=tumour, index=False)


def main():
    process_tumours('PCAWG_CNV_Analysis.xlsx')


if __name__ == "__main__":
    main()