import os
import argparse
import pysam
from collections import Counter
import pandas as pd

def check_file_existance(path):
    abs_path = os.path.abspath(path)
    if os.path.exists(abs_path):
        return abs_path
    else:
        raise argparse.ArgumentTypeError("File {} can not be located".format(path))


def get_pileup_counts(samfile,position,rangelen):
    '''
    Get read counts for a given position or positions if range is > 1
    :param samfile: PySAM connection to BAM file
    :param position: query nucleotide position in the pileup
    :param rangelen: query range length (usually equals 1)
    :return:
        a) dictionary of k-mer counts (where k could be 1 if single nucleotide position is queried)
        b) number of deleted reads
        c) total number of reads at a given position of a pileup
    '''

    kmers = []; reads_del_counter = 0; total_reads=0
    if len(samfile.references) == 1:
        reference_name=samfile.get_reference_name(0)
    else:
        raise Exception("Incompatible BAM file with several references")


    for column in samfile.pileup(reference_name,position,position+rangelen, truncate=True,
                                 max_depth=1000000, min_base_quality=0):
        print("Ref position {} and Nreads {}".format(column.pos+1, column.n))
        total_reads = column.n



        for pileupread in column.pileups:
            # for k-mers of size >1 this needs to be modified to take average of >1 posisions
            # or look at BAM flag for deletions
            if pileupread.is_del == 1:
                reads_del_counter = reads_del_counter+1
            if pileupread.query_position is not None:
                kmers.append(pileupread.alignment.query_sequence[pileupread.query_position:
                                                                 pileupread.query_position+rangelen])



    return(Counter(kmers),reads_del_counter,total_reads)



def argumnets():
    parser = argparse.ArgumentParser("Generate VOCs Waste Water report\n")
    parser.add_argument('-i', '--input', required=True,
                        type=check_file_existance,
						help="BAM file(s)")
    parser.add_argument('-o', '--output_name',
						 required=True, help="Output summary file name "
                        )

    return parser.parse_args()

def get_average_coverage(samfile):
    coverage_per_base = samfile.count_coverage(samfile.get_reference_name(0), quality_threshold=1,
                                               read_callback="nofilter")

    total_coverage = 0
    for array in coverage_per_base:
        total_coverage = total_coverage + sum(array)
    return round(total_coverage / len(coverage_per_base[0]))



def main():
    analysis_stats = {}

    args = argumnets()
    filepath = args.input

    variants_df = pd.read_csv('watchlistSNVs/variants.tsv', sep="\t")
    variants_df = variants_df[(variants_df.VOC == "UK") & (variants_df.Type == "Sub")]

    samplefile = str(os.path.basename(args.input))
    analysis_stats[samplefile] = {"AvgGenomeCoverage": "-",
                                  "UK_VOC_SNVs": {},
                                  "UK_VOC_SNVs_Total": variants_df.shape[0],
                                  "UK_VOC_Count": 0,
                                  "UK_VOCs_absent": []}

    index_file = os.path.abspath(filepath) + '.bai'
    pysam.index(args.input, index_file)
    samfile = pysam.AlignmentFile(filepath, "rb")

    analysis_stats[samplefile]["AvgGenomeCoverage"] = get_average_coverage(samfile)

    for index, row in variants_df.iterrows():

        kmers_counts, reads_del_counter, total_reads = get_pileup_counts(samfile, row.Position - 1, row.Length)

        if total_reads > 0:
            freq_snv = round(kmers_counts[row.Alt] / total_reads * 100, 1)
            analysis_stats[samplefile]['UK_VOC_Count'] = analysis_stats[samplefile]['UK_VOC_Count'] + 1
        else:
            freq_snv = 0
            analysis_stats[samplefile]["UK_VOCs_absent"] = row.AAName

        analysis_stats[samplefile]["UK_VOC_SNVs"][row.AAName] = freq_snv

    print(analysis_stats)

    samfile.close()
    print(pd.DataFrame(analysis_stats[samplefile]["UK_VOC_SNVs"].items(), columns=["SNV", "%Freq"]))


if __name__ == "__main__":
    main()
