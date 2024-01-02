from os import listdir, path
from ppart.io.message import Message
from ppart.io.file_io import AlnIO
from ppart.func.classifier import Classifier
import argparse


def get_opts():
    group = argparse.ArgumentParser()

    group.add_argument("-i", "--input", help="Input aln file", required=True)
    group.add_argument("-b", "--bam", help="Input bam directory", required=True)
    group.add_argument("-o", "--output", help="Output path file", required=True)
    group.add_argument("-t", "--thread", help="Threads, default=10", type=int, default=10)

    return group.parse_args()


def main():
    opts = get_opts()
    in_aln = opts.input
    in_bam_dir = opts.bam
    out_info = opts.output
    threads = opts.thread

    Message.info("Loading Aln file")
    gid = in_aln.split('/')[-1].split('\\')[-1].split('.')[0]
    aln_io = AlnIO(in_aln)
    aln_io.load()

    Message.info("Converting")
    aln_io.to_path()

    Message.info("Getting sample specific variant regions")
    specific_vars = aln_io.get_sample_var_regions()
    fa_db = aln_io.get_ori_seqs()
    classifier = Classifier()
    for smp in specific_vars:
        for bam_file in listdir(in_bam_dir):
            if not bam_file.endswith(".bam"):
                continue
            print(smp, bam_file)
            bam_path = path.join(in_bam_dir, bam_file)
            classifier.classify(smp, fa_db[smp], gid, specific_vars[smp], bam_path)

    Message.info("Writing result")
    aln_io.save_path(out_info)

    Message.info("Finished")
