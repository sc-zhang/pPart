from ppart.io.message import Message
from ppart.io.file_io import AlnIO
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
    aln_io = AlnIO(in_aln)
    aln_io.load()

    Message.info("Converting")
    aln_io.to_path()

    Message.info("Writing result")
    aln_io.save_path(out_info)

    Message.info("Finished")
