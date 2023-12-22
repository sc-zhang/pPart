#!/usr/bin/env python3
from ppart.io.message import Message
from ppart.io.file_io import AlnIO, BamIO
import pysam
import argparse


def get_opts():
    group = argparse.ArgumentParser()

    group.add_argument("-i", "--input", help="Input aln file", required=True)
    group.add_argument("-o", "--output", help="Output path file", required=True)

    return group.parse_args()


def main():
    opts = get_opts()
    in_aln = opts.input
    out_info = opts.output

    Message.info("Loading Aln file")
    aln_io = AlnIO(in_aln)
    aln_io.load()

    Message.info("Converting")
    aln_io.to_path()

    Message.info("Writing result")
    aln_io.save_path(out_info)

    Message.info("Finished")


if __name__ == "__main__":
    main()
