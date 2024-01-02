from ppart.io.message import Message
from ppart.func.bam_operate import BamOperate
import pysam


class Classifier:
    def __init__(self):
        pass

    @staticmethod
    def __binary_search(regions, pos):
        left = 0
        right = len(regions)-1
        while left <= right:
            mid = (right-left)//2+left
            if regions[mid][0] > pos:
                right = mid-1
            elif regions[mid][0] < pos:
                left = mid+1
            else:
                return mid
        if regions[right][0] <= pos <= regions[right][1]:
            return right
        else:
            return -1

    def classify(self, var_regions, bam_file):
        Message.info("Loading bam file")
        bam_op = BamOperate
        with pysam.AlignmentFile(bam_file, 'rb') as fin:
            for record in fin:
                qry_start = record.query_alignment_start
                qry_seq = record.query_sequence
                ref_start = record.reference_start
                ref_seq = record.get_reference_sequence()
                region_idx = self.__binary_search(var_regions, ref_start)
                if region_idx == -1:
                    continue
                var_start = var_regions[region_idx][0]
                var_end = var_regions[region_idx][1]
                ref_aln_seq = ref_seq[var_start-1: var_end]
                qry_aln_seq = ""
                for pos in range(var_start, var_end+1):
                    offset = pos-ref_start
                    base = bam_op.get_base(qry_seq, qry_start, offset, record.cigartuples)
                    qry_aln_seq += base
