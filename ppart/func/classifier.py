from ppart.io.message import Message
from ppart.func.bam_operate import BamOperate
import pysam


class Classifier:
    def __init__(self):
        pass

    @staticmethod
    def __binary_search(regions: list, pos: int):
        left = 0
        right = len(regions) - 1
        while left <= right:
            mid = (right - left) // 2 + left
            if regions[mid][0] > pos:
                right = mid - 1
            elif regions[mid][0] < pos:
                left = mid + 1
            else:
                return mid
        if regions[right][0] <= pos <= regions[right][1]:
            return right
        else:
            return -1

    def get_similarity(self, ref_sample: str, ref_seq: str, gene_id: str, var_regions: list, bam_file: str) -> list:
        bam_op = BamOperate

        similarity = {_[0]: 0 for _ in var_regions}
        with pysam.AlignmentFile(bam_file, 'rb') as fin:
            for record in fin:
                ref_name = record.reference_name

                # only on test data
                if ref_sample.split('.')[0] not in ref_name or gene_id not in ref_name:
                    continue

                qry_start = record.query_alignment_start
                qry_seq = record.query_sequence
                ref_start = record.reference_start

                region_idx = self.__binary_search(var_regions, ref_start)
                if region_idx == -1:
                    continue
                var_start, var_end, var_type = var_regions[region_idx]

                qry_aln_seq = []
                for pos in range(var_start, var_end + 1):
                    offset = pos - ref_start
                    if offset < 0 or offset >= len(qry_seq):
                        base = '-'
                    else:
                        base = bam_op.get_base(qry_seq, qry_start, offset, record.cigartuples)

                    qry_aln_seq.append(base)

                qry_aln_seq = ''.join(qry_aln_seq)
                ref_aln_seq = ref_seq[var_start: var_end + 1]
                cur_similarity = 0
                for _ in range(len(qry_aln_seq)):
                    if qry_aln_seq[_] == ref_aln_seq[_]:
                        cur_similarity += 1
                cur_similarity = cur_similarity*100./len(qry_aln_seq)
                if cur_similarity > similarity[var_start]:
                    similarity[var_start] = cur_similarity

        return [similarity[_] for _ in sorted(similarity)]
