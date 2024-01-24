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

    def get_sites(self, ref_sample: str, gene_id: str, var_regions: list, bam_file: str) -> dict:
        bam_op = BamOperate
        site_db = {}
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
                smp_var_start, smp_var_end, smp_var_type, aln_var_sp, aln_var_ep = var_regions[region_idx]

                qry_aln_seq = []
                for pos in range(smp_var_start, smp_var_end + 1):
                    offset = pos - ref_start
                    if offset < 0 or offset >= len(qry_seq):
                        base = '-'
                    else:
                        base = bam_op.get_base(qry_seq, qry_start, offset, record.cigartuples)

                    qry_aln_seq.append(base)

                qry_aln_seq = ''.join(qry_aln_seq)
                if aln_var_sp not in site_db:
                    site_db[aln_var_sp] = qry_aln_seq

        return site_db

    @staticmethod
    def get_path(site_db, var_idx_db):
        path_list = []
        for site in var_idx_db:
            if site - 1 in site_db:
                seq = site_db[site - 1]
                if seq in var_idx_db[site]:
                    idx = var_idx_db[site][seq]
                else:
                    idx = -1
            else:
                idx = -1
            path_list.append(idx)
        return path_list

    @staticmethod
    def get_similarity(qry_path, ref_path):
        cnt = 0
        for _ in range(len(qry_path)):
            if qry_path[_] == ref_path[_]:
                cnt += 1
        return cnt * 1. / len(qry_path)
