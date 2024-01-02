class BamOperate:
    def __init__(self):
        pass

    @staticmethod
    def get_base(qry_seq, qry_start, ref_offset, cigars):
        qry_offset = 0
        # alignment	meaning	operation
        # M	MATCH	0
        # I	INS	1
        # D	DEL	2
        # N	REF_SKIP	3
        # S	SOFT_CLIP	4
        # H	HARD_CLIP	5
        # P	PAD	6
        # =	EQUAL	7
        # X	DIFF	8
        # B	BACK	9
        qry_only = {1, 4}
        ref_only = {2, 3}
        qry_ref_both = {0, 7, 8}
        for align_type, length in cigars:
            if ref_offset == 0:
                break
            if align_type in qry_only:
                qry_offset += length
            elif align_type in ref_only:
                ref_offset -= length
            elif align_type in qry_ref_both:
                if ref_offset >= length:
                    ref_offset -= length
                    qry_offset += length
                else:
                    qry_offset += ref_offset
                    ref_offset = 0
        return qry_seq[qry_start + qry_offset]
