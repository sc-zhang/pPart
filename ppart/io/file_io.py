from ppart.utils.variant_types import Variants


class AlnIO:
    def __init__(self, in_aln_file):
        self.__aln_file = in_aln_file
        self.__aln_db = {}
        self.__aln_len = 0
        self.__consensus_seq = ""
        self.__diff_pos = []

        # path_db is a dictionary store combinations of variants for all genes,
        # gene id => [var1_idx, var2_idx, var3_idx, ..., varN_idx]
        self.__path_db = {}
        # pos_db is a dictionary store variants information for all variant regions
        # start position of regions (1-based) => variants sequence => index of variants (0-based)
        self.__pos_db = {}
        # regions is a list store variants, (0-based)
        # [start position, end position, "N"(Consensus region)|"V"(Variant region)]
        self.__regions = []

    # Read alignment file and store it to a dictionary
    # key: gene id,
    # value: aligned sequence
    def load(self):
        self.__aln_db = {}
        with open(self.__aln_file, 'r') as fin:
            for line in fin:
                if line[0] == '>':
                    gid = line.strip().split()[0][1:].replace('_R_', '')
                    self.__aln_db[gid] = []
                else:
                    self.__aln_db[gid].append(line.strip().upper())

        for gid in self.__aln_db:
            self.__aln_db[gid] = ''.join(self.__aln_db[gid])
            if self.__aln_len == 0:
                self.__aln_len = len(self.__aln_db[gid])

    # Get consensus sequence and difference positions from alignment sequences
    def __consensus(self):
        for idx in range(self.__aln_len):
            cnt_db = {}
            for gid in self.__aln_db:
                base = self.__aln_db[gid][idx]
                if base not in cnt_db:
                    cnt_db[base] = 0
                cnt_db[base] += 1
            if len(cnt_db) > 1:
                self.__diff_pos.append(idx)
            self.__consensus_seq += sorted(cnt_db, key=lambda x: cnt_db[x], reverse=True)[0]

    # Calculate edit distance between two sequences
    @staticmethod
    def __distance(seq1, seq2):
        dis = 0
        for _ in range(len(seq1)):
            if seq1[_] != seq2[_]:
                dis += 1
        return dis

    # Convert alignments to path
    def to_path(self):
        # Get all difference positions and if two positions are continuous, merge positions to regions.
        self.__consensus()
        merged_diff_pos = []
        for pos in self.__diff_pos:
            if len(merged_diff_pos) == 0:
                merged_diff_pos.append([pos, pos])
            elif pos - 1 == merged_diff_pos[-1][-1]:
                merged_diff_pos[-1][-1] = pos
            else:
                merged_diff_pos.append([pos, pos])

        # mark all consensus regions as "N",
        # mark all variants regions as "V"
        for ps, pe in merged_diff_pos:
            if not self.__regions:
                if ps != 0:
                    self.__regions.append([0, ps - 1, Variants.NORMAL])
                self.__regions.append([ps, pe, Variants.ALLVAR])
            else:
                if self.__regions[-1][1] + 1 != ps:
                    self.__regions.append([self.__regions[-1][1] + 1, ps - 1, Variants.NORMAL])
                self.__regions.append([ps, pe, Variants.ALLVAR])
        if self.__regions[-1][1] + 1 != self.__aln_len:
            self.__regions.append([self.__regions[-1][1] + 1, self.__aln_len - 1, Variants.NORMAL])

        self.__path_db = {gid: [] for gid in self.__aln_db}
        self.__pos_db = {}
        for ps, pe, _ in self.__regions:
            type_db = {}
            type_idx = 0
            for gid in sorted(self.__aln_db):
                type_seq = self.__aln_db[gid][ps: pe + 1]
                if type_seq not in type_db:
                    type_db[type_seq] = type_idx
                    type_idx += 1
                self.__path_db[gid].append(type_db[type_seq])
            self.__pos_db[ps + 1] = type_db

    # Get variant combinations of all genes/samples
    def get_path(self):
        return self.__path_db

    # Get variant index of all variants
    def get_var_idx(self):
        return self.__pos_db

    def save_path(self, out_file):
        with open(out_file, 'w') as fout:
            fout.write("#Sample|POSITIONS\t%s\n" % ('\t'.join(map(str, sorted(self.__pos_db)))))
            for gid in sorted(self.__path_db):
                fout.write("%s\t%s\n" % (gid, '\t'.join(map(str, self.__path_db[gid]))))

            fout.write("#POSITION\tBASE\tINDEX\n")
            for pos in sorted(self.__pos_db):
                for type_seq in sorted(self.__pos_db[pos], key=lambda x: self.__pos_db[pos][x]):
                    fout.write("%d\t%s\t%d\n" % (pos, type_seq, self.__pos_db[pos][type_seq]))

    def get_sample_var_regions(self):
        sample_var_regions = {}
        for smp in self.__aln_db:
            sample_var_regions[smp] = []
            offset = 0
            for sp, ep, var_type in self.__regions:
                if var_type == Variants.NORMAL:
                    continue
                gap_cnt = 0
                for _ in self.__aln_db[smp][sp: ep + 1]:
                    if _ == '-':
                        gap_cnt += 1
                smp_sp = sp - offset
                offset += gap_cnt
                smp_ep = ep - offset
                sample_var_regions[smp].append([smp_sp,
                                                smp_ep,
                                                Variants.SNP if gap_cnt == 0 else Variants.INDEL,
                                                sp,
                                                ep])
        return sample_var_regions

    def get_ori_seqs(self):
        ori_seqs = {}
        for smp in self.__aln_db:
            ori_seqs[smp] = self.__aln_db[smp].upper().replace('-', '')
        return ori_seqs
