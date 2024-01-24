from pathos import multiprocessing
from os import listdir, path
from ppart.io.message import Message
from ppart.io.file_io import AlnIO
from ppart.func.classifier import Classifier
import argparse


def get_opts():
    group = argparse.ArgumentParser()

    group.add_argument("-i", "--input", help="Input aln file", required=True)
    group.add_argument("-b", "--bam", help="Input bam directory", required=True)
    group.add_argument("-o", "--output", help="Output prefix", required=True)
    group.add_argument("-t", "--thread", help="Threads, default=10", type=int, default=10)

    return group.parse_args()


def main():
    opts = get_opts()
    in_aln = opts.input
    in_bam_dir = opts.bam
    out_prefix = opts.output
    threads = opts.thread

    Message.info("Loading Aln file")
    gid = in_aln.split('/')[-1].split('\\')[-1].split('.')[0]
    aln_io = AlnIO(in_aln)
    aln_io.load()

    Message.info("Converting")
    aln_io.to_path()

    Message.info("Getting sample specific variant regions")
    specific_vars = aln_io.get_sample_var_regions()
    classifier = Classifier()

    res = []
    pool = multiprocessing.Pool(processes=threads)
    for ref in specific_vars:
        for bam_file in listdir(in_bam_dir):
            if not bam_file.endswith(".bam"):
                continue
            smp = bam_file.split('.')[0]
            Message.info("\tChecking {} with {}".format(ref, smp))
            bam_path = path.join(in_bam_dir, bam_file)
            r = pool.apply_async(classifier.get_sites, (ref, gid, specific_vars[ref], bam_path,))
            res.append([smp, ref, r])

    pool.close()
    pool.join()

    Message.info("Getting sample paths")
    smp_var_db = {}
    for smp, ref, r in res:
        if smp not in smp_var_db:
            smp_var_db[smp] = {}
        try:
            cur_db = r.get()
            for sp in cur_db:
                if sp not in smp_var_db[smp]:
                    smp_var_db[smp][sp] = cur_db[sp]
        except Exception as e:
            Message.warn("Exception found: {}".format(e))

    smp_path_db = {}
    for smp in smp_var_db:
        smp_path_db[smp] = classifier.get_path(smp_var_db[smp], aln_io.get_var_idx())

    Message.info("Comparing paths")

    res = []
    ref_path_db = aln_io.get_path()
    pool = multiprocessing.Pool(processes=threads)
    for smp in smp_path_db:
        for ref in specific_vars:
            Message.info("\tComparing {} with {}".format(smp, ref))
            r = pool.apply_async(classifier.get_similarity, (smp_path_db[smp], ref_path_db[ref]))
            res.append([smp, ref, r])
    pool.close()
    pool.join()

    best_match_db = {}
    for smp, ref, r in res:
        if smp not in best_match_db:
            best_match_db[smp] = ["", 0]
        try:
            similarity = r.get()
            if similarity > best_match_db[smp][-1]:
                best_match_db[smp][0] = ref
                best_match_db[smp][1] = similarity
        except Exception as e:
            Message.warn("Exception found: {}".format(e))

    Message.info("Writing result")
    aln_io.save_path(out_prefix + ".path")
    with open(out_prefix + ".match", 'w') as fout:
        for smp in sorted(best_match_db):
            fout.write("%s Match %s, similarity: %.2f%%\n" % (smp,
                                                              best_match_db[smp][0],
                                                              best_match_db[smp][1] * 100.))

    Message.info("Finished")
