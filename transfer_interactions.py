import numpy as np
from rich.progress import track
from collections import defaultdict
import logging
import sys

logger = logging.getLogger('download_snapshot')
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    "%Y-%m-%d %H:%M:%S")
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

def buf_count_newlines_gen(fname):
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b

    with open(fname, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count


def run(homolog_file, interaction_database, out_file):
    logger.info(f"Building homolog dictionary for {homolog_file}")
    homologs = defaultdict(list) 
    for line in open(homolog_file):
        target_protein, interactor, evalue, percent_identity = line.strip().split("\t")
        evalue = float(evalue)
        if evalue >= 1e-10:
            continue
        percent_identity = float(percent_identity)
        homologs[interactor].append((target_protein, percent_identity, evalue))

    logger.info("Finished building the homolog dictionary")
    logger.info("Transferring interologs")
    interolog_file = open(out_file, "w")
    for line in open(interaction_database):
        fs = line.strip().split("\t")
        try:
            if len(fs) == 6:
                s1, s2, org_id, det_type, int_type, pubmed = *fs,
            else:
                s1, s2, org_id, det_type, int_type, = *fs,
                pubmed = "-"
        except ValueError:
            print(line)
            sys.exit(1)

        if not s1 in homologs or not s2 in homologs:
            continue
        for p1, perc1, eval1 in homologs[s1]:
            for p2, perc2, eval2 in homologs[s2]:
                if p1 != p2:
                    interactors = sorted([p1, p2])
                    interolog_quality = np.sqrt(perc1 * perc2)
                    if interolog_quality < 75:
                        continue
                    evalue = max(eval1, eval2)
                    ss1 = s1 if interactors[0] == p1 else s2
                    ss2 = s1 if interactors[1] == p1 else s2
                    interolog_file.write(
                        f"{interactors[0]}\t{interactors[1]}\t{ss1}\t{ss2}\t{interolog_quality}\t{evalue}\t{org_id}\t{det_type}\t{int_type}\t{pubmed}\n"
                    )

    interolog_file.close()
    logger.info("Done")



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="formats the raw output of blast into the homolog information for ProteoBOOSTER")
    parser.add_argument("homolog_file", help="processed BLAST file")
    parser.add_argument("interaction_database", help="path to the combined inteactions file")
    parser.add_argument("interolog_file", help="path to write the interologs")
    args = parser.parse_args()
    run(args.homolog_file, args.interaction_database, args.interolog_file)

