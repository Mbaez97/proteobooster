import logging

logger = logging.getLogger('download_snapshot')
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    "%Y-%m-%d %H:%M:%S")
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)


def run(fasta_file, interaction_database, outfile):
    logger.info(f"Extracting accessions from {fasta_file}")
    proteins = []
    with open(fasta_file) as f:
        for line in f:
            if line[0] == ">":
                proteins.append(line.split("|")[1])
    logger.info("Processing interaction database")

    o = open(outfile, "w")
    with open(interaction_database) as f:
        for line in f:
            fs = line.strip().split("\t")
            p1 = fs[0]
            p2 = fs[1]
            if p1 in proteins and p2 in proteins:
                o.write(line)
    o.close()
    logger.info("Done")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="filters an interactions database using a fasta file")
    parser.add_argument("fasta_file", help="processed BLAST file")
    parser.add_argument("interaction_database",
                        help="path to the combined inteactions file")
    parser.add_argument("output", help="path to write the results")
    args = parser.parse_args()
    run(args.fasta_file, args.interaction_database, args.output)
