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


def run(names_file, accession_list, output_file):
    logger.info(f"Gathering accessions from {accession_list}")
    accessions = set()
    with open(accession_list, "r") as acc_f:
        for line in acc_f:
            accessions.add(line.strip())
    valid_lines = []
    logger.info(f"Filtering {names_file}")
    with open(names_file, 'r') as names:
        for line in names:
            accession = line.split("|")[1]
            if accession in accessions:
                valid_lines.append(line)
                # search should be faster from now on
                accessions.remove(accession)

    logger.info(f"Writing {output_file}")
    with open(output_file, 'w') as output:
        output.writelines(valid_lines)
    logger.info("Done")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="filters an interactions database using a fasta file")
    parser.add_argument("names_file",
                        help="file to filter")
    parser.add_argument("accession_list",
                        help="file with lines to look up")
    parser.add_argument("output", help="path to write the results")
    args = parser.parse_args()
    run(args.names_file, args.accession_list, args.output)
