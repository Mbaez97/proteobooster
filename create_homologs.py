from rich.progress import track


def run(blast_file, homologs_file):
    evalues = {}
    percent_identities = {}
    for line in open(blast_file):
        output_fields = line.strip().split("\t")
        # here we're assuming uniprot accession structure
        org_accession = output_fields[0].split("|")[1]
        interactor_accession = output_fields[1].split("|")[1]
        evalue = float(output_fields[10])
        if evalue > 1e-2:
            continue
        # we don't really care if the proteins are the same protein
        if org_accession == interactor_accession:
            continue
        key = f"{org_accession}\t{interactor_accession}"
        if key not in evalues or evalues[key] > evalue:
            evalues[key] = evalue
            percent_identities[key] = output_fields[2]

    with open(homologs_file, "w") as homologs:
        for key in track(evalues.keys(), total=len(evalues)):
            org_prot, other_prot = key.split()
            homologs.write(f"{org_prot}\t{other_prot}\t{evalues[key]}"
                           f"\t{percent_identities[key]}\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="formats the raw output of blast into the"
                    " homolog information for ProteoBOOSTER")
    parser.add_argument("blast_file", help="raw BLAST output")
    parser.add_argument("homolog_file",
                        help="path where the homolog file will be written")
    args = parser.parse_args()
    run(args.blast_file, args.homolog_file)
