import logging
from golib.core.gene_ontology import GeneOntology
from rich.progress import track
import scipy.stats as stats

logger = logging.getLogger('overrepresentation')
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    "%Y-%m-%d %H:%M:%S")
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

def read_complexes(path):
    complexes = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("Clu"):
                fields = line.strip().split(",")
                complexes[int(fields[0])] = fields[7].replace('"', "").split()
    return complexes

def get_proteins_from_fasta_file(path):
    proteins = []
    with open(path) as f:
        for line in f:
            if line[0] == ">":
                proteins.append(line.strip().split("|")[1])
    return proteins

def run(proteome_file, complexes_file, goa_file, obo_file, out_file, 
        pvalue_tau=0.05, min_group_count=1):
    logger.info(f"Parsing proteome fasta file {proteome_file}...")
    background = get_proteins_from_fasta_file(proteome_file)
    total_background = len(background)

    logger.info("Building Ontology in memory...")
    go = GeneOntology(obo=obo_file)
    go.build_ontology()

    logger.info("Loading GO annotations for this proteome...")
    go.load_gaf_file(goa_file, "overrep")
    go.up_propagate_annotations("overrep")
    annotations = go.annotations("overrep")

    logger.info("Processing Complexes file...")
    complexes = read_complexes(complexes_file)
    
    logger.info(f"Found {len(complexes)} complexes, analyzing overrepresentation")
    cond = annotations["Protein"].isin(background)
    goterms = annotations[cond]["GO ID"].unique()
    overrepresented_goterms = []
    num_tested_hypotheses = len(goterms) 
    for complex_id, proteins in complexes.items():
        logger.info(f"Analyzing complex {complex_id}/{len(complexes)} ({complex_id/len(complexes) * 100.0:.2f}%)) ...")
        total_group = len(proteins)
        if total_group < min_group_count:
            continue
        for goterm in track(goterms, description="Analyzing..."):
            cond_background = (annotations["GO ID"] == goterm) & annotations["Protein"].isin(background)
            cond_group = (annotations["GO ID"] == goterm) & annotations["Protein"].isin(proteins)
            background_count = annotations[cond_background].shape[0]
            group_count = annotations[cond_group].shape[0]
            contingency_table = [[group_count, background_count],
                                 [total_group - group_count, total_background - background_count]]
            _, pvalue = stats.fisher_exact(contingency_table, alternative="greater")
            corrected_pvalue = pvalue * num_tested_hypotheses
            if corrected_pvalue < pvalue_tau:
                overrepresented_goterms.append((complex_id, goterm, corrected_pvalue))

    logger.info("Writing overrepresentation file...")
    with open(out_file, "w") as out:
        for complex_id, goterm, pvalue in overrepresented_goterms:
            out.write(f"{complex_id}\t{goterm}\t{pvalue}\n")

    logger.info("Done")
        

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="formats the raw output of blast into the homolog information for ProteoBOOSTER")
    parser.add_argument("proteome_file", help="Path to the proteome fasta file")
    parser.add_argument("complexes_file", help="File with complexes to analyze")
    parser.add_argument("goa_file", help="path to GOA file")
    parser.add_argument("obo_file", help="path to go.obo file")
    parser.add_argument("output_file", help="path to write the results")
    args = parser.parse_args()
    run(args.proteome_file, args.complexes_file, args.goa_file, args.obo_file, args.output_file)


