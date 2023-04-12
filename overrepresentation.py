import logging
from golib.core.gene_ontology import GeneOntology

logger = logging.getLogger('overrepresentation')
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    "%Y-%m-%d %H:%M:%S")
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

def run(complexes, goa, obo, outfile):
    pass

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="formats the raw output of blast into the homolog information for ProteoBOOSTER")
    parser.add_argument("complexes_file", help="File with complexes to analyze")
    parser.add_argument("goa_file", help="path to GOA file")
    parser.add_argument("obo_file", help="path to go.obo file")
    parser.add_argument("output_file", help="path to write the results")
    args = parser.parse_args()
    run(args.complexes_file, args.goa_file, args.obo_file, args.output_file)


