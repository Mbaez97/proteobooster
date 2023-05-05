import logging
from golib.core import gene_ontology
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple

logger = logging.getLogger('prepare_database_files')
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    "%Y-%m-%d %H:%M:%S")
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)


def get_tax_info(taxonomy_info_file) -> pd.DataFrame:

    def extract_kingdom(lineage):
        domain_map = {
            "Archaea": "a",
            "Viruses": "v",
            "Bacteria": "b",
            "Eukaryota": "e",
        }
        for name, abbr in domain_map.items():
            try:
                if name in lineage:
                    return abbr
            except TypeError:
                return "NONE"
        return "NONE"

    taxonomy_info = pd.read_table(taxonomy_info_file)
    taxonomy_info["domain"] = taxonomy_info["Lineage"].apply(extract_kingdom)
    cond = taxonomy_info["domain"] == "NONE"
    taxonomy_info = taxonomy_info[~cond][
        ["Taxon Id", "Scientific name", "domain"]]
    return taxonomy_info


def get_goterm_data(obo_file: str,
                    annotation_file: str,
                    proteins: pd.DataFrame) -> Tuple(pd.DataFrame,
                                                     pd.DataFrame):
    go = gene_ontology.GeneOntology(obo_file)
    go.build_ontology()
    go_data = {k: [] for k in ["go_id", "function", "ontology", "db_id"]}
    sub_ontology_map = {
        "biological_process": "bp",
        "molecular_function": "mf",
        "cellular_component": "cc",
    }
    for db_id, t in enumerate(go._terms, start=1):
        go_data["go_id"].append(int(t.split(":")[-1]))
        go_data["function"].append(go._terms[t].name)
        go_data["ontology"].append(sub_ontology_map[go._terms[t].domain])
        go_data["db_id"].append(db_id)

    go.load_gaf_file(annotation_file, "goa")
    go.up_propagate_annotations("goa")
    annotations = go.annotations("goa")
    annotations = annotations.merge(proteins,
                                    left_on="Protein",
                                    right_on="protein_id")
    annotations["go_id"] = annotations["GO_ID"].apply(
        lambda x: int(x.split(":")[-1])
    )
    annotations = annotations.merge(go_data,
                                    left_on="go_id",
                                    right_on="go_id")[["db_id", "DB_ID"]]
    annotations.columns = ["DB_ID_GO", "DB_ID_protein"]
    return pd.DataFrame(go_data), annotations


def load_interologs(interolog_file: str) -> pd.DataFrame:
    logger.info("Loading interologs...")
    cols = [
        "target1",
        "target2",
        "source1",
        "source1",
        "quality",
        "evalue",
        "org_id",
        "det_type",
        "int_type",
        "pubmed",
    ]
    interologs = pd.read_table(interolog_file, names=cols)

    interologs['target1'],\
        interologs['target2'],\
        interologs['source1'],\
        interologs['source2'] = (
            np.where(
                interologs['source1'] > interologs['source2'],

                (interologs['target2'], interologs['target1'],
                 interologs['source2'], interologs['source1']),

                (interologs['target1'], interologs['target2'],
                 interologs['source1'], interologs['source2'])
            )
        )
    return interologs


def run(alias: str,
        obo_file: str,
        goa_file: str,
        interolog_file: str,
        target_proteins_info_file: str,
        source_proteins_info_file: str,
        output_dir: str) -> None:
    outdir = Path(output_dir)

    interologs = load_interologs(interolog_file)

    logger.info("Extracting proteins from interolog table...")
    proteins = pd.DataFrame(
        np.unique(
            interologs[["target1", "target2", "source1", "source2"]]
            .values.flatten()
        ),
        columns=["protein_id"]
    ).reset_index(names="DB_ID")
    proteins["DB_ID"] += 1

    logger.info("Loading Gene Ontology file...")
    go_df, annotations = get_goterm_data(obo_file)
    go_df[["db_id", "go_id", "function", "ontology"]].to_csv(
        outdir / f"{alias}-goterms.tsv",
        sep="\t",
        index=False
    )
    # TODO(mateo): figure out when/how we should get the protein names
    # from unirprot, fow now we assume this is done manually by the user
    # and fed into this script

    # source_proteins = pd.DataFrame(
    #     np.unique(interologs[["source1", "source2"]].values.flatten()),
    #     columns=["protein_id"]
    # )
    logger.info("Preparing source proteins...")
    source_proteins_info = pd.read_table(source_proteins_info_file)
    source_proteins_info["database"] = source_proteins_info["Reviewed"].apply(
        lambda x: 0 if x == "Reviewed" else 1
    )
    source_proteins_info["has_experimental_interactions"] = 1

    logger.info("Preparing target proteins...")
    target_proteins_info = pd.read_table(target_proteins_info_file)
    target_proteins_info["database"] = source_proteins_info["Reviewed"].apply(
        lambda x: 0 if x == "Reviewed" else 1
    )
    # TODO(mateo): calculate this from the interactions file
    source_proteins_info["has_experimental_interactions"] = 1


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="generates files to load the ProteoBOOSTER database")
    parser.add_argument("proteome_file",
                        help="Path to the proteome fasta file")
    parser.add_argument("complexes_file",
                        help="File with complexes to analyze")
    parser.add_argument("goa_file", help="path to GOA file")
    parser.add_argument("obo_file", help="path to go.obo file")
    parser.add_argument("output_dir",
                        help="path to a directory to write the results")
    args = parser.parse_args()
    run(args.proteome_file, args.complexes_file, args.goa_file,
        args.obo_file, args.output_file)
