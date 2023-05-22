import logging
from golib.core import gene_ontology
from golib.io import obo_parser
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
    taxonomy_info["Taxon Id"] = taxonomy_info["Taxon Id"].astype(str)
    return taxonomy_info


def get_mi_ontology(mi_obo: str) -> pd.DataFrame:
    op = obo_parser.OboParser(mi_obo)
    mi_data = {k: [] for k in ["mi_id", "description"]}
    for stanza in op:
        mi_data["mi_id"].append(stanza.tags["id"][0].value)
        mi_data["description"].append(stanza.tags["name"][0].value)
    mi_ontology = pd.DataFrame(mi_data)
    return mi_ontology.reset_index(names="DB_ID")


def get_goterm_data(obo_file: str,
                    annotation_file: str,
                    proteins: pd.DataFrame) -> Tuple[pd.DataFrame,
                                                     pd.DataFrame]:
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
    go_data = pd.DataFrame(go_data)

    go.load_gaf_file(annotation_file, "goa")
    go.up_propagate_annotations("goa")

    logger.info("Preparing GO data...")
    annotations = go.annotations("goa")
    annotations = annotations.merge(proteins,
                                    left_on="Protein",
                                    right_on="protein_id")
    annotations["go_id"] = annotations["GO ID"].apply(
        lambda x: int(x.split(":")[-1])
    )
    annotations = annotations.merge(go_data,
                                    left_on="go_id",
                                    right_on="go_id")[["db_id", "DB_ID"]]
    annotations.columns = ["DB_ID_GO", "DB_ID_protein"]

    return go_data, annotations


def load_interologs(interolog_file: str) -> pd.DataFrame:
    logger.info("Loading interologs...")
    cols = [
        "target1",
        "target2",
        "source1",
        "source2",
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
    interologs["org_id"] = interologs["org_id"].astype(str)
    return interologs


def load_exp_interactions(interactions_file: str) -> pd.DataFrame:
    """
    Loads the experimental interactions belonging to an organism
    Note
    ----
        This function loads the interactions belonging to the target organism
        only. Interacttions from source organisms will be derived from the
        interologs table, as we don't want to replicate the entirety of
        the source database in each instance.
    """
    logger.info(f"Loading interactions from {interactions_file}...")
    cols = [
        "p1",
        "p2",
        "org_id",
        "det_type",
        "int_type",
        "pubmed"
    ]
    interactions = pd.read_table(interactions_file, names=cols)
    interactions["p1"], interactions["p2"] = (
        np.where(
            interactions["p1"] > interactions["p2"],
            (interactions["p2"], interactions["p1"]),
            (interactions["p1"], interactions["p2"])
        )
    )
    interactions["org_id"] = interactions["org_id"].astype(str)
    return interactions


def parse_fasta_line(line: str) -> Tuple[str, str, str, str, str]:
    db_and_names, tax_and_other = line.strip().split("OX=")
    db_acc, name = db_and_names.strip().split(" ", maxsplit=1)
    db, acc, entry = db_acc[1:].split("|")
    tax = tax_and_other.split(" ", maxsplit=1)[0]
    return db, acc, entry, name, tax


def load_protein_info(protein_info_file: str) -> pd.DataFrame:
    data = {k: [] for k in ["accession", "entry", "name", "db", "taxon"]}
    for line in open(protein_info_file):
        db, acc, entry, name, tax = parse_fasta_line(line)
        data["accession"].append(acc)
        data["db"].append(db)
        data["entry"].append(entry)
        data["name"].append(name)
        data["taxon"].append(tax)
    return pd.DataFrame(data)


def get_complex_data(complexes_file) -> Tuple[pd.DataFrame, pd.DataFrame]:
    complexes = pd.read_csv(complexes_file)
    complex_protein = {k: [] for k in ["protein", "complex_id"]}
    for i, r in complexes.iterrows():
        prots = r["Members"].split()
        for protein in prots:
            complex_protein["protein"].append(protein)
            complex_protein["complex_id"].append(r["Cluster"])
    complex_protein = pd.DataFrame(complex_protein)
    return complexes, complex_protein


def add_protein_metadata(proteins_info_file: str,
                         interologs: pd.DataFrame,
                         interactions: pd.DataFrame,
                         proteins: pd.DataFrame) -> pd.DataFrame:
    interolog_source_proteins = np.unique(
        interologs[["source1", "source2"]].values.flatten())
    interaction_proteins = np.unique(
        interactions[["p1", "p2"]].values.flatten())
    proteins_with_exp_interactions = np.unique(
        np.hstack([
            interolog_source_proteins,
            interaction_proteins
        ])
    )
    proteins_info = load_protein_info(proteins_info_file)
    proteins_info["database"] = proteins_info["db"].apply(
        lambda x: 0 if x == "sp" else 1
    )
    proteins_info["has_experimental_interaction"] = (
        proteins_info["accession"].isin(proteins_with_exp_interactions))
    return proteins.merge(proteins_info,
                          left_on="protein_id",
                          right_on="accession")


def load_homolog_data(homologs_file: str,
                      proteins: pd.DataFrame) -> pd.DataFrame:
    homologs = pd.read_table(homologs_file,
                             names=["target", "source", "evalue", "perc_id"])
    homologs = (homologs
                .merge(proteins[["DB_ID", "accession"]],
                       left_on="target", right_on="accession")
                .drop(columns=["accession"])
                .merge(proteins[["DB_ID", "accession"]],
                       left_on="source", right_on="accession",
                       suffixes=["_target", "_source"])
                .reset_index(names="DB_ID"))
    homologs["DB_ID"] += 1
    return homologs


def get_evidence(interactions: pd.DataFrame) -> pd.DataFrame:
    evidence = interactions[["pubmed", "DB_ID_detection_type", "DB_ID",
                             "DB_ID_interaction_type"]].rename(
                                     columns={
                                         "DB_ID": "interaction_supported"})
    ev_data = {k: [] for k in ["DB_ID", "pubmed", "detection", "interaction",
                               "supported_interaction"]}
    ev_id = 1
    evidence["pubmed"] = evidence["pubmed"].astype(str)
    for _, r in evidence.iterrows():
        if "," in r["pubmed"]:
            pubmeds = r["pubmed"].split(",")
        else:
            pubmeds = [r["pubmed"]]
        for pub in pubmeds:
            if pub == "-":
                continue
            ev_data["DB_ID"].append(ev_id)
            ev_data["pubmed"].append(pub)
            ev_data["detection"].append(r["DB_ID_detection_type"])
            ev_data["supported_interaction"].append(r["interaction_supported"])
            ev_data["interaction"].append(r["DB_ID_interaction_type"])
            ev_id += 1
    return pd.DataFrame(ev_data)


def run(alias: str,
        obo_file: str,
        goa_file: str,
        interolog_file: str,
        interactions_file: str,
        proteins_info_file: str,
        homologs_file: str,
        complexes_file: str,
        overrepresentation_file: str,
        mi_obo: str,
        taxonomy_file: str,
        output_dir: str):
    outdir = Path(output_dir)
    interologs = load_interologs(interolog_file)
    interactions = load_exp_interactions(interactions_file)

    logger.info("Loading all taxonomies...")
    tax_info = get_tax_info(taxonomy_file)

    logger.info("Gathering protein information...")
    proteins = pd.DataFrame(
        np.unique(
            np.hstack([
                interologs[["target1", "target2",
                            "source1", "source2"]].values.flatten(),
                interactions[["p1", "p2"]].values.flatten()
            ])
        ),
        columns=["protein_id"]
    ).reset_index(names="DB_ID")
    proteins["DB_ID"] += 1

    logger.info("Loading Gene Ontology File")
    go_df, annotations = get_goterm_data(obo_file, goa_file, proteins)

    goterm_db_file = outdir / f"{alias}-db-goterms.tsv"
    logger.info(f"Writing GO terms information into {goterm_db_file}...")
    go_df[["db_id", "go_id", "function", "ontology"]].to_csv(
        goterm_db_file,
        sep="\t",
        index=False
    )

    annotations_db_file = outdir / f"{alias}-db-go_annotations.tsv"
    logger.info(f"Writing GO annotations into {annotations_db_file}...")
    annotations[["DB_ID_GO", "DB_ID_protein"]].to_csv(
        annotations_db_file,
        sep="\t",
        index=False
    )

    logger.info(f"Loading information from MI ontology file {mi_obo}...")
    mi_ontology = get_mi_ontology(mi_obo)

    mi_db_file = outdir / f"{alias}-db-mi_ontology.tsv"
    logger.info(f"Writing MI ontology as TSV into {mi_db_file}...")
    mi_ontology.to_csv(mi_db_file, sep="\t", index=False)

    logger.info(f"Loading proteins metadata from {proteins_info_file}...")
    proteins = add_protein_metadata(proteins_info_file,
                                    interologs,
                                    interactions,
                                    proteins)

    logger.info("Filtering taxonomy data...")
    cond = tax_info["Taxon Id"].isin(proteins["taxon"].unique())
    tax_info_keep = (tax_info[cond]
                     .reset_index(drop=True)
                     .reset_index(names="DB_ID"))
    tax_info_keep["DB_ID"] += 1
    org_db_file = outdir / f"{alias}-db-organisms.tsv"
    logger.info(f"Writing organisms file into {org_db_file}...")
    tax_info_keep.to_csv(org_db_file, sep="\t", index=False)

    proteins_db_file = outdir / f"{alias}-db-proteins.tsv"
    logger.info(f"Writing organisms file into {proteins_db_file}...")
    proteins = (proteins
                .merge(tax_info_keep[["DB_ID", "Taxon Id"]],
                       left_on="taxon", right_on="Taxon Id")
                .rename(columns={"DB_ID_x": "DB_ID",
                                 "DB_ID_y": "DB_ID_org"}))[[
                                     "DB_ID", "accession", "entry", "name",
                                     "database", "DB_ID_org",
                                     "has_experimental_interaction"]]
    proteins.to_csv(proteins_db_file, sep="\t", index=False)

    logger.info(f"Loading homologs from {homologs_file}...")
    homologs = load_homolog_data(homologs_file, proteins)
    homologs_db_file = outdir / f"{alias}-db-homologs.tsv"
    logger.info(f"Writing homologs file into {homologs_db_file}...")
    homologs[["DB_ID", "DB_ID_target", "DB_ID_source",
              "evalue", "perc_id"]].to_csv(homologs_db_file,
                                           sep="\t", index=False)

    logger.info("Merging interologs...")
    p_cols = ["accession", "DB_ID"]
    interologs = (interologs
                  .merge(proteins[p_cols],
                         left_on="target1", right_on="accession")
                  .merge(proteins[p_cols],
                         left_on="target2", right_on="accession",
                         suffixes=["_1", "_2"])
                  .merge(proteins[p_cols],
                         left_on="source1", right_on="accession",
                         suffixes=["_source_1", "_source_2"])
                  .merge(proteins[p_cols],
                         left_on="source2", right_on="accession",
                         suffixes=["_source_1", "_source_2"]))
    keep_cols = ["DB_ID_1", "DB_ID_2", "DB_ID_source_1", "DB_ID_source_2",
                 "evalue", "quality", "org_id",
                 "det_type", "int_type", "pubmed"]
    km1 = keep_cols + ["DB_ID"]
    km2 = keep_cols + ["DB_ID_homology_1", "DB_ID_homology_2"]
    km3 = keep_cols + ["DB_ID_homology_1", "DB_ID_homology_2", "DB_ID"]
    keep_final = [
        "DB_ID_1", "DB_ID_2", "DB_ID_source_1", "DB_ID_source_2",
        "evalue", "quality", "DB_ID_org",  "pubmed",
        "DB_ID_homology_1", "DB_ID_homology_2",
        "DB_ID_detection_type", "DB_ID_interaction_type"
    ]
    hom_cols = ["DB_ID", "DB_ID_source", "DB_ID_target"]
    interologs = (interologs[keep_cols]
                  .merge(homologs[hom_cols],
                         left_on=["DB_ID_1", "DB_ID_source_1"],
                         right_on=["DB_ID_target", "DB_ID_source"])[km1]
                  .merge(homologs[hom_cols],
                         left_on=["DB_ID_2", "DB_ID_source_2"],
                         right_on=["DB_ID_target", "DB_ID_source"],
                         suffixes=["_homology_1", "_homology_2"])[km2]
                  .merge(mi_ontology,
                         left_on="det_type",
                         right_on="mi_id")[km3]
                  .merge(mi_ontology,
                         left_on="int_type",
                         right_on="mi_id",
                         suffixes=["_detection_type", "_interaction_type"])
                  .merge(tax_info_keep[["DB_ID", "Taxon Id"]],
                         left_on="org_id",
                         right_on="Taxon Id")
                  .rename(columns={"DB_ID": "DB_ID_org"}))[keep_final]

    logger.info("Merging experimental interactions...")
    interactions = (interactions
                    .merge(proteins[["accession", "DB_ID"]],
                           left_on="p1", right_on="accession")
                    .merge(proteins[["accession", "DB_ID"]],
                           left_on="p2", right_on="accession",
                           suffixes=["_1", "_2"])
                    .merge(mi_ontology,
                           left_on="det_type",
                           right_on="mi_id")
                    .merge(mi_ontology,
                           left_on="int_type",
                           right_on="mi_id",
                           suffixes=["_detection_type", "_interaction_type"])
                    .merge(tax_info_keep[["DB_ID", "Taxon Id"]],
                           left_on="org_id",
                           right_on="Taxon Id")
                    .rename(columns={"DB_ID": "DB_ID_org"})
                    )[["DB_ID_1", "DB_ID_2", "DB_ID_org",
                       "DB_ID_interaction_type", "DB_ID_detection_type",
                       "pubmed"]]
    target_org_id = interactions["DB_ID_org"].unique()[0]

    source_interactions = interologs[[
        "DB_ID_source_1", "DB_ID_source_2",
        "DB_ID_org", "pubmed", "DB_ID_detection_type",
        "DB_ID_interaction_type"]]
    source_interactions = (source_interactions
                           .drop_duplicates()
                           .reset_index(drop=True)
                           .rename(columns={
                               "DB_ID_source_1": "DB_ID_1",
                               "DB_ID_source_2": "DB_ID_2"}))
    interactions = pd.concat([
        interactions, source_interactions]).reset_index(names="DB_ID")
    interactions["DB_ID"] += 1
    interactions["PB_ID"] = interactions["DB_ID"]
    interactions["interaction_type"] = 0
    interactions["pubmed"] = interactions["pubmed"].astype(str)
    interactions_db_file = outdir / f"{alias}-db-interactions.tsv"
    logger.info(f"Writing interactions into {interactions_db_file}...")
    w_cols = [
        "DB_ID", "interaction_type", "PB_ID", "DB_ID_1", "DB_ID_2", "DB_ID_org"
    ]
    interactions[w_cols].to_csv(interactions_db_file, sep="\t", index=False)

    logger.info("Extracting evidence metadata...")
    evidence = get_evidence(interactions)

    evidence_db_file = outdir / f"{alias}-db-evidence.tsv"
    logger.info(f"Writing evidence into {evidence_db_file}...")
    evidence.to_csv(evidence_db_file, sep="\t", index=False)

    predicted_interactions = interologs[
        ["DB_ID_1", "DB_ID_2", "DB_ID_source_1", "DB_ID_source_2", "evalue",
         "quality", "DB_ID_homology_1", "DB_ID_homology_2"]]
    predicted_interactions = predicted_interactions.merge(
        interactions[["DB_ID", "DB_ID_1", "DB_ID_2"]].rename(
            columns={"DB_ID_1": "source_ID_1", "DB_ID_2": "source_ID_2"}),
        left_on=["DB_ID_source_1", "DB_ID_source_2"],
        right_on=["source_ID_1", "source_ID_2"]).rename(
            columns={"DB_ID": "DB_ID_exp_interaction"}).reset_index(
                names="DB_ID")
    predicted_interactions["DB_ID"] += 1
    predicted_interactions["interaction_type"] = 1
    predicted_interactions["DB_ID_organism"] = target_org_id
    predicted_interactions["is_best"] = False

    interolog_db_file = outdir / f"{alias}-db-interologs.tsv"
    logger.info(f"Writing interologs into {interolog_db_file}...")
    predicted_interactions.to_csv(interolog_db_file, sep="\t", index=False)

    logger.info(f"Loading protein complexes from {complexes_file}...")
    complexes, complexes_protein = get_complex_data(complexes_file)
    complexes_protein = complexes_protein.merge(
        proteins[["DB_ID", "accession"]],
        left_on="protein",
        right_on="accession").drop(
                columns=["accession", "protein"]).rename(
                        columns={"DB_ID": "DB_ID_protein"})

    complexes_db_file = outdir / f"{alias}-db-complexes.tsv"
    logger.info(f"Writing complexes into {complexes_db_file}...")
    complexes[["Cluster", "Size", "Density", "Quality", "P-value"]].to_csv(
        complexes_db_file, sep="\t", index=False
    )

    cprots_db_file = outdir / f"{alias}-db-complexes_proteins.tsv"
    logger.info(f"Writing complexes proteins into {cprots_db_file}...")
    complexes_protein.to_csv(cprots_db_file, sep="\t", index=False)

    logger.info("Building complex-interolog relationships...")
    complex_interolog = {k: [] for k in ["complex_id", "interolog_id"]}
    for c_id in complexes_protein["complex_id"].unique():
        cond = complexes_protein["complex_id"] == c_id
        prots = complexes_protein[cond]["DB_ID_protein"].unique()
        cond = (predicted_interactions["DB_ID_1"].isin(prots) &
                predicted_interactions["DB_ID_2"].isin(prots))
        interologs = predicted_interactions[cond]["DB_ID"].unique()
        for i in interologs:
            complex_interolog["complex_id"].append(c_id)
            complex_interolog["interolog_id"].append(i)
    complex_interolog = pd.DataFrame(complex_interolog)

    c_interolog_db_file = outdir / f"{alias}-db-comp_interologs.tsv"
    logger.info(f"Writing complex-interologs into {c_interolog_db_file}...")
    complex_interolog.to_csv(c_interolog_db_file, sep="\t", index=False)

    logger.info("Building complex-interolog relationships...")
    complex_interaction = {k: [] for k in ["complex_id", "interaction_id"]}
    for c_id in complexes_protein["complex_id"].unique():
        cond = complexes_protein["complex_id"] == c_id
        prots = complexes_protein[cond]["DB_ID_protein"].unique()
        cond = (interactions["DB_ID_1"].isin(prots) &
                interactions["DB_ID_2"].isin(prots))
        interologs = interactions[cond]["DB_ID"].unique()
        for i in interologs:
            complex_interaction["complex_id"].append(c_id)
            complex_interaction["interaction_id"].append(i)
    complex_interaction = pd.DataFrame(complex_interaction)

    c_interaction_db_file = outdir / f"{alias}-db-comp_interactions.tsv"
    logger.info(f"Writing complex-interactions into {c_interolog_db_file}...")
    complex_interaction.to_csv(c_interaction_db_file, sep="\t", index=False)

    logger.info(f"Loading functional info from {overrepresentation_file}...")
    overrep = pd.read_table(overrepresentation_file,
                            names=["complex_id", "goterm", "pvalue"])
    overrep["goterm"] = overrep["goterm"].apply(
        lambda x: int(x.split(":")[1]))
    overrep = overrep.merge(
        go_df[["db_id", "go_id"]],
        left_on="goterm",
        right_on="go_id")[["complex_id", "db_id", "pvalue"]]
    overrep_db_file = outdir / f"{alias}-db-overrep.tsv"
    logger.info(f"Writing functional info into {overrep_db_file}...")
    overrep.to_csv(overrep_db_file, sep="\t", index=False)

    logger.info("Done")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="generates files for ProteoBOOSTER-web")
    parser.add_argument("alias", help="asd")
    parser.add_argument("obo_file", help="asd")
    parser.add_argument("goa_file", help="asd")
    parser.add_argument("interolog_file", help="asd")
    parser.add_argument("interactions_file", help="asd")
    parser.add_argument("homologs_file", help="asd")
    parser.add_argument("proteins_info_file", help="asd")
    parser.add_argument("complexes_file", help="asd")
    parser.add_argument("overrepresentation_file", help="asd")
    parser.add_argument("mi_obo", help="asd")
    parser.add_argument("taxonomy_file", help="asd")
    parser.add_argument("output_dir", help="asd")

    args = parser.parse_args()

    run(args.alias,
        args.obo_file,
        args.goa_file,
        args.interolog_file,
        args.interactions_file,
        args.proteins_info_file,
        args.homologs_file,
        args.complexes_file,
        args.overrepresentation_file,
        args.mi_obo,
        args.taxonomy_file,
        args.output_dir)
