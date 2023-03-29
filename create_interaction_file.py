#!/usr/bin/env python
"""
Usage: create_interaction_file.py [snapshot_subdirectory] [interaction_files_directory]

Creates a file containing all the unique interactions found in a snapshot.
"""

import os
import re
import shutil
from collections import defaultdict
import operator
import urllib.request
import requests
import sys
import logging
from dataclasses import dataclass


@dataclass
class OrganismInfo:
    code: str
    taxon: str
    kingdom: str
    name: str
    common_name: str
    synonym: str


class CreateInteractionFile:
    
    def __init__(self, snapshot_subdirectory, interaction_file_subdirectory, cache_domain_file, logger):
        #self.proteomes = defaultdict(set)
        self.snapshot_subdirectory = snapshot_subdirectory
        self.interaction_file_subdirectory = interaction_file_subdirectory
        self.unique_interaction_file_entries = defaultdict(set)
        self.cache_domain_file = cache_domain_file
        self.cache_domain = dict()
        self.logger = logger
        self.valid_interaction_detection_methods = set(
            "MI:0004 MI:0006 MI:0007 MI:0008 MI:0009 MI:0010 MI:0011 MI:0012 MI:0013 MI:0014 "
            "MI:0016 MI:0017 MI:0018 MI:0019 MI:0020 MI:0027 MI:0028 MI:0029 MI:0030 MI:0031 "
            "MI:0034 MI:0038 MI:0040 MI:0041 MI:0042 MI:0043 MI:0045 MI:0047 MI:0048 MI:0049 "
            "MI:0051 MI:0052 MI:0053 MI:0054 MI:0055 MI:0065 MI:0066 MI:0067 MI:0069 MI:0071 "
            "MI:0073 MI:0077 MI:0081 MI:0084 MI:0089 MI:0090 MI:0091 MI:0092 MI:0095 MI:0096 "
            "MI:0097 MI:0098 MI:0099 MI:0104 MI:0107 MI:0108 MI:0111 MI:0112 MI:0114 MI:0115 "
            "MI:0225 MI:0226 MI:0227 MI:0231 MI:0232 MI:0254 MI:0255 MI:0256 MI:0257 MI:0276 "
            "MI:0369 MI:0370 MI:0397 MI:0398 MI:0399 MI:0400 MI:0401 MI:0402 MI:0404 MI:0405 "
            "MI:0406 MI:0410 MI:0411 MI:0412 MI:0413 MI:0415 MI:0416 MI:0417 MI:0419 MI:0420 "
            "MI:0423 MI:0424 MI:0425 MI:0426 MI:0428 MI:0430 MI:0432 MI:0434 MI:0435 MI:0437 "
            "MI:0438 MI:0439 MI:0440 MI:0441 MI:0508 MI:0509 MI:0510 MI:0511 MI:0512 MI:0513 "
            "MI:0514 MI:0515 MI:0516 MI:0588 MI:0602 MI:0603 MI:0604 MI:0605 MI:0606 MI:0655 "
            "MI:0657 MI:0663 MI:0676 MI:0678 MI:0695 MI:0696 MI:0697 MI:0698 MI:0699 MI:0700 "
            "MI:0726 MI:0727 MI:0728 MI:0729 MI:0807 MI:0808 MI:0809 MI:0813 MI:0814 MI:0824 "
            "MI:0825 MI:0826 MI:0827 MI:0841 MI:0858 MI:0859 MI:0870 MI:0872 MI:0879 MI:0880 "
            "MI:0887 MI:0888 MI:0889 MI:0891 MI:0892 MI:0893 MI:0894 MI:0895 MI:0899 MI:0900 "
            "MI:0901 MI:0905 MI:0916 MI:0920 MI:0921 MI:0928 MI:0938 MI:0943 MI:0944 MI:0946 "
            "MI:0947 MI:0949 MI:0953 MI:0963 MI:0964 MI:0965 MI:0966 MI:0968 MI:0969 MI:0972 "
            "MI:0976 MI:0979 MI:0982 MI:0983 MI:0984 MI:0989 MI:0990 MI:0991 MI:0992 MI:0993 "
            "MI:0994 MI:0995 MI:0996 MI:0997 MI:0998 MI:0999 MI:1000 MI:1001 MI:1002 MI:1003 "
            "MI:1004 MI:1005 MI:1006 MI:1007 MI:1008 MI:1009 MI:1010 MI:1011 MI:1016 MI:1017 "
            "MI:1019 MI:1022 MI:1024 MI:1026 MI:1028 MI:1029 MI:1030 MI:1031 MI:1034 MI:1035 "
            "MI:1036 MI:1037 MI:1038 MI:1086 MI:1087 MI:1088 MI:1089 MI:1103 MI:1104 MI:1111 "
            "MI:1112 MI:1113 MI:1137 MI:1138 MI:1142 MI:1145 MI:1147 MI:1183 MI:1184 MI:1187 "
            "MI:1189 MI:1190 MI:1191 MI:1192 MI:1203 MI:1204 MI:1211 MI:1218 MI:1219 MI:1229 "
            "MI:1232 MI:1235 MI:1236 MI:1238 MI:1246 MI:1247 MI:1249 MI:1252 MI:1309 MI:1311 "
            "MI:1312 MI:1313 MI:1314 MI:1320 MI:1321 MI:1325".split())
        self.valid_interaction_types = set(
            "MI:0192 MI:0193 MI:0194 MI:0195 MI:0197 MI:0198 MI:0199 MI:0200 MI:0201 MI:0202 "
            "MI:0203 MI:0204 MI:0206 MI:0207 MI:0209 MI:0210 MI:0211 MI:0212 MI:0213 MI:0214 "
            "MI:0216 MI:0217 MI:0220 MI:0407 MI:0408 MI:0414 MI:0556 MI:0557 MI:0558 MI:0559 "
            "MI:0566 MI:0567 MI:0568 MI:0569 MI:0570 MI:0571 MI:0572 MI:0701 MI:0844 MI:0871 "
            "MI:0881 MI:0882 MI:0883 MI:0902 MI:0910 MI:0914 MI:0915 MI:0945 MI:0971 MI:0985 "
            "MI:0986 MI:0987 MI:1027 MI:1126 MI:1127 MI:1139 MI:1140 MI:1143 MI:1146 MI:1148 "
            "MI:1230 MI:1237 MI:1250 MI:1251 MI:1310 MI:1327".split())
    
    def load_proteomes(self):
        c = 0
        for line in open(os.path.join(self.snapshot_subdirectory, 'proteomes.tab')):
            fields = line.split('\t')
            c += 1
            print('\r', c, len(fields), end='')
            for protein in fields[6].split(','):
                self.proteomes[protein].add(fields[0])

    def create_interaction_file_subdirectory(self):
        """ Creates the interaction directory if it does not exist
        """
        if not os.path.exists(self.interaction_file_subdirectory):
            os.makedirs(self.interaction_file_subdirectory)
    
    def create_biogrid_with_uniprot_accessions(self):
        """ Creates BioGrid file in which UniProtKB accessions are used
            (the default one only uses genes).
        """

        biogrid_file = os.path.join(self.snapshot_subdirectory, "biogrid.tab")
        biogrid_with_uniprot_file = os.path.join(self.interaction_file_subdirectory,
                                                 "biogrid_with_uniprot_accessions.tab")
        if not os.path.exists(biogrid_with_uniprot_file):
            # 1.- we read the mapping filename into a dictionary
            mapping_filename = os.path.join(self.snapshot_subdirectory, "idmapping_selected.tab")
            entrez_to_uniprot_mapping_dictionary = {}
            for line in open(mapping_filename):
                if line.startswith("#"):
                    continue
                fields = line.split('\t')
                entrez_to_uniprot_mapping_dictionary[fields[2]] = fields[0]

            # 2.- we read the biogrid file and we write a translated version (biogrid_with_uniprot_accessions) 
            # in which no gene names are used
            with open(biogrid_with_uniprot_file, "w") as out:
                for num, line in enumerate(open(biogrid_file)):
                    if num == 0:
                        out.write(line)
                        continue
                    fields = line.split('\t')
                    for i in (0, 1):
                        entrez_identifier = fields[i][fields[i].find("entrez gene/locuslink:") + 22:]
                        if entrez_identifier in entrez_to_uniprot_mapping_dictionary:
                            fields[i] = "{0}|uniprotkb:{1}".\
                                format(fields[i], entrez_to_uniprot_mapping_dictionary[entrez_identifier])
                    out.write("\t".join(fields))
    
    def create_interaction_file(self):
        """ Creates a unique interaction file (tab-separated) with the following format
            for each line:
                [PROTEIN1] [PROTEIN2] [NCBI_ID] [DETECTION_TYPE] [INTERACTION_TYPE] [LIST PUBMED IDS]
            where PROTEIN1 and PROTEIN 2 are two (distinct) UniProtKB accession numbers of 
            proteins which belong to the same organism (identified by its NCBI_ID), 
            the DETECTION_TYPE is the detection method used in the experiment
            (term in the MI ontology), the interaction type is the actual type 
            of interaction (term in the MI ontology), and the LIST PUBMED IDS 
            is a CSV list of the PubMed ids that support that particular interaction.

            Note that a given pair of proteins may appear as different interactions (lines in 
            the file), but a given pair of proteins plus a detection type and an interaction
            type should not appear more than once in the interaction file as all the appearances
            should be included as different PubMed ids (different experimental evidences). 
        """
        outfile = os.path.join(self.interaction_file_subdirectory, "interaction_file.tab")
        if not os.path.exists(outfile):
            self.all_prots = {} 
            self.__parse_interaction_database(os.path.join(self.interaction_file_subdirectory,
                                                           "biogrid_with_uniprot_accessions.tab"))
            self.__parse_interaction_database(os.path.join(self.snapshot_subdirectory, "intact.tab"))
            self.__parse_interaction_database(os.path.join(self.snapshot_subdirectory, "mint.tab"))
            self.__parse_interaction_database(os.path.join(self.snapshot_subdirectory, "dip.tab"))
            self.__get_uniprot_ncbi_ids()
            with open(outfile, 'w') as interaction_file:
                for _unique_interaction_file_entry in self.unique_interaction_file_entries:
                    fields = _unique_interaction_file_entry.split()
                    if fields[0] not in self.all_prots or fields[1] not in self.all_prots:
                        continue
                    elif self.all_prots[fields[0]] != self.all_prots[fields[1]] or self.all_prots[fields[0]] == -1 or self.all_prots[fields[1]] == -1:
                        continue
                    unique_interaction_file_entry = "{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1],
                        self.all_prots[fields[0]], fields[2], fields[3]                
                    )
                    pubmeds = ",".join(sorted(self.unique_interaction_file_entries[_unique_interaction_file_entry]))
                    interaction_file.write("{}\t{}\n".format(unique_interaction_file_entry, pubmeds))
        else:
            self.logger.info("interaction file exists, skipping")

    def __get_uniprot_ncbi_ids(self):
        found = set()
        notDigit = set()
        lineSmall = 0
        for line in open(os.path.join(self.snapshot_subdirectory, "idmapping_selected.tab")):
            fields = line.strip().split("\t")
            ## small fix to make it work with the new snapshot
            if len(fields) > 12:
                if fields[12].isdigit():
                    uniprot_id, ncbi_id = fields[0], int(fields[12])
                    if uniprot_id in self.all_prots:
                        found.add(uniprot_id)
                        self.all_prots[uniprot_id] = ncbi_id
                else:
                    notDigit.add(fields[12])
            else:
                lineSmall += 1
        print("Not found: ", len(self.all_prots) - len(found))
        print("Not Digit: ", len(self.all_prots) - len(notDigit))
        print("Line Small: ", lineSmall)
    
    def __parse_interaction_database(self, interaction_database_path):
        with open(interaction_database_path, encoding='utf-8') as interaction_database:
            interaction_database.readline()
            good_interactions = 0
            no_uniprot = 0
            lineno = -1
            try:
                for lineno, line in enumerate(interaction_database):
                    fields = line.split('\t')
                    
                    identifier_start_increment = 0
                    identifier_start = -1
                    while identifier_start_increment != -1:
                        identifier_start += 1
                        identifier_start_increment = fields[6][identifier_start:].upper().find("MI:")
                        identifier_start += identifier_start_increment
                        if fields[6][identifier_start + 3:identifier_start + 7].isdigit():
                            interaction_detection_method = fields[6][identifier_start:identifier_start + 7].upper()
                            break
                    # (a) discard interactions where the detection method is not valid
                    # (e.g. computational)
                    if not interaction_detection_method in self.valid_interaction_detection_methods:
                        continue

                    identifier_start_increment = 0
                    identifier_start = -1
                    while identifier_start_increment != -1:
                        identifier_start += 1
                        identifier_start_increment = fields[11][identifier_start:].upper().find("MI:")
                        identifier_start += identifier_start_increment
                        if fields[11][identifier_start + 3:identifier_start + 7].isdigit():
                            interaction_type = fields[11][identifier_start:identifier_start + 7].upper()
                            break
                    # (b) discard interactions which are not of the valid types
                    # (e.g. genetic interactions)
                    if not interaction_type in self.valid_interaction_types:
                        continue

                    # (c) only entries having UniProtKB accessions are processed
                    if "uniprotkb:" in fields[0] and "uniprotkb:" in fields[1]:
                        uniprot_accessions = []
                        for i in range(2):
                            identifiers = fields[i].split('|')
                            for identifier in identifiers:
                                if "uniprotkb:" in identifier:
                                    accession = identifier[10:]
                                    if "-PRO_" in accession:
                                        accession = accession.split('-PRO_')[0]
                                    uniprot_accessions.append(accession)
                                    break

                        # (d) self-interactions are discarded, non-mappable interactions or interactions
                        # from a different organism
                        if uniprot_accessions[0] == uniprot_accessions[1]:
                            continue

                        # (e) discard interactions which are not within the same proteome
                        # if len(self.proteomes[uniprot_accessions[0]] | self.proteomes[uniprot_accessions[1]]) < 1:
                        #     continue

                        good_interactions += 1
                        uniprot_accessions = sorted(uniprot_accessions)
                        #organism_ncbi_identifier = re.findall(r'\d+', fields[9])[0]
                        publication_identifiers = fields[8].split('|')
                        unique_pubmed_identifiers = set()
                        for publication_identifier in publication_identifiers:
                            if publication_identifier.startswith("pubmed:") and publication_identifier[7:].isdigit():
                                unique_pubmed_identifiers.add(publication_identifier[7:])
                        interaction_file_entry = "{0}\t{1}\t{2}\t{3}".format(
                            uniprot_accessions[0], uniprot_accessions[1],
                            interaction_detection_method, interaction_type)
                        self.unique_interaction_file_entries[interaction_file_entry].update(unique_pubmed_identifiers)
                        self.all_prots[uniprot_accessions[0]] = -1
                        self.all_prots[uniprot_accessions[1]] = -1
                    else:
                        no_uniprot += 1
            except UnicodeDecodeError:
                self.logger.info(f"lineno: {lineno}")
                sys.exit(1)
        self.logger.info("\tprocessed file {}: {} good interactions found, {} without uniprot ids".format(
            interaction_database_path, good_interactions, no_uniprot
        ))


    def parse_specfile(self):

        with open(os.path.join(self.snapshot_subdirectory, "speclist.txt")) as f:
            section = "preamble"
            starts_with_equal_counter = 0
            starts_with_dash_counter = 0
            org = ""
            taxon = ""
            kingdom = ""
            name = ""
            common_name = ""
            synonym = ""
            organisms = {}
            for lineno, line in enumerate(f):
                try:
                    if line.startswith("="):
                        starts_with_equal_counter += 1
                        if starts_with_equal_counter == 2:
                            section = "real organisms"
                            starts_with_dash_counter = 0
                        elif starts_with_equal_counter == 4:
                            section = "virtual organisms"
                    elif section in ["real organisms", "virtual organisms"]:
                        if line[0] in ["-", "_"]:
                            starts_with_dash_counter += 1
                            continue
                        if 2 > starts_with_dash_counter >= 1:
                            if line.startswith(" "):
                                if org:
                                    k, v = line.strip().split("=")
                                    match k:
                                        case "C":
                                            common_name = v
                                        case "S":
                                            synonym = v
                            elif line[0].isalnum():
                                if org:
                                    o = OrganismInfo(org, taxon, kingdom, name, common_name, synonym)
                                    organisms[int(taxon)] = o
                                org, kingdom, taxon, name = line.strip().split(maxsplit=3)
                                taxon = taxon[:-1]
                                name = name[2:]
                                common_name = ""
                                synonym = ""
                    else:
                        continue
                except ValueError:
                    print(f"LN = {lineno}")
                    sys.exit(1)
            return organisms
            
    
    def create_organism_statistics(self):
        org_info = self.parse_specfile()
        organism_interaction_amounts = defaultdict(int)
        for line in open(os.path.join(self.interaction_file_subdirectory, "interaction_file.tab")):
            organism_ncbi_identifier = line.split('\t')[2]
            organism_interaction_amounts[int(organism_ncbi_identifier)] += 1
           
        organism_interaction_amounts = sorted(organism_interaction_amounts.items(), key=operator.itemgetter(0))
        with open(os.path.join(self.interaction_file_subdirectory, "organism_statistics.tab"), 'w') \
                as organism_statistics:
            for organism_interaction_amount in organism_interaction_amounts:
                ncbi_id, quantity = organism_interaction_amount
                organism_domain = org_info[ncbi_id].kingdom
                organism_name = org_info[ncbi_id].name
                organism_statistics.write("{0}\t{1}\t{2}\t{3}\n".format(
                    ncbi_id, organism_name, organism_domain, quantity))

    def create_unique_interactions_file(self):
        """ Reads the interaction_file.tab and produce another 
            (unique_interactions_file.tab) where only 4 fields
            are kept: 2 accessions, ncbi ID and comma-separated
            list of pubmed ids.
        """
        interactions = os.path.join(self.interaction_file_subdirectory, "interaction_file.tab")
        unique_interactions = os.path.join(self.interaction_file_subdirectory, "unique_interaction_file.tab")
        if not os.path.exists(unique_interactions):
            with open(unique_interactions, "w") as out:
                data = defaultdict(set)
                for line in open(interactions):
                    fields = line.strip().split("\t")
                    # we discard interactions without PubMed ids
                    if len(fields) < 6:
                        continue
                    data[(fields[0], fields[1], fields[2])].update(map(int, fields[5].split(",")))
                for uni1, uni2, ncbi_id in sorted(data.keys(), key=lambda x: int(x[2])):
                    pubmeds = ",".join(map(str, sorted(data[(uni1, uni2, ncbi_id)])))
                    out.write("{}\t{}\t{}\t{}\n".format(uni1, uni2, ncbi_id, pubmeds)) 
        else:
            self.logger.info("unique interaction file exists")

    def __parse_sequence_file(self, file_name, valid_sequences, out):
        """ Parses a Fasta UniProtKB-like file identified by `file_name`,
            and extracts the sequences whose accessions match those in
            `valid_sequences` to the stream represented by `out`. The
            set of found accessions is returned.
        """
        found_accessions = set()
        with open(file_name) as file_input:
            line = file_input.readline()
            while line:
                uniprot_accession = line.split("|")[1]
                old_line = line
                line = file_input.readline()
                if uniprot_accession in valid_sequences:
                    out.write(old_line)
                    while not line.startswith(">") and line.rstrip():
                        out.write(line)
                        line = file_input.readline()
                    found_accessions.add(uniprot_accession)
                else:
                    while not line.startswith(">") and line.rstrip():
                        line = file_input.readline()
        return found_accessions

    def filter_interaction_file(self):
        """ Filters the interaction file. Only those entries having 
            UniProtKB accessions which were found in the sequence
            file will be kept. The old (unfiltered) file will be
            renamed to "interaction_file_unfiltered.tab".
        """
        unfiltered = os.path.join(self.interaction_file_subdirectory, "interaction_file_unfiltered.tab")
        filtered = os.path.join(self.interaction_file_subdirectory, "interaction_file.tab")
        if os.path.exists(filtered) and os.path.exists(unfiltered):
            self.logger.info("filtered file exists, skipping computation")
        else:
            os.rename(filtered, unfiltered)

            num_filtered = 0
            unique_pairs = set()
            with open(filtered, "w") as out:
                for line in open(unfiltered):
                    fields = line.strip().split("\t")
                    if fields[0] not in self.not_found_accs and fields[1] not in self.not_found_accs:
                        out.write(line)
                    else:
                        num_filtered += 1
                        unique_pairs.add((fields[0], fields[1]))
            self.logger.info("Number of interactions filtered: {} ({} unique pairs)".format(num_filtered,
                                                                                            len(unique_pairs)))
    
    def create_sequences(self):

        seqfile = os.path.join(self.interaction_file_subdirectory, "sequences.fasta")

        if not os.path.exists(seqfile):
            # 1.- UniProtKB accessions are retrieved from all possible interactions
            unique_uniprot_accessions = set()
            for line in open(os.path.join(self.interaction_file_subdirectory, "interaction_file.tab")):
                fields = line.split('\t')
                for field in fields[0:2]:
                    unique_uniprot_accessions.add(field)

            # 2.- The sequences corresponding to those accessions are extracted 
            # into "sequences.fasta"
            found_accessions = set()
            with open(seqfile, 'w') as sequences:
                self.logger.info("... parsing Swiss-Prot")
                swissprot_file = os.path.join(self.snapshot_subdirectory, "swissprot.fasta")
                found_accessions |= self.__parse_sequence_file(swissprot_file, unique_uniprot_accessions, sequences)
                self.logger.info("... parsing TremBL")
                trembl_file = os.path.join(self.snapshot_subdirectory, "trembl.fasta")
                found_accessions |= self.__parse_sequence_file(trembl_file, unique_uniprot_accessions, sequences)

            # 3.- The list of not-found accessions is dumped into a file (acc_not_found.txt)
            self.logger.info("Writing the file with not-found accessions")
            with open(os.path.join(self.interaction_file_subdirectory, "acc_not_found.txt"), 'w') as out:
                self.not_found_accs = unique_uniprot_accessions - found_accessions
                out.write("{}".format("\n".join(sorted(self.not_found_accs))))
            self.logger.info("... {} proteins were not found".format(len(self.not_found_accs)))
        else:
            self.logger.info("sequences file exists, skipping computation")

def main():
    if len(sys.argv) != 3:
        print("ERROR: Wrong number of arguments")
        print(__doc__)
    elif not os.path.exists(sys.argv[1]) \
        or not os.path.exists(sys.argv[2]):
        print("ERROR: A directory you have given does not exist")
        print(__doc__)
    else:
        snapshot_subdirectory = sys.argv[1]
        interaction_file_subdirectory = os.path.join(sys.argv[2], os.path.basename(
            os.path.normpath(snapshot_subdirectory)))
        cache_domain_file = os.path.join(sys.argv[2], "cache_domain_file")
        logger = logging.getLogger('create_interaction_file')
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M:%S")
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        
        interaction_file_creator = CreateInteractionFile(snapshot_subdirectory, 
                                    interaction_file_subdirectory, cache_domain_file, logger)
        
        # 1.- the interaction directory is created
        logger.info("creating interaction file subdirectory")
        interaction_file_creator.create_interaction_file_subdirectory()
        # 2.- Biogrid file is translated into a version with UniProt accessions
        logger.info("creating BioGRID with UniProt accessions")
        interaction_file_creator.create_biogrid_with_uniprot_accessions()
        # # 2.5 - Load proteomes file
        # logger.info('loading proteomes file')
        # interaction_file_creator.load_proteomes()
        # 3.- Create the interaction file
        logger.info("creating interaction file")
        interaction_file_creator.create_interaction_file()
        # 4.- Sequences participating in the extracted
        # interactions are extracted
        logger.info("creating sequences")
        interaction_file_creator.create_sequences()
        # 5.- The interaction file is filtered to remove 
        # interactions with non valid sequence ids.
        interaction_file_creator.filter_interaction_file()
        # 6.- A file with unique interactions is also produced
        interaction_file_creator.create_unique_interactions_file()
        # 7.- Statistics for every organisms are produced
        logger.info("creating organism statistics")
        interaction_file_creator.create_organism_statistics()

if __name__ == '__main__':
    main()
