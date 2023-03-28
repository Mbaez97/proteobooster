#!/usr/bin/env python
"""
Downloads a snapshot of the data input files (PPI + sequences + other stuff)
which will to be used for prediction.
"""

import argparse
import os
import zipfile
import gzip
import sys
import datetime
import time
import logging
import urllib.request
import shutil
import tarfile
import requests
import json

class DownloadSnapshot:

    def __init__(self, 
                 snapshot_subdirectory, 
                 dip_user,
                 dip_pass,
                 overwrite=False):
        self.snapshot_subdirectory = snapshot_subdirectory
        self.overwrite = overwrite
        self.dip_user = dip_user
        self.dip_pass = dip_pass

    def create_snapshot_subdirectory(self):
        if not os.path.exists(self.snapshot_subdirectory):
            os.makedirs(self.snapshot_subdirectory)

    def download_biogrid(self):
        url = "http://thebiogrid.org/downloads/archives/" +\
                "Latest%20Release/BIOGRID-ALL-LATEST.mitab.zip"
        local_file_name = "biogrid.zip"
        uncompressed = os.path.join(self.snapshot_subdirectory,
                                    "biogrid.tab"
                                    )
        if not os.path.isfile(uncompressed) or self.overwrite:
            full_file_path = self.__download_file(url, local_file_name)
            self.__unzip_file(
                full_file_path,
                "BIOGRID-ALL-",
                full_file_path.replace(
                    ".zip",
                    ".tab"))

    def download_intact(self):
        url = "ftp://ftp.ebi.ac.uk/pub/databases/intact/" +\
                "current/psimitab/intact.zip"
        local_file_name = "intact.zip"
        uncompressed = os.path.join(self.snapshot_subdirectory,
                                    "intact.tab"
                                    )
        if not os.path.isfile(uncompressed) or self.overwrite:
            full_file_path = self.__download_file(url, local_file_name)
            self.__unzip_file(
                full_file_path,
                "intact.txt",
                full_file_path.replace(
                    ".zip",
                    ".tab"))

    def download_mint(self):
        url = "http://mint.bio.uniroma2.it/mitab/MINT_MiTab.txt"
        url = "http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/*"
        self.__download_file(url, "mint.tab")

    def download_dip(self):
        local_uncompressed_data = os.path.join(self.snapshot_subdirectory, "dip.tab")
        if os.path.isfile(local_uncompressed_data) and not self.overwrite:
            return
        url_prefix = f"ftp://{self.dip_user}:{self.dip_pass}@dip.doe-mbi.ucla.edu/"
        response = urllib.request.urlopen(url_prefix)
        most_recent_year = 0
        for file_entry in response.read().splitlines():
            filename = file_entry.split()[-1]
            if filename.isdigit() and filename > most_recent_year:
                most_recent_year = filename
        response = urllib.request.urlopen(
            "{0}{1}/tab25/".format(url_prefix, most_recent_year))
        most_recent_date_in_filename = 0
        for file_entry in response.read().splitlines():
            filename = file_entry.split()[-1]
            if filename.startswith("dip") and filename.endswith(".txt.gz"):
                date_in_filename = filename[3:-7]
                if date_in_filename.isdigit() and\
                        date_in_filename > most_recent_date_in_filename:
                    most_recent_date_in_filename = date_in_filename
        url = "{0}{1}/tab25/dip{2}.txt.gz".format(
            url_prefix,
            most_recent_year,
            most_recent_date_in_filename)
        local_file = "dip.gz"
        local_file_name = self.__download_file(url, local_file)
        self.__gunzip_file(local_file_name, local_uncompressed_data)

    def download_biogrid_to_uniprot_mapping(self):
        url = "http://thebiogrid.org/downloads/archives/"\
                "External%20Database%20Builds/UNIPROT.tab.txt"
        local_file = "biogrid_to_uniprot_mapping.tab"
        self.__download_file(url, local_file)

    def download_uniprot_mapping(self):
        url = "ftp://ftp.uniprot.org/pub/databases/uniprot/"\
                "current_release/knowledgebase/idmapping/"\
                "idmapping_selected.tab.gz"
        local_file_name = "idmapping_selected.tab.gz"
        local_uncompressed_data = os.path.join(
            self.snapshot_subdirectory, "idmapping_selected.tab")
        if not os.path.isfile(local_uncompressed_data) or self.overwrite:
            local_compressed_data = self.__download_file(url, local_file_name)
            self.__gunzip_file(local_compressed_data, local_uncompressed_data)

    def download_swissprot(self):
        url = "ftp://ftp.uniprot.org/pub/databases/uniprot/"\
                "current_release/knowledgebase/complete/"\
                "uniprot_sprot.fasta.gz"
        local_file_name = "swissprot.gz"
        local_uncompressed_data = os.path.join(
            self.snapshot_subdirectory, "swissprot.fasta")
        if not os.path.isfile(local_uncompressed_data) or self.overwrite:
            local_compressed_data = self.__download_file(url, local_file_name)
            self.__gunzip_file(local_compressed_data, local_uncompressed_data)

    def download_trembl(self):
        url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/"\
                "knowledgebase/complete/uniprot_trembl.fasta.gz"
        local_file_name = "trembl.gz"
        local_uncompressed_data = os.path.join(
            self.snapshot_subdirectory, "trembl.fasta")
        if not os.path.isfile(local_uncompressed_data) or self.overwrite:
            local_compressed_data = self.__download_file(url, local_file_name)
            self.__gunzip_file(local_compressed_data, local_uncompressed_data)

    def download_species_list(self):
        url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/"\
                "knowledgebase/complete/docs/speclist.txt"
        local_file_name = "species_list.txt"
        self.__download_file(url, local_file_name)

    def download_proteomes(self):
        out = open(os.path.join(self.snapshot_subdirectory, 'proteomes.tab'), 'w')
        requestURL = "https://www.ebi.ac.uk/proteins/api/proteomes?offset=0&size=-1&is_redundant=false"
        r = requests.get(requestURL, headers={ "Accept" : "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        responseBody = json.loads(r.text)

        columns = [
            'upid', 
            'name',
            'strain',
            'taxonomy',
            'sourceTaxonomy',
            'superregnum',
        ]

        zero_proteomes = set()

        total = len(responseBody)
        for i, p in enumerate(responseBody):
            d = {}
            for c in columns:
                try:
                    d[c] = p[c]
                except KeyError:
                    d[c] = 'none'
            # retrieve proteins for this proteome
            requestURL = 'https://www.ebi.ac.uk/proteins/api/proteomes/proteins/{upid}'.format(upid=p['upid'])
            prots_r = requests.get(requestURL, headers={ "Accept" : "application/json"})
            if not prots_r.ok:
                prots_r.raise_for_status()
                sys.exit()
            prots_responseBody = json.loads(prots_r.text)
            proteins = set()
            for component in prots_responseBody['component']:
                try:
                    for prot in component['protein']:
                        proteins.add(prot['accession'])
                except KeyError:
                    zero_proteomes.add(p['upid'])
            d['prots'] = ','.join(proteins)
            out.write('{upid}\t{name}\t{strain}\t{taxonomy}\t{sourceTaxonomy}\t{superregnum}\t{prots}\n'.format(**d))
            print('\r{i}/{total} ({perc:.2f}%)'.format(i=i, total=total, perc=(i/total)), end='')
        out.close()
        print()

    def download_uniprot_goa(self):
        url = "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/"\
                "goa_uniprot_all.gaf.gz"
        local_file_name = "uniprot_goa.gz"
        local_uncompressed_data = os.path.join(
            self.snapshot_subdirectory, "uniprot_goa.tab")
        if not os.path.isfile(local_uncompressed_data) or self.overwrite:
            local_compressed_data = self.__download_file(url, local_file_name)
            self.__gunzip_file(local_compressed_data, local_uncompressed_data)

    def download_gene_ontology(self):
        url = "http://geneontology.org/ontology/go-basic.obo"
        local_file_name = "go-basic.obo"
        self.__download_file(url, local_file_name)

    def download_taxonomy_dump(self):
        url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        local_file_name = "taxdump.tar.gz"
        local_uncompressed_data = os.path.join(
            self.snapshot_subdirectory, 
            "names.dmp")
        if not os.path.isfile(local_uncompressed_data) or self.overwrite:
            local_compressed_data = self.__download_file(url, local_file_name)
            self.__untargz_target(local_compressed_data, local_uncompressed_data, "names.dmp")

    @staticmethod
    def __unzip_file(
            full_path_zip_file,
            file_name,
            full_path_uncompressed_file,
            remove=True):
        """ Unzips a text file given its `full_path_zip_file` into the
            `full_path_uncompressed_file`.
            Only the `file_name` will be uncompressed. If `remove` is
            True, the compressed file will be removed afterwards.
        """
        with open(full_path_zip_file, 'rb') as fh:
            z = zipfile.ZipFile(fh)
            outpath = os.path.dirname(full_path_uncompressed_file)
            for name in z.namelist():
                if file_name not in name:
                    continue
                z.extract(name, outpath)
                full_new_name = os.path.join(outpath, name)
                if full_new_name != full_path_uncompressed_file:
                    os.rename(full_new_name, full_path_uncompressed_file)
        if remove:
            os.remove(full_path_zip_file)

    @staticmethod
    def __untargz_target(
            full_path_targz_file,
            full_path_uncompressed_file,
            target_file,
            remove=True):
        """ Extracts the =`target_file` from the `full_path_targz_file` 
            tar.gz provided into `full_path_uncompressed_file`. If `remove` is
            True, the compressed file will be removed, as well.
        """
        #from http://stackoverflow.com/a/37753786/943138
        with tarfile.open(full_path_targz_file, "r|*") as tar:
            counter = 0
            for member in tar:
                if member.isfile():
                    filename = os.path.basename(member.name)
                    if filename != target_file: # do your check
                        continue
                    with open(full_path_uncompressed_file, "wb") as output: 
                        shutil.copyfileobj(tar.fileobj, output, member.size)
                    break # got our file
                counter += 1
                if counter % 1000 == 0:
                    tar.members = [] # free ram... yes we have to do this manually
            tar.members = [] # free ram... yes we have to do this manually
        if remove:
            os.remove(full_path_zip_file)

    @staticmethod
    def __gunzip_file(
            full_path_gzip_file,
            full_path_uncompressed_file,
            remove=True):
        """ Gunzips a text file given its `full_path_gzip_file` into the
            `full_path_uncompressed_file`. If `remove` is True, the compressed
            file will be removed, as well.
        """
        with open(full_path_uncompressed_file, "wb") as out:
            with gzip.open(full_path_gzip_file) as f_in:
                shutil.copyfileobj(f_in, out)
        if remove:
            os.remove(full_path_gzip_file)

    def __download_file(self, url, localfile):
        """ Downloads the file pointed by `url` into a file named `localfile`.
            The full path for downloading the file is constructed by prepending
            the snapshot subdirectory to the local file.
            This full path is returned.
        """
        CHUNK = 512 * 1024
        full_path_file = os.path.join(self.snapshot_subdirectory, localfile)
        if os.path.isfile(full_path_file) and not self.overwrite:
            return full_path_file
        temp_file_name = full_path_file + ".tmp"
        response = urllib.request.urlopen(url)
        with open(temp_file_name, 'wb') as fp:
            while True:
                chunk = response.read(CHUNK)
                if not chunk:
                    break
                fp.write(chunk)
        os.rename(temp_file_name, full_path_file)
        return full_path_file

    def download_taxonomy(self):
        """ Downloads the taxonomy data from NCBI
        """
        taxdump = os.path.join(self.snapshot_subdirectory, "taxdump.tar.gz")
        urllib.request.urlretrieve(
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            taxdump)
        taxcat = os.path.join(self.snapshot_subdirectory, "taxcat.tar.gz")
        urllib.request.urlretrieve(
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.tar.gz", taxcat)
        taxo_dir = os.path.join(self.snapshot_subdirectory, "taxonomy")
        if not os.path.exists(taxo_dir):
            os.makedirs(taxo_dir)

        import tarfile
        ttaxdump = tarfile.open(taxdump, 'r:gz')
        ttaxdump.extractall(taxo_dir)
        ttaxcat = tarfile.open(taxcat, 'r:gz')
        ttaxcat.extractall(taxo_dir)
        os.remove(taxdump)
        os.remove(taxcat)


def main(args):
    if not os.path.exists(args.directory):
        print("ERROR: A directory you have given does not exist")
        print(__doc__)
    else:
        if args.snapshot_directory is None:
            snapshot_subdirectory = os.path.join(
                args.directory, datetime.datetime.fromtimestamp(
                    time.time()).strftime("%Y_%m_%d"))
        else:
            snapshot_subdirectory = os.path.join(
                args.directory, args.snapshot_directory
            )

        download_snapshot = DownloadSnapshot(snapshot_subdirectory,
                                             args.dip_user,
                                             args.dip_pass,
                                             overwrite=args.overwrite)
        logger = logging.getLogger('download_snapshot')
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            "%Y-%m-%d %H:%M:%S")
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        logger.info("creating snapshot subdirectory")
        download_snapshot.create_snapshot_subdirectory()
        logger.info("downloading BioGRID")
        download_snapshot.download_biogrid()
        logger.info("downloading IntAct")
        download_snapshot.download_intact()
        logger.info("downloading Mint")
        download_snapshot.download_mint()
        logger.info("downloading DIP")
        download_snapshot.download_dip()
        logger.info("downloading BioGRID to UniProt mapping")
        download_snapshot.download_biogrid_to_uniprot_mapping()
        logger.info("downloading UniProt Mapping")
        download_snapshot.download_uniprot_mapping()
        logger.info("downloading Swiss-Prot")
        download_snapshot.download_swissprot()
        logger.info("downloading TrEMBL")
        download_snapshot.download_trembl()
        logger.info("downloading species list")
        download_snapshot.download_species_list()
        # logger.info('downloading Non-Redundant Proteomes')
        # download_snapshot.download_proteomes()
        logger.info("downloading UniProtKB-GOA file")
        download_snapshot.download_uniprot_goa()
        logger.info("downloading Gene Ontology file")
        download_snapshot.download_gene_ontology()
        download_snapshot.download_taxonomy()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="downloads a snapshot " +
        "with all available protein interaction datasets")
    parser.add_argument(
        "directory",
        help="parent directory where the snapshot" +
        " directory will be placed")
    parser.add_argument(
        "--dip-user",
        default="DIP_USER",
        help="user to download the DIP database")
    parser.add_argument(
        "--dip-pass",
        default="DIP_PASS",
        help="Password for downloading the DIP database")
    parser.add_argument(
        "--snapshot-directory",
        "-s",
        help="directory which " +
        "will be appended to the parent directory to make the full path for " +
        "the snapshot. If not given, a timestamp will be used instead")
    parser.add_argument(
        "--overwrite",
        "-o",
        help="by default, files in the " +
        "snapshot folder will not be overwritten, allowing the continuation " +
        "of a failed download. If this flag is set, files will be " +
        "overwritten",
        action="store_true")
    args = parser.parse_args()
    main(args)

