#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHORS

    Mickael Silva
    github: @

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import os
import re
import json
import shutil
import requests
import itertools
import multiprocessing
from SPARQLWrapper import SPARQLWrapper, JSON
from urllib.parse import urlparse, urlencode, urlsplit, parse_qs

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

UNIPROT_SERVER = SPARQLWrapper("http://sparql.uniprot.org/sparql")

def verify_cpu_usage(cpu_to_use):
    """ Verify the cpu usage for chewBBACA.

        Args:
            cpu_to_use (int): the number of cpu provided to chewBBACA

        Returns:
            cpu_to_use (int): the number of cpu to use after verification

        Example:

            >>> verify_cpu_usage(6)
            6
    """
    total_cpu = multiprocessing.cpu_count()

    # do not allow a value of cpuToUse greater than the number of cores/threads
    if cpu_to_use >= total_cpu:
        print('Warning! You have provided a CPU core count value '
              'that is equal to or exceeds the number of CPU '
              'cores/threads in your machine!')
        print('Setting a different value...')
        # define a value that is safe according to the number of
        # available cores/threads
        if total_cpu > 2:
            cpu_to_use = total_cpu - 2
        elif total_cpu == 2:
            cpu_to_use = 1
        print('CPU cores/threads value set to: {0}'.format(cpu_to_use))

    elif cpu_to_use == (total_cpu - 1):
        print('Warning! You have provided a CPU core count value '
              'that is close to the maximum core count of your '
              'machine ({0}/{1}). This may affect your system '
              'responsiveness.'.format(cpu_to_use, total_cpu))

    return cpu_to_use


def check_ptf(ptf_path, main_dir):
    """
    """

    if os.path.isfile(ptf_path) is False:
        ptf_basename = os.path.basename(ptf_path)
        chewie_dir = main_dir
        ptfs_dir = os.path.join(os.path.dirname(chewie_dir),
                                'prodigal_training_files')

        ptfs = os.listdir(ptfs_dir)
        if ptf_basename in ptfs:
            chosenTrainingFile = os.path.join(ptfs_dir, ptf_basename)
        else:
            chosenTrainingFile = False

    return chosenTrainingFile


def check_prodigal_output_files(path_to_temp, list_of_genomes):
    """ Checks if Prodigal created ORF files
        equal to the number of genome files provided.

        Args:
            path_to_temp (str): the full path to the 'temp'
            directory created by chewBBACA.
            list_of_genomes (list): list containing the full
            path to the input genomes.

        Returns:
            prints a message if Prodigal created all the
            necessary files, otherwise raises a ValueError
    """

    # list ORF files created in the 'temp' directory
    listOfORFCreated = []
    for orffile in os.listdir(path_to_temp):
        if orffile.endswith('_ORF.txt'):
            listOfORFCreated.append(orffile)

    # raise exception if the number of ORF files is not equal
    # to the number of genome files provided
    if len(list_of_genomes) > len(listOfORFCreated):
        message = ('Missing some ORF files from Prodigal run. '
                   'Missing {0} ORF files out of {1} genome files '
                   'provided.'.format(len(list_of_genomes) - len(listOfORFCreated),
                                      len(list_of_genomes)))
        # remove 'temp' directory
        shutil.rmtree(path_to_temp)
        raise ValueError(message)
    else:
        print('Prodigal created all the necessary files.')


def is_fasta(filename):
    """ Checks if a file is a FASTA file.

        Args:
            filename (str): the full path to the FASTA file

        Returns:
            True if FASTA file,
            False otherwise
    """

    with open(filename, 'r') as handle:
        fasta = SeqIO.parse(handle, 'fasta')

        # returns True if FASTA file, False otherwise
        return any(fasta)


def check_input_type(input_path):
    """ Checks if the input is a file or a directory.

        Args:
            folder_or_list (str): the full path to the file or directory

        Returns:
            list_files (str) if folder_or_list is a path to a file,
            list_files (list) if folder_or_list is a path to a directory,
            Raises Exception otherwise
    """

    # check if input argument is a file or a directory
    if os.path.isfile(input_path):
        list_files = input_path

    elif os.path.isdir(input_path):

        fasta_files = []

        for genome in os.listdir(input_path):

            genepath = os.path.join(input_path, genome)

            # do not continue if genepath is a dir
            if os.path.isdir(genepath):
                continue

            # check if file is a FASTA file
            if is_fasta(genepath):
                fasta_files.append(os.path.abspath(genepath))

        # if there are FASTA files
        if fasta_files:
            # store full paths to FASTA files
            with open('listGenomes2Call.txt', 'w') as f:
                for genome in fasta_files:
                    f.write(genome + '\n')
        else:
            raise Exception('There were no FASTA files in the given '
                            'directory. Please provide a directory'
                            'with FASTA files or a file with the '
                            'list of full paths to the FASTA files.')

        list_files = 'listGenomes2Call.txt'

    else:
        raise Exception('Input argument is not a valid directory or '
                        'file with a list of paths. Please provide a '
                        'valid input, either a folder with FASTA files '
                        'or a file with the list of full paths to FASTA '
                        'files (one per line).')

    return list_files


def escape_special_characters(a_string):
    """ Escapes strings to use in regex

        Args:
            a_string (str): string containing characters to escape

        Returns:
            escaped (str): escaped string
    """

    escaped = re.escape(a_string)

    return escaped


def replace_multiple_characters(namefile):
    """ Replaces multiple characters in a string

        Args:
            namefile (str): string containing the name of the contig
            with characters to replace.

        Returns:
            replaced (str): string containing the name of the contig
            without characters to replace.
    """

    replaced = namefile.replace("|", "_")\
                       .replace("_", "-")\
                       .replace("(", "")\
                       .replace(")", "")\
                       .replace("'", "")\
                       .replace("\"", "")\
                       .replace(":", "")

    return replaced


def listdir_fullpath(path):
    """ Gets the full path of the files from a directory

        Args:
            path (str): full path to a directory

        Returns:
            list containing the full path of every file
            contained in the input directory.
    """

    return [os.path.join(path, f) for f in os.listdir(path)]


def flatten_list(list_to_flatten):
    """Flattens one level of a nested list

        Args:
            list_to_flatten (list)

        Returns:
            flattened list

        Example:

            >>> flatten_list([[[1,2],[3,4]]])
            [[1, 2], [3, 4]]

    """

    return list(itertools.chain(*list_to_flatten))


def invert_dictionary(dictionary):
    """ Inverts a dictionary. Keys become values and vice-versa

        Args:
            dictionary (dict)

        Returns:
            inverted (dict): inverted dictionary

        Example:

            >>> inverted_dictionary({key:value})
            {value:key}
    """

    inverted = {value: key for key, value in dictionary.items()}

    return inverted


def threads_for_blast(files_to_blast, cpu_to_apply):
    """ Define the number of threads for BLAST

        Args:
            files_to_blast (list): list containing the full
            path to the files to BLAST
            cpu_to_apply (int): number of cpu to use

        Returns:
            blast_threads (list): list contaning the number
            of threads to use for each file.
            proc (int): Number of processes to use in multiprocessing
    """

    # define number of processes and available cores for each BLAST
    if len(files_to_blast) >= cpu_to_apply:
        blast_threads = [1 for protogenome in files_to_blast]
        proc = cpu_to_apply

    elif cpu_to_apply % len(files_to_blast) == 0:
        blast_threads = [int(cpu_to_apply / len(files_to_blast))
                         for protogenome in files_to_blast]
        proc = len(blast_threads)

    elif cpu_to_apply % len(files_to_blast) == 1:
        blast_threads = [2] + [1 for protogenome in range(0,len(files_to_blast)-1)]
        proc = len(blast_threads)

    elif cpu_to_apply % len(files_to_blast) > 1:
        base_cpu = int(cpu_to_apply / len(files_to_blast))
        blast_threads = [base_cpu
                         for protogenome in range(0, len(files_to_blast))]
        extra_cpu = cpu_to_apply - sum(blast_threads)
        i = 0
        while extra_cpu > 0:
            blast_threads[i] += 1
            extra_cpu -= 1
            i += 1
        proc = len(blast_threads)

    return blast_threads, proc


def isListEmpty(inList):
    """ Checks if a nested list is empty
    """
    if isinstance(inList, list): # Is a list
        return all(map(isListEmpty, inList)) if isinstance(inList, list) else False


def reverseComplement(strDNA):

    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:
        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]


def translateSeq(DNASeq):

    seq = DNASeq
    tableid = 11
    try:
        myseq = Seq(seq)
        protseq = Seq.translate(myseq, table=tableid, cds=True)
    except:
        try:
            seq = reverseComplement(seq)
            myseq = Seq(seq)
            protseq = Seq.translate(myseq, table=tableid, cds=True)
        except:
            try:
                seq = seq[::-1]
                myseq = Seq(seq)
                protseq = Seq.translate(myseq, table=tableid, cds=True)
            except:
                try:
                    seq = seq[::-1]
                    seq = reverseComplement(seq)
                    myseq = Seq(seq)
                    protseq = Seq.translate(myseq, table=tableid, cds=True)
                except Exception as e:
                    # print("translation error")
                    # print(e)
                    raise

    return protseq


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Args:
            dna_str (str): DNA sequence as string type.
            table_id (int): translation table identifier.

        Returns:
            protseq (str): protein sequence created by translating
            the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def is_url(url):
    """ Checks if a url is valid
    
        Args: 
        url (str): the url to be checked

        Returns:
        True if url is valid, False otherwise.
    
    """
    
    try:
        
        result = urlparse(url)
        return all([result.scheme, result.netloc, result.path])
    
    except:
        return False


def make_url(base_url , *res, **params):
    """ Creates a url. 
    
        Args: 
            base_url (str): the base url
            res (str): endpoint(s) to add to the base url
            params (str): addtional parameters (WIP)

        Returns:
            url (str) with the provided parameters.
            Otherwise, returns base_url.

    """
    
    url = base_url
    
    # Check if the url is valid
    if is_url(url):
        
        if url[-1] == "/":
            url = url[:-1]
    
        # Add the endpoints
        for r in res:
    #        url = '{}/{}'.format(url, r)
            url = f'{url}/{r}'
        
        # Add params if they are provided
        if params:
    #        url = '{}?{}'.format(url, urllib.urlencode(params))
            url = f'{url}?{urlencode(params)}'
        
        return url
    
    else:
        return "An invalid URL was provided."


def get_sequence_from_url(url):
    """
    """
    
    seq = parse_qs(urlsplit(url).query)["sequence"][0]
    
    return seq


def login_user_to_NS(server_url, email, password):
    """ Logs a user in Nomenclature Server
    
        Args:
            server_url (str): url of Nomeclature Server API
            email (str): email of the user in NS
            password (str): password of the user in NS
            
        Returns:
            token (str): authorization token to perform requests to NS
    """
    
    auth_params = {}
    auth_params["email"] = email 
    auth_params["password"] = password
    
    auth_headers = {}
    auth_headers["Content-Type"] = "application/json"
    auth_headers["accepts"] = "application/json"
    
    auth_url = make_url(server_url, "auth", "login")
    
    auth_r = requests.post(auth_url, data=json.dumps(auth_params), headers=auth_headers)
    
    auth_result = auth_r.json() 

    token = auth_result["access_token"]
    
    return token
    

def send_data(sparql_query, url_send_local_virtuoso, virtuoso_user, virtuoso_pass):
    """ Sends data to virtuoso.
    
        Args: 
            sparql_query (str): the sparql query to use
            url_send_local_virtuoso (str): the url for the local Virtuoso folder
            virtuoso_user (str): the Virtuoso user
            virtuoso_pass (str): the Virtuoso password

        Returns:
            r (Response) of the Virtuoso server
    """
    
    url = url_send_local_virtuoso
    headers = {'content-type': 'application/sparql-query'}
    r = requests.post(url, data=sparql_query, headers=headers, auth=requests.auth.HTTPBasicAuth(virtuoso_user, virtuoso_pass))

    #sometimes virtuoso returns 405 God knows why ¯\_(ツ)_/¯ retry in 2 sec
    if r.status_code >201:
        time.sleep(2)
        r = requests.post(url, data=sparql_query, headers=headers, auth=requests.auth.HTTPBasicAuth(virtuoso_user, virtuoso_pass))
        
    return r


def send_post(loci_uri, sequence, token, noCDSCheck):
    """ """

    params = {}
    params['sequence'] = sequence

    if not noCDSCheck:
        params['enforceCDS'] = "False"
    
    headers = {'X-API-KEY': token,
               'Content-type': 'application/json',
               'accept': 'application/json'}
    
#    url = loci_uri + "/alleles"
    
    url = make_url(loci_uri, "alleles")

    req_success = False
    sleepfactor = 4
    while not req_success:
        try:
            r = requests.post(url, data=json.dumps(params), headers=headers, timeout=30)
            
            if r.status_code == 418:
                print("Sequence is already attributed to a loci/allele")
            
            elif r.status_code > 201:
                print(r)
                print("failed sending sequence, retrying in seconds "
                      + str(sleepfactor))
                time.sleep(sleepfactor)
                sleepfactor = sleepfactor * 2
            else:
                req_success = True
        except:
            time.sleep(sleepfactor)
            sleepfactor = sleepfactor * 2
            pass

    req_code = int(r.status_code)
    # allele_url=((r.content).decode("utf-8")).replace('"', '').strip()

    return req_code


def send_sequence(token, sequence, loci_uri, noCDSCheck):
    """ """

    req_success = False
    sleepfactor = 4
    while not req_success:

        reqCode = send_post(loci_uri, sequence, token, noCDSCheck)
        if reqCode > 201:
            print("failed, retrying in seconds "+str(sleepfactor))
            time.sleep(sleepfactor)
            sleepfactor = sleepfactor * 2

        else:
            req_success = True
    
    # if reqCode==401:
        # print ("Token is not valid")
    # elif reqCode>201:
        #
        # try:
            #~ allele_url,reqCode=send_post(loci_uri,sequence,token)
        # except:
            #~ print ("Server returned code "+str(reqCode))
            #~ print(loci_uri)
    # else:
        # new_allele_id=str(int(allele_url.split("/")[-1]))
        
    return reqCode


def process_locus(gene, token, loci_url, auxBar, noCDSCheck):
    """ """

    for allele in SeqIO.parse(gene, "fasta", generic_dna):

        sequence = (str(allele.seq)).upper()
        try:
            sequence = translateSeq(sequence)
            sequence = str(sequence)

        except:
            continue

        #reqCode = send_sequence(token, sequence, loci_url, noCDSCheck)
        reqCode = send_post(loci_url, sequence, token, noCDSCheck)
        
#    if reqCode == 418:
#        print(gene)

    if gene in auxBar:
        auxlen = len(auxBar)
        index = auxBar.index(gene)
        print("[" + "=" * index + ">" +
                " " * (auxlen - index) +
                "] Sending alleles " +
                    str(int((index / auxlen) * 100)) + "%")

    return reqCode


def get_data(sparql_query):
    """ Gets data from Virtuoso """
    
    try:
        UNIPROT_SERVER.setQuery(sparql_query)
        UNIPROT_SERVER.setReturnFormat(JSON)
        UNIPROT_SERVER.setTimeout(20)
        result = UNIPROT_SERVER.query().convert()
    except Exception as e:
        time.sleep(5)
        try:
            UNIPROT_SERVER.setQuery(sparql_query)
            UNIPROT_SERVER.setReturnFormat(JSON)
            UNIPROT_SERVER.setTimeout(20)
            result = UNIPROT_SERVER.query().convert()
        except Exception as e:
            result = e
            
    return result
