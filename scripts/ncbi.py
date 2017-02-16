# coding=utf-8
import errno
import os
import signal  # would NOT work properly on WINDOWS. UNIX only
import time
import urllib
from functools import wraps

import requests
from Bio import Entrez

#  TODO: Use OMIM API
#  TODO: Use multiprocessing

from scripts.utils import params


class TimeoutError(Exception):
    pass


def timeout(error_message=os.strerror(errno.ETIME), seconds=params.timeout):
    """
    Function creates a decorator, which applied to other function would break it if that function
    takes too long to execute.

    Note that this will only work on UNIX.
    The basic idea is to use signal handlers to set an alarm for some time interval and raise
    an exception once that timer expires.
    The process for timing out an operations is described in the documentation for signal module.
    :param seconds: specifies "too long" in seconds
    :param error_message: custom error message
    :return:
    """

    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator


def get_search_query_options(database, term):
    '''
    Returns result of eSearch query on NCBI server: Number of results, QueryKey and WebEnv.

    See details in docs: http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    :param database: database name on NCBI server
    :param term: desired search query on NCBI, one can obtain it from 'Search details' window
     on NCBI cite to the right (ex.: http://www.ncbi.nlm.nih.gov/snp )
    :return:
    '''
    # print(term)  # DEBUG
    handle = Entrez.esearch(db=database, retmode="xml", term=term, usehistory='y')
    records = Entrez.read(handle)
    # print(records)  # DEBUG
    count = int(records['Count'])
    query_key = records['QueryKey']
    web_env = records['WebEnv']
    return count, query_key, web_env


def harvest_url(datab, query_key, web_env, retmax=5, retStart=0):
    '''
    Writes results obtained from NCBI (use QueryKey and WebEnv from eSearch) to file specified.

    See details in docs: http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    :param datab: database name on NCBI server
    :param retmax: Total number of records from the input set to be retrieved, up to
    a maximum of 10,000. Optionally, for a large set the value of retstart can be
    iterated while holding retmax constant, thereby downloading the entire set in
    batches of size retmax.
    :param query_key: QueryKey from eSearch
    :param web_env: WebEnv from eSearch
    :param retStart: Sequential index of the first record to be retrieved (default=0,
    corresponding to the first record of the entire set). This parameter can be used
    in conjunction with retmax to download an arbitrary subset of records from the input set.
    :return:
    '''
    handle_fetch = Entrez.efetch(db=datab, rettype='txt', retmode='text', retmax=retmax, query_key=query_key,
                                 WebEnv=web_env, retstart=retStart)
    return handle_fetch.url


def specify_parameters():
    '''
    Promts user for parameters. Default values on 'enter'.
    :return:
    '''
    print("Please specify some parameters. Hit 'Enter' for default values.")
    retMax_def = params.ret_max

    allele, database, isOmim = input_params()

    term_def = '("Homo sapiens"[Organism] OR homo sapiens[All Fields]) AND (snp[Snp_Class] AND ' \
               '(00000.0500[GLOBAL_MAF] : 00000.5000[GLOBAL_MAF]) AND by 1000g[Validation] NOT merged rs[Filter] ' \
               'AND {}[ALLELE])'.format(allele)

    term_def_OMIM = '"Homo sapiens"[Organism] OR homo sapiens[All Fields]) AND (snp[Snp_Class] AND snp_omim[Filter] AND(00000.0500[GLOBAL_MAF] : 00000.5000[GLOBAL_MAF]) AND by 1000g[Validation] NOT merged rs[Filter] AND {}[ALLELE]'.format(allele)

    term = input("Specify desired NCBI eSearch query [default]: ".format(term_def))
    if isOmim:
        # print(term)
        term = term_def_OMIM
    if not term: term = term_def
    print('Term: {}'.format(term))

    while True:
        inp = input(
            "Do you want to continue downloading dataset from {}? (Y/n): ".format(database))
        reset_flag = False if inp.lower() == 'y' else True
        if (inp.lower() == 'y') or inp.lower() == 'n':
            break
    print(reset_flag)
    retMax = retMax_def
    print('retMax: {}'.format(retMax))
    while True:
        email = input("Enter your email (NCBI will use it to spy and spam you down): ".format(
            retMax_def))  # E-mail is necessary for NCBI API queries
        if not email:
            email = 'andriy.lazorenko@gmail.com'
        if email:
            if len(email.split('@')) != 2: continue  # just sanity check
            break
    return database, term, retMax, email, reset_flag, allele, isOmim


def input_params():
    allele_def = "Y"
    is_omim = False
    database_def = 'SNP'
    database = input("Specify NCBI Database [def: '{}']: ".format(database_def))
    if not database: database = database_def
    if database == "OMIM":
        is_omim = True
        database = database_def
    else:
        database = database_def
    print('Database: {}'.format(database))
    allele = input("Select allele [def: '{}']: ".format(allele_def))
    if not allele: allele = allele_def
    print('Allele: {}'.format(allele.upper()))
    return allele.upper(), database, is_omim


def cook_url(line, param, value):
    """
    Replaces value of PARAM in LINE (URL query) with provided VALUE
    :param line:
    :param param:
    :param value:
    :return:
    """
    i = line.index(param)
    begin = i + len(param) + 1
    sli = line[begin:]
    end = begin + len(sli.split('&')[0])
    return line[:begin] + str(value) + line[end:]


def chopper(dataset, database, count, query_key, web_env, retMax, reset_flag=False):
    with open(os.path.join(os.path.dirname(__file__), dataset)) as fh_read:
        num_lines = sum(1 for line in fh_read if line.rstrip())
        print(num_lines)
        with open(os.path.join(os.path.dirname(__file__), dataset), 'a+') as fh:
            i = num_lines // retMax
            residue = num_lines % retMax
            print(i)
            urly = harvest_url(database, query_key=query_key, web_env=web_env, retmax=retMax)
            while True:
                try:
                    t0 = time.time()
                    chop(retMax, i, urly, fh, residue=residue)
                    i += 1
                    t1 = time.time() - t0
                    print('batch #{} ({} of {} lines): {} sec.'.format(i + 1, (i + 1) * retMax, count, t1))
                    if ((i - 1) * retMax) > count:
                        print("Database may be downloaded. Check '{}' file. Have fun.".format(dataset))
                        break
                except TimeoutError as err:
                    print(err, 'Trying again: ret_max = {}'.format(retMax))
                    pass
                except urllib.error.HTTPError as err2:
                    print(err2)
                    print("Database may be downloaded. Check '{}' file. Have fun.".format(dataset))
                    break
        print('{} entries in {}.'.format(len(list(fh_read.readlines())), dataset))


@timeout(os.strerror(errno.ETIMEDOUT))
def chop(retMax, i, url, fh, residue=0):
    retStart = i * retMax + residue
    url = cook_url(url, 'retstart', retStart)
    print(url)
    req = requests.post(url)
    fh.write(req.text)


def main():
    """
    Runs the code and downloads selected data from NCBI, saving it into txt file.
    :return:
    """
    # base_url_esearch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    # base_url_efetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

    database, term, retMax, Entrez.email, reset_flag, allele, is_omim = specify_parameters()

    file_db = 'OMIM' if is_omim else database
    dataset = '../resources/dataset_raw_{}_{}.txt'.format(file_db, allele)  # destination filename

    print("Results will be downloaded in '{}' file. Be patient.".format(dataset))
    dataset_path = os.path.join(os.path.dirname(__file__), dataset)
    try:
        if not os.path.isfile(dataset_path):
            print("File is created...")
            file = open(dataset_path, 'w')
            file.close()
    except FileNotFoundError as err:
        pass

    count, query_key, web_env = get_search_query_options(database, term)
    print('QueryKey: {}, WebEnv: {}'.format(query_key, web_env))
    chopper(dataset, database, count, query_key, web_env, retMax, reset_flag=reset_flag)


if __name__ == "__main__":
    main()
