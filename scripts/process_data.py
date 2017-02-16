import os
from collections import defaultdict
from multiprocessing.pool import ThreadPool

import pandas as pd

from scripts import ncbi
from scripts.allele_map import get_proper_name
from scripts.utils.time_wrap import timer_time

RESOURCES = '../resources/'


def load_data(allele, database):
    st = 'dataset_raw_'
    path = RESOURCES + st + database + '_' + allele + '.txt'
    full_path = os.path.join(os.path.dirname(__file__), path)
    with open(full_path, 'r') as f:
        content = f.readlines()
    return content


@timer_time
def count_frequencies(content, allele):
    for_ret = defaultdict(int)
    allele = get_proper_name(allele)
    # iterable_lines = iter(content.splitlines())
    # for line in iterable_lines:
    for line in content:
        for_ret['total'] += 1
        l = line.lower()
        # print(l)  # DEBUG
        if ('aa\t' + allele) in l:
            for_ret['aa_b'] += 1
        elif ('ac\t' + allele) in l:
            for_ret['ac_b'] += 1
        elif ('ag\t' + allele) in l:
            for_ret['ag_b'] += 1
        elif ('at\t' + allele) in l:
            for_ret['at_b'] += 1
        elif ('ca\t' + allele) in l:
            for_ret['ca_b'] += 1
        elif ('cc\t' + allele) in l:
            for_ret['cc_b'] += 1
        elif ('cg\t' + allele) in l:
            for_ret['cg_b'] += 1
        elif ('ct\t' + allele) in l:
            for_ret['ct_b'] += 1
        elif ('ga\t' + allele) in l:
            for_ret['ga_b'] += 1
        elif ('gc\t' + allele) in l:
            for_ret['gc_b'] += 1
        elif ('gg\t' + allele) in l:
            for_ret['gg_b'] += 1
        elif ('gt\t' + allele) in l:
            for_ret['gt_b'] += 1
        elif ('ta\t' + allele) in l:
            for_ret['ta_b'] += 1
        elif ('tc\t' + allele) in l:
            for_ret['tc_b'] += 1
        elif ('tg\t' + allele) in l:
            for_ret['tg_b'] += 1
        elif ('tt\t' + allele) in l:
            for_ret['tt_b'] += 1

        if (allele + '\taa') in l:
            for_ret['aa_a'] += 1
        elif (allele + '\tac') in l:
            for_ret['ac_a'] += 1
        elif (allele + '\tag') in l:
            for_ret['ag_a'] += 1
        elif (allele + '\tat') in l:
            for_ret['at_a'] += 1
        elif (allele + '\tca') in l:
            for_ret['ca_a'] += 1
        elif (allele + '\tcc') in l:
            for_ret['cc_a'] += 1
        elif (allele + '\tcg') in l:
            for_ret['cg_a'] += 1
        elif (allele + '\tct') in l:
            for_ret['ct_a'] += 1
        elif (allele + '\tga') in l:
            for_ret['ga_a'] += 1
        elif (allele + '\tgc') in l:
            for_ret['gc_a'] += 1
        elif (allele + '\tgg') in l:
            for_ret['gg_a'] += 1
        elif (allele + '\tgt') in l:
            for_ret['gt_a'] += 1
        elif (allele + '\tta') in l:
            for_ret['ta_a'] += 1
        elif (allele + '\ttc') in l:
            for_ret['tc_a'] += 1
        elif (allele + '\ttg') in l:
            for_ret['tg_a'] += 1
        elif (allele + '\ttt') in l:
            for_ret['tt_a'] += 1
    return for_ret


def produce_output_path(file_db, RESOURCES, allele):
    pth = RESOURCES + file_db + '_' + allele + '.csv'
    pth = os.path.join(os.path.dirname(__file__), pth)
    return pth


def main():
    """
    The method prompts user for input of dinucleotide and database to process and writes the output
    to respective file
    :return:
    """
    is_running = True
    while is_running:
        allele, database, is_omim = ncbi.input_params()
        file_db = 'OMIM' if is_omim else database
        data = load_data(allele, file_db)
        results = count_frequencies(data, allele)
        # print(results)  # DEBUG
        res = dict(results)
        df = pd.DataFrame(res, index=[0])
        path = produce_output_path(file_db, RESOURCES, allele)
        df.to_csv(path)
        # print(sorted(res.items(), key=lambda kv: kv[1], reverse=True))  # DEBUG
        inp = input('Do you want to process more data? Y/N')
        is_running = True if inp.lower() == 'y' else False


def process_file(file_name: str):
    if 'raw_SNP' in file_name:
        allele = file_name[-5:-4]
        # print(allele)  # DEBUG
        data = load_data(allele, 'SNP')
        return count_frequencies(data, allele)
    else:
        return None


def count_all_dinucleotides():
    """
    The method automatically collects data on occurrences of all dinucleotides from SNP database
    and writes them to file 'SNP_TOTAL.csv'
    :return:
    """
    pool = ThreadPool(4)
    results = pool.map(process_file, os.listdir(os.path.join(os.path.dirname(__file__), RESOURCES)))
    # print(results)  # DEBUG
    reduce = defaultdict(int)
    for di in results:
        # print(type(di).__name__)  # DEBUG
        if type(di).__name__ == 'defaultdict':
            for key, value in di.items():
                reduce[key] += value
    res = dict(reduce)
    df = pd.DataFrame(res, index=[0])
    path = produce_output_path('SNP', RESOURCES, 'TOTAL')
    df.to_csv(path)
    print('all dinucleotides counted')

main()
# count_all_dinucleotides()
