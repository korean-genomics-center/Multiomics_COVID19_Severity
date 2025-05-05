# %%
import glob
import gzip
import os
import pickle
from functools import partial
from multiprocessing import Pool
from pathlib import Path


# %%
def main(path_marker_table, dir_methylcpgmin, num_threads, outfilename):
    set_markers = get_set_markers(path_marker_table)
    list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
    list_dict_cpgmin_samples = get_list_dict_cpgmin_samples(list_files_methylcpgmin, set_markers, num_threads)
    list_name_sample = get_list_name_sample(list_files_methylcpgmin)
    dict_cpgmin_all_samples = get_dict_dict_cpgmin_all_samples(list_name_sample, list_dict_cpgmin_samples)
    dump_pickle(dict_cpgmin_all_samples, outfilename)

def get_set_markers(path_marker_table):
    list_markers = list()
    with open(path_marker_table, mode="r") as fr:
        _skiprow = fr.readline()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            chrpos = str(record[0])
            startpos = str(record[1])
            marker = (chrpos, startpos)
            list_markers.append(marker)
    set_markers = set(list_markers)

    return set_markers

def get_list_files_methylcpgmin(dir_methylcpgmin):
    list_files_methylcpgmin = glob.glob(f"{dir_methylcpgmin}/**/*pair_merged.methyl_cpg_min.tsv", recursive=True)

    return list_files_methylcpgmin

def get_list_name_sample(list_files_methylcpgmin):
    list_name_sample = list()
    for file_methylcpgmin in list_files_methylcpgmin:
        dir_methylcpgmin = os.path.basename(os.path.dirname(file_methylcpgmin))
        if dir_methylcpgmin == "HealthyControl":
            name_sample = os.path.basename(file_methylcpgmin).split(".")[0]
        else:
            name_sample = os.path.basename(dir_methylcpgmin)
        list_name_sample.append(name_sample)

    return list_name_sample

def get_dict_cpgmin_sample(file_methylcpgmin, markers=set()):
    dict_marker_freqC = dict()
    with open(file_methylcpgmin, mode="r") as fr:
        _skiprow = fr.readline()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            chrpos = record[0]
            startpos = record[1]
            chr_start = tuple([chrpos, startpos])
            freqC = record[-1]
            if chr_start in markers:
                dict_marker_freqC[chr_start] = freqC
    
    return dict_marker_freqC

def get_list_dict_cpgmin_samples(list_files_methylcpgmin, set_markers, num_threads):
    with Pool(processes=num_threads) as pool:
        func = partial(get_dict_cpgmin_sample, markers=set_markers)
        list_dict_cpgmin_samples = pool.map(func, list_files_methylcpgmin)
        pool.close()
        pool.join()

    return list_dict_cpgmin_samples

def get_dict_dict_cpgmin_all_samples(list_name_sample, list_dict_cpgmin_samples):
    dict_cpgmin_all_samples = dict()
    for sample, dict_cpgmin in zip(list_name_sample, list_dict_cpgmin_samples):
        dict_cpgmin_all_samples[sample] = dict_cpgmin
    
    return dict_cpgmin_all_samples

def dump_pickle(data, outfilename):
    os.makedirs(os.path.dirname(outfilename), exist_ok=True)
    with gzip.open(outfilename, 'wb') as f:
        pickle.dump(data, f)

    return None

def load_pickle(loadfilename):
    with gzip.open(loadfilename,'rb') as f:
        data = pickle.load(f)
    
    return data
# %%
if __name__ == "__main__":
    direction = "hypo"
    visit = "first"
    WORKDIR = str(Path(__file__).parents[3])
    path_marker_table = f"{WORKDIR}/methyl_{direction}_severe_mild_commonMarkers_LOO_markerlist.tsv"
    dir_methylcpgmin = f"{WORKDIR}/Results/10_methyl/MethylCpGMin"
    num_threads = 50
    outfilename = f"{WORKDIR}/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/dictionary_marker_freqC_all_samples_20240220.pk.gz"
    outdir = os.path.dirname(outfilename)
    os.makedirs(outdir, exist_ok=True)
    main(path_marker_table, dir_methylcpgmin, num_threads, outfilename)



# %%
