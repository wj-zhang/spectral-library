import simulate_mixture_spectra
import run_pipeline_two_rounds
import run_pipeline_SpectraST
import run_pipeline_MSPLIT
import os
import subprocess
import run_pipeline_reSpect
import pandas as pd


def get_filenames(library_name, out_folder):
    sim_our_library_file = out_folder + '/' + library_name + '_simulation_library.ms2'
    sim_st_library_file = out_folder + '/' + library_name + '_simulation_library_st.ms2'
    sim_query_file = out_folder + '/' + library_name + '_simulation_query.mgf'
    id_file = out_folder + '/' + library_name + '_simulation_id.csv'
    infor_file = out_folder + '/' + library_name + '_simulation_infor.txt'
    return sim_our_library_file, sim_st_library_file, sim_query_file, id_file, infor_file


def generate_simulation_data(library_folder, library_name, out_folder, option, random_seed=0):
    library_file = library_folder + '/' + library_name + '_target.ms2'
    sim_our_library_file, sim_st_library_file, sim_query_file, id_file, infor_file = get_filenames(library_name, out_folder)
    infor = simulate_mixture_spectra.generate(library_file, sim_our_library_file, sim_st_library_file,
                                              sim_query_file, id_file, infor_file, option,
                                              random_state=random_seed, is_shuffle=False, remove_peaks_ratio=0.5)


def run_search_Ours(library_name, library_out_folder, out_folder):
    sim_our_library_file, _, sim_query_file, _, _ = get_filenames(library_name, library_out_folder)
    raw_id_file = out_folder + '/' + library_name + '_simulation_raw_id.csv'
    id_file = out_folder + '/' + library_name + '_simulation_id.csv'
    decoy_file = out_folder + '/' + library_name + '_all.ms2'

    print('generate decoy library ...')
    run_pipeline_two_rounds.run_generate_decoy_library(sim_our_library_file, decoy_file, ratio=1, combined=True)

    print('seperate target and decoy ...')
    library_decoy_file = out_folder + '/' + library_name + '_decoy.ms2'
    library_target_file = out_folder + '/' + library_name + '_target.ms2'
    run_pipeline_two_rounds.run_separate_target_decoy(decoy_file, library_target_file, library_decoy_file)

    raw_target_id_file = out_folder + '/' + library_name + '_target_simulation_raw_id.csv'
    raw_decoy_id_file = out_folder + '/' + library_name + '_decoy_simulation_raw_id.csv'
    print('library search on target spectra...')
    run_pipeline_two_rounds.run_main_library_search(sim_query_file, sim_our_library_file, raw_target_id_file, mode='SIMULATION_SEARCH')#

    print('library search on decoy spectra...')
    run_pipeline_two_rounds.run_main_library_search(sim_query_file, library_decoy_file, raw_decoy_id_file, mode='SIMULATION_SEARCH')#

    print('merge target and decoy ids ...')
    run_pipeline_two_rounds.run_merge_target_decoy_ids(raw_target_id_file, raw_decoy_id_file, raw_id_file)

    print('cut off threshold of weights ...')
    run_pipeline_two_rounds.run_cut_off_ids(raw_id_file, id_file, threshold=0.15)


def run_search_SpectraST(library_name, library_out_folder, out_folder):
    _, sim_st_library_file, sim_query_file, _, _ = get_filenames(library_name, library_out_folder)
    library_file = out_folder + '/' + library_name + '_simulation'
    raw_id_file = out_folder + '/' + library_name + '_simulation_raw_id.xls'
    id_file = out_folder + '/' + library_name + '_simulation_id.csv'

    print('build library from ms2 file ...')
    run_pipeline_SpectraST.run_build_library(sim_st_library_file, library_file)

    print('generate decoy library ...')
    decoy_file = out_folder + '/' + library_name + '_simulation_decoy'
    run_pipeline_SpectraST.run_generate_decoy_library(library_file, decoy_file, ratio=1)

    print('library search ...')
    run_pipeline_SpectraST.run_library_search(library_file, sim_query_file, raw_id_file, mass_tol=3)
    print('formate id file to csv file ...')
    run_pipeline_SpectraST.run_fdr(raw_id_file, id_file, fdr=0.01)


def run_search_MSPLIT(library_name, library_out_folder, out_folder):
    _, _, sim_query_file, _, _ = get_filenames(library_name, library_out_folder)
    library_file = out_folder.replace('MSPLIT', 'SpectraST') + '/' + library_name + '_simulation_decoy.sptxt'
    svm_tmp_dir = out_folder + '/' + library_name + '_simulation_svm_tmp'
    raw_id_file = out_folder + '/' + library_name + '_simulation_raw_id.txt'
    fdr_id_file = out_folder + '/' + library_name + '_simulation_raw_id_fdr.txt'
    id_file = out_folder + '/' + library_name + '_simulation_id.csv'

    print('build library ...')
    run_pipeline_MSPLIT.run_build_library(library_file)

    print('library search ...')
    run_pipeline_MSPLIT.run_library_search(library_file, sim_query_file, raw_id_file, mass_tol=3)
    print('determine the number of ids ...')
    run_pipeline_MSPLIT.run_fdr(raw_id_file, fdr_id_file, svm_tmp_dir, fdr=0.01)
    print('format id file to csv file ...')
    run_pipeline_MSPLIT.run_format_to_csv_file(fdr_id_file, id_file)


def run_search_reSpect(library_name, library_out_folder, out_folder):
    _, sim_st_library_file, sim_query_file, _, _ = get_filenames(library_name, library_out_folder)
    sim_query_file = sim_query_file.replace('.mgf', '_respect.mgf')
    library_file = out_folder + '/' + library_name + '_simulation_decoy'

    if not os.path.isfile(library_file + '.splib'):
        print('build library from ms2 file ...')
        run_pipeline_SpectraST.run_build_library(sim_st_library_file, library_file)

    ms_mzml_file = out_folder + '/' + library_name + '_simulation_query_respect' + '.mzml'
    if not os.path.isfile(ms_mzml_file):
        print('convert raw file to mzml file')
        subprocess.check_call(["msconvert ", sim_query_file, '-o', out_folder + '/'])

    pepxml_id_file = out_folder + '/' + library_name + '_raw_id_1.pep.xml'
    print('library search of first round...')
    run_pipeline_SpectraST.run_library_search(library_file, ms_mzml_file, pepxml_id_file, mass_tol=1.1)

    isow = 3.1
    accumulated_tag = ''
    for round in range(2, 6):
        print('run PeptideProphet and iProphet...')
        prophet_pepxml_file = pepxml_id_file[:-8] + '_interact.pep.xml'
        run_pipeline_reSpect.run_peptide_i_prophet(pepxml_id_file, prophet_pepxml_file)
        print('run respect...')
        iprophet_pepxml_file = pepxml_id_file[:-8] + '_interact.ipro.pep.xml'
        run_pipeline_reSpect.run_respect(iprophet_pepxml_file)

        accumulated_tag += '_rs'

        ms_mzml_file = out_folder + '/' + library_name + '_simulation_query_respect' + '' + '.mzml'
        dot_idx = ms_mzml_file.find('.', 4)
        ms_mzml_file = ms_mzml_file[:dot_idx] + accumulated_tag + ms_mzml_file[dot_idx:]
        tmp_ms_mzml_file = '_' + library_name + '_simulation_query_respect' + '' + accumulated_tag + '.mzML'
        subprocess.check_call(["msconvert", ms_mzml_file, "-o", out_folder + "/", "--outfile", tmp_ms_mzml_file])
        os.remove(ms_mzml_file)
        os.rename(out_folder + '/' + tmp_ms_mzml_file, ms_mzml_file)

        pepxml_id_file = out_folder + '/' + library_name + '_raw_id_' + str(round) + '.pep.xml'
        run_pipeline_SpectraST.run_library_search(library_file, ms_mzml_file, pepxml_id_file, mass_tol=isow)

    print('merge id files...')
    ids = []
    for round in range(1, 6):
        pepxml_id_file = out_folder + '/' + library_name + '_raw_id_' + str(round) + '.pep.xml'
        id_file = out_folder + '/' + library_name + '_raw_id_' + str(round) + '.csv'
        run_pipeline_reSpect.run_pepxml_fdr(pepxml_id_file, id_file, fdr=0.01, simulation=True)
        ids.append(pd.read_csv(id_file))
    combined_csv = pd.concat(ids)

    combined_csv = combined_csv.groupby('Index')
    id_file = out_folder + '/' + library_name + '_simulation_id.csv'
    ids_list = []
    for index, id_df in combined_csv:
        peptide = []
        score = []
        protein = []
        charge = []
        id_df = id_df.copy()
        id_df.drop_duplicates(subset=['Peptide', 'Charge'], keep='first', inplace=True)
        for idx, row in id_df.iterrows():
            peptide.append(row['Peptide'])
            score.append(row['Score'])
            protein.append(row['Protein'])
            charge.append(row['Charge'])

        peptide = ';'.join(peptide)
        score = ';'.join([str(s) for s in score])
        protein = ';'.join(protein)
        charge = ';'.join([str(c) for c in charge])

        ids_list.append([index, peptide, charge, False, score, protein])
    ids = pd.DataFrame(ids_list, columns=['Index', 'Peptide', 'Charge', 'Decoy', 'Score', 'Protein'])
    ids.to_csv(id_file, header=True, index=False)


if __name__ == '__main__':

    real_library_folder = '../tmp/CharmeRT/real/Ours'
    pepxml_name = 'HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1'

    simulation_library_folder = '../tmp/CharmeRT/simulation/data_with_noise'

    out_folder = {'Ours': '../tmp/CharmeRT/simulation/Ours',
                  'SpectraST': '../tmp/CharmeRT/simulation/SpectraST',
                  'MSPLIT': '../tmp/CharmeRT/simulation/MSPLIT',
                  'reSpect': '../tmp/CharmeRT/simulation/reSpect',
                  }

    option = {'ratio': {1: 0.3, 2: 0.2, 3: 0.2, 4: 0.1, 5: 0.1},
              'weights': {2: [[1, 1], [1, 0.2]],
                          3: [[1, 1, 1], [1, 0.6, 0.2]],
                          4: [[1, 1, 1, 1], [1, 0.6, 0.4, 0.2]],
                          5: [[1, 1, 1, 1, 1], [1, 0.8, 0.6, 0.4, 0.2]]},
              'mz range': 2.5
              }

    methods = [
        'Ours',
        'SpectraST',
        'MSPLIT',
        'reSpect',
    ]

    print('\nGenerate simulation data...')
    generate_simulation_data(real_library_folder, pepxml_name, simulation_library_folder, option)

    for method in methods:
        print('\nLibrary search from simulation by ' + method)
        out_folder_str = 'out_folder[\'' + method + '\'])'
        exec('run_search_' + method + '(pepxml_name, simulation_library_folder, ' + out_folder_str)
