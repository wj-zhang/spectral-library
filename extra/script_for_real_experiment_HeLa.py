import run_pipeline_two_rounds
import run_pipeline_SpectraST
import run_pipeline_MSPLIT
import time
import os.path
import subprocess
import run_pipeline_reSpect
import pandas as pd


def run_build_library_Ours(library_folder, pepxml_name, out_folder, ratios):
    library_corrected_suffix = '_corrected'
    pepxml_file = library_folder + '/' + pepxml_name + '.pep.xml'

    library_file = out_folder + '/' + pepxml_name + '_library.ms2'
    library_file_remove_ids = out_folder + '/' + pepxml_name + '_library_remove_ids.ms2'
    library_target_file = out_folder + '/' + pepxml_name + '_target.ms2'

    print('correct mass shifts of library ...')
    corrected_pepxml_file = run_pipeline_two_rounds.run_correct_mass_shifts_from_pepxml(pepxml_file,
                                                                                        out_folder,
                                                                                        library_corrected_suffix,
                                                                                        remove_duplicates=False)
    print('build library from pepxml ...')
    run_pipeline_two_rounds.run_build_library_from_pepxml(corrected_pepxml_file, library_file, remove_duplicates=False)
    print('remove conflicting identifications ...')
    run_pipeline_two_rounds.run_remove_conflicting_identifications(library_file, library_file_remove_ids, precursor_ppm=None)

    for ratio in ratios:
        print('generate decoy library ...')
        decoy_file = out_folder + '/' + pepxml_name + '_' + str(ratio) + '_all.ms2'
        run_pipeline_two_rounds.run_generate_decoy_library(library_file_remove_ids, decoy_file, ratio=ratio, combined=True)
        print('seperate target and decoy ...')
        library_decoy_file = out_folder + '/' + pepxml_name + '_' + str(ratio) + '_decoy.ms2'
        run_pipeline_two_rounds.run_separate_target_decoy(decoy_file, library_target_file, library_decoy_file)


def run_search_Ours(data_folder, data_name, library_folder, library_name, out_folder, ratio):
    single_ms_file = data_folder + '/' + data_name + '_non_chimera.mgf'
    clone_ms_file = data_folder + '/' + data_name + '_chimera.mgf'

    single_corrected_ms_file = out_folder + '/' + data_name + '_' + str(ratio) + '_non_chimera_corrected.mgf'
    combined_ms_file = out_folder + '/' + data_name + '_' + str(ratio) + '_non_chimera_corrected.mgf'

    decoy_file = out_folder + '/' + library_name + '_' + str(ratio) + '_all.ms2'
    library_target_file = out_folder + '/' + library_name + '_target.ms2'
    library_decoy_file = out_folder + '/' + library_name + '_' + str(ratio) + '_decoy.ms2'

    id_target_file = out_folder + '/' + data_name + '_' + str(ratio) + '_raw_id_target.csv'
    id_decoy_file = out_folder + '/' + data_name + '_' + str(ratio) + '_raw_id_decoy.csv'
    raw_id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_raw_id.csv'
    id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_id.csv'
    peaks_id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_PEAKS_id.csv'

    if os.path.isfile(id_file):
        return

    print('pre-library search ...')
    spectra = run_pipeline_two_rounds.run_pre_library_search(single_ms_file, decoy_file)
    print('correct mass shifts of query ...')
    run_pipeline_two_rounds.run_correct_mass_shifts((single_ms_file, spectra), single_corrected_ms_file,
                                                    mode='query')
    print('combine clone spectra ...')
    run_pipeline_two_rounds.run_combine_clone_from_two_mgfs(single_corrected_ms_file, clone_ms_file,
                                                            combined_ms_file)

    start = time.time()
    print('main library search (target) ...')
    run_pipeline_two_rounds.run_main_library_search(combined_ms_file, library_target_file, id_target_file)
    print('main library search (decoy) ...')
    run_pipeline_two_rounds.run_main_library_search(combined_ms_file, library_decoy_file, id_decoy_file)
    print('merge target and decoy ids ...')
    run_pipeline_two_rounds.run_merge_target_decoy_ids(id_target_file, id_decoy_file, raw_id_file)

    print('cut off ids ...')
    run_pipeline_two_rounds.run_cut_off_ids(raw_id_file, id_file, threshold=0.15)

    print('format to PEAKS ids ...')
    run_pipeline_two_rounds.run_format_to_PEAKS_id(id_file, peaks_id_file, clone_ms_file)

    end = time.time()
    print('running time is ' + str(end-start))


def run_build_library_SpectraST(library_folder, pepxml_name, out_folder, ratios):
    PEAKS_pepxml_file = library_folder + '/' + pepxml_name + '.pep.xml'
    pepxml_file = library_folder + '/' + pepxml_name + '.st.pep.xml'

    raw_library_file = out_folder + '/' + pepxml_name
    consensus_library_file = out_folder + '/' + pepxml_name + '_cons'
    quality_library_file = out_folder + '/' + pepxml_name + '_cons_qua'


    print('format pepxml file for SpectraST ...')
    run_pipeline_SpectraST.run_format_pepxml(PEAKS_pepxml_file, pepxml_file)
    print('build raw library from pepxml file ...')
    run_pipeline_SpectraST.run_build_library(pepxml_file, raw_library_file)

    print('build consensus library ...')
    run_pipeline_SpectraST.run_build_consensus_library(raw_library_file, consensus_library_file)
    print('perform quality control of library ...')
    run_pipeline_SpectraST.run_quality_control_library(consensus_library_file, quality_library_file)
    for ratio in ratios:
        print('generate decoy library ...')
        decoy_file = out_folder + '/' + pepxml_name + '_' + str(ratio) + '_cons_qua_dec'
        run_pipeline_SpectraST.run_generate_decoy_library(quality_library_file, decoy_file, ratio=ratio)


def run_search_SpectraST(data_folder, data_name, library_folder, library_name, out_folder, ratio):
    ms_file = data_folder + '/' + data_name + '_non_chimera.mgf'
    decoy_file = out_folder + '/' + library_name + '_' + str(ratio) + '_cons_qua_dec'

    raw_id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_raw_id.xls'
    id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_id.csv'

    # if os.path.isfile(id_file):
    #     return

    print('library search ...')
    idx = data_name.find('isow')
    if idx >= 0:
        isow = data_name[idx + 4]
    else:
        isow = 3
    print('isow: ' + str(isow))
    run_pipeline_SpectraST.run_library_search(decoy_file, ms_file, raw_id_file, mass_tol=isow)
    print('fdr ...')
    run_pipeline_SpectraST.run_fdr(raw_id_file, id_file, fdr=0.01)


def run_search_reSpect(data_folder, data_name, library_folder, library_name, out_folder, ratio):

    ms_file = data_folder + '/' + data_name + '.raw'

    decoy_file = out_folder.replace('reSpect', 'SpectraST') + '/' + library_name + '_' + str(ratio) + '_cons_qua_dec'

    ms_mzml_file = out_folder + '/' + data_name + '.mzml'
    if not os.path.isfile(ms_mzml_file):
        print('convert raw file to mzml file')
        subprocess.check_call(["msconvert ", ms_file, '--filter', "peakPicking true 1-",  '--filter', "zeroSamples removeExtra", '-o', out_folder + '/'])

    pepxml_id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_raw_id_1.pep.xml'
    print('library search of first round...')
    run_pipeline_SpectraST.run_library_search(decoy_file, ms_mzml_file, pepxml_id_file, mass_tol=1.1)

    idx = data_name.find('isow')
    if idx >= 0:
        isow = str(int(data_name[idx + 4]) + 0.1)
    else:
        isow = 3.1
    accumulated_tag = ''
    for round in range(2, 5):
        print('run PeptideProphet and iProphet...')
        prophet_pepxml_file = pepxml_id_file[:-8] + '_interact.pep.xml'
        run_pipeline_reSpect.run_peptide_i_prophet(pepxml_id_file, prophet_pepxml_file)
        print('run respect...')
        iprophet_pepxml_file = pepxml_id_file[:-8] + '_interact.ipro.pep.xml'
        run_pipeline_reSpect.run_respect(iprophet_pepxml_file)

        accumulated_tag += '_rs'

        ms_mzml_file = out_folder + '/' + data_name + '' + accumulated_tag + '.mzml'
        tmp_ms_mzml_file = '_' + data_name + '' + accumulated_tag + '.mzML'
        subprocess.check_call(["msconvert", ms_mzml_file, "-o", out_folder + "/", "--outfile", tmp_ms_mzml_file])
        os.remove(ms_mzml_file)
        os.rename(out_folder + '/' + tmp_ms_mzml_file, ms_mzml_file)

        pepxml_id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_raw_id_' + str(round) + '.pep.xml'
        run_pipeline_SpectraST.run_library_search(decoy_file, ms_mzml_file, pepxml_id_file, mass_tol=isow)

    print('merge id files...')
    ids = []
    for round in range(1, 5):
        pepxml_id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_raw_id_' + str(round) + '.pep.xml'
        id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_raw_id_' + str(round) + '.csv'
        run_pipeline_reSpect.run_pepxml_fdr(pepxml_id_file, id_file)
        ids.append(pd.read_csv(id_file))
    combined_csv = pd.concat(ids)

    combined_csv = combined_csv.groupby('Index')
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

    id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_id.csv'
    ids = pd.DataFrame(ids_list, columns=['Index', 'Peptide', 'Charge', 'Decoy', 'Score', 'Protein'])
    ids.to_csv(id_file, header=True, index=False)


def run_build_library_MSPLIT(library_folder, pepxml_name, out_folder, ratios):
    for ratio in ratios:
        decoy_file = out_folder.replace('MSPLIT', 'SpectraST') + '/' + pepxml_name + '_' + str(ratio) + '_cons_qua_dec.sptxt'
        run_pipeline_MSPLIT.run_build_library(decoy_file)


def run_search_MSPLIT(data_folder, data_name, library_folder, library_name, out_folder, ratio):
    ms_file = data_folder + '/' + data_name + '_non_chimera.mgf'
    decoy_file = out_folder.replace('MSPLIT', 'SpectraST') + '/' + library_name + '_' + str(ratio) + '_cons_qua_dec.sptxt'

    raw_id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_raw_id.txt'
    fdr_id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_id.txt'
    id_file = out_folder + '/' + data_name + '_' + str(ratio) + '_id.csv'

    svm_tmp_dir = out_folder + '/' + data_name + '_' + str(ratio) + '_svm_tmp'

    print('library search ...')
    idx = data_name.find('isow')
    if idx >= 0:
        isow = data_name[idx + 4]
    else:
        isow = 3
    print('isow: ' + str(isow))
    run_pipeline_MSPLIT.run_library_search(decoy_file, ms_file, raw_id_file, mass_tol=isow)
    print('fdr ...')
    run_pipeline_MSPLIT.run_fdr(raw_id_file, fdr_id_file, svm_tmp_dir, 0.01)
    print('format id file to csv file ...')
    run_pipeline_MSPLIT.run_format_to_csv_file(fdr_id_file, id_file)


if __name__ == '__main__':
    library_folder = '../data/CharmeRT/library'
    data_folder = '../data/CharmeRT/search'

    pepxml_name = 'HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1'

    out_folder = {'Ours': '../tmp/CharmeRT/real/Ours',
                  'SpectraST': '../tmp/CharmeRT/real/SpectraST',
                  'MSPLIT': '../tmp/CharmeRT/real/MSPLIT',
                  'reSpect': '../tmp/CharmeRT/real/reSpect',
                  }

    data_names = [
        'HeLa_1ug_isow8_1hgradient_datasetE_rep1',
        'HeLa_1ug_isow8_1hgradient_datasetE_rep2',
        'HeLa_1ug_isow8_1hgradient_datasetE_rep3',

        'HeLa_1ug_isow8_3hgradient_datasetF_rep1',
        'HeLa_1ug_isow8_3hgradient_datasetF_rep2',
        'HeLa_1ug_isow8_3hgradient_datasetF_rep3',

        'HeLa_1ug_isow4_1hgradient_datasetC_rep1',
        'HeLa_1ug_isow4_1hgradient_datasetC_rep2',
        'HeLa_1ug_isow4_1hgradient_datasetC_rep3',

        'HeLa_1ug_isow4_3hgradient_datasetD_rep1',
        'HeLa_1ug_isow4_3hgradient_datasetD_rep2',
        'HeLa_1ug_isow4_3hgradient_datasetD_rep3',

        'HeLa_1ug_isow2_1hgradient_datasetA_rep1',
        'HeLa_1ug_isow2_1hgradient_datasetA_rep2',
        'HeLa_1ug_isow2_1hgradient_datasetA_rep3',
                  ]

    methods = [
        'Ours',
        'SpectraST',
        'MSPLIT',
        'reSpect',
    ]

    ratios = [
        1,
    ]

    for method in methods:
        print('\nBuild library for ' + method)
        out_folder_str = 'out_folder[\'' + method + '\'], ratios)'
        exec('run_build_library_' + method + '(library_folder, pepxml_name, ' + out_folder_str)

    for ratio in ratios:
        for data_name in data_names:
            print('\n\nRun data: ' + data_name)
            for method in methods:
                print('\nLibrary search by ' + method)
                out_folder_str = 'out_folder[\'' + method + '\'],' + str(ratio) + ')'
                exec('run_search_' + method + '(data_folder, data_name, library_folder, pepxml_name, ' + out_folder_str)
