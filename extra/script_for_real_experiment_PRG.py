import script_for_real_experiment_HeLa


if __name__ == '__main__':
    library_folder = '../data/PRG2017/library'
    data_folder = '../data/PRG2017/search'

    pepxml_name = '10222015_ABRF_1000ng_68X_chimera_0.1'

    out_folder = {'Ours': '../tmp/PRG2017/real/Ours',
                  'SpectraST': '../tmp/PRG2017/real/SpectraST',
                  'MSPLIT': '../tmp/PRG2017/real/MSPLIT',
                  'reSpect': '../tmp/PRG2017/real/reSpect',
                  }

    data_names = [
        '10222015_ABRF3_1000ng_68X',
        '10222015_ABRF4_1000ng_68X',
                  ]

    methods = [
        'Ours',
        'SpectraST',
        'MSPLIT',
        'reSpect',
    ]

    ratios = [1]

    for method in methods:
        print('\nBuild library for ' + method)
        out_folder_str = 'out_folder[\'' + method + '\'], ratios)'
        exec('script_for_real_experiment.run_build_library_' + method + '(library_folder, pepxml_name, ' + out_folder_str)

    for ratio in ratios:
        for data_name in data_names:
            print('\n\nRun data: ' + data_name)
            for method in methods:
                print('\nLibrary search by ' + method)
                out_folder_str = 'out_folder[\'' + method + '\'],' + str(ratio) + ')'
                exec('script_for_real_experiment.run_search_' + method + '(data_folder, data_name, library_folder, pepxml_name, ' + out_folder_str)
