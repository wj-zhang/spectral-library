from pyteomics import mzml


ms_file = '../data/10222015_ABRF3_1000ng_68X.mzML'
scan_id = 5992

ms_file = '../data/10222015_ABRF1_1000ng_68X.mzML'
scan_id = 14974

with mzml.read(ms_file) as ms_reader:
    for ms in ms_reader:
        if ms['ms level'] == 2:
            scan = 'scan=' + str(scan_id)
            if scan in ms['id']:
                masses = [str(x) for x in ms['m/z array']]
                masses = ', '.join(masses)
                print('[' + masses + ']')
                intensities = [str(x) for x in ms['intensity array']]
                intensities = ', '.join(intensities)
                print('[' + intensities + ']')

                break

print('ok')
