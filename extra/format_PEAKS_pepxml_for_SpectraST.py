import re


def run(in_file, out_file):
    with open(in_file, 'r') as reader:
        content = reader.read()
        content = re.sub(r'spectrum="(.*)\.raw" start_scan="(.*)" end_scan="(.*)" '
                         r'precursor_neutral_mass="(.*)" assumed_charge="(.*)" index="(.*)"',
                         r'spectrum="\1.\2.\3.\5" start_scan="\2" end_scan="\3" '
                         r'precursor_neutral_mass="\4" assumed_charge="\5" index="\6"',
                         content)
        content = re.sub(r'</search_hit>',
                         r'  <analysis_result analysis="peptideprophet"> <peptideprophet_result probability="1"/> '
                         r'</analysis_result>\n        </search_hit>',
                         content)

    with open(out_file, 'w') as writer:
        writer.write(content)
