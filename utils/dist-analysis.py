from os import listdir, rename
from os.path import isfile, join
import re

def merge(inpath, outpath):
    with open(inpath, 'r') as infile:
        with open(outpath,'a') as outfile:
            for line in infile:
                outfile.write(line)

def organize(dirpath):
    '''
    within dirpath, export data in r_ files to a_ files and rename r_ files to z_ files
    '''
    files = [f for f in listdir(dirpath) if isfile(join(dirpath, f))]
    pattern = re.compile('^r[0-9]+_(.*-.txt$)')
    for file in files:
        filematch = pattern.match(file)
        if filematch is not None:
            afile = 'a_' + filematch[1]
            zfile = 'z' + file[1:]
            merge(file,afile)
            rename(file,zfile)


# above are subroutine to organize the output dist files by the parallel program which is name after r[no. of rank]_ to avoid discrepancy
# the organized file is started with a_ while files already archived (data copied to a_) is started with z_

# below is the part to extract distribution of observables across disorder configuration and output the result to a short report txt file


picksize = re.compile('^a_.*-L([0-9]*)-V.*')
pickint = re.compile('^a_.*-V([0-9.]*)-.txt$')
pick1 = re.compile('^a_([0-9.]*)-.*')
pick2 = re.compile('^a_[0-9.]*-([0-9.]*)-.*')


def extract_data(path):
    with open(path,'r') as datafile:
        data = datafile.read()
    data = data.split()
    data = [float(d) for d in data]
    return data

def count_in_bins(data, nobin):
    space = 1/nobin
    count = [0 for i in range(nobin)]
    for s in data:
        try:
            count[int(s/space)] += 1
        except IndexError: # avoid the case where entanglement entropy density be exactly 1.
            count[-1] += 1
    return count

def summary(dirpath, regfilter, measure, nobin, result_file='0summary.txt'):
    '''
    function to generate short report from massive txt data

    :param dirpath: string of the work directory, eg. './'
    :param regfilter: re.compile object, as a selection condition of files we will analyse, eg. pick1
    :param measure: int, the number of measured observables we want to analyse, eg. 1 for max length ratio, 2 for entanglement entropy
    :param nobin: int, how many bins we divide to chracterize the data distribution
    :param result_file: string, file name of the output report
    '''
    results = []
    files = [f for f in listdir(dirpath) if isfile(join(dirpath, f)) and regfilter.match(f) is not None]
    for file in files:
        data = extract_data(file)
        count = count_in_bins(data[measure::3], nobin)
        result = [pick1.match(file)[1],pick2.match(file)[1],picksize.match(file)[1],pickint.match(file)[1]]
        result.extend(count)
        results.append(result)
    with open(result_file,'w') as output:
        for result in results:
            for item in result:
                output.write(str(item))
                output.write(' ')
            output.write('\n') 

# config this dict to make the script suitable for the analysis you want to conduct
config_dict = {'dirpath': './', 'regfilter': pick1, 'measure': 1, 'nobin': 10, 'result_file': '0summary.txt'}

if __name__ == '__main__':
    organize(config_dict['dirpath'])
    summary(**config_dict)
