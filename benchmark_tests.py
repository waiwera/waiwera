# run all benchmark tests and produce a summary report
from __future__ import print_function

tests_path = "test/benchmark"

import argparse
import datetime
import os
import sys
import glob
import subprocess
import xml.dom.minidom
import json
import platform
from docutils.core import publish_file

start_time = str(datetime.datetime.now())

parser = argparse.ArgumentParser()
parser.add_argument("-np", help = "number of processes")
args = parser.parse_args()
if args.np: num_procs = args.np
else: num_procs = '1'

orig_path = os.getcwd()
summary = {"tests": []}
passed = True
testcount, passcount = 0, 0

for path, dirs, files in os.walk(tests_path):
    os.chdir(path)
    scripts = glob.glob("test_*.py")
    for script in scripts:
        print(os.path.join(path, script), '...', end = "")
        sys.stdout.flush()
        subprocess.call(['python', script, '-np', num_procs], stdout = open(os.devnull, 'wb'))
        if os.path.isdir('output'):
            os.chdir('output')
            outpaths = [p.strip('/') for p in glob.glob('*/')]
            for outpath in outpaths:
                outfilename = 'SysTest-' + outpath + '.xml'
                xmldoc = xml.dom.minidom.parse(os.path.join(outpath, outfilename))
                node = xmldoc.childNodes[0]
                testname = node.attributes['name'].value
                teststatus = node.attributes['status'].value
                testreport = os.path.join(path, 'output', outpath, testname + '-report.html')
                testreport = os.path.relpath(testreport, tests_path)
                summary_item = {'name': testname, 'status': teststatus, 'report': testreport}
                summary['tests'].append(summary_item)
                testpassed = teststatus == 'Pass'
                passed = passed and testpassed
                if testpassed: passcount += 1
                testcount += 1
                print('', teststatus)
    os.chdir(orig_path)

end_time = str(datetime.datetime.now())

if passed:
    print('Result: Pass')
    ret = 0
else:
    print('Result: %d/%d tests passed (%d%%)' %
          (passcount, testcount, int(100. * passcount / float(testcount))))
    ret = 1

summary_rst_filename = os.path.join(tests_path, 'test_summary.rst')
rstfile = open(summary_rst_filename, 'w')

rstfile.write("Benchmark test summary report\n")
rstfile.write("*****************************\n\n")

rstfile.write("Overall Result: %s\n" % ("Pass" if passed else "Fail"))
rstfile.write("====================\n\n")
summary['pass'] = passed

rstfile.write("Description\n")
rstfile.write("===========\n\n")
rstfile.write("Waiwera benchmark tests\n\n")

rstfile.write("Specification\n")
rstfile.write("=============\n\n")
rstfile.write(" * nproc: %s\n\n" % num_procs)

rstfile.write("Provenance\n")
rstfile.write("==========\n\n")
hostname = platform.uname()[1]
rstfile.write(" * node: %s\n" % hostname)
summary['host'] = hostname
platform_name = platform.system()
release_name = platform.release()
summary['platform'] = {'name': platform_name, 'release': release_name}
rstfile.write(" * platform: %s (%s)\n" % (platform_name, release_name))
rstfile.write(" * start time: %s\n" % start_time)
rstfile.write(" * end time: %s\n" % end_time)
rstfile.write("\n")
summary['start_time'] = start_time
summary['end_time'] = end_time

rstfile.write("Tests\n")
rstfile.write("=====\n\n")

for test in summary['tests']:
    rstfile.write(" * `%s <%s>`_: %s\n" % (test['name'], test['report'], test['status']))
rstfile.write("\n")

rstfile.close()

summary_html_filename = os.path.join(tests_path, 'test_summary.html')
html = publish_file(source_path = summary_rst_filename,
                    destination_path = summary_html_filename,
                    writer_name = "html")

summary_json_filename = os.path.join(tests_path, 'test_summary.json')
json.dump(summary, file(summary_json_filename, 'w'))

sys.exit(ret)
