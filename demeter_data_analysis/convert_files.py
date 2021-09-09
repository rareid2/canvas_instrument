import os
from os import listdir
from os.path import isfile, join

mypath = '/home/rileyannereid/workspace/canvas_demeter/data/all_data'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

os.chdir(mypath)
for f in onlyfiles:
    if '.ps' in f:
        cmd = 'ps2pdf ' + f
        os.system(cmd)

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
os.chdir(mypath)

for f in onlyfiles:
    if '.pdf' in f:
        cmd = 'pdftoppm ' + f + ' ' + f[:-3] + ' -png'
        os.system(cmd)

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
converted_files = '/home/rileyannereid/workspace/canvas_demeter/data/burst_pngs/'

for f in onlyfiles:
    if '.png' in f:
        cmd = 'mv ' + f + ' ' + converted_files + f
        os.system(cmd)

onlyfiles = [f for f in listdir(converted_files) if isfile(join(converted_files, f))]
os.chdir(converted_files)
for f in onlyfiles:
    if '-1.png' in f:
        cmd = 'mv ' + f + ' ' + f[:-6] + 'png'
        os.system(cmd)