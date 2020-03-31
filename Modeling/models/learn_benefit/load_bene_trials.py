import os
import json
import numpy as np

fldir = os.path.dirname(__file__)
trialdir = os.path.join(fldir, "learnbenejson")
datflnm = os.path.join(fldir, "BB_BeneData.csv")

# Read trial info
def read_materials_stack(stack):
    dist, shape, size, mat = stack
    matdict = {'W':0, 'B':0, 'I':0}
    matdict[mat] += size
    return (-dist, matdict)

def read_materials_beam(beam):
    blist = []
    for stack in beam['Beam']:
        blist.append(read_materials_stack(stack))
    return blist

learnbene_trials = dict()
learnbene_training_trials = dict()
for flnm in os.listdir(trialdir):
    if flnm[-5:] == '.json' and flnm[:4] == 'Test' and flnm[-6:] != 'R.json':  # NOTE: EXCLUDES INTRO / LEARNIN TRIALS
        fl = open(os.path.join(trialdir, flnm), 'rU')
        if flnm[-6:] == 'N.json':
            trnm = flnm[:-7]
        else:
            trnm = flnm[:-5]
        jsn = json.load(fl)
        learnbene_trials[trnm] = read_materials_beam(jsn)
        fl.close()
    if flnm[-5:] == '.json' and flnm[:5] == 'Learn' and flnm[-6:] != 'R.json':  # NOTE: EXCLUDES INTRO / LEARNIN TRIALS
        fl = open(os.path.join(trialdir, flnm), 'rU')
        if flnm[-6:] == 'N.json':
            trnm = flnm[:-7]
        else:
            trnm = flnm[:-5]
        jsn = json.load(fl)
        learnbene_training_trials[trnm] = read_materials_beam(jsn)
        fl.close()

datfl = open(datflnm,'rU')
next(datfl)
learnbene_empirical = dict()
learnbene_rts = dict()
learnbene_emp_order = dict()
for ln in datfl:
    l = ln.strip().split(',')
    wnm = l[0]
    tnm = l[2]
    tnmprts = tnm.split('_')
    tnm_base = ('_').join(tnmprts[:-1])
    resp = l[7]
    time = float(l[13])
    order = int(l[12])
    if wnm not in learnbene_empirical.keys():
        learnbene_empirical[wnm] = dict()
    learnbene_empirical[wnm][tnm_base] = resp
    if wnm not in learnbene_emp_order.keys():
        learnbene_emp_order[wnm] = dict()
    learnbene_emp_order[wnm][tnm_base] = order
    if wnm not in learnbene_rts.keys():
        learnbene_rts[wnm] = []
    if tnmprts[0] == 'Testing':
        learnbene_rts[wnm].append(time)
datfl.close()

ks = list(learnbene_empirical.keys())
for w in ks:
    if len(learnbene_empirical[w]) != 180 or np.median(learnbene_rts[w]) < 1000:
        del learnbene_empirical[w]

learnbene_WIDS = list(learnbene_empirical.keys())
