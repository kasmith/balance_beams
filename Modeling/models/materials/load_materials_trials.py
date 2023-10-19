import os, json

fldir = os.path.dirname(__file__)
trialdir = os.path.join(fldir, "materialsjson")
datflnm = os.path.join(fldir, "BB_MaterialsData.csv")

# Read trial info
def read_materials_stack(stack):
    side, dist, mats = stack
    matdict = {'W':0, 'B':0, 'I':0}
    for m in mats:
        matdict[m] += 1
    if side == 'Rght':
        dist *= -1
    return (dist, matdict)

def read_materials_beam(beam):
    blist = []
    for stack in beam:
        blist.append(read_materials_stack(stack))
    return blist

materials_trials = dict()
for flnm in os.listdir(trialdir):
    if flnm[-5:] == '.json' and flnm[:5] != 'Intro' and flnm[-6:] != 'R.json':  # NOTE: EXCLUDES INTRO TRIALS
        fl = open(os.path.join(trialdir, flnm), 'r')
        if flnm[-6:] == 'N.json':
            trnm = flnm[:-7]
        else:
            trnm = flnm[:-5]
        jsn = json.load(fl)
        materials_trials[trnm] = read_materials_beam(jsn)
        fl.close()

# Read empirical data
datfl = open(datflnm,'r')
next(datfl)
materials_empirical = dict()
materials_emp_order = dict()
for ln in datfl:
    l = ln.strip().split(',')
    wnm = l[0]
    tnm = l[2]
    ord = int(l[9])
    tnmprts = tnm.split('_')
    resp = l[4]
    if wnm not in materials_empirical.keys():
        materials_empirical[wnm] = dict()
    if wnm not in materials_emp_order.keys():
        materials_emp_order[wnm] = dict()
    materials_empirical[wnm][tnm] = resp
    materials_emp_order[wnm][tnm] = ord
datfl.close()

for w in list(materials_empirical.keys()):
    if len(materials_empirical[w]) != 144:
        del materials_empirical[w]

materials_WIDS = materials_empirical.keys()
