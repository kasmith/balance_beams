import os, json

fldir = os.path.dirname(__file__)
trialdir = os.path.join(fldir, "ferrettijson")
datflnm = os.path.join(fldir, "BB_FerrettiData.csv")

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

ferretti_trials = dict()
for flnm in os.listdir(trialdir):
    if flnm[-5:] == '.json' and flnm[:5] != 'Intro' and flnm[-6:] != 'R.json':  # NOTE: EXCLUDES INTRO TRIALS
        fl = open(os.path.join(trialdir, flnm), 'rU')
        if flnm[-6:] == 'N.json':
            trnm = flnm[:-7]
        else:
            trnm = flnm[:-5]
        jsn = json.load(fl)
        ferretti_trials[trnm] = read_materials_beam(jsn)
        fl.close()

# Read empirical data
datfl = open(datflnm,'rU')
next(datfl)
ferretti_empirical = dict()
for ln in datfl:
    l = ln.strip().split(',')
    wnm = l[0]
    tnm = l[2]
    tnmprts = tnm.split('_')
    resp = l[4]
    if wnm not in ferretti_empirical.keys():
        ferretti_empirical[wnm] = dict()
    ferretti_empirical[wnm][tnm] = resp
datfl.close()

ks = list(ferretti_empirical.keys())
for w in ks:
    if len(ferretti_empirical[w]) != 144:
        del ferretti_empirical[w]

ferretti_WIDS = list(ferretti_empirical.keys())
