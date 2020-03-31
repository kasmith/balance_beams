import os, json

fldir = os.path.dirname(__file__)
trialdir = os.path.join(fldir, "CombineJSON")
datflnm = os.path.join(fldir, "BB_CombData.csv")

combined_trials = dict()
for flnm in os.listdir(trialdir):
    if flnm[-5:] == '.json' and flnm[:5] != 'Intro' and flnm[-6:] != 'R.json':  # NOTE: EXCLUDES INTRO TRIALS
        fl = open(os.path.join(trialdir, flnm), 'rU')
        if flnm[-6:] == 'N.json':
            trnm = flnm[:-7]
        else:
            trnm = flnm[:-5]
        jsn = json.load(fl)
        if jsn['StrutSize'] == 2.0:
            print ('uh oh')
        combined_trials[trnm] = jsn
        fl.close()


# Read empirical data
datfl = open(datflnm,'rU')
next(datfl)
combined_empirical = dict()
for ln in datfl:
    l = ln.strip().split(',')
    wnm = l[0].strip('"')
    tnm = l[2].strip('"') + '_' + l[9].strip('"')
    resp = l[4].strip('"')
    if wnm not in combined_empirical.keys():
        combined_empirical[wnm] = dict()
    combined_empirical[wnm][tnm] = resp
datfl.close()

combined_WIDS = list(combined_empirical.keys())

for w in combined_WIDS:
    if len(combined_empirical[w].keys()) != 192:
        del combined_empirical[w]
