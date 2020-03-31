import os, json

fldir = os.path.dirname(__file__)
trialdir = os.path.join(fldir, "BalanceJSON")
datflnm = os.path.join(fldir, "BB_BalanceData.csv")

balance_trials = dict()
for flnm in os.listdir(trialdir):
    if flnm[-5:] == '.json' and flnm[:5] != 'Intro' and flnm[-6:] != 'R.json':  # NOTE: EXCLUDES INTRO TRIALS
        fl = open(os.path.join(trialdir, flnm), 'rU')
        if flnm[-6:] == 'N.json':
            trnm = flnm[:-7]
        else:
            trnm = flnm[:-5]
        jsn = json.load(fl)
        balance_trials[trnm] = jsn
        fl.close()

# Read empirical data
datfl = open(datflnm,'rU')
next(datfl)
balance_empirical = dict()
balance_emp_order = dict()
for ln in datfl:
    l = ln.strip().split(',')
    wnm = l[0].strip('"')
    tnm = l[2].strip('"')
    ord = int(l[13].strip('"'))
    tnmprts = tnm.split('_')
    resp = l[10].strip('"')
    if wnm not in balance_empirical.keys():
        balance_empirical[wnm] = dict()
    if wnm not in balance_emp_order.keys():
        balance_emp_order[wnm] = dict()
    balance_empirical[wnm][tnm] = resp
    balance_emp_order[wnm][tnm] = ord
datfl.close()


balance_WIDS = list(balance_empirical.keys())
