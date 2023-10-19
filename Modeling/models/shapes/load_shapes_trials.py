import os, json

fldir = os.path.dirname(__file__)
trialdir = os.path.join(fldir, "shapesjson")
datflnm = os.path.join(fldir, "BB_ShapesData.csv")

# Read trial info
def read_shapes_stack(stack):
    if len(stack) == 1:
        stack = stack[0]
        if stack[0] == 'C':
            return ('C' + stack[2], stack[1])
        else:
            return (stack[0], stack[1])
    else:
        if (stack[0][0] != 'B'): raise Exception('Something went wrong')
        return ('B', len(stack))

def read_shapes_side(side):
    sidelist = []
    for i in range(len(side)):
        stack = side[i]
        if len(stack) != 0:
            rs = read_shapes_stack(stack)
            sidelist.append((rs[0], rs[1], i + 1))
    return sidelist

shapes_trials = dict()
for flnm in os.listdir(trialdir):
    if flnm[-6:] == 'N.json' and flnm[:5] != 'Intro':  # NOTE: EXCLUDES INTRO TRIALS
        fl = open(os.path.join(trialdir, flnm), 'r')
        trnm = flnm[:-5]
        jsn = json.load(fl)
        shapes_trials[trnm] = read_shapes_side(jsn['Left']) + [(b, m, -d) for b,m,d in read_shapes_side(jsn['Right'])]
        fl.close()

# Read empirical data
datfl = open(datflnm,'r')
next(datfl)
shapes_empirical = dict()
shapes_emp_order = dict()
for ln in datfl:
    l = ln.strip().split(',')
    wnm = l[0]
    tnm = l[1]
    ord = int(l[6])
    tnmprts = tnm.split('_')
    if tnmprts[3] == 'N':
        resp = l[2]
    else:
        tnm = tnmprts[0]+'_'+tnmprts[1]+'_'+tnmprts[2]+'_N'
        if l[2] == 'L':
            resp = 'R'
        elif l[2] == 'R':
            resp = 'L'
        else:
            resp = 'B'
    if wnm not in shapes_empirical.keys():
        shapes_empirical[wnm] = dict()
    if wnm not in shapes_emp_order.keys():
        shapes_emp_order[wnm] = dict()
    shapes_empirical[wnm][tnm] = resp
    shapes_emp_order[wnm][tnm] = ord
datfl.close()

for w in list(shapes_empirical.keys()):
    if len(shapes_empirical[w]) != 154:
        del shapes_empirical[w]

shapes_WIDS = shapes_empirical.keys()
