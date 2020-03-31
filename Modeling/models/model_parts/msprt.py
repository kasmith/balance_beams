


def load_dirichlet_samples(flnm):
    dsfl = open(flnm, 'rU')
    dsfl.next()
    dirichsamps = dict()
    for ln in dsfl:
        l = ln.strip('\n').split(',')
        N = int(l[0])
        if N not in dirichsamps: dirichsamps[N] = dict()
        dirichsamps[N][(int(l[1]), int(l[2]), int(l[3]))] = (float(l[4]), float(l[5]), float(l[6]))
    dsfl.close()
    return dirichsamps


def pStopAndGo(prev, probs, N, T, dsamps):
    ddict = dsamps[N]
    pstop = [0.,0.,0.]
    remaining = dict()
    for pcts, pprb in prev.items():
        # Add one sample to each prior
        for i in range(3):
            ncts = [pcts[j] for j in range(3)]
            ncts[i] += 1
            ncts = tuple(ncts)
            nprb = pprb * probs[i]
            belief = ddict[ncts]
            if belief[0] > T: pstop[0] += nprb
            elif belief[1] > T: pstop[1] += nprb
            elif belief[2] > T: pstop[2] += nprb
            else:
                if ncts not in remaining:
                    remaining[ncts] = 0
                remaining[ncts] += nprb
    return pstop, remaining


def MSPRT(probs, T, dsamps, stopthresh = .99, normalize = False):
    assert T > .5 and T < 1, 'Inappropriate SPRT threshold ' + str(T)
    prem = 1
    pstop = [0.,0.,0.]
    expnstop = 0.
    n = 1
    maxn = max(dsamps.keys())
    rem = {(0,0,0): 1.}
    while n <= maxn and sum(pstop) < stopthresh:
        thisstp, rem = pStopAndGo(rem, probs, n, T, dsamps)
        expnstop += sum(thisstp) * n
        prem -= sum(thisstp)
        pstop = [prev+now for prev,now in zip(pstop,thisstp)]
        n += 1
    if normalize:
        totp = sum(pstop)
        expnstop += n*prem
        pstop = [p/totp for p in pstop]
    return (pstop, expnstop)