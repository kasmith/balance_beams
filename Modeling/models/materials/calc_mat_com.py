from load_materials_trials import materials_trials

mdict = {'W':1., 'I': 9., 'B':2.4}

ofl = open('AltMatTrs.csv','w')
ofl.write('TrialBase,Dist_L_Same,Dist_R_Same,Mass_L_Same,Mass_R_Same,Torque_L_Same,Torque_R_Same,Type_Same,Falls_Same,')
ofl.write('Dist_L_Heavy,Dist_R_Heavy,Mass_L_Heavy,Mass_R_Heavy,Torque_L_Heavy,Torque_R_Heavy,Type_Heavy,Falls_Heavy\n')

for nm in materials_trials.keys():
    tr = materials_trials[nm]

    # Assume everything is made of the same materials
    dist_l = 0
    dist_r = 0
    mass_l = 0
    mass_r = 0
    trq_l = 0
    trq_r = 0
    for d, mats in tr:
        bls = mats['I'] + mats['B'] + mats['W']
        if d > 0:
            if d > dist_l:
                dist_l = d
            mass_l += bls
            trq_l += bls*d
        else:
            d = -d
            if d > dist_r:
                dist_r = d
            mass_r += bls
            trq_r += bls*d

    if trq_l == trq_r:
        falls = "B"
    elif trq_l > trq_r:
        falls = "L"
    else:
        falls = "R"

    if falls == 'B':
        if mass_l == mass_r:
            tp = 'Bal'
        else:
            tp = 'CB'
    elif mass_l == mass_r:
        tp = 'Dist'
    elif dist_l == dist_r:
        tp = 'Weight'
    elif (dist_l > dist_r and mass_l < mass_r):
        if falls == 'L':
            tp = 'CD'
        else:
            tp = 'CW'
    elif (dist_l < dist_r and mass_l > mass_r):
        if falls == 'L':
            tp = 'CW'
        else:
            tp = 'CD'
    else:
        if (dist_l > dist_r and falls == 'L') or (dist_l < dist_r and falls == 'R'):
            tp = "Obvious"
        else:
            tp = "Non-obvious"

    ofl.write(nm + ',' + str(dist_l) + ',' + str(dist_r) + ',' + str(mass_l) + ',' + str(mass_r) + ',')
    ofl.write(str(trq_l) + ',' + str(trq_r) + ',' + tp + ',' + falls + ',')

    # Assume only the heaviest material exists
    dist_l = 0
    dist_r = 0
    mass_l = 0
    mass_r = 0
    trq_l = 0
    trq_r = 0

    has_i = False
    has_b = False
    for _, mats in tr:
        if mats['I'] > 0:
            has_i = True
        if mats['B'] > 0:
            has_b = True

    for d, mats in tr:
        if has_i:
            bls = mats['I']
        elif has_b:
            bls = mats['B']
        else:
            bls = mats['W']

        if bls > 0:
            if d > 0:
                if d > dist_l:
                    dist_l = d
                mass_l += bls
                trq_l += bls * d
            else:
                d = -d
                if d > dist_r:
                    dist_r = d
                mass_r += bls
                trq_r += bls * d

    if trq_l == trq_r:
        falls = "B"
    elif trq_l > trq_r:
        falls = "L"
    else:
        falls = "R"

    if falls == 'B':
        if mass_l == mass_r:
            tp = 'Bal'
        else:
            tp = 'CB'
    elif mass_l == mass_r:
        tp = 'Dist'
    elif dist_l == dist_r:
        tp = 'Weight'
    elif (dist_l > dist_r and mass_l < mass_r):
        if falls == 'L':
            tp = 'CD'
        else:
            tp = 'CW'
    elif (dist_l < dist_r and mass_l > mass_r):
        if falls == 'L':
            tp = 'CW'
        else:
            tp = 'CD'
    else:
        tp = "Obvious"

    ofl.write(str(dist_l) + ',' + str(dist_r) + ',' + str(mass_l) + ',' + str(mass_r) + ',')
    ofl.write(str(trq_l) + ',' + str(trq_r) + ',' + tp + ',' + falls + '\n')

ofl.close()