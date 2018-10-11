import CoolProp, pandas, json, os, numpy as np
CP = CoolProp.CoolProp

def groups(ofname):
    df = pandas.read_excel('GC-VTPR.xlsx', sheetname='Groups')
    entries = []
    for i,row in df.iterrows():
        entry = {
          "Q_k": row['Qk'], 
          "R_k": row['Rk'], 
          "maingroup_name": row["main group name"], 
          "mgi": row['main group index'], 
          "sgi": row['sub group index'], 
          "subgroup_name": row["sub group name"]
        }
        entries.append(entry)

    with open(ofname, 'w') as fp:
        json.dump(entries, fp, indent = 2, sort_keys = True)

def interaction_parameters(ofname):
    df = pandas.read_excel('GC-VTPR.xlsx', sheetname='InteractionParameters')
    df = df.fillna(0.0)
    entries = []
    for i,row in df.iterrows():
        entry = {
            "a_ij": row['aij / K'],
            "a_ji": row['aji / K'],
            "b_ij": row['bij'],
            "b_ji": row['bji'],
            "c_ij": row['cij / K-1'], 
            "c_ji": row['cji / K-1'], 
            "mgi1": row['i'], 
            "mgi2": row['j']
        }
        entries.append(entry)

    with open(ofname, 'w') as fp:
        json.dump(entries, fp, indent = 2, sort_keys = True)

def fluids(ofname):
    df = pandas.read_excel('GC-VTPR.xlsx', sheetname='Components')

    df = df.fillna(0.0)
    entries = []
    for _index, row in df.iterrows():        

        group_string = row['increments [count * sub group number]']
        count_sgi = group_string.split(';')
        count_sgi = [l.split('*') for l in count_sgi]

        alpha = None
        if not np.isnan(row['L']):
            alpha = {'type': 'Twu', 'c': [row['L'],row['M'],row['N']], "BibTeX": "Guennec-FPE-2016-tcPR"}

        groups = []
        for count, sgi in count_sgi:
            groups.append({"count": int(count.strip()), "sgi": int(sgi.strip())})

        entry = {
            "inchikey": "?????????????", 
            "name": row['english name'],
            "userid": "", 
            "molemass": 0.02, # XXXXXXXXXXXXX
            "acentric": row['acentric'], 
            "pc": row['p_crit (Pa)'], 
            "groups": groups, 
            "registry_number": row['CAS-nr.'], 
            "Tc": row['T_crit (K)'],
            "BibTeX": "Guennec-FPE-2016-tcPR"
        }
        if alpha is not None:
            entry['alpha'] = alpha
        entries.append(entry)

    with open(ofname, 'w') as fp:
        json.dump(entries, fp, indent = 2, sort_keys = True)

if __name__=='__main__':
    assert(os.path.exists('GC-VTPR'))
    groups('GC-VTPR/group_data.json')
    interaction_parameters('GC-VTPR/interaction_parameters.json')
    fluids('GC-VTPR/decompositions.json')