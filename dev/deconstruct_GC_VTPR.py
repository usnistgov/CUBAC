import CoolProp, pandas, json, os
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
    for fluid in CoolProp.__fluids__ + ['1-Propanol']:
        if fluid in CoolProp.__fluids__:
            CAS = CP.get_fluid_param_string(fluid, "CAS")
        else:
            CAS = df.loc[df['english name'] == fluid, 'CAS-nr.'].iloc[0]
        
        if CAS in df['CAS-nr.'].tolist():
            pass
        elif '.ppf' in CAS.lower():
            print "\t\t\t\tSkipping", fluid
            continue
        else:
            print "Unable to match", fluid, CAS
            continue

        group_string = df.loc[df['CAS-nr.'] == CAS,'increments [count * sub group number]'].iloc[0]
        count_sgi = group_string.split(';')
        count_sgi = [l.split('*') for l in count_sgi]

        groups = []
        for count, sgi in count_sgi:
            groups.append({"count": int(count.strip()), "sgi": int(sgi.strip())})

        entry = {
            "inchikey": "?????????????", 
            "name": fluid, 
            "userid": "", 
            "acentric": CP.PropsSI(fluid,'acentric'), 
            "pc": CP.PropsSI(fluid,'pcrit'), 
            "groups": groups, 
            "registry_number": CAS, 
            "Tc": CP.PropsSI(fluid,'Tcrit')
        }
        entries.append(entry)

    with open(ofname, 'w') as fp:
        json.dump(entries, fp, indent = 2, sort_keys = True)

if __name__=='__main__':
    assert(os.path.exists('GC-VTPR'))
    groups('GC-VTPR/group_data.json')
    interaction_parameters('GC-VTPR/interaction_parameters.json')
    fluids('GC-VTPR/decompositions.json')