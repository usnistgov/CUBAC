
import pandas, json
df = pandas.read_excel('Kang.xlsx','Groups')
groups = []
for ind, row in df.iterrows():
    group = {'sgi': row['sub group index'],
             'subgroup_name': row['sub group name'],
             'mgi': row['main group index'],
             'maingroup_name': row['main group name'],
             'R_k': row['Rk'],
             'Q_k': row['Qk']}
    groups.append(group)
with open('Kang_group_data.json','w') as fp:
    json.dump(groups, fp, indent = 1, sort_keys = True)

df = pandas.read_excel('Kang.xlsx','InteractionParameters')

groups = []
first = True
for ind, row in df.iterrows():
    if first:
        mgi1, mgi2 = row['i'], row['j']
        aij, bij, cijx1000= row['aij'],row['bij'],row['1000cij']
        cij = cijx1000/1000.0
        first = False
    else:
        aji, bji, cjix1000= row['aij'],row['bij'],row['1000cij']
        cji = cjix1000/1000.0

        def prep(n):
            return round(n,10)
        group = {'mgi1': int(mgi1),
                 'mgi2': int(mgi2),
                 'a_ij': prep(aij),
                 'a_ji': prep(aji),
                 'b_ij': prep(bij),
                 'b_ji': prep(bji),
                 'c_ij': prep(cij),
                 'c_ji': prep(cji)
                 }
        groups.append(group)
        first = True
with open('Kang_interaction_parameters.json','w') as fp:
    json.dump(groups, fp, indent = 1, sort_keys = True)
