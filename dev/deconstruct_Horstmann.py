
import pandas, json
df = pandas.read_excel('../Horstmann.xlsx','Groups')
groups = []
for ind, row in df.iterrows():
    group = {'sgi': row['sub group index'],
             'subgroup_name': row['sub group name'],
             'mgi': row['main group index'],
             'maingroup_name': row['main group name'],
             'R_k': row['Rk'],
             'Q_k': row['Qk']}
    groups.append(group)
with open('Horstmann_group_data.json','w') as fp:
    json.dump(groups, fp, indent = 1, sort_keys = True)


df = pandas.read_excel('../Horstmann.xlsx','InteractionParameters')
df = df.fillna(0.0)
def clean(s):
    if isinstance(s, (float, int)):
        return s
    else:
        if s.endswith('.'):
            return s.rstrip('.')

groups = []
for ind, row in df.iterrows():
    group = {'mgi1': row['i'],
             'mgi2': row['j'],
             'a_ij': float(row['aij / K']),
             'a_ji': float(clean(row['aji / K'])),
             'b_ij': row['bij'],
             'b_ji': row['bji'],
             'c_ij': row['cij / K-1'],
             'c_ji': row['cji / K-1']
             }
    groups.append(group)
with open('Horstmann_interaction_parameters.json','w') as fp:
    json.dump(groups, fp, indent = 1, sort_keys = True)
