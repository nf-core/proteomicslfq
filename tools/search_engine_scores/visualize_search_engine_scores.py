from pyopenms import *
from typing import List
import pandas as pd
import plotly.express as px
import plotly

prots = []  # type: List[ProteinIdentification]
peps = []  # type: List[PeptideIdentification]
IdXMLFile().load("/Users/pfeuffer/Downloads/debugPLFQ/BSAConsensus/consensus_eval/UPS1_50000amol_R1_comet_idx_feat_perc.xml", prots, peps)

cols = ["comet", "comet_xcorr", "comet_pep", "msgf", "msgf_raw", "msgf_pep", "comet_seq", "msgf_seq", 'target_decoy_comet', 'target_decoy_msgf']

d = {}
for pep in peps:
    best_hit = pep.getHits()[0]  # type: PeptideHit
    d[pep.getMetaValue("spectrum_reference")] = {
        'comet': best_hit.getMetaValue("COMET:lnExpect"),
        'comet_xcorr': best_hit.getMetaValue("MS:1002252"),
        'comet_pep': best_hit.getScore(),
        'msgf': None,
        'msgf_raw': None,
        'msgf_pep': None,
        'comet_seq':  best_hit.getSequence().toString(),
        'msgf_seq': None,
        'target_decoy_comet': best_hit.getMetaValue("target_decoy"),
        'target_decoy_msgf': None,
    }


prots = []  # type: List[ProteinIdentification]
peps = []  # type: List[PeptideIdentification]
IdXMLFile().load("/Users/pfeuffer/Downloads/debugPLFQ/BSAConsensus/consensus_eval/UPS1_50000amol_R1_msgf_idx_feat_perc.xml", prots, peps)

for pep in peps:
    best_hit = pep.getHits()[0]  # type: PeptideHit
    spec_ref = pep.getMetaValue("spectrum_reference")
    if spec_ref in d:
        mydict = d[spec_ref]
        mydict["msgf"] = np.log10(best_hit.getMetaValue("MS:1002052"))
        mydict["msgf_seq"] = best_hit.getSequence().toString()
        mydict["msgf_raw"] = best_hit.getMetaValue("MS:1002049")
        mydict["msgf_pep"] = best_hit.getScore()
        mydict["target_decoy_msgf"] = best_hit.getMetaValue("target_decoy")

    else:
        d[spec_ref] = {
            'comet': None,
            'comet_pep': None,
            'msgf': best_hit.getMetaValue("MS:1002052"),  # MSGF SpecEval
            'msgf_raw': best_hit.getMetaValue("MS:1002049"),
            'msgf_pep': best_hit.getScore(),
            'comet_seq':  None,
            'msgf_seq': best_hit.getSequence().toString(),
            'target_decoy_comet': None,
            'target_decoy_msgf': best_hit.getMetaValue("target_decoy")
        }


df = pd.DataFrame.from_dict(d, orient="index", columns=cols)
df['same_seq'] = df.apply(lambda x: x['comet_seq'] == x['msgf_seq'], axis=1)


def getTDstatus(x):
    if x['same_seq']:
        return x['target_decoy_comet']
    else:
        if x['target_decoy_comet'] is not None and x['target_decoy_msgf'] is not None:
            if x['target_decoy_comet'] == x['target_decoy_msgf']:
                return x['target_decoy_comet']
            else:
                return 'mixed'
        elif x['target_decoy_comet'] is not None:
            return x['target_decoy_comet']
        else:
            return x['target_decoy_msgf']


df['target_decoy'] = df.apply(getTDstatus, axis=1)
df['pep_max'] = df.apply(lambda x: max(x["comet_pep"],x["msgf_pep"]), axis=1)

print(df['comet_seq'].head(10))
print(df['msgf_seq'].head(10))

df.to_csv("/Users/pfeuffer/Downloads/debugPLFQ/BSAConsensus/consensus_eval/pandas.csv")
#tidy_df = df.melt(id_vars=[])
#print(tidy_df.head())
fig = px.scatter(df, x='comet_xcorr', y='msgf_raw', color='same_seq', symbol='target_decoy',
                 hover_data=["target_decoy_comet", "comet_seq", "msgf_seq"],
                 marginal_x="violin",
                 marginal_y="violin", )
fig.show()
plotly.offline.plot(fig, "/Users/pfeuffer/Downloads/debugPLFQ/BSAConsensus/consensus_eval/plotly.html")

fig = px.violin(df, x='target_decoy', y='pep_max', color='target_decoy')
fig.show()


