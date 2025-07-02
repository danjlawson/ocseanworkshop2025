from lingpy import *
import os
ocseanwl = Wordlist('OCSEAN_processed_joineddata.tsv')
ocseanlex = LexStat(ocseanwl,check=True)
ocseanlex.get_scorer()
ocseanlex.cluster(method='lexstat', threshold=0.55, cluster_method='infomap', ref='cogid')
ocseanlex.output('tsv', filename=os.path.join('data','OCSEAN_processed_cognatedetected'))
