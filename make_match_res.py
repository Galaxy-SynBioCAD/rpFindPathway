import rpToolServe_findpath as rpFindPathway
import rpTool as rpReader
import rpToolCache
import glob
import json

selenzyme_results_path= glob.glob('results/*/rpSelenzyme.tar.xz')
all_res = {}

for selen_path in glob.glob('results/*/rpSelenzyme.tar.xz'):
    meas_path = 'validation/'+str(selen_path.split('/')[-2])+'.sbml'
    all_res[selen_path.split('/')[-2]] = rpFindPathway.runFindPathway_hdd(meas_path, selen_path, 'rp_pathway')

with open('match_results.json', 'w') as fp:
    json.dump(all_res, fp)
