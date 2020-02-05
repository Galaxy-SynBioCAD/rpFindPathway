import tarfile
import tempfile
import glob

import rpSBML
import rpTool

#######################################################################
##################### Detect ##########################################
#######################################################################


def runFindPathway_hdd(measured_rpsbml_path, inputTar, pathway_id='rp_pathway'):
    dict_global = {}
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        tar = tarfile.open(inputTar, 'r:xz')
        tar.extractall(path=tmpOutputFolder)
        tar.close()
        measured_rpsbml = rpSBML.rpSBML('measured')
        measured_rpsbml.readSBML(measured_rpsbml_path)
        for sbml_path in glob.glob(tmpOutputFolder+'/*'):
            fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.rpsbml', '').replace('.xml', '')
            rpsbml = rpSBML.rpSBML(fileName)
            rpsbml.readSBML(sbml_path)
            found, reactions_score, reactions_std, reactions_ec_score, reactions_ec_std, dict_result = rpTool.compareRPpathways(measured_rpsbml, sim_rpsbml, scrict_length=True, pathway_id='rp_pathway')
            if found:
                dict_global[fileName] = {'reactions_score': reactions_score, 'reactions_std': reactions_std, 'reactions_ec_score': reactions_ec_score, 'reactions_ec_std': reactions_ec_std, 'dict_result': dict_result}
    return dict_global
