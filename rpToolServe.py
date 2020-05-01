import tarfile
import tempfile
import glob
import logging

import rpSBML
import rpTool

'''
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

logging.disable(logging.NOTSET)
my_logger = logging.getLogger('MyLogger')
my_logger.setLevel(logging.DEBUG)

'''

#logging.disable(logging.INFO)
#logging.disable(logging.WARNING)



def runFindPathway_hdd(measured_rpsbml_path, inputTar, strict_length=False, pathway_id='rp_pathway'):
    dict_global = {}
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        tar = tarfile.open(inputTar, 'r')
        tar.extractall(path=tmpOutputFolder)
        tar.close()
        measured_rpsbml = rpSBML.rpSBML('measured')
        measured_rpsbml.readSBML(measured_rpsbml_path)
        for sbml_path in glob.glob(tmpOutputFolder+'/*'):
            fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.rpsbml', '').replace('.xml', '')
            rpsbml = rpSBML.rpSBML(fileName)
            rpsbml.readSBML(sbml_path)
            found, score, dict_result = rpTool.compareRPpathways(measured_rpsbml, rpsbml, strict_length, pathway_id)
            #if score>0.0:
            dict_global[fileName] = {'score': score, 'found': found, 'dict_result': dict_result}
            #dict_global[fileName] = {'reactions_score': reactions_score, 'reactions_std': reactions_std, 'reactions_ec_score': reactions_ec_score, 'reactions_ec_std': reactions_ec_std, 'dict_result': dict_result}
    return dict_global
