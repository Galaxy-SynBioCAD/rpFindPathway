#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpCofactors REST service

"""
import argparse
import tarfile
import tempfile
import logging
import os
import sys
import json
import csv
import argparse

sys.path.insert(0, '/home/')
import rpFindPathwayServe

logging.basicConfig(
    #level=logging.DEBUG,
    #level=logging.WARNING,
    level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


'''
Example species input: 
    {'input_type': {'input_format': 'tar'}, 'search': {'db_name': 'chebi', 'id': '38407', 'inchi': 'InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+', 'search_type': 'species'}, 'output_type': {'output_format': 'csv'}, 'adv': {'pathway_id': 'rp_pathway'}}
Example reaction input:
    {'input_type': {'input_format': 'tar'}, 'search': {'ec': [{'id': '1.1.1.1'}], 'reactants': [{'db_name': 'chebi', 'id': '17333', 'inchi': 'InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+'}], 'products': [{'db_name': 'chebi', 'id': '38407', 'inchi': 'InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+'}], 'search_type': 'reaction'}, 'output_type': {'output_format': 'csv'}, 'adv': {'pathway_id': 'rp_pathway'}}
Example pathway input:
    {'input_type': {'input_format': 'tar'}, 'search': {'reactions': [{'ec': [{'id': '1.1.1.1'}], 'reactants': [{'name': 'one', 'db_name': 'mnx', 'id': 'MNXM3', 'inchi': 'InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+'}], 'products': [{'name': 'two', 'db_name': 'mnx', 'id': 'MNXM4', 'inchi': 'InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+'}]}, {'ec': [{'id': '2.2.2.2'}], 'reactants': [{'name': 'three', 'db_name': 'mnx', 'id': 'MNXM4', 'inchi': 'InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+'}], 'products': [{'name': 'four', 'db_name': 'mnx', 'id': 'MNXM5', 'inchi': 'InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+'}]}], 'search_type': 'pathway'}, 'output_type': {'output_format': 'csv'}, 'adv': {'pathway_id': 'rp_pathway'}}
'''



##
# (measured_rpsbml_path, inputTar, output, output_format, pathway_id='rp_pathway')
# TODO: change this to JSON instead of JSON
if __name__ == "__main__":
    input_json_path = sys.argv[1]
    input_path = sys.argv[2]
    output_path = sys.argv[3]
    dict_input = json.load(open(input_json_path, 'r'))
    global_match = {}
    logging.debug(input_path)
    logging.debug(output_path)
    logging.debug(dict_input)
    ############## based on the input format run the find pathway algorithm ######
    with tempfile.TemporaryDirectory() as tmpInputFolder:
        inputTar = None
        if dict_input['input_type']['input_format']=='tar':
            inputTar = input_path 
        elif dict_input['input_type']['input_format']=='sbml':
            #make the tar.xz 
            inputTar = tmpInputFolder+'/tmp_input.tar.xz'
            with tarfile.open(inputTar, mode='w:gz') as tf:
                info = tarfile.TarInfo('single.rpsbml.xml') #need to change the name since galaxy creates .dat files
                info.size = os.path.getsize(input_path)
                tf.addfile(tarinfo=info, fileobj=open(input_path, 'rb'))
        else:
            logging.debug('Cannot identify the input/output format: '+str(dict_input['input_type']['input_format']))
            exit(1)
        if dict_input['search']['search_type']=='species':
            measured_rpsbml = rpFindPathwayServe.makeSpecies(dict_input['search'], dict_input['adv']['pathway_id'], dict_input['adv']['species_group_id'])
            global_match = rpFindPathwayServe.findSpecies(measured_rpsbml,
                                                         inputTar,
                                                         dict_input['adv']['pathway_id'],
                                                         dict_input['adv']['species_group_id'])

        elif dict_input['search']['search_type']=='reaction':
            measured_rpsbml = rpFindPathwayServe.makeReaction(dict_input['search'], dict_input['adv']['pathway_id'], dict_input['adv']['species_group_id'])
            global_match = rpFindPathwayServe.findReaction(measured_rpsbml,
                                                          inputTar,
                                                          dict_input['adv']['pathway_id'],
                                                          dict_input['adv']['species_group_id'])
        elif dict_input['search']['search_type']=='pathway':
            measured_rpsbml = rpFindPathwayServe.makePathway(dict_input['search'],
                                                             dict_input['adv']['pathway_id'],
                                                             dict_input['adv']['species_group_id'])
            if dict_input['search']['ordered']=='True' or dict_input['search']['ordered']=='true' or dict_input['search']['ordered']=='T' or dict_input['search']['ordered']==True:
                global_match = rpFindPathwayServe.findOrderedPathway(measured_rpsbml,
                                                                    inputTar,
                                                                    dict_input['adv']['pathway_id'],
                                                                    dict_input['adv']['species_group_id'])
            elif dict_input['search']['ordered']=='False' or dict_input['search']['ordered']=='false' or dict_input['search']['ordered']=='F' or dict_input['search']['ordered']==False:
                global_match = rpFindPathwayServe.findReactions(measured_rpsbml,
                                                                inputTar,
                                                                dict_input['adv']['pathway_id'])
            else:
                logging.error('Cannot detect the ordered input: '+str(dict_input['ordered']['ordered']))
                exit(1)
        else:
            logging.error('Cannot detect the search type: '+str(dict_input['search']['search_type']))
            exit(1)
    #output the results
    #'output_type': {'output_format': 'csv'}
    if dict_input['output_type']['output_format']=='csv':
        if dict_input['search']['search_type']=='species':
            with open(output_path, 'w') as fp:
                fw = csv.writer(fp, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                fw.writerow(['pathway_id', 'species_id', 'score'])
                for rpsbml_id in global_match:
                    assert len(global_match[rpsbml_id])==1
                    for match_spe in global_match[rpsbml_id][list(global_match[rpsbml_id].keys())[0]]:
                        fw.writerow([rpsbml_id, match_spe, global_match[rpsbml_id][list(global_match[rpsbml_id].keys())[0]][match_spe]])
        elif dict_input['search']['search_type']=='reaction':
            with open(output_path, 'w') as fp:
                fw = csv.writer(fp, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                fw.writerow(['pathway_id', 'reaction_id', 'score'])
                for rpsbml_id in global_match:
                    if global_match[rpsbml_id]['id']:
                        for match_reac in global_match[rpsbml_id]['id']:
                            fw.writerow([rpsbml_id, match_reac, global_match[rpsbml_id]['score']])
        elif dict_input['search']['search_type']=='pathway':
            with open(output_path, 'w') as fp:
                fw = csv.writer(fp, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                fw.writerow(['pathway_id', 'score'])
                for rpsbml_id in global_match:
                    fw.writerow([rpsbml_id, global_match[rpsbml_id][0]]) 
        else:
            logging.error('Cannot detect the search type: '+str(dict_input['search']['search_type']))
            exit(1)
    elif dict_input['output_type']['output_format']=='json':
        with open(output_path, 'w') as fp:
            json.dump(global_match, fp)
    else:
        logging.error('Cannot interpret the output_format: '+str(dict_input['output_type']['output_format']))
        exit(1)
