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

import argparse

sys.path.insert(0, '/home/')
import rpFindPathwayServe


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
    logging.error(input_path)
    logging.error(output_path)
    logging.error(dict_input)
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
            logging.error('Cannot identify the input/output format: '+str(dict_input['input_type']['input_format']))
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
            if dict_input['ordered']=='True' or dict_input['ordered']=='true' or dict_input['ordered']=='T':
                global_match = rpFindPathwayServe.findOrderedPathway(measured_rpsbml,
                                                                    inputTar,
                                                                    dict_input['adv']['pathway_id'],
                                                                    dict_input['adv']['species_group_id'])
            elif dict_input['ordered']=='False' or dict_input['ordered']=='false' or dict_input['ordered']=='F':
                global_match = findReactions(measured_rpsbml,
                                            inputTar,
                                            dict_input['adv']['pathway_id'])
            else:
                logging.error('Cannot detect the ordered input: '+str(dict_input['ordered']))
                exit(1)
        else:
            logging.error('Cannot detect the search type: '+str(dict_input['search']['search_type']))
            exit(1)
    #output the results
    with open(output_path, 'w') as fp:
        json.dump(global_match, fp)
