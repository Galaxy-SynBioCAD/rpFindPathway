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
import rpToolServe

##
# (measured_rpsbml_path, inputTar, output, output_format, pathway_id='rp_pathway')
# TODO: change this to JSON instead of JSON
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Comapre one or more rpSBML')
    parser.add_argument('-input_rpsbml', type=str)
    parser.add_argument('-input_target', type=str)
    parser.add_argument('-target_format', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    params = parser.parse_args()
    ############## based on the input format run the find pathway algorithm ######
    if params.target_format=='tar':
        dict_result = rpToolServe.runFindPathway_hdd(params.input_rpsbml, params.input_target, params.pathway_id, params.species_group_id)
    elif params.target_format=='sbml':
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            inputTar = tmpOutputFolder+'/tmp_input.tar.xz'
            with tarfile.open(inputTar, mode='w:gz') as tf:
                info = tarfile.TarInfo('single.rpsbml.xml') #need to change the name since galaxy creates .dat files
                info.size = os.path.getsize(params.input_target)
                tf.addfile(tarinfo=info, fileobj=open(params.input_target, 'rb'))
            dict_result = rpToolServe.runFindPathway_hdd(params.input_rpsbml, inputTar, params.pathway_id, params.species_group_id)
    elif params.target_format=='json':
        #in this case you need to convert, using rpReader, the pathway to 
    else:
        logging.error('Cannot identify the input/output format: '+str(params.target_format))
        exit(1)
    with open(params.output, 'w') as fp:
        json.dump(dict_result, fp)
        exit(1)
