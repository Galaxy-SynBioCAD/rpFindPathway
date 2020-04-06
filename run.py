#!/usr/bin/env python3
"""
Created on March 17 2020

@author: Melchior du Lac
@description: rpFindPathway

"""
import argparse
import tempfile
import os
import logging
import shutil
import docker


##
#
#
def main(input_rpsbml,
         input_target,
         target_format,
         output,
         pathway_id):
    docker_client = docker.from_env()
    image_str = 'brsynth/rpfindpathways-standalone:dev'
    try:
        image = docker_client.images.get(image_str)
    except docker.errors.ImageNotFound:
        logging.warning('Could not find the image, trying to pull it')
        try:
            docker_client.images.pull(image_str)
            image = docker_client.images.get(image_str)
        except docker.errors.ImageNotFound:
            logging.error('Cannot pull image: '+str(image_str))
            exit(1)
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        shutil.copy(input_rpsbml, tmpOutputFolder+'/input_rpsbml.dat')
        shutil.copy(input_target, tmpOutputFolder+'/input_target.dat')
        command = ['python',
                   '/home/tool_rpFindPathway.py',
                   '-input_rpsbml',
                   '/home/tmp_output/input_rpsbml.dat',
                   '-input_target',
                   '/home/tmp_output/input_target.dat',
                   '-target_format',
                   str(target_format),
                   '-output',
                   '/home/tmp_output/output.dat',
                   '-pathway_id',
                   str(pathway_id)]
        container = docker_client.containers.run(image_str,
                                                 command,
                                                 detach=True,
                                                 stderr=True,
                                                 volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        container.wait()
        err = container.logs(stdout=False, stderr=True)
        print(err)
        shutil.copy(tmpOutputFolder+'/output.dat', output)
        container.remove()



##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Convert the results of RP2 and rp2paths to SBML files')
    parser.add_argument('-input_rpsbml', type=str)
    parser.add_argument('-input_target', type=str)
    parser.add_argument('-target_format', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    params = parser.parse_args()
    main(params.input_rpsbml,
         params.input_target,
         params.target_format,
         params.output,
         params.pathway_id)
