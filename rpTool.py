import rpSBML
import numpy as np
import tempfile


#TODO: add the MNXR to the validation table and compare (with low score) to see
# if you can recover from associated reaction rule MNXR
# contained within the measured pathway and vise versa
#TODO: make a scoring function taking into account how precise the match is:
#   - EC match: if full match or 3rd match --> full = 1, 3rd = 0.5, no = 0
#   - If the species are same (mean) --> 1-(abs(total_species-found_species)/total_species)
#   - The difference between the measured number of steps and the predicted (normalised) --> 
#     1-(abs(measured-predicted)/measured)
# Result is mean of the three
def compareRPpathways(rpsbml, measured_sbml):
    #return all the species annotations of the RP pathways
    meas_rp = None
    rp_rp = None
    try:
        meas_rp = measured_sbml.readRPspecies()
        found_meas_rp = rpsbml.readRPspecies()
        for meas_step_id in meas_rp:
            meas_rp[meas_step_id]['annotation'] = measured_sbml.model.getReaction(meas_step_id).getAnnotation()
            for spe_name in meas_rp[meas_step_id]['reactants']:
                meas_rp[meas_step_id]['reactants'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
            for spe_name in meas_rp[meas_step_id]['products']:
                meas_rp[meas_step_id]['products'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
        rp_rp = rpsbml.readRPspecies()
        found_meas_rp['length_score'] = 0.0
        found_meas_rp['global_score'] = 0.0
        found_meas_rp['ec_score'] = 0.0
        found_meas_rp['species_score'] = 0.0
        for rp_step_id in rp_rp:
            rp_rp[rp_step_id]['annotation'] = rpsbml.model.getReaction(rp_step_id).getAnnotation()
            found_meas_rp[rp_step_id]['found'] = False
            found_meas_rp[rp_step_id]['rp_step_id'] = None
            for spe_name in rp_rp[rp_step_id]['reactants']:
                rp_rp[rp_step_id]['reactants'][spe_name] = rpsbml.model.getSpecies(spe_name).getAnnotation()
                found_meas_rp[rp_step_id]['reactants'][spe_name] = False
            for spe_name in rp_rp[rp_step_id]['products']:
                rp_rp[rp_step_id]['products'][spe_name] = rpsbml.model.getSpecies(spe_name).getAnnotation()
                found_meas_rp[rp_step_id]['products'][spe_name] = False
        try:
            #IMPORTANT - must not consider the sink
            del rp_rp['targetSink']
        except KeyError:
            logging.error('The generated RP pathway does not have a targetSink')
    except AttributeError:
        logging.error('TODO: debug, for some reason some are passed as None here')
        return False, found_meas_rp
    #def compareRPpathways(self, measured_sbml):
    #compare the number of steps in the pathway
    ################ Length #####################
    if len(meas_rp)<len(rp_rp): #add one for the targetSink
        logging.warning('The predicted pathway is longer than the measured one')
        return False, found_meas_rp
    found_meas_rp['length_score'] = 1.0-(np.abs(len(meas_rp)-len(rp_rp))/len(meas_rp))
    #compare the simulated with the measured. i.e. if simulated is contained within the measured
    rp_reac_num = list(rp_rp.keys())
    rp_reac_num = [int(i.replace('RP', '')) for i in rp_reac_num]
    if len(rp_reac_num)>1:
        rp_reac_num.sort()
    for rp_num in rp_reac_num:
        ############## compare using the reactions ###################
        #for meas_step_id in list(meas_rp.keys()).sort(key=lambda x: int(x[-1]), reverse=True):
        #WARNING: this is instead of a loop through all reactions. We assume that 
        #the predicted pathway is the same as the other
        ###################################
        meas_step_id = 'M'+str(rp_num)
        rp_step_id = 'RP'+str(rp_num)
        ###################################
        source_dict = readMIRIAMAnnotation(rp_rp[rp_step_id]['annotation'])
        target_dict = readMIRIAMAnnotation(meas_rp[meas_step_id]['annotation'])
        try:
            found_meas_rp[rp_step_id]['ec_source'] = source_dict['ec-code']
        except KeyError:
            found_meas_rp[rp_step_id]['ec_source'] = None
        try:
            found_meas_rp[rp_step_id]['ec_target'] = target_dict['ec-code']
        except KeyError:
            found_meas_rp[rp_step_id]['ec_target'] = None
        #test perfect matches (ec or others)
        '''# could have the same uniprot ID and not the same stuff
        if compareMIRIAMAnnotations(rp_rp[rp_step_id]['annotation'], meas_rp[meas_step_id]['annotation']): 
            found_meas_rp[rp_step_id]['rp_step_id'] = meas_step_id
            found_meas_rp[rp_step_id]['found'] = True
            found_meas_rp['ec_score'] += 1.0
        '''
        #test to the third EC number
        if 'ec-code' in source_dict and 'ec-code' in target_dict:
            if len(set(source_dict['ec-code']).intersection(set(target_dict['ec-code'])))>0:
                found_meas_rp[rp_step_id]['rp_step_id'] = meas_step_id
                found_meas_rp[rp_step_id]['found'] = True
                found_meas_rp['ec_score'] += 1.0
            elif [i for i, j in zip(source_dict['ec-code'], target_dict['ec-code']) if i.split('.')[:-1]==j.split('.')[:-1]]:
                found_meas_rp[rp_step_id]['rp_step_id'] = meas_step_id
                found_meas_rp[rp_step_id]['found'] = True
                found_meas_rp['ec_score'] += 0.5
        ############## compare using the species ###################
        #for rp_step_id in rp_rp:
        #for rp_step_id in list(rp_rp.keys()).sort(key=lambda x: int(x[-1])):
        # We test to see if the meas reaction elements all exist in rp reaction and not the opposite
        #because the measured pathways may not contain all the elements
        #WARNING: using the above set species since we assume that the order of the simulations should be the same
        ########## reactants ##########
        for rp_spe_id in rp_rp[rp_step_id]['reactants']:
            for meas_spe_id in meas_rp[meas_step_id]['reactants']:
                if compareMIRIAMAnnotations(meas_rp[meas_step_id]['reactants'][meas_spe_id], rp_rp[rp_step_id]['reactants'][rp_spe_id]):
                    found_meas_rp[rp_step_id]['reactants'][meas_spe_id] = True
                    break
                else:
                    if compareIBISBAAnnotations(meas_rp[meas_step_id]['reactants'][meas_spe_id], rp_rp[rp_step_id]['reactants'][rp_spe_id]):
                        found_meas_rp[rp_step_id]['reactants'][meas_spe_id] = True
                        break
        ########### products ###########
        for rp_spe_id in rp_rp[rp_step_id]['products']:
            for meas_spe_id in meas_rp[meas_step_id]['products']:
                if compareMIRIAMAnnotations(meas_rp[meas_step_id]['products'][meas_spe_id], rp_rp[rp_step_id]['products'][rp_spe_id]):
                    found_meas_rp[rp_step_id]['products'][meas_spe_id] = True
                    break
                else:
                    if compareIBISBAAnnotations(meas_rp[meas_step_id]['products'][meas_spe_id], rp_rp[rp_step_id]['products'][rp_spe_id]):
                        found_meas_rp[rp_step_id]['products'][meas_spe_id] = True
                        break
        ######### test to see the difference
        pro_found = [found_meas_rp[rp_step_id]['products'][i] for i in found_meas_rp[rp_step_id]['products']]
        rea_found = [found_meas_rp[rp_step_id]['reactants'][i] for i in found_meas_rp[rp_step_id]['reactants']]
        if pro_found and rea_found:
            if all(pro_found) and all(rea_found):
                found_meas_rp[rp_step_id]['found'] = True
                found_meas_rp[rp_step_id]['rp_step_id'] = rp_step_id
    #calculate the score
    spe_found = [True if not found_meas_rp[i]['rp_step_id']==None else False for i in rp_rp]
    found_meas_rp['species_score'] = 1.0-(np.abs(len(spe_found)-spe_found.count(True))/len(spe_found))
    found_meas_rp['ec_score'] = found_meas_rp['ec_score']/len(rp_rp)
    ################# Now see if all steps have been found ############
    if all(found_meas_rp[i]['found'] for i in found_meas_rp):
        found_meas_rp['measured_model_id'] = measured_sbml.model.getId()
        found_meas_rp['rp_model_id'] = rpsbml.model.getId()
        found_meas_rp['global_score'] = (found_meas_rp['ec_score']+found_meas_rp['species_score']+found_meas_rp['length_score'])/3.0
        return True, found_meas_rp
    else:
        return False, found_meas_rp


#######################################################################
##################### Detect ##########################################
#######################################################################


def detectPathways_mem(measured_sbml_paths, global_results, global_runtime_stats, status_file, path_id, path_to_res):
    sbml_paths = readrpSBMLtar(path_to_res+'/rpglobalscore_sbml.tar.xz')
    ######## detect wich one is the measured pathway ######
    at_least_one = False
    for rp_rpsbml in sbml_paths:
        #found, measReac_rpReac = sbml_paths[rp_rpsbml].compareRPpathways(measured_rpreader.sbml_paths['measured_'+str(path_id)])
        found, measReac_rpReac = compareRPpathways(sbml_paths[rp_rpsbml], measured_sbml_paths['measured_'+str(path_id)])
        global_results[path_id][rp_rpsbml] = measReac_rpReac
        if found:
            if not at_least_one:
                global_runtime_stats['found'] += 1
            at_least_one = True
            logging.info('Found!')
            logging.info('rp_rpsbml: '+str(rp_rpsbml))
            logging.info('measReac_rpReac: '+str(measReac_rpReac))
            status_file.write('rp_rpsbml: '+str(rp_rpsbml)+'\n')
            status_file.write('measReac_rpReac: '+str(measReac_rpReac)+'\n')
            status_file.write('#########################\n')
    if not at_least_one:
        status_file.write('Did not find any pathways\n')
        global_runtime_stats['notfound'] += 1
    status_file.close()
    sbml_paths = None


def detectPathways_hdd(global_results, global_runtime_stats, status_file, measured_model, inputTar, path_id):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        tar = tarfile.open(inputTar, 'r:xz')
        tar.extractall(path=tmpOutputFolder)
        tar.close()
        measured_rpsbml = rpSBML.rpSBML('measured_'+str(path_id))
        measured_rpsbml.readSBML(measured_model)
        for sbml_path in glob.glob(tmpOutputFolder+'/*'):
            fileName = sbml_path.split('/')[-1].replace('.sbml', '')
            rpsbml = rpSBML.rpSBML('Measured_'+str(path_id))
            rpsbml.readSBML(sbml_path)
            found, measReac_rpReac = compareRPpathways(rpsbml, measured_rpsbml)
            global_results[path_id][fileName] = measReac_rpReac
            if found:
                if not at_least_one:
                    global_runtime_stats['found'] += 1
                at_least_one = True
                logging.info('Found!')
                logging.info('rp_rpsbml: '+str(rp_rpsbml))
                logging.info('measReac_rpReac: '+str(measReac_rpReac))
                status_file.write('rp_rpsbml: '+str(rp_rpsbml)+'\n')
                status_file.write('measReac_rpReac: '+str(measReac_rpReac)+'\n')
                status_file.write('#########################\n')

