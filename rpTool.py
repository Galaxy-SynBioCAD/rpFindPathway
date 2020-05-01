import numpy as np
import tempfile
import logging
import pandas as pd
from sklearn.metrics import jaccard_score

#To undo the logging info
#logging.disable(logging.NOTSET)

## Compute the Jaccard index similarity coefficient score between two MIRIAM dicts
# 1.0 is a perfect match and 0.0 is a complete miss
# We assume that the meas has the "complete" information and simulated has the incomplete info
#TODO: interchange meas and sim since actually the measured has inclomple info and sim has more info
def jaccardMIRIAM(meas_miriam, sim_miriam):
    values = list(set([ x for y in list(meas_miriam.values())+list(sim_miriam.values()) for x in y]))
    meas_data = {}
    sim_data = {}
    for key in set(list(meas_miriam.keys())+list(sim_miriam.keys())):
        tmp_meas_row = []
        tmp_sim_row = []
        for value in values:
            if key in meas_miriam:
                if value in meas_miriam[key]:
                    tmp_meas_row.append(1)
                else:
                    tmp_meas_row.append(0)
            else:
                tmp_meas_row.append(0)
            if key in sim_miriam:
                if value in sim_miriam[key]:
                    tmp_sim_row.append(1)
                else:
                    tmp_sim_row.append(0)
            else:
                tmp_sim_row.append(0)
        meas_data[key] = tmp_meas_row
        sim_data[key] = tmp_sim_row
    meas_dataf = pd.DataFrame(meas_data, index=values)
    sim_dataf = pd.DataFrame(sim_data, index=values)
    #return meas_dataf, sim_dataf, jaccard_score(meas_dataf, sim_dataf, average='weighted')
    return jaccard_score(meas_dataf, sim_dataf, average='weighted')




## Match all the measured chemical species to the simulated chemical species between two SBML 
#
# TODO: for all the measured species compare with the simualted one. Then find the measured and simulated species that match the best and exclude the 
# simulated species from potentially matching with another
#
def compareSpecies(measured_rpsbml, sim_rpsbml):
    ############## compare species ###################
    meas_species_match = {}
    species_match = {}
    for measured_species in measured_rpsbml.model.getListOfSpecies():
        logging.info('--- Trying to match chemical species: '+str(measured_species.getId())+' ---')
        meas_species_match[measured_species.getId()] = {}
        species_match[measured_species.getId()] = {'id': None, 'score': 0.0, 'found': False}
        #TODO: need to exclude from the match if a simulated chemical species is already matched with a higher score to another measured species
        for sim_species in sim_rpsbml.model.getListOfSpecies():
            meas_species_match[measured_species.getId()][sim_species.getId()] = {'score': 0.0, 'found': False}
            measured_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(measured_species.getAnnotation())
            sim_rpsbml_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(sim_species.getAnnotation())
            measured_miriam_annot = sim_rpsbml.readMIRIAMAnnotation(measured_species.getAnnotation())
            sim_miriam_annot = sim_rpsbml.readMIRIAMAnnotation(sim_species.getAnnotation())
            #find according to xref
            #logging.info('=====================================')
            #logging.info('Measured MIRIAM: '+str(sim_rpsbml.readMIRIAMAnnotation(measured_species.getAnnotation())))
            #logging.info('Simulated MIRIAM: '+str(sim_rpsbml.readMIRIAMAnnotation(sim_species.getAnnotation())))
            #### MIRIAM ####
            if sim_rpsbml.compareMIRIAMAnnotations(measured_species.getAnnotation(), sim_species.getAnnotation()):
                logging.info('--> Matched MIRIAM: '+str(sim_species.getId()))
                #meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.4
                meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.4*jaccardMIRIAM(sim_miriam_annot, measured_miriam_annot)
                meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
            #logging.info('=====================================')
            ##### InChIKey ##########
            #find according to the inchikey -- allow partial matches
            if 'inchikey' in measured_brsynth_annot and 'inchikey' in sim_rpsbml_brsynth_annot:
                measured_inchikey_split = measured_brsynth_annot['inchikey'].split('-')
                sim_rpsbml_inchikey_split = sim_rpsbml_brsynth_annot['inchikey'].split('-')
                if measured_inchikey_split[0]==sim_rpsbml_inchikey_split[0]:
                    logging.info('Matched first layer InChIkey: ('+str(measured_brsynth_annot['inchikey'])+' -- '+str(sim_rpsbml_brsynth_annot['inchikey'])+')')
                    meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                    if measured_inchikey_split[1]==sim_rpsbml_inchikey_split[1]:
                        logging.info('Matched second layer InChIkey: ('+str(measured_brsynth_annot['inchikey'])+' -- '+str(sim_rpsbml_brsynth_annot['inchikey'])+')')
                        meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                        meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
                        if measured_inchikey_split[2]==sim_rpsbml_inchikey_split[2]:
                            logging.info('Matched third layer InChIkey: ('+str(measured_brsynth_annot['inchikey'])+' -- '+str(sim_rpsbml_brsynth_annot['inchikey'])+')')
                            meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                            meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
            if species_match[measured_species.getId()]['score']<meas_species_match[measured_species.getId()][sim_species.getId()]['score'] and meas_species_match[measured_species.getId()][sim_species.getId()]['found']:
                logging.info('Measured species '+str(measured_species.getId())+' ('+str(species_match[measured_species.getId()]['score'])+') match has found a better score '+str(sim_species.getId())+' ('+str(meas_species_match[measured_species.getId()][sim_species.getId()]['score'])+')')
                species_match[measured_species.getId()]['id'] = sim_species.getId()
                species_match[measured_species.getId()]['score'] = round(meas_species_match[measured_species.getId()][sim_species.getId()]['score'], 4)
                species_match[measured_species.getId()]['found'] = meas_species_match[measured_species.getId()][sim_species.getId()]['found']
    logging.info('species_match:')
    logging.info(species_match)
    logging.info('-----------------------')
    return species_match


##
# Compare that all the measured species of a reactions are found within sim species to match with a reaction.
# We assume that there cannot be two reactions that have the same species and reactants. This is maintained by SBML
# Compare also by EC number, from the third ec to the full EC
# TODO: need to remove from the list reactions simulated reactions that have matched
def compareReactions(measured_rpsbml, sim_rpsbml, species_match, pathway_id='rp_pathway'):
    ############## compare the reactions #######################
    #construct sim reactions with species
    logging.info('------ Comparing reactions --------')
    #match the reactants and products conversion to sim species
    measured_reactions_match = {}
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        logging.info('Species match of measured reaction: '+str(measured_reaction_id))
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        measured_reaction_miriam = measured_rpsbml.readMIRIAMAnnotation(measured_reaction.getAnnotation())
        ################ construct the dict transforming the species #######
        tmp_measured_reactions_match = {}
        for sim_reaction_id in sim_rpsbml.readRPpathwayIDs(pathway_id):
            logging.info('\t=========== '+str(sim_reaction_id)+' ==========')
            logging.info('\t+++++++ Species match +++++++')
            tmp_measured_reactions_match[sim_reaction_id] = {'reactants': {},
                                                             'reactants_score': 0.0,
                                                             'products': {},
                                                             'products_score': 0.0,
                                                             'species_score': 0.0,
                                                             'species_std': 0.0,
                                                             'species_reaction': None,
                                                             'ec_score': 0.0,
                                                             'ec_reaction': None,
                                                             'score': 0.0,
                                                             'found': False}
            sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
            sim_reactants_id = [reactant.species for reactant in sim_reaction.getListOfReactants()]
            sim_products_id = [product.species for product in sim_reaction.getListOfProducts()]
            ############ species ############
            logging.info('\tspecies_match: '+str(species_match.keys()))
            logging.info('\tsim_reactants_id: '+str(sim_reactants_id))
            logging.info('\tmeasured_reactants_id: '+str([i.species for i in measured_reaction.getListOfReactants()]))
            logging.info('\tsim_products_id: '+str(sim_products_id))
            logging.info('\tmeasured_products_id: '+str([i.species for i in measured_reaction.getListOfProducts()]))
            for reactant in measured_reaction.getListOfReactants():
                if reactant.species in species_match:
                    if species_match[reactant.species]['id'] in sim_reactants_id:
                        tmp_measured_reactions_match[sim_reaction_id]['reactants'][reactant.species] = species_match[reactant.species]
                        logging.info('\t\tMatched measured reactant species: '+str(reactant.species)+' with simulated species: '+str(species_match[reactant.species]['id']))
                    else:
                        tmp_measured_reactions_match[sim_reaction_id]['reactants'][reactant.species] = {'id': None, 'score': 0.0, 'found': False}
                        logging.info('\t\tCould not find the folowing measured reactant in the currrent reaction: '+str(reactant.species))
                else:
                    tmp_measured_reactions_match[sim_reaction_id]['reactants'][reactant.species] = {'id': None, 'score': 0.0, 'found': False}
                    logging.info('\t\tCould not find the following measured reactant in the matched species: '+str(reactant.species))
            for product in measured_reaction.getListOfProducts():
                if product.species in species_match:
                    if species_match[product.species]['id'] in sim_products_id:
                        tmp_measured_reactions_match[sim_reaction_id]['products'][product.species] = species_match[product.species]
                        logging.info('\t\tMatched measured product species: '+str(product.species)+' with simulated species: '+str(species_match[product.species]['id']))
                    else:
                        tmp_measured_reactions_match[sim_reaction_id]['products'][product.species] = {'id': None, 'score': 0.0, 'found': False}
                        logging.info('\t\tCould not find the folowing measured product in the currrent reaction: '+str(product.species))
                else:
                    tmp_measured_reactions_match[sim_reaction_id]['products'][product.species] = {'id': None, 'score': 0.0, 'found': False}
                    logging.info('\t\tCould not find the following measured product in the matched species: '+str(product.species))
            reactants_score = [tmp_measured_reactions_match[sim_reaction_id]['reactants'][i]['score'] for i in tmp_measured_reactions_match[sim_reaction_id]['reactants']]
            tmp_measured_reactions_match[sim_reaction_id]['reactants_score'] = np.mean(reactants_score)
            products_score = [tmp_measured_reactions_match[sim_reaction_id]['products'][i]['score'] for i in tmp_measured_reactions_match[sim_reaction_id]['products']]
            tmp_measured_reactions_match[sim_reaction_id]['products_score'] = np.mean(products_score)
            ### calculate pathway species score
            tmp_measured_reactions_match[sim_reaction_id]['species_score'] = np.mean(reactants_score+products_score)
            tmp_measured_reactions_match[sim_reaction_id]['species_std'] = np.std(reactants_score+products_score)
            tmp_measured_reactions_match[sim_reaction_id]['species_reaction'] = sim_reaction_id
            #tmp_measured_reactions_match[sim_reaction_id]['found'] = True
            #break #only if we assume that one match is all that can happen TODO: calculate all matches and take the highest scoring one
            #continue #if you want the match to be more continuous
            ########## EC number ############
            #Warning we only match a single reaction at a time -- assume that there cannot be more than one to match at a given time
            logging.info('\t+++++ EC match +++++++')
            if 'ec-code' in measured_reaction_miriam:
                sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
                sim_reaction_miriam = sim_rpsbml.readMIRIAMAnnotation(sim_reaction.getAnnotation())
                if 'ec-code' in sim_reaction_miriam:
                    #we only need one match here
                    measured_ec = [i for i in measured_reaction_miriam['ec-code']]
                    sim_ec = [i for i in sim_reaction_miriam['ec-code']]
                    #perfect match - one can have multiple ec score per reaction
                    logging.info('\tMeasured EC: '+str(measured_ec))
                    logging.info('\tSimulated EC: '+str(sim_ec))
                    #change the EC match to compare each layer
                    if any(i in measured_ec for i in sim_ec) and tmp_measured_reactions_match[sim_reaction_id]['ec_score']<1.0:
                        tmp_measured_reactions_match[sim_reaction_id]['found'] = True
                        tmp_measured_reactions_match[sim_reaction_id]['ec_reaction'] = sim_reaction_id
                        tmp_measured_reactions_match[sim_reaction_id]['ec_score'] = 1.0
                        logging.info('\t--> EC match (perfect) between measured reaction: '+str(measured_reaction_id)+' and simulated reaction: '+str(sim_reaction_id))
                        #break #only if you assume that one match is all that is possible
                        #continue if you want the match to be more continuous
                    ### partial match ####
                    #use this and ignore any "-"
                    measured_frac_ec = [i.split('.')[:-1] for i in measured_reaction_miriam['ec-code']]
                    measured_frac_ec = ['.'.join(i) for i in measured_frac_ec]
                    sim_frac_ec = [i.split('.')[:-1] for i in sim_reaction_miriam['ec-code']]
                    sim_frac_ec = ['.'.join(i) for i in sim_frac_ec]
                    #up to the third ec
                    if any(i in measured_frac_ec for i in sim_frac_ec) and tmp_measured_reactions_match[sim_reaction_id]['ec_score']<0.5:
                        tmp_measured_reactions_match[sim_reaction_id]['found'] = True
                        tmp_measured_reactions_match[sim_reaction_id]['ec_reaction'] = sim_reaction_id
                        tmp_measured_reactions_match[sim_reaction_id]['ec_score'] = 0.5
                        logging.info('\t--> EC match (third) between measured reaction: '+str(measured_reaction_id)+' and simulated reaction: '+str(sim_reaction_id))
                        #break #only if you assume that one match is all that is possible
                        #continue #if you want to have more continuous
            #WRNING: Here 80% for species match and 20% for ec match
            tmp_measured_reactions_match[sim_reaction_id]['score'] = np.average([tmp_measured_reactions_match[sim_reaction_id]['species_score'], tmp_measured_reactions_match[sim_reaction_id]['ec_score']], weights=[0.8, 0.2])
        #Select the best scoring pathway among all -- if not WARNING (and consider choosing the pathway that is numbered)
        #TODO: check that the measured pathway is in retro
        logging.info('\ttmp_measured_reactions_match: ')
        logging.info('\t'+str(tmp_measured_reactions_match))
        ordered_keys = [k for k,v in sorted(tmp_measured_reactions_match.items(), key=lambda item:item[1]['score'], reverse=True)]
        logging.info('\t'+str([{i: tmp_measured_reactions_match[i]['score']} for i in ordered_keys]))
        if len(ordered_keys)>=2:
            if not tmp_measured_reactions_match[ordered_keys[0]]['score']==0.0 and not tmp_measured_reactions_match[ordered_keys[1]]['score']==0.0:
                if tmp_measured_reactions_match[ordered_keys[0]]['score']==tmp_measured_reactions_match[ordered_keys[1]]['score']:
                    logging.warning('Multiple simulated reactions have the same score match with '+str(measured_reaction_id))
                    logging.warning([{i: tmp_measured_reactions_match[i]['score']} for i in ordered_keys if tmp_measured_reactions_match[i]['score']>0.0])
                    #hack - take the one that is the same step number
                    logging.info('HACK: taking RP'+measured_reaction_id[-1:])
                    measured_reactions_match[measured_reaction_id] = tmp_measured_reactions_match['RP'+measured_reaction_id[-1:]] 
                    continue
        measured_reactions_match[measured_reaction_id] = tmp_measured_reactions_match[ordered_keys[0]]
    #### compile a reaction score based on the ec and species scores
    logging.info(measured_reactions_match)
    logging.info('-------------------------------')
    return measured_reactions_match


## Compare a measured to sim rpSBML pathway 
#
# Works by trying to find a measured pathway contained within a simulated one. Does not perform perfect match! Since 
# simulated ones contain full cofactors while the measured ones have impefect information
def compareRPpathways(measured_rpsbml, sim_rpsbml, strict_length=False, pathway_id='rp_pathway'):
    logging.info('##################### '+str(sim_rpsbml.model.getId())+' ######################')
    penalty_length = 1.0
    if strict_length:
        if not len(measured_rpsbml.readRPpathwayIDs(pathway_id))==len(sim_rpsbml.readRPpathwayIDs(pathway_id)):
            return False, 0.0, {}
    else:
        #calculate the penatly score for the legnth of the pathways not being the same length
        meas_path_length = len(measured_rpsbml.readRPpathwayIDs(pathway_id))
        sim_path_length = len(sim_rpsbml.readRPpathwayIDs(pathway_id))
        #add a penatly to the length only if the simulated pathway is longer than the measured one
        #if its smaller then we will not be able to retreive all reactions and the scor will not be good in any case
        #if meas_path_length<sim_path_length:
        penalty_length = 1.0-np.abs(meas_path_length-sim_path_length)/meas_path_length
    logging.info('species_match')
    species_match = compareSpecies(measured_rpsbml, sim_rpsbml)
    logging.info('reactions_match')
    reactions_match = compareReactions(measured_rpsbml, sim_rpsbml, species_match, pathway_id)
    global_score = []
    global_found = []
    for measured_reaction_id in reactions_match:
        #make sure that EC and reaction match are the same
        if reactions_match[measured_reaction_id]['ec_reaction'] and reactions_match[measured_reaction_id]['species_reaction']:
            try:
                assert reactions_match[measured_reaction_id]['ec_reaction']==reactions_match[measured_reaction_id]['species_reaction']
            except AssertionError:
                logging.error(measured_rpsbml.model.getId())
                logging.error(sim_rpsbml.model.getId())
                logging.error(reactions_match[measured_reaction_id])
                logging.error(reactions_match[measured_reaction_id]['ec_reaction'])
                logging.error(reactions_match[measured_reaction_id]['species_reaction'])
        #80% for the species and 20% for EC
        global_score.append(np.average([reactions_match[measured_reaction_id]['species_score'], reactions_match[measured_reaction_id]['ec_score']], weights=[0.8, 0.2]))
        global_found.append(reactions_match[measured_reaction_id]['found'])
    if all(global_found):
        return True, np.mean(global_score)*penalty_length, reactions_match
    else:
        return False, np.mean(global_score)*penalty_length, reactions_match
