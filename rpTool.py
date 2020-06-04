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


## Function to find the unique species
#
# pd_matrix is organised such that the rows are the simulated species and the columns are the measured ones
#
def findUniqueRowColumn(pd_matrix):
    logging.info(pd_matrix)
    to_ret = {}
    ######################## filter by the global top values ################
    logging.info('################ Filter best #############')
    #transform to np.array
    x = pd_matrix.values
    #resolve the rouding issues to find the max
    x = np.around(x, decimals=5)
    #first round involves finding the highest values and if found set to 0.0 the rows and columns (if unique)
    top = np.where(x==np.max(x))
    #as long as its unique keep looping
    if np.count_nonzero(x)==0:
        return to_ret
    while len(top[0])==1 and len(top[1])==1: 
        if np.count_nonzero(x)==0:
            return to_ret
        pd_entry = pd_matrix.iloc[[top[0][0]],[top[1][0]]]
        row_name = str(pd_entry.index[0])
        col_name = str(pd_entry.columns[0])
        if col_name in to_ret:
            logging.error('Overwriting (1): '+str(col_name))
            logging.error(x)
        to_ret[col_name] = [row_name]
        #delete the rows and the columns 
        logging.info('==================')
        logging.info('Column: '+str(col_name))
        logging.info('Row: '+str(row_name))
        pd_matrix.loc[:, col_name] = 0.0
        pd_matrix.loc[row_name, :] = 0.0
        x = pd_matrix.values
        x = np.around(x, decimals=5)
        top = np.where(x==np.max(x))
        logging.info(pd_matrix)
        logging.info(top)
        logging.info('==================')
    #################### filter by columns (measured) top values ##############
    logging.info('################ Filter by column best ############')
    x = pd_matrix.values
    x = np.around(x, decimals=5)
    if np.count_nonzero(x)==0:
        return to_ret
    reloop = True
    while reloop:
        if np.count_nonzero(x)==0:
            return to_ret
        reloop = False
        for col in range(len(x[0])):
            if np.count_nonzero(x[:,col])==0:
                continue
            top_row = np.where(x[:,col]==np.max(x[:,col]))[0]
            if len(top_row)==1:
                top_row = top_row[0]
                #if top_row==0.0:
                #    continue
                #check to see if any other measured pathways have the same or larger score (accross)
                row = list(x[top_row, :])
                #remove current score consideration
                row.pop(col)
                if max(row)>=x[top_row, col]:
                    logging.info('For col '+str(col)+' there are either better or equal values: '+str(row))
                    logging.info(x)
                    continue
                #if you perform any changes on the rows and columns, then you can perform the loop again
                reloop = True
                pd_entry = pd_matrix.iloc[[top_row],[col]]
                logging.info('==================')
                row_name = pd_entry.index[0]
                col_name = pd_entry.columns[0]
                logging.info('Column: '+str(col_name))
                logging.info('Row: '+str(row_name))
                if col_name in to_ret:
                    logging.error('Overwriting (2): '+str(col_name))
                    logging.error(pd_matrix.values)
                to_ret[col_name] = [row_name]
                #delete the rows and the columns 
                pd_matrix.loc[:, col_name] = 0.0
                pd_matrix.loc[row_name, :] = 0.0
                x = pd_matrix.values
                x = np.around(x, decimals=5)
                logging.info(pd_matrix)
                logging.info('==================')
    ################## laslty if there are multiple values that are not 0.0 then account for that ######
    logging.info('################# get the rest ##########')
    x = pd_matrix.values
    x = np.around(x, decimals=5)
    if np.count_nonzero(x)==0:
        return to_ret
    for col in range(len(x[0])):
        if not np.count_nonzero(x[:,col])==0:
            top_rows = np.where(x[:,col]==np.max(x[:,col]))[0]
            if len(top_rows)==1:
                top_row = top_rows[0]
                pd_entry = pd_matrix.iloc[[top_row],[col]]
                row_name = pd_entry.index[0]
                col_name = pd_entry.columns[0]
                if not col_name in to_ret:
                    to_ret[col_name] = [row_name]
                else:
                    logging.warning('At this point should never have only one: '+str(x[:,col]))
                    logging.warning(x)
            else:
                for top_row in top_rows:
                    pd_entry = pd_matrix.iloc[[top_row],[col]]
                    row_name = pd_entry.index[0]
                    col_name = pd_entry.columns[0]
                    if not col_name in to_ret:
                        to_ret[col_name] = []
                    to_ret[col_name].append(row_name)
    logging.info(pd_matrix)
    logging.info('###################')
    return to_ret

## Match all the measured chemical species to the simulated chemical species between two SBML 
#
# TODO: for all the measured species compare with the simualted one. Then find the measured and simulated species that match the best and exclude the 
# simulated species from potentially matching with another
#
def compareSpecies(measured_rpsbml, sim_rpsbml):
    ############## compare species ###################
    meas_sim = {}
    sim_meas = {}
    species_match = {}
    for measured_species in measured_rpsbml.model.getListOfSpecies():
        logging.info('--- Trying to match chemical species: '+str(measured_species.getId())+' ---')
        meas_sim[measured_species.getId()] = {}
        species_match[measured_species.getId()] = {'id': None, 'score': 0.0, 'found': False}
        #TODO: need to exclude from the match if a simulated chemical species is already matched with a higher score to another measured species
        for sim_species in sim_rpsbml.model.getListOfSpecies():
            meas_sim[measured_species.getId()][sim_species.getId()] = {'score': 0.0, 'found': False}
            if not sim_species.getId() in sim_meas:
                sim_meas[sim_species.getId()] = {}
            sim_meas[sim_species.getId()][measured_species.getId()] = {'score': 0.0, 'found': False}
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
                #meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.4
                meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.2+0.2*jaccardMIRIAM(sim_miriam_annot, measured_miriam_annot)
                meas_sim[measured_species.getId()][sim_species.getId()]['found'] = True
            #logging.info('=====================================')
            ##### InChIKey ##########
            #find according to the inchikey -- allow partial matches
            if 'inchikey' in measured_brsynth_annot and 'inchikey' in sim_rpsbml_brsynth_annot:
                measured_inchikey_split = measured_brsynth_annot['inchikey'].split('-')
                sim_rpsbml_inchikey_split = sim_rpsbml_brsynth_annot['inchikey'].split('-')
                if measured_inchikey_split[0]==sim_rpsbml_inchikey_split[0]:
                    logging.info('Matched first layer InChIkey: ('+str(measured_brsynth_annot['inchikey'])+' -- '+str(sim_rpsbml_brsynth_annot['inchikey'])+')')
                    meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                    if measured_inchikey_split[1]==sim_rpsbml_inchikey_split[1]:
                        logging.info('Matched second layer InChIkey: ('+str(measured_brsynth_annot['inchikey'])+' -- '+str(sim_rpsbml_brsynth_annot['inchikey'])+')')
                        meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                        meas_sim[measured_species.getId()][sim_species.getId()]['found'] = True
                        if measured_inchikey_split[2]==sim_rpsbml_inchikey_split[2]:
                            logging.info('Matched third layer InChIkey: ('+str(measured_brsynth_annot['inchikey'])+' -- '+str(sim_rpsbml_brsynth_annot['inchikey'])+')')
                            meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                            meas_sim[measured_species.getId()][sim_species.getId()]['found'] = True
            sim_meas[sim_species.getId()][measured_species.getId()]['score'] = meas_sim[measured_species.getId()][sim_species.getId()]['score']
            sim_meas[sim_species.getId()][measured_species.getId()]['found'] = meas_sim[measured_species.getId()][sim_species.getId()]['found']
    #build the matrix to send
    meas_sim_mat = {}
    for i in meas_sim:
        meas_sim_mat[i] = {}
        for y in meas_sim[i]:
            meas_sim_mat[i][y] = meas_sim[i][y]['score']
    unique = findUniqueRowColumn(pd.DataFrame(meas_sim_mat))
    logging.info('findUniqueRowColumn:')
    logging.info(unique)
    for meas in meas_sim:
        if meas in unique:
            species_match[meas]['id'] = unique[meas]
            species_match[meas]['score'] = round(meas_sim[meas][unique[meas][0]]['score'], 5)
            species_match[meas]['found'] = meas_sim[meas][unique[meas][0]]['found']
    logging.info('#########################')
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
    tmp_reaction_match = {}
    meas_sim = {}
    sim_meas = {}
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        logging.info('Species match of measured reaction: '+str(measured_reaction_id))
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        measured_reaction_miriam = measured_rpsbml.readMIRIAMAnnotation(measured_reaction.getAnnotation())
        ################ construct the dict transforming the species #######
        meas_sim[measured_reaction_id] = {}
        tmp_reaction_match[measured_reaction_id] = {}
        for sim_reaction_id in sim_rpsbml.readRPpathwayIDs(pathway_id):
            if not sim_reaction_id in sim_meas:
                sim_meas[sim_reaction_id] = {}
            sim_meas[sim_reaction_id][measured_reaction_id] = {}
            meas_sim[measured_reaction_id][sim_reaction_id] = {}
            logging.info('\t=========== '+str(sim_reaction_id)+' ==========')
            logging.info('\t+++++++ Species match +++++++')
            tmp_reaction_match[measured_reaction_id][sim_reaction_id] = {'reactants': {},
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
            logging.info('\tspecies_match: '+str(species_match))
            logging.info('\tspecies_match: '+str(species_match.keys()))
            logging.info('\tsim_reactants_id: '+str(sim_reactants_id))
            logging.info('\tmeasured_reactants_id: '+str([i.species for i in measured_reaction.getListOfReactants()]))
            logging.info('\tsim_products_id: '+str(sim_products_id))
            logging.info('\tmeasured_products_id: '+str([i.species for i in measured_reaction.getListOfProducts()]))
            #ensure that the match is 1:1
            #1)Here we assume that a reaction cannot have twice the same species
            cannotBeSpecies = []
            #if there is a match then we loop again since removing it from the list of potential matches would be appropriate
            keep_going = True
            while keep_going:
                logging.info('\t\t----------------------------')
                keep_going = False
                for reactant in measured_reaction.getListOfReactants():
                    logging.info('\t\tReactant: '+str(reactant.species))
                    #if a species match has been found AND if such a match has been found
                    founReaIDs = [tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants'][i]['id'] for i in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants'] if not tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants'][i]['id']==None]
                    logging.info('\t\tfounReaIDs: '+str(founReaIDs))
                    if reactant.species in species_match and not species_match[reactant.species]['id']==None and not reactant.species in founReaIDs:
                        #return all the similat entries
                        speMatch = list(set(species_match[reactant.species]['id'])&set(sim_reactants_id))
                        speMatch = list(set(speMatch)-set(cannotBeSpecies))
                        logging.info('\t\tspeMatch: '+str(speMatch))
                        if len(speMatch)==1:
                            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants'][reactant.species] = {'id': speMatch[0], 'score': species_match[reactant.species]['score'], 'found': True}
                            cannotBeSpecies.append(speMatch[0])
                            keep_going = True
                            logging.info('\t\tMatched measured reactant species: '+str(reactant.species)+' with simulated species: '+str(speMatch[0])) 
                        elif not reactant.species in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants']:
                            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants'][reactant.species] = {'id': None, 'score': 0.0, 'found': False}
                            #logging.info('\t\tCould not find the folowing measured reactant in the currrent reaction: '+str(reactant.species))
                    elif not reactant.species in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants']:
                        tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants'][reactant.species] = {'id': None, 'score': 0.0, 'found': False}
                        #logging.info('\t\tCould not find the following measured reactant in the matched species: '+str(reactant.species))
                for product in measured_reaction.getListOfProducts():
                    logging.info('\t\tProduct: '+str(product.species))
                    foundProIDs = [tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products'][i]['id'] for i in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products'] if not tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products'][i]['id']==None]
                    logging.info('\t\tfoundProIDs: '+str(foundProIDs))
                    if product.species in species_match and not species_match[product.species]['id']==None and not product.species in foundProIDs:
                        #return all the similat entries
                        speMatch = list(set(species_match[product.species]['id'])&set(sim_products_id))
                        speMatch = list(set(speMatch)-set(cannotBeSpecies))
                        logging.info('\t\tspeMatch: '+str(speMatch))
                        if len(speMatch)==1:
                            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products'][product.species] = {'id': speMatch[0], 'score': species_match[product.species]['score'], 'found': True}
                            cannotBeSpecies.append(speMatch[0])
                            keep_going = True
                            logging.info('\t\tMatched measured product species: '+str(product.species)+' with simulated species: '+str(speMatch[0]))    
                        elif not product.species in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products']:
                            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products'][product.species] = {'id': None, 'score': 0.0, 'found': False}
                            #logging.info('\t\tCould not find the following measured product in the matched species: '+str(product.species))
                    elif not product.species in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products']:
                        tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products'][product.species] = {'id': None, 'score': 0.0, 'found': False}
                        #logging.info('\t\tCould not find the following measured product in the matched species: '+str(product.species))
                logging.info('\t\tcannotBeSpecies: '+str(cannotBeSpecies))
            reactants_score = [tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants'][i]['score'] for i in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants']]
            reactants_found = [tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants'][i]['found'] for i in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants']]
            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['reactants_score'] = np.mean(reactants_score)
            products_score = [tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products'][i]['score'] for i in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products']]
            products_found = [tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products'][i]['found'] for i in tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products']]
            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['products_score'] = np.mean(products_score)
            ### calculate pathway species score
            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['species_score'] = np.mean(reactants_score+products_score)
            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['species_std'] = np.std(reactants_score+products_score)
            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['species_reaction'] = sim_reaction_id
            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['found'] = all(reactants_found+products_found)
            #tmp_reaction_match[measured_reaction_id][sim_reaction_id]['found'] = True
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
                    #perfect match - one can have multiple ec score per reaction - keep the score for the highest matching one
                    logging.info('\t~~~~~~~~~~~~~~~~~~~~')
                    logging.info('\tMeasured EC: '+str(measured_ec))
                    logging.info('\tSimulated EC: '+str(sim_ec)) 
                    measured_frac_ec = [[y for y in ec.split('.') if not y=='-'] for ec in measured_reaction_miriam['ec-code']]
                    sim_frac_ec = [[y for y in ec.split('.') if not y=='-'] for ec in sim_reaction_miriam['ec-code']]
                    #complete the ec numbers with None to be length of 4
                    for i in range(len(measured_frac_ec)):
                        for y in range(len(measured_frac_ec[i]),4):
                            measured_frac_ec[i].append(None)
                    for i in range(len(sim_frac_ec)):
                        for y in range(len(sim_frac_ec[i]),4):
                            sim_frac_ec[i].append(None)
                    logging.info('\t'+str(measured_frac_ec))
                    logging.info('\t'+str(sim_frac_ec))
                    best_ec_compare = {'meas_ec': [], 'sim_ec': [], 'score': 0.0, 'found': False}
                    for ec_m in measured_frac_ec:
                        for ec_s in sim_frac_ec:
                            tmp_score = 0.0
                            for i in range(4):
                                if not ec_m[i]==None and not ec_s[i]==None:
                                    if ec_m[i]==ec_s[i]:
                                        tmp_score += 0.25
                                        if i==2:
                                            best_ec_compare['found'] = True
                                    else:
                                        break
                            if tmp_score>best_ec_compare['score']:
                                best_ec_compare['meas_ec'] = ec_m
                                best_ec_compare['sim_ec'] = ec_s
                                best_ec_compare['score'] = tmp_score
                    logging.info('\t'+str(best_ec_compare))
                    if best_ec_compare['found']:
                        tmp_reaction_match[measured_reaction_id][sim_reaction_id]['found'] = True
                    tmp_reaction_match[measured_reaction_id][sim_reaction_id]['ec_reaction'] = sim_reaction_id
                    tmp_reaction_match[measured_reaction_id][sim_reaction_id]['ec_score'] = best_ec_compare['score']
                    logging.info('\t~~~~~~~~~~~~~~~~~~~~')
            #WRNING: Here 80% for species match and 20% for ec match
            tmp_reaction_match[measured_reaction_id][sim_reaction_id]['score'] = np.average([tmp_reaction_match[measured_reaction_id][sim_reaction_id]['species_score'], tmp_reaction_match[measured_reaction_id][sim_reaction_id]['ec_score']], weights=[0.8, 0.2])
            sim_meas[sim_reaction_id][measured_reaction_id] = tmp_reaction_match[measured_reaction_id][sim_reaction_id]['score']
            meas_sim[measured_reaction_id][sim_reaction_id] = tmp_reaction_match[measured_reaction_id][sim_reaction_id]['score']
    ### matrix compare #####
    unique = findUniqueRowColumn(pd.DataFrame(meas_sim))
    logging.info('findUniqueRowColumn')
    logging.info(unique)
    reaction_match = {} 
    for meas in meas_sim:
        reaction_match[meas] = {'id': None, 'score': 0.0, 'found': False}
        if meas in unique:
            #if len(unique[meas])>1:
            #    logging.warning('Multiple values may match, choosing the first arbitrarily')
            reaction_match[meas]['id'] = unique[meas]
            reaction_match[meas]['score'] = round(tmp_reaction_match[meas][unique[meas][0]]['score'], 5)
            reaction_match[meas]['found'] = tmp_reaction_match[meas][unique[meas][0]]['found']
    #### compile a reaction score based on the ec and species scores
    logging.info(tmp_reaction_match)
    logging.info(reaction_match)
    logging.info('-------------------------------')
    return reaction_match, tmp_reaction_match


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
        if meas_path_length>sim_path_length:
            penalty_length = 1.0-np.abs(meas_path_length-sim_path_length)/meas_path_length
        elif meas_path_length<=sim_path_length:
            penalty_length = 1.0-np.abs(meas_path_length-sim_path_length)/sim_path_length
    species_match = compareSpecies(measured_rpsbml, sim_rpsbml)
    reaction_match, all_rection_match_info = compareReactions(measured_rpsbml, sim_rpsbml, species_match, pathway_id)
    logging.info(penalty_length)
    logging.info([reaction_match[i]['score'] for i in reaction_match])
    logging.info([reaction_match[i]['found'] for i in reaction_match])
    global_score = np.mean([reaction_match[i]['score'] for i in reaction_match]) 
    global_found = [reaction_match[i]['found'] for i in reaction_match]
    if all(global_found):
        return True, np.mean(global_score)*penalty_length, all_rection_match_info
    else:
       return False, np.mean(global_score)*penalty_length, all_rection_match_info
