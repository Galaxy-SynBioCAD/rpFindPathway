import numpy as np
import tempfile
import logging
import pandas as pd
from sklearn.metrics import jaccard_score
import rpGraph

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

############ GRAPH WAY (ORDERED) ##########

## Match all the measured chemical species to the simulated chemical species between two SBML 
#
# TODO: for all the measured species compare with the simualted one. Then find the measured and simulated species that match the best and exclude the 
# simulated species from potentially matching with another
#
def compareSpecies_graph(measured_rpsbml, sim_rpsbml):
    ############## compare species ###################
    meas_sim = {}
    sim_meas = {}
    species_match = {}
    for measured_species in measured_rpsbml.model.getListOfSpecies():
        logging.info('--- Trying to match chemical species: '+str(measured_species.getId())+' ---')
        meas_sim[measured_species.getId()] = {}
        species_match[measured_species.getId()] = {}
        #species_match[measured_species.getId()] = {'id': None, 'score': 0.0, 'found': False}
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
            #### MIRIAM ####
            if sim_rpsbml.compareMIRIAMAnnotations(measured_species.getAnnotation(), sim_species.getAnnotation()):
                logging.info('--> Matched MIRIAM: '+str(sim_species.getId()))
                #meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.4
                meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.2+0.2*jaccardMIRIAM(sim_miriam_annot, measured_miriam_annot)
                meas_sim[measured_species.getId()][sim_species.getId()]['found'] = True
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
            species_match[meas] = {}
            for unique_spe in unique[meas]:
                species_match[meas][unique_spe] = round(meas_sim[meas][unique[meas][0]]['score'], 5)
        else:
            logging.warning('Cannot find a species match for the measured species: '+str(meas))
    logging.info('#########################')
    logging.info('species_match:')
    logging.info(species_match)
    logging.info('-----------------------')
    return species_match


## Compare individual reactions and see if the measured pathway is contained within the simulated one
#
# Note that we assure that the match is 1:1 between species using the species match
#
def compareReaction_graph(species_match, meas_reac, sim_reac):
    scores = []
    all_match = True
    ########### reactants #######
    ignore_reactants = []
    for meas_reactant in meas_reac.getListOfReactants():
        if meas_reactant.species in species_match:
            spe_found = False
            for sim_reactant in sim_reac.getListOfReactants():
                if sim_reactant.species in species_match[meas_reactant.species] and not sim_reactant.species in ignore_reactants:
                    scores.append(species_match[meas_reactant.species][sim_reactant.species])
                    ignore_reactants.append(sim_reactant.species)
                    spe_found = True
                    break
            if not spe_found:
                scores.append(0.0)
                all_match = False
        else:
            logging.warning('Cannot find the measured species '+str(meas_reactant.species)+' in the the matched species: '+str(species_match))
            scores.append(0.0)
            all_match = False
    #products
    ignore_products = []
    for meas_product in meas_reac.getListOfProducts():
        if meas_product.species in species_match:
            pro_found = False
            for sim_product in sim_reac.getListOfProducts():
                if sim_product.species in species_match[meas_product.species] and not sim_product.species in ignore_products:
                    scores.append(species_match[meas_product.species][sim_product.species])
                    ignore_products.append(sim_product.species)
                    pro_found = True
                    break
            if not pro_found:
                scores.append(0.0)
                all_match = False
        else:
            logging.warning('Cannot find the measured species '+str(meas_product.species)+' in the the matched species: '+str(species_match))
            scores.append(0.0)
            all_match = False
    return np.mean(scores), all_match



########## EC number ############
def compareEC(meas_reac_miriam, sim_reac_miriam):
    #Warning we only match a single reaction at a time -- assume that there cannot be more than one to match at a given time
    if 'ec-code' in meas_reac_miriam and 'ec-code' in sim_reac_miriam:
        measured_frac_ec = [[y for y in ec.split('.') if not y=='-'] for ec in meas_reac_miriam['ec-code']]
        sim_frac_ec = [[y for y in ec.split('.') if not y=='-'] for ec in sim_reac_miriam['ec-code']]
        #complete the ec numbers with None to be length of 4
        for i in range(len(measured_frac_ec)):
            for y in range(len(measured_frac_ec[i]), 4):
                measured_frac_ec[i].append(None)
        for i in range(len(sim_frac_ec)):
            for y in range(len(sim_frac_ec[i]), 4):
                sim_frac_ec[i].append(None)
        logging.info('Measured: ')
        logging.info(measured_frac_ec)
        logging.info('Simulated: ')
        logging.info(sim_frac_ec)
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
        return best_ec_compare['score']
    else:
        logging.warning('One of the two reactions does not have any EC entries.\nMeasured: '+str(meas_reac_miriam)+' \nSimulated: '+str(sim_reac_miriam))
        return 0.0

##
#
# Note: remember that we are trying to find the measured species equivalent
#
def compareOrderedReactions(measured_rpsbml,
                            sim_rpsbml,
                            pathway_id='rp_pathway',
                            species_group_id='central_species'):
    measured_rpgraph = rpGraph.rpGraph(measured_rpsbml, pathway_id, species_group_id)
    sim_rpgraph = rpGraph.rpGraph(sim_rpsbml, pathway_id, species_group_id)
    measured_ordered_reac = measured_rpgraph.orderedRetroReactions()
    sim_ordered_reac = sim_rpgraph.orderedRetroReactions()
    species_match = compareSpecies_graph(measured_rpsbml, sim_rpsbml)
    scores = []
    if len(measured_ordered_reac)>len(sim_ordered_reac):
        for i in range(len(sim_ordered_reac)):
            logging.info('measured_ordered_reac['+str(i)+']: '+str(measured_ordered_reac))
            logging.info('sim_ordered_reac['+str(i)+']: '+str(sim_ordered_reac))
            spe_score, is_full_match = compareReaction_graph(species_match, 
                                                             measured_rpsbml.model.getReaction(measured_ordered_reac[i]), 
                                                             sim_rpsbml.model.getReaction(sim_ordered_reac[i]))
            ec_score = compareEC(measured_rpsbml.readMIRIAMAnnotation(measured_rpsbml.model.getReaction(measured_ordered_reac[i]).getAnnotation()),
                                 sim_rpsbml.readMIRIAMAnnotation(sim_rpsbml.model.getReaction(sim_ordered_reac[i]).getAnnotation()))
            scores.append(np.average([spe_score, ec_score], weights=[0.8, 0.2]))
        return np.mean(scores)*( 1.0-np.abs(len(measured_ordered_reac)-len(sim_ordered_reac))/len(measured_ordered_reac) ), False
    elif len(measured_ordered_reac)<len(sim_ordered_reac):
        for i in range(len(measured_ordered_reac)):
            logging.info('measured_ordered_reac['+str(i)+']: '+str(measured_ordered_reac))
            logging.info('sim_ordered_reac['+str(i)+']: '+str(sim_ordered_reac))
            spe_score, is_full_match = compareReaction_graph(species_match, 
                                                             measured_rpsbml.model.getReaction(measured_ordered_reac[i]), 
                                                             sim_rpsbml.model.getReaction(sim_ordered_reac[i]))
            ec_score = compareEC(measured_rpsbml.readMIRIAMAnnotation(measured_rpsbml.model.getReaction(measured_ordered_reac[i]).getAnnotation()),
                                 sim_rpsbml.readMIRIAMAnnotation(sim_rpsbml.model.getReaction(sim_ordered_reac[i]).getAnnotation()))
            scores.append(np.average([spe_score, ec_score], weights=[0.8, 0.2]))
        return np.mean(scores)*( 1.0-np.abs(len(measured_ordered_reac)-len(sim_ordered_reac))/len(sim_ordered_reac) ), False
    #if the pathways are of the same length is the only time when the match may be perfect
    elif len(measured_ordered_reac)==len(sim_ordered_reac):
        perfect_match = True
        for i in range(len(sim_ordered_reac)):
            logging.info('measured_ordered_reac['+str(i)+']: '+str(measured_ordered_reac))
            logging.info('sim_ordered_reac['+str(i)+']: '+str(sim_ordered_reac))
            spe_score, is_full_match = compareReaction_graph(species_match, 
                                                             measured_rpsbml.model.getReaction(measured_ordered_reac[i]), 
                                                             sim_rpsbml.model.getReaction(sim_ordered_reac[i]))
            ec_score = compareEC(measured_rpsbml.readMIRIAMAnnotation(measured_rpsbml.model.getReaction(measured_ordered_reac[i]).getAnnotation()),
                                 sim_rpsbml.readMIRIAMAnnotation(sim_rpsbml.model.getReaction(sim_ordered_reac[i]).getAnnotation()))
            scores.append(np.average([spe_score, ec_score], weights=[0.8, 0.2]))
            if ec_score==0.0 and not is_full_match:
                prefect_match = False
        return np.mean(scores), perfect_match


