import rpSBML
import numpy as np
import tempfile
import logging

## Compare a measured to sim rpSBML pathway 
#
# Works by trying to find a measured pathway contained within a simulated one. Does not perform perfect match! Since 
# simulated ones contain full cofactors while the measured ones have impefect information
def compareRPpathways(sim_rpsbml, measured_rpsbml, pathway_id='rp_pathway'):
    ############## compare species ###################
    meas_species_match = {}
    for measured_species in measured_rpsbml.model.getListOfSpecies():
        meas_species_match[measured_species.getId()] = []
        for sim_species in sim_rpsbml.model.getListOfSpecies():
            match_score = 0.0
            measured_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(measured_species.getAnnotation())
            sim_rpsbml_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(sim_species.getAnnotation())
            #find according to xref
            if sim_rpsbml.compareMIRIAMAnnotations(measured_species.getAnnotation(), sim_species.getAnnotation()):
                match_score += 0.3
            #find according to the SMILES
            if measured_brsynth_annot['smiles'] and sim_rpsbml_brsynth_annot['smiles']:
                if measured_brsynth_annot['smiles']==sim_rpsbml_brsynth_annot['smiles']:
                    match_score += 0.2
            #find accorsing to the inchi
            if measured_brsynth_annot['inchi'] and sim_rpsbml_brsynth_annot['inchi']:
                if measured_brsynth_annot['inchi']==sim_rpsbml_brsynth_annot['inchi']:
                    match_score += 0.2
            #find according to the inchikey -- allow partial matches
            if measured_brsynth_annot['inchikey'] and sim_rpsbml_brsynth_annot['inchikey']:
                measured_inchikey_split = measured_brsynth_annot['inchikey'].split('-')
                sim_rpsbml_inchikey_split = sim_rpsbml_brsynth_annot['inchikey'].split('-')
                if measured_inchikey_split[0]==sim_rpsbml_inchikey_split[0]:
                    match_score += 0.1
                    if measured_inchikey_split[1]==sim_rpsbml_inchikey_split[1]:
                        match_score += 0.1
                        if measured_inchikey_split[2]==sim_rpsbml_inchikey_split[2]:
                            match_score += 0.1
            if match_score>0.0:
                #only one match
                try:
                    if meas_species_match[measured_species.getId()][1]<match_score:
                        meas_species_match[measured_species.getId()] = [sim_species.getId(), match_score]
                except IndexError:
                    meas_species_match[measured_species.getId()] = [sim_species.getId(), match_score]
    #TODO: need to deal if there are more than one match
    ############## compare the reactions #######################
    sim_reactions = {}
    for sim_reaction_id in sim_rpsbml.readRPpathwayIDs(pathway_id):
        sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
        sim_reactions[sim_reaction_id] = {'reactants': [], 'products': []}
        for sim_reactant in sim_reaction.getListOfReactants():
            sim_reactions[sim_reaction_id]['reactants'].append(sim_reactant.species)
        for sim_product in sim_reaction.getListOfProducts():
            sim_reactions[sim_reaction_id]['products'].append(sim_product.species)
    measured_reactions_match = {}
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        ################ construct the dict transforming the species #######
        measured_reactions_match[measured_reaction.getId()] = {
                'species_match': {
                    'reactants': [],
                    'reactants_score': [],
                    'products': [],
                    'products_score': [],
                    'score': 0.0,
                    'found': False,
                    'sim_reaction': []},
                'ec_match': {
                    'found': False,
                    'score': 0.0},
                'confidence_score': 0.0,
                'found': False}
        for reactant in measured_reaction.getListOfReactants():
            try:
                measured_reactions_match[measured_reaction.getId()]['species_match']['reactants'].append(
                        meas_species_match[reactant.species][0])
                measured_reactions_match[measured_reaction.getId()]['species_match']['reactants_score'].append(
                        meas_species_match[reactant.species][1])
            except (KeyError, IndexError) as e:
                measured_reactions_match[measured_reaction.getId()]['species_match']['reactants'].append(reactant.species)
                measured_reactions_match[measured_reaction.getId()]['species_match']['reactants_score'].append(0.0)
        for product in measured_reaction.getListOfProducts():
            try:
                measured_reactions_match[measured_reaction.getId()]['species_match']['products'].append(
                        meas_species_match[product.species][0])
                measured_reactions_match[measured_reaction.getId()]['species_match']['products_score'].append(
                        meas_species_match[product.species][1])
            except (KeyError, IndexError) as e:
                measured_reactions_match[measured_reaction.getId()]['species_match']['products'].append(product.species)
                measured_reactions_match[measured_reaction.getId()]['species_match']['products_score'].append(0.0)
    for measured_reaction_id in measured_reactions_match:
        for sim_reaction_id in sim_reactions:
            #perfect species match
            #if not set(measured_reactions_match[measured_reaction_id]['species_match']['reactants'])-set(sim_reactions[sim_reaction_id]['reactants']) and not set(measured_reactions_match[measured_reaction_id]['species_match']['products'])-set(sim_reactions[sim_reaction_id]['products'])):
            #is contained within the 
            if set(measured_reactions_match[measured_reaction_id]['species_match']['reactants']).issubset(sim_reactions[sim_reaction_id]['reactants']) and set(measured_reactions_match[measured_reaction_id]['species_match']['products']).issubset(sim_reactions[sim_reaction_id]['products']):
                measured_reactions_match[measured_reaction_id]['species_match']['sim_reaction'].append(sim_reaction_id)
                #measured_reactions_match[measured_reaction_id]['species_match']['score'] = np.mean(measured_reactions_match[measured_reaction_id]['species_match']['reactants_score']+measured_reactions_match[measured_reaction_id]['species_match']['products_score'])
                #break
    #TODO: rewrite this
    #if all([True if i>0.0 else False for i in measured_reactions_match[measured_reaction_id]['species_match']['reactants_score']+measured_reactions_match[measured_reaction_id]['species_match']['products_score']]):
    all_bool = []
    for i in measured_reactions_match[measured_reaction_id]['species_match']['reactants_score']:
        if i>0.0:
            all_bool.append(True)
        else:
            all_bool.append(False)
    for i in measured_reactions_match[measured_reaction_id]['species_match']['reactants_score']:
        if i>0.0:
            all_bool.append(True)
        else:
            all_bool.append(False)
    if all(all_bool):
        measured_reactions_match[measured_reaction_id]['species_match']['found'] = True
    #if all([measured_reactions_match[measured_reaction_id]['species_match']['reactants_score']
    #    measured_reactions_match[measured_reaction_id]['species_match']['found'] = True
    for measured_reaction_id in measured_reactions_match:
        measured_reactions_match[measured_reaction_id]['species_match']['score'] = np.mean(measured_reactions_match[measured_reaction_id]['species_match']['reactants_score']+measured_reactions_match[measured_reaction_id]['species_match']['products_score'])
    ######## compare by EC number ############
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        measured_reaction_miriam = sim_rpsbml.readMIRIAMAnnotation(measured_reaction.getAnnotation())
        if 'ec-code' in measured_reaction_miriam:
            for sim_reaction_id in sim_rpsbml.readRPpathwayIDs(pathway_id):
                sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
                sim_reaction_miriam = sim_rpsbml.readMIRIAMAnnotation(sim_reaction.getAnnotation())
                if 'ec-code' in sim_reaction_miriam:
                    measured_frac_ec = [i.split('.')[:-1] for i in measured_reaction_miriam['ec-code']]
                    measured_frac_ec = ['.'.join(i) for i in measured_frac_ec]
                    sim_frac_ec = [i.split('.')[:-1] for i in sim_reaction_miriam['ec-code']]
                    sim_frac_ec = ['.'.join(i) for i in sim_frac_ec]
                    ##### try to compare them ####
                    if set(measured_reaction_miriam['ec-code']).issubset(set(sim_reaction_miriam['ec-code'])):
                        #recover only the ones that are in the measured
                        measured_reactions_match[measured_reaction.getId()]['ec_match']['score'] = 1.0
                        measured_reactions_match[measured_reaction.getId()]['ec_match']['found'] = True
                    elif set(measured_frac_ec).issubset(set(sim_frac_ec)):
                        measured_reactions_match[measured_reaction.getId()]['ec_match']['score'] = 0.5
                        measured_reactions_match[measured_reaction.getId()]['ec_match']['found'] = True
        measured_reactions_match[measured_reaction.getId()]['confidence_score'] = np.mean([measured_reactions_match[measured_reaction.getId()]['species_match']['score']]+[measured_reactions_match[measured_reaction.getId()]['ec_match']['score']])
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        measured_reactions_match[measured_reaction_id]['found'] = any([measured_reactions_match[measured_reaction_id]['species_match']['found'], measured_reactions_match[measured_reaction_id]['ec_match']['found']])
    #TODO
    ############ compare by length of pathway ######
    #TODO: make sure that 
    # Reaction level stuff
    reactions_id = list(measured_reactions_match.keys())
    conf_score = np.mean([measured_reactions_match[i]['species_match']['score'] for i in reactions_id]+[measured_reactions_match[i]['ec_match']['score'] for i in reactions_id])
    #measured_reactions_match['found'] = any([measured_reactions_match[i]['species_match']['found'] for i in reactions_id]+[measured_reactions_match[i]['ec_match']['found'] for i in reactions_id])
    isFound = all([measured_reactions_match[i]['found'] for i in measured_reactions_match])
    measured_reactions_match['confidence_score'] = conf_score
    measured_reactions_match['found'] = isFound
    return measured_reactions_match['found'], measured_reactions_match['confidence_score'], measured_reactions_match
