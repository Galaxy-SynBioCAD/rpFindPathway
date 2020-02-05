import numpy as np
import tempfile
import logging


##
#
#
def compareSpecies(measured_rpsbml, sim_rpsbml):
    ############## compare species ###################
    meas_species_match = {}
    species_match = {}
    for measured_species in measured_rpsbml.model.getListOfSpecies():
        meas_species_match[measured_species.getId()] = {}
        species_match[measured_species.getId()] = {'id': None, 'score': -1.0, 'found': False}
        for sim_species in sim_rpsbml.model.getListOfSpecies():
            meas_species_match[measured_species.getId()][sim_species.getId()] = {'score': 0.0, 'found': False}
            measured_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(measured_species.getAnnotation())
            sim_rpsbml_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(sim_species.getAnnotation())
            #find according to xref
            if sim_rpsbml.compareMIRIAMAnnotations(measured_species.getAnnotation(), sim_species.getAnnotation()):
                meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.3
                meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
            #find according to the SMILES
            if measured_brsynth_annot['smiles'] and sim_rpsbml_brsynth_annot['smiles']:
                if measured_brsynth_annot['smiles']==sim_rpsbml_brsynth_annot['smiles']:
                    meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                    meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
            #find accorsing to the inchi
            if measured_brsynth_annot['inchi'] and sim_rpsbml_brsynth_annot['inchi']:
                if measured_brsynth_annot['inchi']==sim_rpsbml_brsynth_annot['inchi']:
                    meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                    meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
            #find according to the inchikey -- allow partial matches
            if measured_brsynth_annot['inchikey'] and sim_rpsbml_brsynth_annot['inchikey']:
                measured_inchikey_split = measured_brsynth_annot['inchikey'].split('-')
                sim_rpsbml_inchikey_split = sim_rpsbml_brsynth_annot['inchikey'].split('-')
                if measured_inchikey_split[0]==sim_rpsbml_inchikey_split[0]:
                    meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.1
                    if measured_inchikey_split[1]==sim_rpsbml_inchikey_split[1]:
                        meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.1
                        meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
                        if measured_inchikey_split[2]==sim_rpsbml_inchikey_split[2]:
                            meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.1
                            meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
            if species_match[measured_species.getId()]['score']<meas_species_match[measured_species.getId()][sim_species.getId()]['score'] and meas_species_match[measured_species.getId()][sim_species.getId()]['found']:
                species_match[measured_species.getId()]['id'] = sim_species.getId()
                species_match[measured_species.getId()]['score'] = round(meas_species_match[measured_species.getId()][sim_species.getId()]['score'], 4)
                species_match[measured_species.getId()]['found'] = meas_species_match[measured_species.getId()][sim_species.getId()]['found']
    return species_match


##
# Compare that all the measured species of a reactions are found within sim species to match with a reaction
#
def compareReactions(measured_rpsbml, sim_rpsbml, species_match, pathway_id='rp_pathway'):
    #TODO: need to deal if there are more than one match
    ############## compare the reactions #######################
    #construct sim reactions with species
    sim_reactions = {}
    for sim_reaction_id in sim_rpsbml.readRPpathwayIDs(pathway_id):
        sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
        sim_reactions[sim_reaction_id] = {'reactants': [], 'products': []}
        for sim_reactant in sim_reaction.getListOfReactants():
            sim_reactions[sim_reaction_id]['reactants'].append(sim_reactant.species)
        for sim_product in sim_reaction.getListOfProducts():
            sim_reactions[sim_reaction_id]['products'].append(sim_product.species)
    #match the reactants and products conversion to sim species
    measured_reactions_match = {}
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        ################ construct the dict transforming the species #######
        measured_reactions_match[measured_reaction.getId()] = {
                    'reactants': {},
                    'reactants_score': 0.0,
                    'products': {},
                    'products_score': 0.0,
                    'score': 0.0,
                    'found': False,
                    'reaction': None}
        for reactant in measured_reaction.getListOfReactants():
            try:
                measured_reactions_match[measured_reaction.getId()]['reactants'][reactant.species] = species_match[reactant.species]
            except KeyError:
                measured_reactions_match[measured_reaction.getId()]['reactants'][reactant.species] = {'id': None, 'score': 0.0, 'found': False}
        reactants_score = [measured_reactions_match[measured_reaction.getId()]['reactants'][i]['score'] for i in measured_reactions_match[measured_reaction.getId()]['reactants']]
        measured_reactions_match[measured_reaction.getId()]['reactants_score'] = np.mean(reactants_score)
        #reactants_found = [measured_reactions_match[measured_reaction.getId()]['reactants'][i]['found'] for i in measured_reactions_match[measured_reaction.getId()]['reactants']]
        for product in measured_reaction.getListOfProducts():
            try:
                measured_reactions_match[measured_reaction.getId()]['products'][product.species] = species_match[product.species]
            except KeyError:
                measured_reactions_match[measured_reaction.getId()]['products'][product.species] = {'id': None, 'score': 0.0, 'found': False}
        products_score = [measured_reactions_match[measured_reaction.getId()]['products'][i]['score'] for i in measured_reactions_match[measured_reaction.getId()]['products']]
        measured_reactions_match[measured_reaction.getId()]['products_score'] = np.mean(products_score)
        #products_found = [measured_reactions_match[measured_reaction.getId()]['products'][i]['found'] for i in measured_reactions_match[measured_reaction.getId()]['products']]
        measured_reactions_match[measured_reaction.getId()]['score'] = np.mean(reactants_score+products_score)
    #Try to find the sik reaction
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        for sim_reaction_id in sim_rpsbml.readRPpathwayIDs(pathway_id):
            sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
            #sim reaction species
            sim_reactants = [i.species for i in sim_reaction.getListOfReactants()]
            sim_products = [i.species for i in sim_reaction.getListOfProducts()]
            #measured reaction pecies
            meas_sim_reactants = [measured_reactions_match[measured_reaction_id]['reactants'][i]['id'] for i in measured_reactions_match[measured_reaction_id]['reactants']]
            meas_sim_products = [measured_reactions_match[measured_reaction_id]['products'][i]['id'] for i in measured_reactions_match[measured_reaction_id]['products']]
            if not set(sim_reactants)-set(meas_sim_reactants) and not set(sim_products)-set(meas_sim_products):
                measured_reactions_match[measured_reaction_id]['reaction'] = sim_reaction_id
                measured_reactions_match[measured_reaction_id]['found'] = True
                continue
    return measured_reactions_match

##
#
#
def compareEC(measured_rpsbml, sim_rpsbml, pathway_id='rp_pathway'):
    ######## compare by EC number ############
    ec_reactions_match = {}
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        measured_reaction_miriam = sim_rpsbml.readMIRIAMAnnotation(measured_reaction.getAnnotation())
        ec_reactions_match[measured_reaction_id] = {'score': 0.0, 'found': False}
        if 'ec-code' in measured_reaction_miriam:
            for sim_reaction_id in sim_rpsbml.readRPpathwayIDs(pathway_id):
                sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
                sim_reaction_miriam = sim_rpsbml.readMIRIAMAnnotation(sim_reaction.getAnnotation())
                if 'ec-code' in sim_reaction_miriam:
                    #we only need one match here
                    measured_frac_ec = [i for i in measured_reaction_miriam['ec-code']]
                    sim_frac_ec = [i for i in sim_reaction_miriam['ec-code']]
                    if any(i in measured_frac_ec for i in sim_frac_ec):
                        ec_reactions_match[measured_reaction_id] = {'score': 1.0, 'found': True}
                        continue
                    measured_frac_ec = [i.split('.')[:-1] for i in measured_reaction_miriam['ec-code']]
                    measured_frac_ec = ['.'.join(i) for i in measured_frac_ec]
                    sim_frac_ec = [i.split('.')[:-1] for i in sim_reaction_miriam['ec-code']]
                    sim_frac_ec = ['.'.join(i) for i in sim_frac_ec]
                    if any(i in measured_frac_ec for i in sim_frac_ec):
                        ec_reactions_match[measured_reaction_id] = {'score': 0.5, 'found': True}
                        continue

## Compare a measured to sim rpSBML pathway 
#
# Works by trying to find a measured pathway contained within a simulated one. Does not perform perfect match! Since 
# simulated ones contain full cofactors while the measured ones have impefect information
def compareRPpathways(sim_rpsbml, measured_rpsbml, pathway_id='rp_pathway'):
    species_match = rpTool.compareSpecies(measured_rpsbml, sim_rpsbml)
    reactions_match = rpTool.compareSpecies(measured_rpsbml, sim_rpsbml, species_match, pathway_id)
    reactions_score = np.mean([reactions_match[i]['score'] for i in reactions_match])
    is_reaction_match = all([reactions_match[i]['found'] for i in reactions_match])
    reactions_ec = compareEC(measured_rpsbml, sim_rpsbml, pathway_id)
    reactions_ec_score = np.mean([reactions_ec[i]['score'] for i in reactions_ec])
    is_reactions_ec_match = all([reactions_ec[i]['found'] for i in reactions_ec])
    

    '''
    measured_reactions_match[measured_reaction.getId()]['confidence_score'] = np.mean([measured_reactions_match[measured_reaction.getId()]['species_match']['score']]+[measured_reactions_match[measured_reaction.getId()]['ec_match']['score']])
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        measured_reactions_match[measured_reaction_id]['found'] = any([measured_reactions_match[measured_reaction_id]['species_match']['found'], measured_reactions_match[measured_reaction_id]['ec_match']['found']])
    '''
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
