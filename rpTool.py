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
                meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.4
                meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
            """ Use only the InChI key as a measure
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
            """
            #find according to the inchikey -- allow partial matches
            if measured_brsynth_annot['inchikey'] and sim_rpsbml_brsynth_annot['inchikey']:
                measured_inchikey_split = measured_brsynth_annot['inchikey'].split('-')
                sim_rpsbml_inchikey_split = sim_rpsbml_brsynth_annot['inchikey'].split('-')
                if measured_inchikey_split[0]==sim_rpsbml_inchikey_split[0]:
                    meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                    if measured_inchikey_split[1]==sim_rpsbml_inchikey_split[1]:
                        meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                        meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
                        if measured_inchikey_split[2]==sim_rpsbml_inchikey_split[2]:
                            meas_species_match[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                            meas_species_match[measured_species.getId()][sim_species.getId()]['found'] = True
            if species_match[measured_species.getId()]['score']<meas_species_match[measured_species.getId()][sim_species.getId()]['score'] and meas_species_match[measured_species.getId()][sim_species.getId()]['found']:
                species_match[measured_species.getId()]['id'] = sim_species.getId()
                species_match[measured_species.getId()]['score'] = round(meas_species_match[measured_species.getId()][sim_species.getId()]['score'], 4)
                species_match[measured_species.getId()]['found'] = meas_species_match[measured_species.getId()][sim_species.getId()]['found']
    return species_match


##
# Compare that all the measured species of a reactions are found within sim species to match with a reaction.
# We assume that there cannot be two reactions that have the same species and reactants. This is maintained by SBML
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
                    'std': 0.0,
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
        measured_reactions_match[measured_reaction.getId()]['score'] = np.std(reactants_score+products_score)
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
            if set(meas_sim_reactants).issubset(set(sim_reactants)) and set(meas_sim_products).issubset(set(sim_products)):
                measured_reactions_match[measured_reaction_id]['reaction'] = sim_reaction_id
                measured_reactions_match[measured_reaction_id]['found'] = True
                continue
    return measured_reactions_match

## Match the ec numbers of the reactions in two models
#
# Returns True if at least one ec number in measured pathway exists in a reaction of 
#
def compareEC(measured_rpsbml, sim_rpsbml, pathway_id='rp_pathway'):
    ######## compare by EC number ############
    ec_reactions_match = {}
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(pathway_id):
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        measured_reaction_miriam = sim_rpsbml.readMIRIAMAnnotation(measured_reaction.getAnnotation())
        ec_reactions_match[measured_reaction_id] = {'score': 0.0, 'found': False, 'id': None}
        if 'ec-code' in measured_reaction_miriam:
            for sim_reaction_id in sim_rpsbml.readRPpathwayIDs(pathway_id):
                sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
                sim_reaction_miriam = sim_rpsbml.readMIRIAMAnnotation(sim_reaction.getAnnotation())
                if 'ec-code' in sim_reaction_miriam:
                    #we only need one match here
                    measured_ec = [i for i in measured_reaction_miriam['ec-code']]
                    sim_ec = [i for i in sim_reaction_miriam['ec-code']]
                    if any(i in measured_ec for i in sim_ec) and ec_reactions_match[measured_reaction_id]<1.0:
                        ec_reactions_match[measured_reaction_id] = {'score': 1.0, 'found': True, 'id': sim_reaction_id}
                        continue
                    measured_frac_ec = [i.split('.')[:-1] for i in measured_reaction_miriam['ec-code']]
                    measured_frac_ec = ['.'.join(i) for i in measured_frac_ec]
                    sim_frac_ec = [i.split('.')[:-1] for i in sim_reaction_miriam['ec-code']]
                    sim_frac_ec = ['.'.join(i) for i in sim_frac_ec]
                    if any(i in measured_frac_ec for i in sim_frac_ec) and ec_reactions_match[measured_reaction_id]<0.5:
                        ec_reactions_match[measured_reaction_id] = {'score': 0.5, 'found': True, 'id': sim_reaction_id}
                        continue
    return ec_reactions_match

## Compare a measured to sim rpSBML pathway 
#
# Works by trying to find a measured pathway contained within a simulated one. Does not perform perfect match! Since 
# simulated ones contain full cofactors while the measured ones have impefect information
def compareRPpathways(measured_rpsbml, sim_rpsbml, scrict_length=True, pathway_id='rp_pathway'):
    if scrict_length:
        if not len(measured_rpsbml.readRPpathwayIDs(pathway_id))==len(sim_rpsbml.readRPpathwayIDs(pathway_id)):
            return False, 0.0, 0.0, 0.0, 0.0
    species_match = compareSpecies(measured_rpsbml, sim_rpsbml)
    reactions_match = compareReactions(measured_rpsbml, sim_rpsbml, species_match, pathway_id)
    reactions_score = np.mean([reactions_match[i]['score'] for i in reactions_match])
    reactions_std = np.std([reactions_match[i]['score'] for i in reactions_match])
    is_reaction_match = all([reactions_match[i]['found'] for i in reactions_match])
    reactions_ec = compareEC(measured_rpsbml, sim_rpsbml, pathway_id)
    reactions_ec_score = np.mean([reactions_ec[i]['score'] for i in reactions_ec])
    reactions_ec_std = np.std([reactions_ec[i]['score'] for i in reactions_ec])
    is_reactions_ec_match = all([reactions_ec[i]['found'] for i in reactions_ec])
    if is_reaction_match or is_reactions_ec_match:
        return True, reactions_score, reactions_std, reactions_ec_score, reactions_ec_std, {'species_match': species_match, 'reactions_match': reactions_match, 'ec_match': reactions_ec}
    else:
        return False, reactions_score, reactions_std, reactions_ec_score, reactions_ec_std, {'species_match': species_match, 'reactions_match': reactions_match, 'ec_match': reactions_ec}

