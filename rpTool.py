import rpSBML
import numpy as np
import tempfile


def compareRPpathways_rewrite(sim_rpsbml, measured_rpsbml):
    ############## compare using the species ###################
    print('Species')
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
    print(meas_species_match)
    print('#############################')
    print('Simulated')
    #TODO: need to deal if there are more than one match
    ############## compare the reactions #######################
    sim_reactions = {}
    for sim_reaction_id in sim_rpsbml.readRPpathwayIDs():
        sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
        sim_reactions[sim_reaction_id] = {'reactants': [], 'products': []}
        for sim_reactant in sim_reaction.getListOfReactants():
            sim_reactions[sim_reaction_id]['reactants'].append(sim_reactant.species)
        for sim_product in sim_reaction.getListOfProducts():
            sim_reactions[sim_reaction_id]['products'].append(sim_product.species)
    print(sim_reactions)
    print('#############################')
    print('Measured')
    measured_reactions_match = {}
    sum_score = []
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(): 
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        ################ construct the dict transforming the species #######
        measured_reactions_match[measured_reaction.getId()] = {
                'species_match': {
                    'reactants': [],
                    'reactants_score': [],
                    'products': [], 
                    'products_score': [],
                    'score': 0.0,
                    'sim_reaction': None},
                'ec_match': {
                    'measured': None, 
                    'sim': None, 
                    'score': 0.0},
                'confidence_score': 0.0
                }
        for reactant in measured_reaction.getListOfReactants():
            try:
                measured_reactions_match[measured_reaction.getId()]['species_match']['reactants'].append(
                        meas_species_match[reactant.species][0])
                measured_reactions_match[measured_reaction.getId()]['species_match']['reactants_score'].append(
                        meas_species_match[reactant.species][1])
                sum_score.append(meas_species_match[reactant.species][1])
            except (KeyError, IndexError) as e:
                measured_reactions_match[measured_reaction.getId()]['species_match']['reactants'].append(reactant.species)
                measured_reactions_match[measured_reaction.getId()]['species_match']['reactants_score'].append(0.0)
                sum_score.append(0.0)
        for product in measured_reaction.getListOfProducts():
            try:
                measured_reactions_match[measured_reaction.getId()]['species_match']['products'].append(
                        meas_species_match[product.species][0])
                measured_reactions_match[measured_reaction.getId()]['species_match']['products_score'].append(
                        meas_species_match[product.species][1])
                sum_score.append(meas_species_match[product.species][1])
            except (KeyError, IndexError) as e:
                measured_reactions_match[measured_reaction.getId()]['species_match']['products'].append(product.species)
                measured_reactions_match[measured_reaction.getId()]['species_match']['products_score'].append(0.0)
                sum_score.append(0.0)
    for measured_reaction_id in measured_reactions_match:
        for sim_reaction_id in sim_reactions:
            if set(measured_reactions_match[measured_reaction_id]['species_match']['reactants']).issubset(sim_reactions[sim_reaction_id]['species_match']['reactants']) and set(measured_reactions_match[measured_reaction_id]['species_match']['products']).issubset(sim_reactions[sim_reaction_id]['species_match']['products']):
                measured_reactions_match[measured_reaction_id]['species_match']['sim_reaction'] = sim_reaction_id

                break
    measured_reactions_match[measured_reaction.getId()]['species_match']['score'] = np.mean(sum_score)
    print(measured_reactions_match)
    print('#############################')
    print('EC')
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs():
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        measured_reaction_miriam = rpsbml.readMIRIAMAnnotation(measured_reaction.getAnnotation())
        try:
            measured_reactions_match[measured_reaction_id]['ec_match']['measured'] = measured_reaction_miriam['ec-code']
            for sim_reaction_id in sim_rpsbml.readRPpathwayIDs():
                sim_reaction = sim_rpsbml.model.getReaction(sim_reaction_id)
                sim_reaction_miriam = rpsbml.readMIRIAMAnnotation(sim_reaction.getAnnotation())
                try:
                    measured_reactions_match[measured_reaction_id]['ec_match']['sim'] = sim_reaction_miriam['ec-code']
                ##### try to compare them ####
                if set(measured_reaction_miriam['ec-code']).issubset(set(sim_reaction_miriam['ec-code'])):
                    #recover only the ones that are in the measured
                    found_meas_rp[rp_step_id]['rp_step_id'] = meas_step_id
                    found_meas_rp[rp_step_id]['found'] = True
                    found_meas_rp['ec_score'] += 1.0
                elif [i for i, j in zip(source_dict['ec-code'], target_dict['ec-code']) if i.split('.')[:-1]==j.split('.')[:-1]]:
                    found_meas_rp[rp_step_id]['rp_step_id'] = meas_step_id
                    found_meas_rp[rp_step_id]['found'] = True
                    found_meas_rp['ec_score'] += 0.5
                except KeyError:
                    measured_reactions_match[measured_reaction_id]['ec_match']['sim'] = None

        except KeyError:
            measured_reactions_match[measured_reaction_id]['ec_match']['measured'] = None
            



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
            

    

    '''
    for measured_reaction_id in measured_reactions_match:
        meas_reactants = [measured_reactions_match[measured_reaction_id]['reactants']['id']]
        for sim_reaction_id in sim_reactions:
            measured_reactions_match[measured_reaction_id]
            


        for measured_reactant in measured_reaction.getListOfReactants():
            try:
                measured_reactions_match[measured_reaction.getId()]['reactants'][measured_reactant.species]['species'] = list(meas_species_match[measured_reactant.species].keys())[0]
                measured_reactions_match[measured_reaction.getId()]['reactants'][measured_reactant.species]['score'] = meas_species_match[measured_reactant.species][list(meas_species_match[measured_reactant.species].keys())[0]]
            except (KeyError, IndexError) as e:
                pass
        for measured_product in measured_reaction.getListOfProducts():
            try:
                measured_reactions_match[measured_reaction.getId()]['products'][measured_product.species]['species'] = list(meas_species_match[measured_product.species].keys())[0]
                measured_reactions_match[measured_reaction.getId()]['products'][measured_product.species]['score'] = meas_species_match[measured_product.species][list(meas_species_match[measured_product.species].keys())[0]]
            except (KeyError, IndexError) as e:
                pass
        #calculate 
        measured_reactions_match[measured_reaction.getId()]
        if all([True if measured_reactions_match[measured_reaction.getId()][g]['species'] else False for g in measured_reactions_match[measured_reaction.getId()]['reactants']]) and all([True if measured_reactions_match[measured_reaction.getId()][g]['species'] else False for g in measured_reactions_match[measured_reaction.getId()]['reactants']]):
            measured_reactions_match[measured_reaction.getId()]['sim_reaction'] = 
    print(measured_reactions_match)
    '''



    '''
    ############## compare the reactions
    #1) list through all the measured species and compare with sim_rpsbml ones by:
    #       -> 
    measured_reactions_match = {}
    for measured_reaction_id in measured_rpsbml.readRPpathwayIDs(): 
        measured_reaction = measured_rpsbml.model.getReaction(measured_reaction_id)
        print('########### '+str(measured_reaction.getId())+' ##############')
        measured_reactions_match[measured_reaction.getId()] = {}
        for sim_rpsbml_reaction_id in sim_rpsbml.readRPpathwayIDs():
            sim_rpsbml_reaction = sim_rpsbml.model.getReaction(sim_rpsbml_reaction_id)
            print('\t########## '+str(sim_rpsbml_reaction.getId())+' ###########')
            measured_reactions_match[measured_reaction.getId()][sim_rpsbml_reaction.getId()] = {}
            species_match = {'reactants': {'species': {}, 'num': len(measured_reaction.getListOfReactants())}, 'products': {'species': {}, 'num': len(measured_reaction.getListOfProducts())}, 'score': None}
            #reactant
            for measured_reactant in measured_reaction.getListOfReactants():
                species_match['reactants']['species'][measured_reactant.species] = {}
                for sim_rpsbml_reactant in sim_rpsbml_reaction.getListOfReactants():
                    measured_species = measured_rpsbml.model.getSpecies(measured_reactant.species)
                    sim_rpsbml_species = sim_rpsbml.model.getSpecies(sim_rpsbml_reactant.species)
                    match_score = 0.0
                    measured_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(measured_species.getAnnotation())
                    sim_rpsbml_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(sim_rpsbml_species.getAnnotation())
                    #find according to xref
                    if sim_rpsbml.compareMIRIAMAnnotations(measured_species.getAnnotation(), sim_rpsbml_species.getAnnotation()):
                        match_score += 0.2
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
                            match_score += 0.2
                            if measured_inchikey_split[1]==sim_rpsbml_inchikey_split[1]:
                                match_score += 0.2
                                if measured_inchikey_split[2]==sim_rpsbml_inchikey_split[2]:
                                    match_score += 0.2
                    if match_score>0.0:
                        species_match['reactants']['species'][measured_reactant.species][sim_rpsbml_reactant.species] = match_score
            #product
            for measured_product in measured_reaction.getListOfProducts():
                species_match['products']['species'][measured_product.species] = {}
                for sim_rpsbml_product in sim_rpsbml_reaction.getListOfProducts():
                    measured_species = measured_rpsbml.model.getSpecies(measured_product.species)
                    sim_rpsbml_species = sim_rpsbml.model.getSpecies(sim_rpsbml_product.species)
                    match_score = 0.0
                    measured_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(measured_species.getAnnotation())
                    sim_rpsbml_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(sim_rpsbml_species.getAnnotation())
                    #find according to xref
                    if sim_rpsbml.compareMIRIAMAnnotations(measured_species.getAnnotation(), sim_rpsbml_species.getAnnotation()):
                        match_score += 0.2
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
                            match_score += 0.2
                            if measured_inchikey_split[1]==sim_rpsbml_inchikey_split[1]:
                                match_score += 0.2
                                if measured_inchikey_split[2]==sim_rpsbml_inchikey_split[2]:
                                    match_score += 0.2
                    if match_score>0.0:
                        species_match['products']['species'][measured_product.species][sim_rpsbml_product.species] = match_score
            print('\t'+str(species_match))
            #calculate mean of the scores

            #measured_reactions_match[measured_reaction.getId()][[sim_rpsbml_reaction.getId()] = species_match
        '''


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
        #1) list through all the measured species and compare with rpsbml ones by:
        #       -> 
        measured_reactions_match = {}
        for measured_reaction in measured_sbml.model.getListOfReactions(): 
            measured_reactions_match[measured_reaction.getId()] = {}
            for rpsbml_reaction in rpsbml.model.getListOfReactions():
                measured_reactions_match[measured_reaction.getId()][rpsbml_reaction.getId()] = {}
                species_match = {'reactants': {}, 'products': {}, 'score': None}
                #reactant
                for measured_reactant in measured_reaction.getListOfReactants():
                    species_match['reactants'][measured_reactant.species] = {}
                    for rpsbml_reactant in rpsbml_reaction.getListOfReactants():
                        measured_species = measured_reaction.getSpecies(measured_product.species)
                        rpsbml_species = rpsbml_reaction.getSpecies(rpsbml_product.species)
                        match_score = 0.0
                        measured_brsynth_annot = self.readBRSYNTHAnnotation(measured_species.getAnnotation())
                        rpsbml_brsynth_annot = self.readBRSYNTHAnnotation(rpsbml_species.getAnnotation())
                        #find according to xref
                        if rpsbml.compareMIRIAMAnnotations(measured_species.getAnnotation(), rpsbml_species.getAnnotation()):
                            match_score += 0.2
                        #find according to the SMILES
                        if measured_brsynth_annot['smiles'] and rpsbml_brsynth_annot['smiles']:
                            if measured_brsynth_annot['smiles']==rpsbml_brsynth_annot['smiles']:
                                match_score += 0.2
                        #find accorsing to the inchi
                        if measured_brsynth_annot['inchi'] and rpsbml_brsynth_annot['inchi']:
                            if measured_brsynth_annot['inchi']==rpsbml_brsynth_annot['inchi']:
                                match_score += 0.2
                        #find according to the inchikey -- allow partial matches
                        if measured_brsynth_annot['inchikey'] and rpsbml_brsynth_annot['inchikey']:
                            measured_inchikey_split = measured_brsynth_annot['inchikey'].split('-')
                            rpsbml_inchikey_split = rpsbml_brsynth_annot['inchikey'].split('-')
                            if measured_inchikey_split[0]==rpsbml_inchikey_split[0]:
                                match_score += 0.2
                                if measured_inchikey_split[1]==rpsbml_inchikey_split[1]:
                                    match_score += 0.2
                                    if measured_inchikey_split[2]==rpsbml_inchikey_split[2]:
                                        match_score += 0.2
                        if match_score>0.0:
                            species_match['reactants'][measured_reactant.species][rpsbml_reactant.species] = match_score
                #product
                for measured_product in measured_reaction.getListOfProducts():
                    species_match['products'][measured_product.species] = {}
                    for rpsbml_product in rpsbml_reaction.getListOfProducts():
                        measured_species = measured_reaction.getSpecies(measured_product.species)
                        rpsbml_species = rpsbml_reaction.getSpecies(rpsbml_product.species)
                        match_score = 0.0
                        measured_brsynth_annot = self.readBRSYNTHAnnotation(measured_species.getAnnotation())
                        rpsbml_brsynth_annot = self.readBRSYNTHAnnotation(rpsbml_species.getAnnotation())
                        #find according to xref
                        if rpsbml.compareMIRIAMAnnotations(measured_species.getAnnotation(), rpsbml_species.getAnnotation()):
                            match_score += 0.2
                        #find according to the SMILES
                        if measured_brsynth_annot['smiles'] and rpsbml_brsynth_annot['smiles']:
                            if measured_brsynth_annot['smiles']==rpsbml_brsynth_annot['smiles']:
                                match_score += 0.2
                        #find accorsing to the inchi
                        if measured_brsynth_annot['inchi'] and rpsbml_brsynth_annot['inchi']:
                            if measured_brsynth_annot['inchi']==rpsbml_brsynth_annot['inchi']:
                                match_score += 0.2
                        #find according to the inchikey -- allow partial matches
                        if measured_brsynth_annot['inchikey'] and rpsbml_brsynth_annot['inchikey']:
                            measured_inchikey_split = measured_brsynth_annot['inchikey'].split('-')
                            rpsbml_inchikey_split = rpsbml_brsynth_annot['inchikey'].split('-')
                            if measured_inchikey_split[0]==rpsbml_inchikey_split[0]:
                                match_score += 0.2
                                if measured_inchikey_split[1]==rpsbml_inchikey_split[1]:
                                    match_score += 0.2
                                    if measured_inchikey_split[2]==rpsbml_inchikey_split[2]:
                                        match_score += 0.2
                        if match_score>0.0:
                            species_match['products'][measured_product.species][rpsbml_product.species] = match_score
                #calculate mean of the scores

                #measured_reactions_match[measured_reaction.getId()][[rpsbml_reaction.getId()] = species_match



                #for measured_products in measured_reaction.getListOfProducts():
                #    for rpsbml_products in rpsbml_reaction.getListOfProducts():
                        



        measured_reaction_ids = [i.getId() for i in measured_sbml.model.getListOfReactions()] 
        rpsbml_reaction_ids = [i.getId() for i in rpsbml.model.getListOfReactions()] 
        for measured_reaction in measured_sbml.model.getListOfReactions(): 
            measured_reaction_speciesID = [] 
            for reactant in measured_reaction.getListOfReactants(): 
                measured_reaction_speciesID.append(reactant.species) 
            measured_reaction_productsID = [] 
            for product in measured_reaction.getListOfProducts(): 
                measured_reaction_productsID.append(product.species) 
            for rpsbml_reaction in rpsbml.model.getListOfReactions(): 
                #loop through target model reactions 
                rpsbml_reaction_speciesID = [i.species for i in rpsbml_reaction.getListOfReactants()] 
                rpsbml_reaction_productsID = [i.species for i in rpsbml_reaction.getListOfProducts()] 

                #perfect match
                #if not set(measured_reaction_speciesID)-set(rpsbml_reaction_speciesID) and not set(measured_reaction_productsID)-set(rpsbml_reaction_productsID): 

                #calculate 
                #abs(len(set(measured_reaction_speciesID)-set(rpsbml_reaction_speciesID))-len(measured_reaction_speciesID))


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

