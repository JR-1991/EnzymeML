'''
EnzymeML (c) University of Manchester 2018

EnzymeML is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-arguments
import re
import uuid

from libsbml import SBMLDocument, BIOLOGICAL_QUALIFIER, BQB_IS, CVTerm


def get_document():
    '''Get SBMLDocument.'''
    document = SBMLDocument()

    # Create model with default units:
    model = document.createModel()
    model.setExtentUnits('mole')
    model.setTimeUnits('second')

    # Create default compartment of size 1 litre:
    compartment = model.createCompartment()
    compartment.setId('c')
    compartment.setConstant(True)
    compartment.setSize(1)
    compartment.setSpatialDimensions(3)
    compartment.setUnits('litre')

    return document, model, compartment


def add_reaction(model, name, reversible=True):
    '''Add reaction.'''
    reaction = model.createReaction()
    reaction.setId(get_id(uuid.uuid4()))
    reaction.setName(name)
    reaction.setSBOTerm('SBO:0000176')
    reaction.setReversible(reversible)
    reaction.setFast(False)
    return reaction


def add_substrate(model, reaction, species_id, comp_id, name=None):
    '''Add substrate.'''
    return _add_reaction_participant(model, reaction, species_id, name,
                                     comp_id, True)


def add_product(model, reaction, species_id, comp_id, name):
    '''Add product.'''
    return _add_reaction_participant(model, reaction, species_id, name,
                                     comp_id, False)


def add_enzyme(model, species_id, comp_id, name=None, uniprot_id=None):
    '''Add enzyme.'''
    species = model.createSpecies()
    species.setId(species_id)
    species.setSBOTerm('SBO:0000252')
    species.setCompartment(comp_id)
    species.setConstant(True)
    species.setBoundaryCondition(False)
    species.setHasOnlySubstanceUnits(True)

    if name:
        species.setName(name)

    if uniprot_id:
        add_annotation(species, 'http://identifiers.org/uniprot/' + uniprot_id)

    return species


def add_non_participant(model, species_id, comp_id, sbo_term=0):
    '''Add non-participating species.'''
    species = model.createSpecies()
    species.setId(species_id)
    species.setCompartment(comp_id)
    species.setConstant(True)
    species.setBoundaryCondition(True)
    species.setHasOnlySubstanceUnits(True)

    if sbo_term:
        species.setSBOTerm(sbo_term)

    return species


def add_parameter(kinetic_law, value, units, name, sbo_term=0):
    '''Add parameter.'''
    parameter = kinetic_law.createLocalParameter()
    parameter.setValue(float(value))
    parameter.setUnits(units)
    parameter.setId(get_id(uuid.uuid4()))
    parameter.setName(name)

    if sbo_term:
        parameter.setSBOTerm(sbo_term)


def get_id(id_in):
    '''Format id.'''
    return '_' + re.sub(r'\W+', '_', str(id_in))


def add_annotation(obj, resource, qualifier_type=BIOLOGICAL_QUALIFIER,
                   qualifier_sub_type=BQB_IS):
    '''Add an annotation.'''
    cv_term = CVTerm()
    cv_term.setQualifierType(qualifier_type)

    if qualifier_type is BIOLOGICAL_QUALIFIER:
        cv_term.setBiologicalQualifierType(qualifier_sub_type)

    cv_term.addResource(resource)

    obj.setMetaId('_meta' + obj.getId())
    obj.addCVTerm(cv_term)


def _add_reaction_participant(model, reaction, species_id, name, comp_id,
                              is_substrate):
    '''Add reaction participant.'''
    species = model.createSpecies()
    species.setId(species_id)
    species.setSBOTerm('SBO:0000247')
    species.setCompartment(comp_id)
    species.setHasOnlySubstanceUnits(True)
    species.setConstant(False)
    species.setBoundaryCondition(False)

    if name:
        species.setName(name)

    spec_ref = reaction.createReactant() if is_substrate \
        else reaction.createProduct()

    spec_ref.setSpecies(species_id)
    spec_ref.setConstant(False)

    return species, spec_ref
