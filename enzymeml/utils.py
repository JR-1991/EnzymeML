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


def add_substrate(model, reaction, species_id, comp_id,
                  name=None, stoichiometry=1):
    '''Add substrate.'''
    return _add_substrate_product(model, reaction, species_id, name,
                                  comp_id, stoichiometry, True)


def add_product(model, reaction, species_id, comp_id,
                name=None, stoichiometry=1):
    '''Add product.'''
    return _add_substrate_product(model, reaction, species_id, name,
                                  comp_id, stoichiometry, False)


def add_enzyme(model, reaction, species_id, comp_id,
               name=None, uniprot_id=None):
    '''Add enzyme.'''
    species = _add_species(model, species_id, name, comp_id, 252,
                           constant=True, boundary_condition=False)

    if uniprot_id:
        add_annotation(species, 'http://identifiers.org/uniprot/' + uniprot_id)

    spec_ref = reaction.createModifier()
    spec_ref.setSpecies(species_id)
    spec_ref.setSBOTerm(460)

    return species


def add_non_participant(model, species_id, comp_id,
                        name=None, sbo_term=0):
    '''Add non-participating species.'''
    return _add_species(model, species_id, name, comp_id, sbo_term,
                        constant=True, boundary_condition=True)


def add_parameter(kinetic_law, value, units, name, sbo_term=0):
    '''Add parameter.'''
    parameter = kinetic_law.createLocalParameter()
    parameter.setValue(float(value))
    parameter.setUnits(units)
    parameter.setId(get_id(uuid.uuid4()))
    parameter.setName(name)

    if sbo_term:
        parameter.setSBOTerm(sbo_term)


def set_notes(elem, notes):
    '''Set notes.'''
    elem.setNotes('<body xmlns=\'http://www.w3.org/1999/xhtml\'>' +
                  '<pre>' + notes + '</pre></body>')


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


def _add_substrate_product(model, reaction, species_id, name, comp_id,
                           stoichiometry, is_substrate):
    '''Add reaction participant.'''
    species = _add_species(model, species_id, name, comp_id, 247,
                           constant=False, boundary_condition=False)

    spec_ref = reaction.createReactant() if is_substrate \
        else reaction.createProduct()

    spec_ref.setSpecies(species_id)
    spec_ref.setStoichiometry(stoichiometry)
    spec_ref.setConstant(False)
    spec_ref.setSBOTerm(10 if is_substrate else 11)

    return species, spec_ref


def _add_species(model, species_id, name, comp_id, sbo_term,
                 constant, boundary_condition):
    '''Add species.'''
    species = model.createSpecies()
    species.setId(species_id)
    species.setCompartment(comp_id)
    species.setHasOnlySubstanceUnits(True)
    species.setConstant(constant)
    species.setBoundaryCondition(boundary_condition)

    if name:
        species.setName(name)

    if sbo_term:
        species.setSBOTerm(sbo_term)

    return species
