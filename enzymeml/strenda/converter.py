'''
EnzymeML (c) University of Manchester 2018

EnzymeML is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-instance-attributes
import sys
import xml.sax

from libsbml import UNIT_KIND_MOLE, UNIT_KIND_SECOND, \
    writeSBMLToFile, writeSBMLToString

from enzymeml import utils


class StrendaHandler(xml.sax.ContentHandler):
    '''STRENDA-DB handler.'''

    def __init__(self):
        xml.sax.ContentHandler.__init__(self)

        self.__document, self.__model, self.__compartment = \
            utils.get_document()

        self.__start = False
        self.__element_name = None
        self.__species = None
        self.__reaction = None
        self.__spec_ref = None
        self.__kinetic_law = None
        self.__parent = None

    def startElement(self, name, attrs):
        self.__start = True
        self.__element_name = name

        if name == 'experiment':
            self.__document.setId(attrs['strendaId'])

        elif name == 'assayConditions':
            self.__parent = name

        elif name == 'protein':
            self.__add_protein(attrs)

        elif name == 'dataset':
            self.__add_reaction(attrs)

        elif name == 'smallCompound':
            self.__add_small_compound(name, attrs)

        elif name == 'macromolecule':
            self.__add_macromolecule(name, attrs)

        elif name == 'kineticParameter':
            self.__parent = name

        elif name == 'value':
            self.__add_value(attrs)

    def endElement(self, _):
        self.__start = False

        # if name == 'dataset':
        #    self.__reaction = None

    def characters(self, content):
        if self.__start and content.strip() and content != 'null':
            if self.__element_name == 'name' and self.__species:
                self.__species.setName(content)
            elif self.__element_name == 'cid':
                utils.add_annotation(
                    self.__species,
                    'http://identifiers.org/pubchem.compound/' + content)
            elif self.__element_name == 'chebiId':
                utils.add_annotation(self.__species,
                                     'http://identifiers.org/chebi/CHEBI:' +
                                     content)
            elif self.__element_name == 'inchi':
                utils.add_annotation(self.__species,
                                     'http://identifiers.org/inchi/' + content)
            elif self.__element_name == 'stoichiometry':
                self.__spec_ref.setStoichiometry(float(content))

    def write_sbml_to_file(self, filename):
        '''Write SBML to file.'''
        writeSBMLToFile(self.__document, filename)

    def write_sbml_to_string(self):
        '''Write SBML to string.'''
        return writeSBMLToString(self.__document)

    def __add_small_compound(self, name, attrs):
        '''Add small compound.'''
        species_id = utils.get_id(attrs['refId'])

        if attrs['role'] == 'Substrate':
            self.__species, self.__spec_ref = \
                utils.add_substrate(self.__model, self.__reaction,
                                    species_id,
                                    self.__compartment.getId(),
                                    name)
        elif attrs['role'] == 'Product':
            self.__species, self.__spec_ref = \
                utils.add_product(self.__model, self.__reaction,
                                  species_id,
                                  self.__compartment.getId(),
                                  name)
        else:
            self.__species = self.__model.createSpecies()
            self.__species.setId(species_id)
            self.__species.setCompartment(self.__compartment.getId())
            self.__species.setConstant(True)
            self.__species.setBoundaryCondition(True)
            self.__species.setHasOnlySubstanceUnits(True)
            self.__species.setSBOTerm('SBO:0000247')

        self.__parent = name

    def __add_protein(self, attrs):
        '''Add protein.'''
        self.__species = utils.add_enzyme(self.__model,
                                          utils.get_id(attrs['uniprotKbAC']),
                                          self.__compartment.getId(),
                                          uniprot_id=attrs['uniprotKbAC'])

    def __add_macromolecule(self, name, attrs):
        '''Add macromolecule.'''
        species_id = utils.get_id(attrs['refId'])

        self.__species = self.__model.createSpecies()
        self.__species.setId(species_id)
        self.__species.setCompartment(self.__compartment.getId())
        self.__species.setConstant(True)
        self.__species.setBoundaryCondition(True)
        self.__species.setHasOnlySubstanceUnits(True)

        if attrs['moleculeClass'] == 'Protein':
            self.__species.setSBOTerm('SBO:0000252')

        self.__parent = name

    def __add_reaction(self, attrs):
        '''Add reaction.'''
        self.__reaction = utils.add_reaction(self.__model, attrs['name'])
        self.__kinetic_law = self.__reaction.createKineticLaw()

    def __add_value(self, attrs):
        '''Add value.'''
        if attrs['type'] == 'Concentration':
            conc, units = self.__get_value_units(float(attrs['value']),
                                                 attrs['unit'])
            self.__species.setInitialConcentration(conc)
            self.__species.setUnits(units)
        elif attrs['type'] == 'ConcentrationRange':
            start_conc, start_units = \
                self.__get_value_units(float(attrs['startValue']),
                                       attrs['unit'])
            end_conc, end_units = \
                self.__get_value_units(float(attrs['endValue']),
                                       attrs['unit'])

            self.__species.setInitialConcentration(start_conc)
            self.__species.setUnits(start_units)
            self.__species.setConstant(False)
        elif self.__parent == 'kineticParameter':
            self.__add_parameter(attrs)

    def __add_parameter(self, attrs):
        '''Add parameter.'''
        value, units = self.__get_value_units(float(attrs['value']),
                                              attrs['unit'])

        if attrs['name'] == 'kcat':
            sbo_term = 25
        elif attrs['name'] == 'km':
            sbo_term = 373
        else:
            sbo_term = 0

        utils.add_parameter(self.__kinetic_law, value, units, attrs['name'],
                            sbo_term)

    def __get_value_units(self, value, units):
        '''Get value and units.'''
        if units == 'mM':
            value = value / 10 ** 3
            units = 'mole'

        elif units == 'microM':
            value = value / 10 ** 6
            units = 'mole'

        elif units == 'nM':
            value = value / 10 ** 9
            units = 'mole'

        elif units == 'units-ml1':
            value = value * 10 ** 3
            units = 'item'

        elif units == 's-1':
            unit_def_id = 's_1'

            if not self.__model.getUnitDefinition(unit_def_id):
                unit_def = self.__model.createUnitDefinition()
                unit_def.setId(unit_def_id)
                unit_def.setName(unit_def.getName())
                unit = unit_def.createUnit()
                unit.setScale(1)
                unit.setMultiplier(1)
                unit.setExponent(-1)
                unit.setKind(UNIT_KIND_SECOND)

            units = unit_def_id

        elif units == 'M-1S-1':
            unit_def_id = 'M_1S_1'

            if not self.__model.getUnitDefinition(unit_def_id):
                unit_def = self.__model.createUnitDefinition()
                unit_def.setId(unit_def_id)
                unit_def.setName(unit_def.getName())
                unit = unit_def.createUnit()
                unit.setScale(1)
                unit.setMultiplier(1)
                unit.setExponent(1)
                unit.setKind(UNIT_KIND_MOLE)
                unit = unit_def.createUnit()
                unit.setScale(1)
                unit.setMultiplier(1)
                unit.setExponent(-1)
                unit.setKind(UNIT_KIND_SECOND)

            units = unit_def_id

        return value, units


def convert(in_filename, out_filename='strenda_sbml.xml'):
    '''Convert file.'''
    parser = xml.sax.make_parser()

    handler = StrendaHandler()
    parser.setContentHandler(handler)

    with open(in_filename, 'r') as fle:
        parser.parse(fle)

    handler.write_sbml_to_file(out_filename)


def main(args):
    '''main method.'''
    convert(*args)


if __name__ == '__main__':
    main(sys.argv[1:])
