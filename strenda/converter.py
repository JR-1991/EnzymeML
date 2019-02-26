'''
EnzymeML (c) University of Manchester 2018

EnzymeML is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-branches
# pylint: disable=too-many-instance-attributes
import re
import sys
import uuid
import xml.sax

from libsbml import BIOLOGICAL_QUALIFIER, BQB_IS, CVTerm, SBMLDocument, \
    UNIT_KIND_DIMENSIONLESS, UNIT_KIND_LITRE, UNIT_KIND_MOLE, UNIT_KIND_SECOND, \
    writeSBMLToFile, writeSBMLToString
    


class StrendaHandler(xml.sax.ContentHandler):
    '''STRENDA-DB handler.'''

    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        
        self.__document = SBMLDocument()
        
        self.__model = self.__document.createModel()
        self.__model.setExtentUnits('mole')
        self.__model.setTimeUnits('second')
        
        self.__compartment = self.__model.createCompartment()
        self.__compartment.setId('c')
        self.__compartment.setConstant(True)
        self.__compartment.setSize(1)
        self.__compartment.setSpatialDimensions(3)
        self.__compartment.setUnits('litre')

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
            self.__species = self.__model.createSpecies()
            self.__species.setId(_get_id(attrs['uniprotKbAC']))
            self.__species.setSBOTerm('SBO:0000252')
            self.__species.setCompartment(self.__compartment.getId())
            self.__species.setHasOnlySubstanceUnits(True)
            self.__species.setConstant(True)
            self.__species.setBoundaryCondition(False)

            _add_annotation(self.__species,
                            'http://identifiers.org/uniprot/' +
                            attrs['uniprotKbAC'])

        elif name == 'dataset':
            self.__reaction = self.__model.createReaction()
            self.__reaction.setId(_get_id(uuid.uuid4()))
            self.__reaction.setName(attrs['name'])
            self.__reaction.setSBOTerm('SBO:0000176')
            self.__reaction.setReversible(False)
            self.__reaction.setFast(False)
            self.__kinetic_law = self.__reaction.createKineticLaw()

        elif name == 'smallCompound':
            species_id = _get_id(attrs['refId'])

            self.__species = self.__model.createSpecies()
            self.__species.setId(species_id)
            self.__species.setSBOTerm('SBO:0000247')
            self.__species.setCompartment(self.__compartment.getId())
            self.__species.setHasOnlySubstanceUnits(True)

            if attrs['role'] == 'Substrate':
                self.__spec_ref = self.__reaction.createReactant()
                self.__spec_ref.setSpecies(species_id)
                self.__spec_ref.setConstant(False)
                self.__species.setConstant(False)
                self.__species.setBoundaryCondition(False)
            elif attrs['role'] == 'Product':
                self.__spec_ref = self.__reaction.createProduct()
                self.__spec_ref.setSpecies(species_id)
                self.__spec_ref.setConstant(False)
                self.__species.setConstant(False)
                self.__species.setBoundaryCondition(False)
            else:
                self.__species.setConstant(True)
                self.__species.setBoundaryCondition(True)

            self.__parent = name

        elif name == 'macromolecule':
            species_id = _get_id(attrs['refId'])

            self.__species = self.__model.createSpecies()
            self.__species.setId(species_id)
            self.__species.setCompartment(self.__compartment.getId())
            self.__species.setConstant(True)
            self.__species.setBoundaryCondition(True)
            self.__species.setHasOnlySubstanceUnits(True)

            if attrs['moleculeClass'] == 'Protein':
                self.__species.setSBOTerm('SBO:0000252')

            self.__parent = name

        elif name == 'kineticParameter':
            self.__parent = name

        elif name == 'value':
            if attrs['type'] == 'Concentration':
                conc, units = self.__get_value_units(float(attrs['value']),
                                                     attrs['unit'])
                self.__species.setInitialConcentration(conc)
                self.__species.setUnits(units)
            elif attrs['type'] == 'ConcentrationRange':
                start_conc, start_units = self.__get_value_units(float(attrs['startValue']),
                                                     attrs['unit'])
                end_conc, end_units = self.__get_value_units(float(attrs['endValue']),
                                                     attrs['unit'])
                self.__species.setInitialConcentration(start_conc)
                self.__species.setUnits(start_units)
                self.__species.setConstant(False)
            elif self.__parent == 'kineticParameter':
                self.__add_parameter(attrs)

    def endElement(self, _):
        self.__start = False

        # if name == 'dataset':
        #    self.__reaction = None

    def characters(self, content):
        if self.__start and content.strip() and content != 'null':
            if self.__element_name == 'name' and self.__species:
                self.__species.setName(content)
            elif self.__element_name == 'cid':
                _add_annotation(self.__species,
                                'http://identifiers.org/pubchem.compound/' +
                                content)
            elif self.__element_name == 'chebiId':
                _add_annotation(self.__species,
                                'http://identifiers.org/chebi/CHEBI:' +
                                content)
            elif self.__element_name == 'inchi':
                _add_annotation(self.__species,
                                'http://identifiers.org/inchi/' +
                                content)
            elif self.__element_name == 'stoichiometry':
                self.__spec_ref.setStoichiometry(float(content))

    def write_sbml_to_file(self, filename):
        '''Write SBML to file.'''
        writeSBMLToFile(self.__document, filename)

    def write_sbml_to_string(self):
        '''Write SBML to string.'''
        return writeSBMLToString(self.__document)

    def __add_parameter(self, attrs):
        '''Add parameter.'''
        value, units = self.__get_value_units(float(attrs['value']),
                                              attrs['unit'])
        parameter = self.__kinetic_law.createLocalParameter()
        parameter.setValue(float(value))
        parameter.setUnits(units)
        parameter.setId(_get_id(uuid.uuid4()))
        parameter.setName(attrs['name'])
        
        #TODO: Add SBO Terms

    def __get_value_units(self, value, units):
        '''Get value and units.'''
        if units == 'mM':
            return value / 10 ** 3, 'mole'

        if units == 'microM':
            return value / 10 ** 6, 'mole'

        if units == 'nM':
            return value / 10 ** 9, 'mole'

        if units == 's-1':
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

            return value, unit_def_id

        if units == 'M-1S-1':
            unit_def_id = 'M_1S_1'

            if not self.__model.getUnitDefinition(unit_def_id):
                unit_def = self.__model.createUnitDefinition()
                unit_def.setId(unit_def_id)
                unit_def.setName(unit_def.getName())
                unit = unit_def.createUnit()
                unit.setScale(1)
                unit.setMultiplier(1)
                unit.setExponent(-1)
                unit.setKind(UNIT_KIND_MOLE)
                unit = unit_def.createUnit()
                unit.setScale(1)
                unit.setMultiplier(1)
                unit.setExponent(-1)
                unit.setKind(UNIT_KIND_SECOND)

            return value, unit_def_id
        
            if units == 'units-ml1':
                unit_def_id = 'units_ml1'
    
                if not self.__model.getUnitDefinition(unit_def_id):
                    unit_def = self.__model.createUnitDefinition()
                    unit_def.setId(unit_def_id)
                    unit_def.setName(unit_def.getName())
                    unit = unit_def.createUnit()
                    unit.setScale(1)
                    unit.setMultiplier(1)
                    unit.setExponent(-1)
                    unit.setKind(UNIT_KIND_DIMENSIONLESS)
                    unit = unit_def.createUnit()
                    unit.setScale(-3)
                    unit.setMultiplier(1)
                    unit.setExponent(-1)
                    unit.setKind(UNIT_KIND_LITRE)

            return value, unit_def_id
        return value, units


def _get_id(id_in):
    '''Format id.'''
    return '_' + re.sub(r'\W+', '_', str(id_in))


def _add_annotation(obj, resource, qualifier_type=BIOLOGICAL_QUALIFIER,
                    qualifier_sub_type=BQB_IS):
    '''Add an annotation.'''
    cv_term = CVTerm()
    cv_term.setQualifierType(qualifier_type)

    if qualifier_type is BIOLOGICAL_QUALIFIER:
        cv_term.setBiologicalQualifierType(qualifier_sub_type)

    cv_term.addResource(resource)

    obj.setMetaId('_meta' + obj.getId())
    obj.addCVTerm(cv_term)


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
