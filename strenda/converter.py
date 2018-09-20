'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import re
import sys
import uuid
import xml.sax

from libsbml import BIOLOGICAL_QUALIFIER, BQB_IS, CVTerm, SBMLDocument, \
    writeSBMLToFile, writeSBMLToString


class StrendaHandler(xml.sax.ContentHandler):
    '''STRENDA-DB handler.'''

    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.__document = SBMLDocument()
        self.__model = self.__document.createModel()

        self.__start = False
        self.__element_name = None
        self.__species = None
        self.__reaction = None
        self.__spec_ref = None
        self.__parent = None

    def startElement(self, name, attrs):
        self.__start = True
        self.__element_name = name

        if name == 'experiment':
            self.__document.setId(attrs['strendaId'])
        elif name == 'protein':
            self.__species = self.__model.createSpecies()
            self.__species.setId(_get_id(attrs['uniprotKbAC']))

            _add_annotation(self.__species,
                            'http://identifiers.org/uniprot/' +
                            attrs['uniprotKbAC'])

            self.__parent = name
        elif name == 'dataset':
            self.__reaction = self.__model.createReaction()
            self.__reaction.setId(_get_id(uuid.uuid4()))
            self.__reaction.setName(attrs['name'])

        elif name == 'smallCompound':
            species_id = _get_id(attrs['refId'])

            self.__species = self.__model.createSpecies()
            self.__species.setId(species_id)

            if attrs['role'] == 'Substrate':
                self.__spec_ref = self.__reaction.createReactant()
                self.__spec_ref.setSpecies(species_id)
            elif attrs['role'] == 'Product':
                self.__spec_ref = self.__reaction.createProduct()
                self.__spec_ref.setSpecies(species_id)

            self.__parent = name

        elif name == 'value':
            if attrs['type'] == 'Concentration':
                self.__species.setInitialConcentration(float(attrs['value']))
                # attrs['unit']

    def endElement(self, name):
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

    def writeSBMLToFile(self, filename):
        '''Write SBML to file.'''
        writeSBMLToFile(self.__document, filename)

    def writeSBMLToString(self):
        '''Write SBML to string.'''
        return writeSBMLToString(self.__document)


def _get_id(id_in):
    '''Format id.'''
    return re.sub('\W+', '_', str(id_in))


def _add_annotation(obj, resource):
    '''Add an annotation.'''
    cv_term = CVTerm()
    cv_term.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv_term.setBiologicalQualifierType(BQB_IS)
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

    handler.writeSBMLToFile(out_filename)


def main(args):
    '''main method.'''
    convert(*args)


if __name__ == '__main__':
    main(sys.argv[1:])
