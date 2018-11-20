'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=fixme
# pylint: disable=invalid-name
# pylint: disable=unused-import
import sys
import zipfile
import zlib


def pack(filenames, enzml_filename, compression=zipfile.ZIP_DEFLATED):
    '''Pack files to EnzymeML.'''
    zf = zipfile.ZipFile(enzml_filename, mode='w')

    # TODO: add manifest

    try:
        for filename in filenames:
            zf.write(filename, compress_type=compression)
    finally:
        zf.close()


def unpack(enzml_filename):
    '''Unpack files from EnzymeML.'''
    zf = zipfile.ZipFile(enzml_filename)

    for info in zf.infolist():
        print(zf.read(info.filename))


def main(args):
    '''main method.'''
    pack(args[:-1], args[-1])
    print(unpack(args[-1]))


if __name__ == '__main__':
    main(sys.argv[1:])
