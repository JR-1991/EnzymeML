'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import sys

import pandas as pd


def convert(filename):
    '''Convert Pandas to hdf5.'''
    df = pd.read_csv(filename, sep='\t')

    # Remove hashtag:
    df.columns = [col.replace('# ', ' ') for col in df.columns]

    # Remove empty columns:
    df.dropna(axis=1, how='all', inplace=True)

    # Write to csv:
    df.to_csv('out.csv', index=False)

    # Write to hdf5:
    df.to_hdf('out.h5', key='df', mode='w')


def main(args):
    '''main method.'''
    convert(*args)


if __name__ == '__main__':
    main(sys.argv[1:])
