'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import re
import sys

import h5py

import numpy as np
import pandas as pd


def convert(in_filename, out_filename='out.hdf5'):
    '''Convert Pandas to hdf5.'''
    df = pd.read_csv(in_filename, sep='\t', dtype=np.float64)
    df.name = 'timeseries'

    # Remove hashtag:
    df.columns = [re.sub(r'#\s+', '', col) for col in df.columns]

    # Set index to be time:
    df.set_index('Time', inplace=True)

    # Remove empty columns:
    df.dropna(axis=1, how='all', inplace=True)

    # Write to hdf5:
    df.to_hdf(out_filename, key=df.name, mode='w')

    return out_filename


def read_hdf(filename):
    '''Read hdf5 file.'''
    fle = h5py.File(filename, 'r')

    for key, value in fle.items():
        print(key, value, value.attrs.keys())

        for key, val in value.items():
            print(key, val, val.attrs.keys())


def main(args):
    '''main method.'''
    out_filename = convert(*args)
    read_hdf(out_filename)

    # Write to csv:
    pd.read_hdf(out_filename).to_csv('out.csv')


if __name__ == '__main__':
    main(sys.argv[1:])
