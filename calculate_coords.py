'''
Script to calculate coords in a peptide backbone
using the nerf method, given a file containing
number of residues in the backbone and pairs
of phi psi angles

Rohit Bhattacharya
rohit.bhattachar@gmail.com
'''

# imports
from nerf import nerf
import numpy as np
import argparse


def parse_args():
    '''
    Parse arguments
    '''

    info = 'Calculate coords given torsions file'
    parser = argparse.ArgumentParser(description=info)

    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Torsions input file')

    args = parser.parse_args()
    return vars(args)


def place_residues(torsions):
    '''
    Do the actual placement of residues
    based on nerf method.
    '''

    # place the 0th carbonyl C at origin
    n_i_zero = np.array([0.0, 0.0, 0.0])

    # using CH1E-NH1 bond length = 1.458 ang
    # place C_alpha on pos x axis
    c_a_zero = np.array([1.458, 0.0, 0.0])

    # place C_carbonyl on xy plane using
    # CH1E-C bond length = 1.525 ang
    # C-CH1E-NH1 bond angle = 111.2 deg
    l = 1.525
    theta = 111.2 * np.pi/180  # in radians
    c_o_zero = c_a_zero + np.array([-l*np.cos(theta), l*np.sin(theta), 0])

    # start a listing of atoms starting with the zeroth
    # residue represented as a tuple of N, C_alpha, C_carbonyl
    atom_coords = [(n_i_zero, c_a_zero, c_o_zero)]

    # go through each set of torsion angles
    for i in range(1, len(torsions)):

        # calculate coords of all atoms of
        # the ith residue based on coords
        # from i-1th residue
        n_i, c_a, c_o = atom_coords[i-1]
        phi, psi, omega = torsions[i]
        n_i_plus_1 = nerf(n_i, c_a, c_o, l=1.329, theta=116.2, chi=psi)
        c_a_i_plus_1 = nerf(c_a, c_o, n_i_plus_1,
                            l=1.458, theta=121.7, chi=omega)
        c_o_i_plus_1 = nerf(c_o, n_i_plus_1, c_a_i_plus_1,
                            l=1.525, theta=111.2, chi=phi)
        atom_coords.append((n_i_plus_1, c_a_i_plus_1, c_o_i_plus_1))

    # return the calculated coords
    return atom_coords


def write_to_pdb(residue_coords):
    '''
    Write out calculated coords
    to PDB
    '''

    atom_names = ['N', 'CA', 'C']*(len(residue_coords))
    out_file = open('torsions.pdb', 'w')
    atom_coords = []
    for res_tuple in residue_coords:
        for atom_coord in res_tuple:
            atom_coords.append([round(c, 3) for c in atom_coord])

    res_number = 0
    for i, coord in enumerate(atom_coords):

        if i % 3 == 0:
            res_number += 1
        out_file.write('ATOM' +
                       ' '*2 + ' '*(5-len(str(i+1))) + ' ' +
                       str(i+1) + ' '*(4-len(atom_names[i])) +
                       atom_names[i] + ' ' + 'ALA' +
                       ' ' + 'A' + ' '*(4-len(str(res_number))) +
                       str(res_number) + ' '*4 +
                       ' '*(8-len(str(coord[0]))) + str(coord[0]) +
                       ' '*(8-len(str(coord[1]))) + str(coord[1]) +
                       ' '*(8-len(str(coord[2]))) + str(coord[2]) +
                       '\n')

    out_file.close()


def main():

    opts = parse_args()

    # make sure all the inputs
    # look okay
    try:
        in_file = open(opts['input'])
        n_residues = int(in_file.readline().strip())

    except IOError:
        print 'Input file', opts['input'], 'does not exist'

    except ValueError:
        print 'Num residues must be an integer'

    torsions = []
    for line in in_file:
        line = line.strip().split(',')

        # must be a pair of phi psi angles
        assert len(line) == 2, 'Must have a pair of comma separated phi psi angles'
        try:
            phi = float(line[0])
            psi = float(line[1])
            omega = 180
        except ValueError:
            print 'Angles must be valid float'
        torsions.append((phi, psi, omega))

    # finally make sure we have same number of phi psi
    # pairs as number of residues
    try:
        assert len(torsions) == n_residues, 'Must have same number of torsions as residues'

    except AssertionError as e:
        print e

    # actually calculate the coordinates
    calculated_coords = place_residues(torsions)
    write_to_pdb(calculated_coords)

if __name__ == '__main__':
    main()
