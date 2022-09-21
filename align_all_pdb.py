import sys, os, math
from os import path
from os import listdir
from os.path import isfile, join
import Bio.PDB


def main(name, argv):
    if (len(argv) <> 3):
        print_usage(name)
        return
    start_res = int(argv[1])
    end_res = int(argv[2])
    mypath = argv[0]
    onlyfiles = [a for a in listdir(mypath) if isfile(join(mypath, a))]
    for fi in onlyfiles:
        for fj in onlyfiles:
            print 'fi=', fi, 'fj=', fj
            Align_two(mypath, fi, fj, start_res, end_res)
    return ()


def Align_two(path, reference, sample, start_id, end_id):

    ## Select residues to align
    atoms_to_be_aligned = range(start_id, end_id + 1)

    ## Get structures
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    ref_structure = pdb_parser.get_structure("reference", path + reference)
    sample_structure = pdb_parser.get_structure("sample", path + sample)

    ## Use the first model in the pdb-files for alignment
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]

    ## Make a list of the atoms to align
    ref_atoms = []
    sample_atoms = []
    ref_res_id = []
    sample_res_id = []

    for ref_chain in ref_model:
        for ref_res in ref_chain:
            if ref_res.get_id()[1] in atoms_to_be_aligned:
                ref_res_id.append(int(ref_res.get_id()[1]))

    ## Same for sample structure
    for sample_chain in sample_model:
        for sample_res in sample_chain:
            if sample_res.get_id()[1] in atoms_to_be_aligned:
                sample_res_id.append(int(sample_res.get_id()[1]))

    ## Filter out elements that are not in both of lists: sample_atoms and ref_atoms
    i = 0
    j = 0
    print len(ref_res_id)
    print len(sample_res_id)
    while (i < len(ref_res_id)) or (j < len(sample_res_id)):
        if ref_res_id[i] == sample_res_id[j]:
            i += 1
            j += 1
        elif ref_res_id[i] > sample_res_id[j]:
            del sample_res_id[j]
        else:
            del ref_res_id[i]

    ## Go through all atoms to find CA
    for ref_chain in ref_model:
        for ref_res in ref_chain:
            if int(ref_res.get_id()[1]) in ref_res_id:
                ref_atoms.append(ref_res['CA'])

    ## Same for sample structure
    for sample_chain in sample_model:
        for sample_res in sample_chain:
            if int(sample_res.get_id()[1]) in sample_res_id:
                sample_atoms.append(sample_res['CA'])
                print sample_res.get_id()[1]
    print "___________"
    print len(ref_res_id)
    print len(sample_res_id)
    print "\n"

    ## Initiate superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())

    # Write RMSD to file:
    # print super_imposer.rms
    # with open("rmsd_matrix.txt", "a") as ouf:
    #	ouf.write(str(super_imposer.rms))
    #	ouf.write(' ')

    ## Save aligned structure in reference.pdb
    io = Bio.PDB.PDBIO()
    io.set_structure(sample_structure)
    io.save('aligned_pdbs/' + reference + "_aligned_" + sample + ".pdb")
    io.set_structure(ref_structure)
    io.save('reference_structures/' + reference + '_ref_' + sample + '.pdb')


def print_usage(name):
    print "Usage : " + name + "<Folder> <start_id> <end_id>"


if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])