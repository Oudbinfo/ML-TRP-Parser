import json
from Bio.PDB import Select


def get_res_index(res_num: str, chain_residues: list):
    """Returns the residue index of a residue"""
    res_index = [i for i, chain_res in enumerate(chain_residues) if str(chain_res.id[1]) == str(res_num)]
    return res_index[0]


def get_chain_range(start, end, chain_residues: list):
    """Returns the range of subsection in a chain"""
    chain_start = get_res_index(start, chain_residues)
    chain_end = get_res_index(end, chain_residues)
    chain_range = chain_residues[chain_start:chain_end + 1]
    return chain_range


def get_res_frames(res_range, length, step):
    """Returns a list of fragmented frames of a residue range"""
    frames = []
    start = 0
    for i in range(len(res_range)):
        frame = res_range[start:start + length]
        if len(frame) >= round(length * 0.75):
            frames.append(frame)
            start += step
    return frames


def get_structure(res_range, res_chain, structure, out_path, io_handler):
    """Trims and saves the structure in the range of start and end residues on a chain"""

    # Select residues within the specified range and chain
    class ResSelect(Select):
        def accept_residue(self, res):
            if res in res_range and res.parent.id == res_chain:
                return True
            else:
                return False

    # Create a PDBIO object, set the structure, and save the selected residues to the output path
    io_handler.set_structure(structure)
    io_handler.save(out_path, ResSelect())


