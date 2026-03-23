#####################
# Import Statements #
#####################
import pandas as pd
from pymol import cmd, util, stored
import os


##############
# Script I/O #
##############
output_dir = os.path.abspath(r"D:\OneDrive\Projects\Deinococcus\PprA\2015-20xx Robert Szabla\Papers\PprA\deposited_scripts\PprA_interface_analyzer")
output_filename = "pprA_residue_contacts.csv"
model_sources = {
    '6bdu' : 'https://modelarchive.org/api/projects/ma-hl5ox?type=basic__original_pdb',
    '6mc6' : 'https://modelarchive.org/api/projects/ma-dm7xl?type=basic__original_pdb',
    #'6neo' : 'https://modelarchive.org/api/projects/ma-bfx40?type=basic__original_pdb',
    #'9om8' : 'https://modelarchive.org/api/projects/ma-y9ezg?type=basic__original_pdb',
    #'9or6' : 'https://modelarchive.org/api/projects/ma-uw4cd?type=basic__original_pdb',
    #'9yi3' : 'https://modelarchive.org/api/projects/ma-1s1jd?type=basic__original_pdb',
    #'9yl4' : 'https://modelarchive.org/api/projects/ma-wefg9?type=basic__original_pdb',
    #'9yup' : 'https://modelarchive.org/api/projects/ma-jbqd3?type=basic__original_pdb',
    #'6mc8' : 'https://modelarchive.org/api/projects/ma-1xc76?type=basic__original_pdb',
    '6o5l' : 'https://modelarchive.org/api/projects/ma-157xn?type=basic__original_pdb',
    #'9y6e' : 'https://modelarchive.org/api/projects/ma-jtmui?type=basic__original_pdb', # uncomment when publicly available
    #'6a27' : 'https://files.rcsb.org/download/6A27.pdb',
    #'6a28' : 'https://files.rcsb.org/download/6A28.pdb',
    #'6a29' : 'https://files.rcsb.org/download/6A29.pdb',
}


######################
# Processing Options #
######################
# These two variables controls whether the find contacts routine should be run. 
# If neither inter- or intra- modes are activated, just display the dimer pairs
find_intramolecular = True
find_intermolecular = True

# Recommended to keep these values as default:
remove_hydrogens=False

distance_cutoffs = {
    'pi-pi'       : 5.5,
    'ionic'       : 4.7,
    'H-bond'      : 3.6,
    'pi-cat'      : 6.0,
    'hydrophobic' : 4.5,
}

interaction_definitions = {
    'ionic'          : [('((resn ARG + resn LYS + resn HIS) and sc. and elem N)', '((resn ASP + resn GLU) and sc. and elem O)'), 
                        ('((resn ASP + resn GLU) and sc. and elem O)', '((resn ARG + resn LYS + resn HIS) and sc. and elem N)')],
                     
    'pi-pi'          : [('(name ARO)', '(name ARO)')],
                     
    'H-bond'         : [('(elem H and (bound_to (elem O+N)))', '(elem O+N)'), 
                        ('(elem O+N)', '(elem H and (bound_to (elem O+N)))')],
                     
    'pi-cat'         : [('((resn ARG + resn LYS + resn HIS) and sc. and elem N)', '(name ARO)'),
                        ('(name ARO)', '((resn ARG + resn LYS + resn HIS) and sc. and elem N)')],
                              
    'hydrophobic'    : [('((resn PRO + resn TRP + resn MET + resn MSE + resn PHE + resn ILE + resn LEU + resn VAL + resn ALA + resn TYR + resn LYS + resn ARG) and sc. and elem C)', '((resn PRO + resn TRP + resn MET + resn MSE + resn PHE + resn ILE + resn LEU + resn VAL + resn ALA + resn TYR + resn LYS + resn ARG) and sc. and elem C)')]
}


##########################
# PprA Interface library #
##########################
interface_pairs = {
    'P' : {
            '6o5l_P_AA' : ['6o5l', '/6o5l_sym_00000000//A', '/6o5l_sym_03000000//A'], 
            '6bdu_P_AB' : ['6bdu', '/6bdu_sym_00000000//A', '/6bdu_sym_02010000//B'],
            '6mc6_P_AB' : ['6mc6', '/6mc6_sym_00000000//A', '/6mc6_sym_02-10000//B'],
            '6neo_P_BB' : ['6neo', '/6neo_sym_00000000//B', '/6neo_sym_09000000//B'],
            '9y6e_P_DD' : ['9y6e', '/9y6e_sym_00000000//D', '/9y6e_sym_0500-100//D'],
            '6a29_P_FH' : ['6a29', '/6a29_sym_00000000//F', '/6a29_sym_01010101//H'],
          },
    'F' : {
            '6o5l_F_AA' : ['6o5l', '/6o5l_sym_00000000//A', '/6o5l_sym_06000000//A'],  
            '6bdu_F_AA' : ['6bdu', '/6bdu_sym_00000000//A', '/6bdu_sym_0200-100//A'],
            '6mc6_F_AA' : ['6mc6', '/6mc6_sym_00000000//A', '/6mc6_sym_0200-100//A'],
            '9om8_F_AB' : ['9om8', '/9om8_sym_00000000//A', '/9om8_sym_030000-1//B'],  
            '9om8_F_CC' : ['9om8', '/9om8_sym_00000000//C', '/9om8_sym_01000000//C'],
            '9om8_F_DD' : ['9om8', '/9om8_sym_00000000//D', '/9om8_sym_010000-1//D'],
            '9or6_F_BB' : ['9or6', '/9or6_sym_00000000//B', '/9or6_sym_01-100-1//B'],
            '9yi3_F_AA' : ['9yi3', '/9yi3_sym_00000000//A', '/9yi3_sym_01000000//A'],
            '9yi3_F_BB' : ['9yi3', '/9yi3_sym_00000000//B', '/9yi3_sym_01000100//B'],
            '9yl4_F_AF' : ['9yl4', '/9yl4_sym_00000000//A', '/9yl4_sym_01-10001//F'],
            '9yl4_F_BC' : ['9yl4', '/9yl4_sym_00000000//B', '/9yl4_sym_00000000//C'],
            '9yl4_F_DE' : ['9yl4', '/9yl4_sym_00000000//D', '/9yl4_sym_00000000//E'],
            '9yl4_F_GL' : ['9yl4', '/9yl4_sym_00000000//G', '/9yl4_sym_01000001//L'],
            '9yl4_F_HI' : ['9yl4', '/9yl4_sym_00000000//H', '/9yl4_sym_00000000//I'],
            '9yl4_F_JK' : ['9yl4', '/9yl4_sym_00000000//J', '/9yl4_sym_00000000//K'],
            '9yup_F_AC' : ['9yup', '/9yup_sym_00000000//A', '/9yup_sym_00000000//C'],
            '6a29_F_AB' : ['6a29', '/6a29_sym_00000000//A', '/6a29_sym_00000000//B'],
            '6a29_F_CD' : ['6a29', '/6a29_sym_00000000//C', '/6a29_sym_00000000//D'],
            '6a29_F_EF' : ['6a29', '/6a29_sym_00000000//E', '/6a29_sym_00000000//F'],
            '6a29_F_GH' : ['6a29', '/6a29_sym_00000000//G', '/6a29_sym_00000000//H'],
          },  
    'S' : {
            '6o5l_S_AA' : ['6o5l', '/6o5l_sym_00000000//A', '/6o5l_sym_10000000//A'],
            '6bdu_S_AB' : ['6bdu', '/6bdu_sym_00000000//A', '/6bdu_sym_00000000//B'],
            '6mc6_S_AB' : ['6mc6', '/6mc6_sym_00000000//A', '/6mc6_sym_00000000//B'],
            '6mc8_S_AB' : ['6mc8', '/6mc8_sym_00000000//A', '/6mc8_sym_00000000//B'],
            '6neo_S_BB' : ['6neo', '/6neo_sym_00000000//B', '/6neo_sym_06000001//B'],
            '9om8_S_AD' : ['9om8', '/9om8_sym_00000000//A', '/9om8_sym_00000100//D'],
            '9om8_S_BC' : ['9om8', '/9om8_sym_00000000//B', '/9om8_sym_00000100//C'],
            '9or6_S_AB' : ['9or6', '/9or6_sym_00000000//A', '/9or6_sym_00000000//B'],
            '9yi3_S_AB' : ['9yi3', '/9yi3_sym_00000000//A', '/9yi3_sym_00000000//B'],
            '9yl4_S_AB' : ['9yl4', '/9yl4_sym_00000000//A', '/9yl4_sym_00000000//B'],
            '9yl4_S_CD' : ['9yl4', '/9yl4_sym_00000000//C', '/9yl4_sym_00000000//D'],
            '9yl4_S_EF' : ['9yl4', '/9yl4_sym_00000000//E', '/9yl4_sym_00000000//F'],
            '9yl4_S_GH' : ['9yl4', '/9yl4_sym_00000000//G', '/9yl4_sym_00000000//H'],
            '9yl4_S_IJ' : ['9yl4', '/9yl4_sym_00000000//I', '/9yl4_sym_00000000//J'],
            '9yl4_S_KL' : ['9yl4', '/9yl4_sym_00000000//K', '/9yl4_sym_00000000//L'],
            '9yup_S_AB' : ['9yup', '/9yup_sym_00000000//A', '/9yup_sym_00000000//B'],
            '9yup_S_CD' : ['9yup', '/9yup_sym_00000000//C', '/9yup_sym_00000000//D'],
            '9y6e_S_AB' : ['9y6e', '/9y6e_sym_00000000//A', '/9y6e_sym_00000000//B'],
            '9y6e_S_CD' : ['9y6e', '/9y6e_sym_00000000//C', '/9y6e_sym_00000000//D'],
            '6a27_S_AB' : ['6a27', '/6a27_sym_00000000//A', '/6a27_sym_00000000//B'],
            '6a28_S_AB' : ['6a28', '/6a28_sym_00000000//A', '/6a28_sym_00000000//B'],
          },
    'M' : {
            '6o5l_M_A' : ['6o5l', '/6o5l_sym_00000000//A', '/6o5l_sym_00000000//A'],
            '6bdu_M_A' : ['6bdu', '/6bdu_sym_00000000//A', '/6bdu_sym_00000000//A'],
            '6bdu_M_B' : ['6bdu', '/6bdu_sym_00000000//B', '/6bdu_sym_00000000//B'],
            '6mc8_M_A' : ['6mc8', '/6mc8_sym_00000000//A', '/6mc8_sym_00000000//A'],
            '6mc8_M_B' : ['6mc8', '/6mc8_sym_00000000//B', '/6mc8_sym_00000000//B'],
            '6mc6_M_A' : ['6mc6', '/6mc6_sym_00000000//A', '/6mc6_sym_00000000//A'],
            '6mc6_M_B' : ['6mc6', '/6mc6_sym_00000000//B', '/6mc6_sym_00000000//B'],
            '6neo_M_B' : ['6neo', '/6neo_sym_00000000//B', '/6neo_sym_00000000//B'],
            '9om8_M_A' : ['9om8', '/9om8_sym_00000000//A', '/9om8_sym_00000000//A'],
            '9om8_M_B' : ['9om8', '/9om8_sym_00000000//B', '/9om8_sym_00000000//B'],
            '9om8_M_C' : ['9om8', '/9om8_sym_00000000//C', '/9om8_sym_00000000//C'],
            '9om8_M_D' : ['9om8', '/9om8_sym_00000000//D', '/9om8_sym_00000000//D'],
            '9or6_M_A' : ['9or6', '/9or6_sym_00000000//A', '/9or6_sym_00000000//A'],
            '9or6_M_B' : ['9or6', '/9or6_sym_00000000//B', '/9or6_sym_00000000//B'],
            '9yi3_M_A' : ['9yi3', '/9yi3_sym_00000000//A', '/9yi3_sym_00000000//A'],
            '9yi3_M_B' : ['9yi3', '/9yi3_sym_00000000//B', '/9yi3_sym_00000000//B'],
            '9yl4_M_A' : ['9yl4', '/9yl4_sym_00000000//A', '/9yl4_sym_00000000//A'],
            '9yl4_M_B' : ['9yl4', '/9yl4_sym_00000000//B', '/9yl4_sym_00000000//B'],
            '9yl4_M_C' : ['9yl4', '/9yl4_sym_00000000//C', '/9yl4_sym_00000000//C'],
            '9yl4_M_D' : ['9yl4', '/9yl4_sym_00000000//D', '/9yl4_sym_00000000//D'],
            '9yl4_M_E' : ['9yl4', '/9yl4_sym_00000000//E', '/9yl4_sym_00000000//E'],
            '9yl4_M_F' : ['9yl4', '/9yl4_sym_00000000//F', '/9yl4_sym_00000000//F'],
            '9yl4_M_G' : ['9yl4', '/9yl4_sym_00000000//G', '/9yl4_sym_00000000//G'],
            '9yl4_M_H' : ['9yl4', '/9yl4_sym_00000000//H', '/9yl4_sym_00000000//H'],
            '9yl4_M_I' : ['9yl4', '/9yl4_sym_00000000//I', '/9yl4_sym_00000000//I'],
            '9yl4_M_J' : ['9yl4', '/9yl4_sym_00000000//J', '/9yl4_sym_00000000//J'],
            '9yl4_M_K' : ['9yl4', '/9yl4_sym_00000000//K', '/9yl4_sym_00000000//K'],
            '9yl4_M_L' : ['9yl4', '/9yl4_sym_00000000//L', '/9yl4_sym_00000000//L'],
            '9yup_M_A' : ['9yup', '/9yup_sym_00000000//A', '/9yup_sym_00000000//A'],
            '9yup_M_B' : ['9yup', '/9yup_sym_00000000//B', '/9yup_sym_00000000//B'],
            '9yup_M_C' : ['9yup', '/9yup_sym_00000000//C', '/9yup_sym_00000000//C'],
            '9yup_M_D' : ['9yup', '/9yup_sym_00000000//D', '/9yup_sym_00000000//D'],
            '9y6e_M_A' : ['9y6e', '/9y6e_sym_00000000//A', '/9y6e_sym_00000000//A'],
            '9y6e_M_B' : ['9y6e', '/9y6e_sym_00000000//B', '/9y6e_sym_00000000//B'],
            '9y6e_M_C' : ['9y6e', '/9y6e_sym_00000000//C', '/9y6e_sym_00000000//C'],
            '9y6e_M_D' : ['9y6e', '/9y6e_sym_00000000//D', '/9y6e_sym_00000000//D'],
            '6a27_M_A' : ['6a27', '/6a27_sym_00000000//A', '/6a27_sym_00000000//A'],
            '6a27_M_B' : ['6a27', '/6a27_sym_00000000//B', '/6a27_sym_00000000//B'],
            '6a28_M_A' : ['6a28', '/6a28_sym_00000000//A', '/6a28_sym_00000000//A'],
            '6a28_M_B' : ['6a28', '/6a28_sym_00000000//B', '/6a28_sym_00000000//B'],
            '6a29_M_A' : ['6a29', '/6a29_sym_00000000//A', '/6a29_sym_00000000//A'],
            '6a29_M_B' : ['6a29', '/6a29_sym_00000000//B', '/6a29_sym_00000000//B'],
            '6a29_M_C' : ['6a29', '/6a29_sym_00000000//C', '/6a29_sym_00000000//C'],
            '6a29_M_D' : ['6a29', '/6a29_sym_00000000//D', '/6a29_sym_00000000//D'],
            '6a29_M_E' : ['6a29', '/6a29_sym_00000000//E', '/6a29_sym_00000000//E'],
            '6a29_M_F' : ['6a29', '/6a29_sym_00000000//F', '/6a29_sym_00000000//F'],
            '6a29_M_G' : ['6a29', '/6a29_sym_00000000//G', '/6a29_sym_00000000//G'],
            '6a29_M_H' : ['6a29', '/6a29_sym_00000000//H', '/6a29_sym_00000000//H'],
          }
}


######################
# Rendering Settings #
######################
# Align all dimer pairs to these reference structures:
ref_dimer = {
    'P' : '6o5l_P_AA',
    'F' : '6o5l_F_AA',
    'S' : '6o5l_S_AA',
    'M' : '6o5l_M_A',
}
# Use either a single chain or both chains when aligning:
align_both_chains = True
# Domain colors to be applied after domain selections are created:
domain_colors = {
    'P_domain' : 'deeppurple', 
    'S_domain' : 'brightorange', 
    'F_domain' : 'forest', 
    'linker'   : 'gray50',
}
interaction_colors = {
    'ionic'       : 'red',
    'pi-pi'       : 'red',
    'H-bond'      : 'slate',        
    'pi-cat'      : 'slate',
    'hydrophobic' : 'wheat',
}


###########################
# PprA Domain Definitions #
###########################
# Map PDB codes to Deinococcus species for proper domain selection boundaries.
ppra_homologs = {
    '6o5l' : 'per',
    '6bdu' : 'rad',
    '6mc8' : 'des',
    '6mc6' : 'rad',
    '6neo' : 'rad',
    '6a27' : 'rad',
    '6a28' : 'rad',
    '6a29' : 'rad',
    '9om8' : 'rad',
    '9or6' : 'rad',
    '9yi3' : 'rad',
    '9yl4' : 'rad',
    '9yup' : 'rad',
    '9y6e' : 'rad',
}
# PprA domain definitions for different homologs (by residue index)
domain_boundaries = {
    'rad' : {
             'linker'   : '(polymer.protein and (resi -8))',
             'P_domain' : '(polymer.protein and (resi 9-111))',
             'S_domain' : '(polymer.protein and (resi 112-172 + resi 233-268))',
             'F_domain' : '(polymer.protein and (resi 173-232 + resi 269-284))',
            },
    'geo' : {
             'linker'   : '(polymer.protein and (resi -24))',
             'P_domain' : '(polymer.protein and (resi 25-131))',
             'S_domain' : '(polymer.protein and (resi 132-192 + resi 253-288))',
             'F_domain' : '(polymer.protein and (resi 193-252 + resi 289-302))',
            },
    'per' : {
             'linker'   : '(polymer.protein and (resi -18))',
             'P_domain' : '(polymer.protein and (resi 19-123))',
             'S_domain' : '(polymer.protein and (resi 124-184 + resi 245-280))',
             'F_domain' : '(polymer.protein and (resi 185-244 + resi 281-294))',
            },
    'des' : {
             'linker'   : '(polymer.protein and (resi -16))',
             'P_domain' : '(polymer.protein and (resi 17-123))',
             'S_domain' : '(polymer.protein and (resi 124-184 + resi 245-280))',
             'F_domain' : '(polymer.protein and (resi 185-244 + resi 281-294))',
            },
}



########################
# Function Definitions #
########################

# Function to make atom selections grouped by domains.
def select_domains(selection:str, homolog:str, prefix:str=None, group_name:str=None):
    print(f"Preparing domain selections for {selection}...")
    group_member_list = []
    for domain in domain_boundaries[homolog].keys():
        selection_expression = f"(({selection}) and ({domain_boundaries[homolog][domain]}))"
        if prefix:
            selection_name = f"{prefix}_{domain}_sele"
        else:
            selection_name = f"{domain}_sele"
        if cmd.count_atoms(selection_expression):
            cmd.select(selection_name, selection_expression, merge=1)
            group_member_list.append(selection_name)

    # group selections together
    if group_name:
        group_members = ' '.join(group_member_list)
        cmd.group(group_name, group_members)
    print("Done!")
    return



# Extract relevant dimer from crystallographic symmetry
def dimer_from_symmetry(dimer_triplet, dimer_name, cutoff=5.0):
    base_symmate = dimer_triplet[0]
    chain_A = dimer_triplet[1]
    chain_B = dimer_triplet[2]
    print(f'Extracting {dimer_name} from crystallographic symmetry...')
    cmd.copy(f"{base_symmate}_sym_00000000", base_symmate)
    chains_present = False
    while (not chains_present) and (cutoff < 1000):
        try:
            chains_present = (cmd.select('sym_A', chain_A) and cmd.select('sym_B', chain_B))  #Evaluates to true if both chains are present
            print(f'Found {chain_A} and {chain_B}!')
        except:
            print(f'Could not find {chain_A} and {chain_B}... expanding search to {cutoff+50}A')
            cmd.symexp(base_symmate+'_sym_', base_symmate, base_symmate, cutoff)
            cutoff+=50
    
    tmp_A = f'{dimer_name}_A'
    tmp_B = f'{dimer_name}_B'
    cmd.copy_to(tmp_A, chain_A, rename='')
    cmd.copy_to(tmp_B, chain_B, rename='')
    cmd.alter(tmp_A, f"segi='A'")
    cmd.alter(tmp_B, f"segi='B'")
    cmd.create(dimer_name, f'{tmp_A} + {tmp_B}')
    cmd.delete(base_symmate+'_sym_*')
    cmd.delete(f'{tmp_A} + {tmp_B} + sym_A + sym_B')
    return



def invert_interaction(interaction):
    interaction_inv = {
        'segi1':interaction['segi2'], 'chain1':interaction['chain2'], 'resi1':interaction['resi2'], 'resn1':interaction['resn2'],
        'segi2':interaction['segi1'], 'chain2':interaction['chain1'], 'resi2':interaction['resi1'], 'resn2':interaction['resn1'],
        'contact_type' : interaction['contact_type'],}
    return interaction_inv



def find_aromatics(model):
    print(f'\tFinding aromatic centers for {model}...')
    stored.aromatic_list = {}
    aromatic_resn = '(resn PHE + resn TYR + resn TRP + resn HIS)'
    cmd.iterate(f'{model} and n. CA and {aromatic_resn}' , "stored.aromatic_list[(segi, chain, resi, resn)]=''" )
    for resi in list(stored.aromatic_list.keys()):
        cmd.pseudoatom(model, f'( /{model}/{resi[0]}/{resi[1]}/{resi[2]} and sc. and not (elem H + elem O + n. CB))', segi=resi[0], chain=resi[1], resi=resi[2], resn=resi[3], name='ARO')



# function to convert from cmd.find_pairs() output, to a useable contact_list format
def contact_list_from_atom_pairs(atom_pair_list, contact_type=''):
    # Convert atom indices in atom_pair_list to (segi, chain, resi, resn) identifiers
    stored.index_conversions = {}
    atom_indices_int = []
    for atom_pair in atom_pair_list:
        index_A = atom_pair[0][1]
        index_B = atom_pair[1][1]
        atom_indices_int.append(index_A)
        atom_indices_int.append(index_B)
    atom_indices_str = [str(num) for num in atom_indices_int]
    model = atom_pair_list[0][0][0]
    atom_selection = f"{model} and (index {'+'.join(atom_indices_str)})"
    cmd.iterate(atom_selection, "stored.index_conversions[index] = (segi, chain, resi, resn)")

    # construct interaction_list
    contact_list = []
    for atom_pair in atom_pair_list:
        index_A = atom_pair[0][1]
        index_B = atom_pair[1][1]
        contact = {
            'segi1'  : stored.index_conversions[index_A][0], 
            'chain1' : stored.index_conversions[index_A][1], 
            'resi1'  : stored.index_conversions[index_A][2], 
            'resn1'  : stored.index_conversions[index_A][3],
            'segi2'  : stored.index_conversions[index_B][0], 
            'chain2' : stored.index_conversions[index_B][1], 
            'resi2'  : stored.index_conversions[index_B][2], 
            'resn2'  : stored.index_conversions[index_B][3],
            'contact_type' : contact_type,
        }
        if contact not in contact_list:
            contact_list.append(contact)
    return contact_list




def find_vdW_contacts(model, selection_1, selection_2):
    # Hydrophobic interactions
    print(f'\tFinding van-der-Waals interactions for {model} between {selection_1} and {selection_2}...')
    if selection_1 == selection_2:
        dist_obj = f'{model}_intra_vdW'
    else:
        dist_obj = f'{model}_inter_vdW'
    interaction_list=[]
    #stored.resi1_list = {}
    #cmd.iterate(f'{selection_1} and n. CA and {nonpolar_resn}' , "stored.resi1_list[(segi, chain, resi, resn)]=''" )
    #for resi1 in list(stored.resi1_list.keys()):
    #    #print(f'\n\tChecking {resi1[0]}{resi1[1]}-{resi1[2]}', end='')
    #    resi1_sele = f'(/{model}/{resi1[0]}/{resi1[1]}/{resi1[2]} and sc. and elem C)'
    #    stored.resi2_list = {}
    #    resi2_sele = f'(({selection_2} and sc. and elem C and {nonpolar_resn} and not (/{model}/{resi1[0]}/{resi1[1]}/{resi1[2]})) within {vdw_cutoff} of {resi1_sele})'
    #    cmd.iterate(f'(byres ({resi2_sele})) and n. CA' , "stored.resi2_list[(segi, chain, resi, resn)]=''" )
    #    if (len(stored.resi2_list.keys()) > 0):# and (resi1[0],resi1[1],resi1[2]) != (resi2[0],resi2[1],resi1[2]):
    #        cmd.distance(dist_obj, resi1_sele, resi2_sele, label=0, cutoff=vdw_cutoff)
    #        for resi2 in list(stored.resi2_list.keys()):
    #            interaction = {
    #            'segi1':resi1[0], 'chain1':resi1[1], 'resi1':resi1[2], 'resn1':resi1[3],
    #            'segi2':resi2[0], 'chain2':resi2[1], 'resi2':resi2[2], 'resn2':resi2[3],
    #            'contact_type' : 'vdW',}
    #            interaction_inv = invert_interaction(interaction)
    #            if (interaction not in interaction_list) and (interaction_inv not in interaction_list):
    #                interaction_list.append(interaction)
    return interaction_list



def remove_items(test_list, item): 
    # using list comprehension to perform the task 
    res = [i for i in test_list if i != item] 
    return res 


def remove_duplicate_h_bonds(h_bond_contacts, salt_bridge_contacts):
    # Remove H-bonds that are also salt bridges:
    print(f'\tRemoving duplicate H-bonds...')
    for salt_bridge in salt_bridge_contacts:
        h_bond = salt_bridge.copy()
        h_bond['contact_type'] = 'H-Bond'
        h_bond_inv = invert_interaction(h_bond)
        if h_bond in h_bond_contacts:
            print(f'\t\t{h_bond}*\n\t\t{salt_bridge}')
            h_bond_contacts = remove_items(h_bond_contacts, h_bond)
        if h_bond_inv in h_bond_contacts:
            print(f'\t\t{h_bond_inv}*\n\t\t{salt_bridge}')
            h_bond_contacts = remove_items(h_bond_contacts, h_bond_inv)
    return h_bond_contacts


def find_interactions(model, sele_1, sele_2, column_constants={}):
    print(f'Finding intermolecular interactions for {model}...')

    contact_list=[]

    for contact_type, definitions in interaction_definitions.items():
        cutoff = distance_cutoffs[contact_type]
        print(f'\tFinding {contact_type} interactions between {sele_1} and {sele_2} within {cutoff}A...')
        
        if sele_1 == sele_2: dist_obj = f'{model}_intra_{contact_type}'
        else:                dist_obj = f'{model}_inter_{contact_type}'
        found_contacts = False
        
        for subsele_1, subsele_2, in definitions:
            atom_sele_1 = f"(({sele_1}) and ({subsele_1}))"
            atom_sele_2 = f"(({sele_2}) and ({subsele_2}))"
            atom_pair_list = cmd.find_pairs(atom_sele_1, atom_sele_2, cutoff=cutoff, mode=0)
            cmd.distance(dist_obj, atom_sele_1, atom_sele_2, cutoff=cutoff, mode=0, label=0)

            if atom_pair_list:
                contact_list = contact_list + contact_list_from_atom_pairs(atom_pair_list, contact_type=contact_type)
                found_contacts = True
        
        cmd.color(interaction_colors[contact_type], dist_obj)
        if not found_contacts:
            cmd.delete(dist_obj)


    #salt_bridge_contacts = find_salt_bridge_contacts(model, selection_1, selection_2)
    #pi_pi_contacts       =       find_pi_pi_contacts(model, selection_1, selection_2)
    ##h_bond_contacts      =      find_h_bond_contacts(model, selection_1, selection_2)
    ##h_bond_contacts = remove_duplicate_h_bonds(h_bond_contacts, salt_bridge_contacts)
    ##pi_cat_contacts      =      find_pi_cat_contacts(model, selection_1, selection_2)
    ##vdW_contacts         =      find_vdW_contacts(model, selection_1, selection_2)


    contacts = pd.DataFrame(contact_list)
    for name,value in column_constants.items():
        contacts[name] = value

    return contacts



def main():
    contacts = pd.DataFrame(columns=[
        'interface',
        'homolog',
        'pdb',
        'sym1',
        'segi1',
        'chain1', 
        'resi1', 
        'resn1',
        'sym2', 
        'segi2',
        'chain2', 
        'resi2', 
        'resn2', 
        'contact_type',
    ])
    
    for pdb_code, url in model_sources.items():
        print(f"Loading {pdb_code} model from {url} ...")
        cmd.load(url, object=pdb_code, format="pdb")
        cmd.remove(f"{pdb_code} and not polymer.protein")
        if remove_hydrogens:
            cmd.remove(f"{pdb_code} and elem H")
        else:
            if cmd.count_atoms(f"{pdb_code} and elem H") == 0:
                cmd.h_add(pdb_code)
        find_aromatics(pdb_code)
        cmd.group('crystal_structures', pdb_code)
    

    for interface in interface_pairs.keys():
        for dimer_name, dimer_triplet in interface_pairs[interface].items():
            pdb_code = dimer_triplet[0]
            if (pdb_code in model_sources.keys()) and (interface in ref_dimer.keys()):
                homolog=ppra_homologs[pdb_code]
                sym1 = interface_pairs[interface][dimer_name][1].split('/')[1].split('_sym_')[1]
                sym2 = interface_pairs[interface][dimer_name][2].split('/')[1].split('_sym_')[1]
                segi1='A'
                segi2='B'
                chain1=interface_pairs[interface][dimer_name][1].split('/')[3]
                chain2=interface_pairs[interface][dimer_name][2].split('/')[3]
                sele1=f"/{dimer_name}/{segi1}/{chain1}"
                sele2=f"/{dimer_name}/{segi2}/{chain2}"
                ''
                dimer_from_symmetry(dimer_triplet, dimer_name)
                if interface == 'M':
                    cmd.remove(f"/{dimer_name}/B")
                    segi2=segi1
                    sele2=sele1
                    
                if ref_dimer[interface] in cmd.get_object_list():
                    if align_both_chains:
                        mobile = f'/{dimer_name}'
                        target = f'/{ref_dimer[interface]}'
                    else:
                        mobile = f'/{dimer_name}/A'
                        target = f'/{ref_dimer[interface]}/A'
                else:
                    mobile = f'/{dimer_name}/A'
                    target = f"bychain (first ({next(iter(model_sources))}))"
                if not mobile == target:
                    cmd.super(mobile, target)
                
                column_constants = {
                    'interface':interface,
                    'homolog':homolog,
                    'pdb':pdb_code,
                    'sym1':sym1,
                    'segi1':segi1,
                    'chain1':chain1,
                    'sym2':sym2, 
                    'segi2':segi2,
                    'chain2':chain2, 
                }
                if ((interface == 'M') and find_intramolecular) or ((interface != 'M') and find_intermolecular):
                    single_dimer_contacts = find_interactions(model=dimer_name, sele_1=sele1, sele_2=sele2, column_constants=column_constants)
                    contacts = pd.concat([contacts, single_dimer_contacts], ignore_index=True)
                    csv_file = os.path.join(output_dir, output_filename)
                    contacts.to_csv(csv_file)
                    #group_list = [model, model+'_inter*', model+'_intra*', model+'_surface']
                    cmd.group(f"{dimer_name}_group", f"{dimer_name}*")
            if ((interface == 'M') and find_intramolecular) or ((interface != 'M') and find_intermolecular):
                cmd.group(f"{interface}_interface", f"{dimer_name}_group", action='add')
            if not ((interface == 'M') and find_intramolecular) or ((interface != 'M') and find_intermolecular):
                cmd.group(f"{interface}_interface", dimer_name, action='add')
    
    for pdb_code in model_sources.keys():
        homolog = ppra_homologs[pdb_code]
        select_domains(pdb_code, homolog=homolog, prefix='xtal', group_name='crystal_structures')
        select_domains(f"{pdb_code}_*", homolog=homolog, group_name='domain_selections')
        for domain_name, domain_sele in domain_boundaries[homolog].items():
            color = domain_colors[domain_name]
            atom_selection = f"{pdb_code}* and {domain_sele}"
            cmd.color(color, atom_selection)
        #cmd.spectrum('segi')
    util.cnc()
    cmd.show('lines')
    cmd.hide("(all and hydro and (elem C extend 1))") # hide nonpolar hydrogens
    cmd.disable('crystal_structures')
    cmd.orient("*_interface")
    #cmd.center(ref_dimer)
    return contacts


contacts = main()

print(contacts)