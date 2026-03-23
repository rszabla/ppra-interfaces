#####################
# Import Statements #
#####################
import cairo
import pandas as pd
import math
import os
import argparse



##############
# Script I/O #
##############
# Input and output are now controlled via command line arguments in main()


##########################
# Figure Render Settings #
##########################
w = 144 #288                       # width of figure - units in px @ 96 px/inch
outline_width = 0.75          # Width of rectangle lines - units in px
interaction_width = 0.75      # Width of interaction lines - units in px
font_size = 5                 # Size of font in px
transparent_bg = True         # Transparent background
arc_angle = 100               # Angle of circle that is spanned by intra-molecular arc (in degrees)
border_width = outline_width  # 
chain_height = 10             #
interchain_space = 15         #
text_pad = outline_width      # Text box padding 

# Domain colors
domain_colors = {
    'P':[0.6, 0.1, 0.6, 1], # deeppurple
    'S':[1.0, 0.7, 0.2, 1], # brightorange
    'F':[0.2, 0.6, 0.2, 1], # forest
}
# Interaction colors
contact_rgba = {
    'hydrophobic' : [0.5, 0.5, 1.0, 0.60], # bluewhite
    'pi-cation'   : [0.1, 0.1, 1.0, 0.80], # br0
    'H-Bond'      : [0.1, 0.1, 1.0, 0.80], # br0
    'pi-pi'       : [1.0, 0.1, 0.1, 0.80], # br9
    'ionic'       : [1.0, 0.1, 0.1, 0.80], # br9
}
# Interaction line dashes
contact_dashes = {
    'hydrophobic' : [1.0,  0.0], # 
    'pi-cation'   : [1.0,  0.0], # 
    'H-Bond'      : [1.0,  0.2], # 
    'pi-pi'       : [1.0,  0.2], # 
    'ionic'       : [1.0,  0.0], # 
}

black_rgba = [0,0,0,1]
white_rgba = [1,1,1,1]

document_unit = 4
# 1) USER User unit, a value in the current coordinate system. If used in the root element for the initial coordinate systems it corresponds to pixels
# 2) EM The size of the element's font
# 3) EX The x-height of the element's font
# 4) PX Pixels (1px = 1/96th of 1in)
# 5) IN Inches (1in = 2.54cm = 96px)
# 6) CM Centimeters (1cm = 96px/2.54)
# 7) MM Millimeters (1mm = 1/10th of 1cm)
# 8) PT Points (1pt = 1/72th of 1in)
# 9) PC Picas (1pc = 1/6th of 1in)
# 10) PERCENT Percent, a value that is some fraction of another reference value



#################################
# PprA Sequence and Domain Info #
#################################
# PprA homolog that will be used as a reference for residue numbering
ref_homolog = 'rad'
prot_length = 284

homolog_alignment = {
    'rad' : 'MARAKAKD----------------QTDGIYAAFDTLMSTAGVDSQIAALAASEADAGTLDAALTQSLQEAQGRWGLGLHHLRHEARLTDD----GDIEILTDGRPSARVSEGFGALAQAYAPMQALDERGLSQWAALGEGYRAPGDLPLAQLKVLIEHARDFETDWSAGRGETFQRVWRKGDTLFVEVARPASAEAALSDAAWDVIASIKDRAFQRELMRRSEKDGMLGALLGARHAGAKANLAQLPEAHFTVQAFVQTLSGAAARNAEEYRAALKTAAAALEEYQGVTTRQLSEVLRHGLRES',
    'geo' : 'MTKTKQKDRNALQESPRPNVSGPVGSEDVLKSFDALMATADVDSQIHALAESGADEETLGRELTLALQLAQDRWGLGLLHLRHDAALARTPEGTPDVVLRADGAVVARLSDGPAAIARSYASMQALGAEGLSEWGVLPDGHRVTLKGGSGQLRVLVEDARDFETHWTAERGGVWSRTWRQGETLVVEVHRPASPATVLADAAWDVITSIKDRNFQRELMERSNSVGMLGALLGARHSGAGSALDRLPSAHFTVRSAVIRESGVSARSLERWKAMLREGMEQLEALQKTVTRELAEVLSHGLR--',
    'per' : 'MTKASKKSDA-----P-RTREAPTPREDALRGFDALMATAGVESTIVKHAASGADSQTLNDELTRSLQLAHDRWGLGLLHLRHEARLDRGE--DTDVILLVDGREVARLSQGAAAISATYETMRAQNADDLSDWGVLPEGHRVTLKAGNNQMRVLVEDARDFETHWSSERGGAFVRTWRQGETLAVEVHRPASPGTALADAAWDAIMSIKDRNFQRELMERSNSVGMLGALLGARHKDAGRALERLPEAHFAVRSTVVRMTGGAQREFDQWRSMVREGLDQLDELQKTTTRHLTEILRHGLK--',
    'des' : 'MTKTTRRKAT-----PT---EPVTSSVNPLARFAELVATAGLQSDVQALADSGADDTTLEAQLTQELRLAHDRWGLGLLHLQHSARLIHTDGVPSDIALLVDGAPRAQLSDGARAIAGTYASMQAPGPEGRSEWGILPEGHRVTLRPGLGQLRVLIEDARDFETHWTPGAAQTWTRTWRQGETLAVEVHRPATPATALADAAWDVITSIKDRTFQRELMERSNQVGMLGALLGARHSGAGDALNQLPEAHFAVSSAVVRETGREGREVDRWKAMQREATETLDELQKAATRRLAAVLSGGLR--',
}
# Describes the expected secondary structure of each residue - This only affects the look of the cartoon representation
ss_definitions = {
          #  MARAKAKD----------------QTDGIYAAFDTLMSTAGVDSQIAALAASEADAGTLDAALTQSLQEAQGRWGLGLHHLRHEARLTDD----GDIEILTDGRPSARVSEGFGALAQAYAPMQALDERGLSQWAALGEGYRAPGDLPLAQLKVLIEHARDFETDWSAGRGETFQRVWRKGDTLFVEVARPASAEAALSDAAWDVIASIKDRAFQRELMRRSEKDGMLGALLGARHAGAKANLAQLPEAHFTVQAFVQTLSGAAARNAEEYRAALKTAAAALEEYQGVTTRQLSEVLRHGLRES
    'rad' : 'LLLLLLLL----------------LLLLHHHHHHHHHHHHLLLLLHHHHHHHLLLHHHHHHHHHHHHHHHHHHHLLLLLLLLSSSSSSLL----SSSSSSSLLSSSSSHHHHHHHHHHHHHHHLLLLLLLLLLLLLLLLLLLLLLLLLHHHHHHHHHHLLLSSSLLSSSSLLSSSSSSSSLLSSSSSSSLLLLHHHHHHHHHHHHHHHHLLHHHHHHHHHHHHHHHHHHHHHLLLHHHHHHHHHLLHHHHLLLLSSSSSSSLHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHLL',
          #  MTKTKQKDRNALQESPRPNVSGPVGSEDVLKSFDALMATADVDSQIHALAESGADEETLGRELTLALQLAQDRWGLGLLHLRHDAALARTPEGTPDVVLRADGAVVARLSDGPAAIARSYASMQALGAEGLSEWGVLPDGHRVTLKGGSGQLRVLVEDARDFETHWTAERGGVWSRTWRQGETLVVEVHRPASPATVLADAAWDVITSIKDRNFQRELMERSNSVGMLGALLGARHSGAGSALDRLPSAHFTVRSAVIRESGVSARSLERWKAMLREGMEQLEALQKTVTRELAEVLSHGLR--
    'geo' : 'LLLLLLLLLLLLLLLLLLLLLLLLLLLLHHHHHHHHHHHHLLLLLHHHHHHHLLLHHHHHHHHHHHHHHHHHHHLLLLLLLLSSSSSSLLLLLLSSSSSSSLLSSSSSHHHHHHHHHHHHHHHLLLLLLLLLLLLLLLLLLLLLLLLLHHHHHHHHHHLLLSSSLLSSSSLLSSSSSSSSLLSSSSSSSLLLLHHHHHHHHHHHHHHHHLLHHHHHHHHHHHHHHHHHHHHHLLLHHHHHHHHHLLHHHHLLLLSSSSSSSLHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH--',
          #  MTKASKKSDA-----P-RTREAPTPREDALRGFDALMATAGVESTIVKHAASGADSQTLNDELTRSLQLAHDRWGLGLLHLRHEARLDRGE--DTDVILLVDGREVARLSQGAAAISATYETMRAQNADDLSDWGVLPEGHRVTLKAGNNQMRVLVEDARDFETHWSSERGGAFVRTWRQGETLAVEVHRPASPGTALADAAWDAIMSIKDRNFQRELMERSNSVGMLGALLGARHKDAGRALERLPEAHFAVRSTVVRMTGGAQREFDQWRSMVREGLDQLDELQKTTTRHLTEILRHGLK--
    'per' : 'LLLLLLLLLL-----L-LLLLLLLLLLLHHHHHHHHHHHHLLLLLHHHHHHHLLLHHHHHHHHHHHHHHHHHHHLLLLLLLLSSSSSSLLL--LSSSSSSSLLSSSSSHHHHHHHHHHHHHHHLLLLLLLLLLLLLLLLLLLLLLLLLHHHHHHHHHHLLLSSSLLSSSSLLSSSSSSSSLLSSSSSSSLLLLHHHHHHHHHHHHHHHHLLHHHHHHHHHHHHHHHHHHHHHHHHLLLLLLHHHHLHHHHLLLLSSSSSSSLHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH--',
          #  MTKTTRRKAT-----PT---EPVTSSVNPLARFAELVATAGLQSDVQALADSGADDTTLEAQLTQELRLAHDRWGLGLLHLQHSARLIHTDGVPSDIALLVDGAPRAQLSDGARAIAGTYASMQAPGPEGRSEWGILPEGHRVTLRPGLGQLRVLIEDARDFETHWTPGAAQTWTRTWRQGETLAVEVHRPATPATALADAAWDVITSIKDRTFQRELMERSNQVGMLGALLGARHSGAGDALNQLPEAHFAVSSAVVRETGREGREVDRWKAMQREATETLDELQKAATRRLAAVLSGGLR--
    'des' : 'LLLLLLLLLL-----LL---LLLLLLLLHHHHHHHHHHHHLLLLHHHHHHHHLLLHHHHHHHHHHHHHHHHHHHLLLLLLLLSSSSSSSLLLLLSSSSSSSLLSSSSSHHHHHHHHHHHHHHHLLLLLLLLLLLLLLLLLLLLLLLLLHHHHHHHHHHLLLSSSLLSSSSLLSSSSSSSSLLSSSSSSSLLLLHHHHHHHHHHHHHHHHLLHHHHHHHHHHHHHHHHHHHHHHHLLLHHHHHHHHLHHHHLLLLSSSSSSSLHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH--',
}
# Protein domain boundaries
domain_boundaries = {
    'rad' : {
             'P' : [(  1, 111)],
             'S' : [(112, 172), (233, 268)],
             'F' : [(173, 232), (269, 284)],
            },
    'geo' : {
             'P' : [( 25, 131)],
             'S' : [(132, 192), (253, 288)],
             'F' : [(193, 252), (289, 302)],
            },
    'per' : {
             'P' : [( 19, 123)],
             'S' : [(124, 184), (245, 280)],
             'F' : [(185, 244), (281, 294)],
            },
    'des' : {
             'P' : [( 17, 123)],
             'S' : [(124, 184), (245, 280)],
             'F' : [(185, 244), (281, 294)],
            },
}
visible_ranges = {
    ('6o5l', 'A') : [(21, 291)],
    ('6bdu', 'A') : [(12, 280)],
    ('6bdu', 'B') : [(12, 280)],
    ('6mc8', 'A') : [(20, 294)],
    ('6mc8', 'B') : [(21, 294)],
    ('6mc6', 'A') : [(12, 281)],
    ('6mc6', 'B') : [(11, 280)],
    ('6neo', 'B') : [(14, 213), (222, 279)],
    ('6a27', 'A') : [(12, 172), (226, 266)],
    ('6a27', 'B') : [(12, 171), (227, 266)],
    ('6a28', 'A') : [(10, 189), (221, 278)],
    ('6a28', 'B') : [(10, 185), (226, 278)],
    ('6a29', 'A') : [(10, 284)],
    ('6a29', 'B') : [(10, 284)],
    ('6a29', 'C') : [(10, 284)],
    ('6a29', 'D') : [(10, 284)],
    ('6a29', 'E') : [( 7, 284)],
    ('6a29', 'F') : [(10, 284)],
    ('6a29', 'G') : [(10, 284)],
    ('6a29', 'H') : [(10, 284)],
    ('9om8', 'A') : [(11, 215), (223, 281)],
    ('9om8', 'B') : [(11, 213), (222, 279)],
    ('9om8', 'C') : [(11, 204), (224, 278)],
    ('9om8', 'D') : [(11, 188), (223, 273)],
    ('9or6', 'A') : [(11, 281)],
    ('9or6', 'B') : [(11, 214), (223, 278)],
    ('9yi3', 'A') : [( 10, 188), (193, 213), (222, 277)],
    ('9yi3', 'B') : [( 12, 186), (226, 276)],
    ('9yl4', 'A') : [( 9, 284)],
    ('9yl4', 'B') : [( 9, 284)],
    ('9yl4', 'C') : [( 9, 284)],
    ('9yl4', 'D') : [( 9, 284)],
    ('9yl4', 'E') : [( 9, 284)],
    ('9yl4', 'F') : [( 9, 284)],
    ('9yl4', 'G') : [(10, 284)],
    ('9yl4', 'H') : [( 9, 284)],
    ('9yl4', 'I') : [(10, 284)],
    ('9yl4', 'J') : [( 9, 284)],
    ('9yl4', 'K') : [(10, 280)],
    ('9yl4', 'L') : [( 9, 284)],
    ('9yup', 'A') : [( 9, 284)],
    ('9yup', 'B') : [( 9, 212), (224, 281)],
    ('9yup', 'C') : [(12, 284)],
    ('9yup', 'D') : [(11, 211), (224, 279)],
    ('9y6e', 'A') : [(12, 280)],
    ('9y6e', 'B') : [(12, 280)],
    ('9y6e', 'C') : [(12, 280)],
    ('9y6e', 'D') : [(12, 280)],
}
#visible_ranges=[[1,284]]

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


########################
# Function Definitions #
########################

def convert_index(input_resi: int, homolog_1: str, homolog_2: str, homolog_alignment: dict) -> int:
    """
    DESCRIPTION:
        Function to convert a residue index from one homolog to another while handling gaps correctly.

    ARGUMENTS:
        `input_resi`:    (required `int`)
                        The residue index in the input homolog (1-based indexing).
        `homolog_1`:     (required `str`)
                        The name of the input homolog.
        `homolog_2`:     (required `str`)
                        The name of the output homolog.
        `homolog_alignment`: (required `dict`)
                             Dictionary containing aligned homolog sequences with homolog names as keys.

    RETURNS:
        `int`: The corresponding residue index in homolog_2 (1-based indexing), 
               or the closest residue on homolog_2 if the input residue aligns with a gap in homolog_2 
               or `-1` if the index is out of range.
    """   
    # Retrieve the aligned sequences
    seq1 = homolog_alignment[homolog_1]
    seq2 = homolog_alignment[homolog_2]
    # Initialize counters for residue indices (1-based)
    resi1_count = 0
    resi2_count = 0
    # Iterate through the aligned sequences
    for i in range(len(seq1)):
        if seq1[i] != '-':
            resi1_count += 1
        if seq2[i] != '-':
            resi2_count += 1
        # Check if we reached the input residue index in homolog_1
        if resi1_count == input_resi:
            if seq2[i] != '-':
                return resi2_count
            # Find the nearest residue in homolog_2
            left, right = i - 1, i + 1
            while left >= 0 or right < len(seq2):
                if left >= 0 and seq2[left] != '-':
                    return sum(1 for j in range(left + 1) if seq2[j] != '-')
                if right < len(seq2) and seq2[right] != '-':
                    return sum(1 for j in range(right + 1) if seq2[j] != '-')
                left -= 1
                right += 1
    return -1  # Return -1 if input_resi is out of range


def get_ss_segments(ss_str_seq: str, ignore_gaps=True) -> dict:
    ss_segments = {}
    start_idx = None
    current_ss = None
    residue_idx = 1  # Start indexing from 1
    
    for i, ss in enumerate(ss_str_seq):
        if ignore_gaps and ss == '-':
            continue
        elif not ignore_gaps and ss == '-':
            ss = 'UNK'  # Treat gaps as an unknown secondary structure
        
        if current_ss is None:
            current_ss = ss
            start_idx = residue_idx
        elif ss != current_ss:
            ss_segments[(start_idx, residue_idx - 1)] = current_ss
            current_ss = ss
            start_idx = residue_idx
        
        residue_idx += 1
    
    if current_ss is not None:
        ss_segments[(start_idx, residue_idx - 1)] = current_ss
    
    return ss_segments



def load_interactions(csv_file):
    col_dtype = {
        'interface'   : str,
        'homolog'     : str,
        'pdb'         : str,
        'sym1'        : str,
        'segi1'       : str,
        'chain1'      : str, 
        'resi1'       : int, 
        'resn1'       : str,
        'sym2'        : str, 
        'segi2'       : str,
        'chain2'      : str, 
        'resi2'       : int, 
        'resn2'       : str, 
        'contact_type' : str,
    }
    contacts = pd.read_csv(csv_file, index_col=False, dtype=col_dtype)
    # drop unnamed column 0
    contacts = contacts.loc[:, ~contacts.columns.str.contains('^Unnamed')]
    return contacts



def process_columns(contacts, homolog_alignment):
    """
    DESCRIPTION:
        Function to convert residue indices to the residue numbering of a reference homolog.

    ARGUMENTS:
        `contacts`: (required `pd.DataFrame`)
                   Dataframe containing residue interaction data.
        `homolog_alignment`: (required `dict`)
                              Dictionary containing aligned homolog sequences.

    RETURNS:
        `pd.DataFrame`: Updated dataframe with converted residue indices.
    """
    for idx, row in contacts.iterrows():
        homolog = row['homolog']
        if homolog != 'rad':
            contacts.at[idx, 'resi1'] = convert_index(row['resi1'], homolog, ref_homolog, homolog_alignment)
            contacts.at[idx, 'resi2'] = convert_index(row['resi2'], homolog, ref_homolog, homolog_alignment)
    return contacts



def list_chain_instances(contacts:pd.DataFrame, identifier='chain', mode='both'):
    # Extract unique (pdb, chain1) combinations
    chain1_instances = list(contacts[['pdb', f'{identifier}1']].drop_duplicates().itertuples(index=False, name=None))
    # Extract unique (pdb, chain2) combinations
    chain2_instances = list(contacts[['pdb', f'{identifier}2']].drop_duplicates().itertuples(index=False, name=None))
    
    if mode == 'chain1':
        return chain1_instances
    if mode == 'chain2':
        return chain2_instances
    if mode == 'both':
        return chain1_instances + chain2_instances
    else:
        raise ValueError("unknown value for list_chain_instances mode. Values can only be 'chain1', 'chain2' or 'both'")



# white background
def draw_background(ctx:cairo.Context, rgba_color, w,h):
    ctx.set_source_rgba(*rgba_color)
    ctx.rectangle(0, 0, w, h)
    ctx.stroke_preserve()
    ctx.fill()
    ctx.set_source_rgba(*black_rgba)  # black outline color



def compute_alpha_stack(stack_height, target_alpha=1.0, threshold=0.99):
    """
    Function to calculate the color transparency that is required to achieve a target alpha value
    when the color is stacked on top of itself.
    """
    target_alpha = target_alpha*threshold
    return 1 - (1 - target_alpha) ** (1 / stack_height)



# Fill in domain color only for sections of that are visible in the crystal structure
# as defined in visible_ranges dict
def draw_domain(ctx:cairo.Context, x, y, width, height, domain, show_num=False, chain_instances=[]):
    for domain_segment in domain_boundaries[ref_homolog][domain]:
        x_domain = x + width*((domain_segment[0]-1)/prot_length)
        width_domain = width * ((1+ domain_segment[1] - domain_segment[0])/prot_length)
        transparency =  compute_alpha_stack(len(chain_instances), target_alpha=1.0)#1 - 0.075**(1/len(chain_instances)) #1/len(chain_instances)
        for chain_instance in chain_instances:
            vrange_list = visible_ranges[chain_instance]
            for vrange in vrange_list:
                if (vrange[0] < domain_segment[1]) and (vrange[1] > domain_segment[0]):
                    if vrange[0] < domain_segment[0]:
                        vrange = [domain_segment[0], vrange[1]]
                    if vrange[1] > domain_segment[1]:
                        vrange = [vrange[0], domain_segment[1]]
                    #print(f'\tVisible: {vrange[0]}-{vrange[1]}')
                    x_vrange = x + width*((vrange[0]-1)/prot_length)
                    width_vrange = width * ((1 + vrange[1] - vrange[0])/prot_length)
                    ctx.rectangle(x_vrange, y, width_vrange, height)
                    rgba = [domain_colors[domain][0], domain_colors[domain][1], domain_colors[domain][2], transparency]
                    ctx.set_source_rgba(*rgba)
                    ctx.fill()

        if show_num:
            text = str(domain_segment[1])
            x_text = x_domain + width_domain - ctx.text_extents(text)[4] - text_pad
            y_text = y + height/2 + ctx.text_extents(text)[3]/2
            ctx.move_to(x_text, y_text)
            ctx.set_source_rgba(*black_rgba)
            ctx.show_text(text)
        ctx.rectangle(x_domain, y, width_domain, height)
        ctx.set_source_rgba(*black_rgba)
        ctx.set_line_width(outline_width)
        ctx.stroke()


def roundrect(ctx: cairo.Context, x, y, width, height, r):
    ctx.move_to(x, y+r)
    ctx.arc(x+r, y+r, r, math.pi, 3*math.pi/2)
    ctx.arc(x+width-r, y+r, r, 3*math.pi/2, 0)
    ctx.arc(x+width-r, y+height-r, r, 0, math.pi/2)
    ctx.arc(x+r, y+height-r, r, math.pi/2, math.pi)
    ctx.close_path()



def draw_arrow(ctx: cairo.Context, x_start, y_start, arrow_length, arrow_width, tip_width):
    body_length = arrow_length - tip_width  # Length of the arrow body
    body_thickness = arrow_width / 2  # Thickness of the arrow body
    half_width = arrow_width / 2
    body_top = y_start + (arrow_width - body_thickness) / 2
    body_bottom = body_top + body_thickness
    
    ctx.move_to(x_start, body_top)  # Start at top-left of the body
    ctx.line_to(x_start + body_length, body_top)  # Draw body to the right
    ctx.line_to(x_start + body_length, y_start)  # Move up to arrow tip base
    ctx.line_to(x_start + arrow_length, y_start + half_width)  # Tip of arrow
    ctx.line_to(x_start + body_length, y_start + arrow_width)  # Move down to arrow tip base
    ctx.line_to(x_start + body_length, body_bottom)  # Move back to body bottom
    ctx.line_to(x_start, body_bottom)  # Draw body back to start
    ctx.close_path()
    


def draw_sec_struc(ctx:cairo.Context, x, y, width, height):
    ss_segments = get_ss_segments(ss_definitions[ref_homolog], ignore_gaps=True)
    for segment_range, ss in ss_segments.items():
        x_segment = x + width*((segment_range[0]-1)/prot_length)
        width_segment = width * ((1+ segment_range[1] - segment_range[0])/prot_length)
        if ss == 'H':
            y_segment = y + height/4
            height_segment = height/2
            helix_freq = width*3/prot_length
            curved_width = width*1.5/prot_length
            #draw_spring(ctx, x_segment, y_segment, width_segment, height_segment, helix_freq)
            roundrect(ctx, x_segment, y_segment, width_segment, height_segment, curved_width)
            ctx.set_source_rgba(0,0,0,0.5)
            ctx.set_line_width(outline_width)
            ctx.stroke()
        if ss == 'S':
            y_segment = y + height/4
            height_segment = height/2
            tip_width = width*2/prot_length
            draw_arrow(ctx, x_segment, y_segment, width_segment, height_segment, tip_width)
            ctx.set_source_rgba(0,0,0,0.5)
            ctx.set_line_width(outline_width)
            ctx.stroke()
            #ctx.fill()   
    return



# Draw a rectangle representing the protein chain - colored by domain
def draw_chain(ctx:cairo.Context, x, y, width, height, chain_instances=[]):
    # Draw P, S and F domains
    draw_domain(ctx, x, y, width, height, 'P', chain_instances=chain_instances, show_num=False)
    draw_domain(ctx, x, y, width, height, 'S', chain_instances=chain_instances, show_num=False)
    draw_domain(ctx, x, y, width, height, 'F', chain_instances=chain_instances, show_num=False)
    draw_sec_struc(ctx, x, y, width, height)


def draw_intra_contacts(intra_contacts, output_file):

    def draw_intra_contact(row):
        """
        Function called by pd.DataFrame.apply to draw intramolecular residue contacts
        Function is called once for each row in the interactions dataframe
        Each intramolecular interaction is represented by an arced line connecting two residue indices within the same rectangle
        """
        resi1 = row['resi1']
        resi2 = row['resi2']
        contact_type = row['contact_type']
        transparency = compute_alpha_stack(top_interaction_count, target_alpha=contact_rgba[contact_type][3]) #contact_rgba[contact_type][3] - 0.075**(1/top_interaction_count)
        color = [contact_rgba[contact_type][0], contact_rgba[contact_type][1], contact_rgba[contact_type][2], transparency]
        dash = contact_dashes[contact_type]
        x1 = x0 + (resi1 / prot_length) * chain_width
        x2 = x0 + (resi2 / prot_length) * chain_width
        
        # Calculate distance between two points
        d = abs(x2 - x1)
        # Calculate the radius of the circular arc connecting the points
        r = d/(2*math.sin(math.radians(arc_angle/2)))
        # Calculate length of h (line from center of circle to centerpoint of d)
        arc_depth = math.sqrt( (r**2) - (d/2)**2  )
        # Calculate x position of arc center
        x_c = d/2 + min([x1, x2])
        # Calculate y position of arc center and start/end angles of the arc (in radians)
        angle_start = (180 + (180-arc_angle)/2) * (math.pi/180)
        angle_end = (360 - (180-arc_angle)/2) * (math.pi/180)
        y_c = y0 - outline_width/2 + arc_depth
        
        ctx.set_source_rgba(*color)
        ctx.set_line_width(interaction_width)
        ctx.set_dash(dash)  # set line pattern
        ctx.arc(x_c, y_c, r, angle_start, angle_end)
        ctx.stroke()
        ctx.set_line_width(outline_width)
        ctx.set_dash([1,0])  # set solid line pattern
        ctx.set_source_rgba(*black_rgba)
        return

    # Calculate chain width (in px)
    chain_width = w - (2*border_width)
    # Find largest intra-chain distance
    d_max = max((abs(intra_contacts['resi1'] - intra_contacts['resi2'])/prot_length) * chain_width)
    # Calculate largest arc:
    r = d_max/(2*math.sin(math.radians(arc_angle/2)))
    arc_depth = math.sqrt( (r**2) - (d_max/2)**2  )
    arc_height = r - arc_depth
    # calculate height of figure
    h = 2*border_width + arc_height + chain_height

    # Set chain parameters
    x0 = border_width
    y0 = h - chain_height - border_width

    # Create a new surface (canvas) to draw on
    surface = cairo.SVGSurface(output_file, w, h)
    # Setting SVG unit
    surface.set_document_unit(document_unit)
    # Create a new context (pen) to draw with
    ctx = cairo.Context(surface)
    # Set pen width
    ctx.set_line_width(outline_width)
    # Set Font
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    ctx.set_font_size(font_size)

    # Find the most frequent interaction across crystal structures
    top_interaction = intra_contacts.groupby(['resi1', 'resi2','contact_type']).size().idxmax()
    top_interaction_count = intra_contacts.groupby(['resi1', 'resi2','contact_type']).size().max()

    # draw background
    if not transparent_bg:
        draw_background(ctx, white_rgba,w,h)

    # Draw Intra-molecular contacts
    for contact_type in list(contact_rgba.keys()):
        intra_contacts[(intra_contacts['contact_type'] == contact_type)].apply( draw_intra_contact, axis=1)
    
    # Draw chain
    chain_instances = list_chain_instances(intra_contacts, identifier='chain', mode='both')
    draw_chain(ctx, x0, y0, chain_width, chain_height, chain_instances=chain_instances)



def draw_inter_contacts(inter_contacts, output_file):
        
    def draw_inter_contact(row):
        """
        Function called by pd.DataFrame.apply to draw intermolecular residue contacts
        Function is called once for each row in the interactions dataframe
        Each intermolecular interaction is represented by a straight line connecting the two rectangles at positions corresponding to the residue indices
        """
        if row['segi1'] == 'A':
            a = row['resi1']
            b = row['resi2']
        else:
            a = row['resi2']
            b = row['resi1']
        a_x = a_x0 + (a / prot_length) * chain_width
        a_y = a_y0 + chain_height #+ outline_width/2
        b_x = b_x0 + (b / prot_length) * chain_width
        b_y = b_y0 #- outline_width/2
        contact_type = row['contact_type']
        transparency = compute_alpha_stack(top_interaction_count, target_alpha=contact_rgba[contact_type][3]) #contact_rgba[contact_type][3] - 0.075**(1/top_interaction_count)
        color = [contact_rgba[contact_type][0], contact_rgba[contact_type][1], contact_rgba[contact_type][2], transparency]
        dash = contact_dashes[contact_type]
        ctx.set_source_rgba(*color)
        ctx.set_line_width(interaction_width)
        ctx.set_dash(dash)  # dashed line pattern
        ctx.move_to(a_x, a_y)
        ctx.line_to(b_x, b_y)
        ctx.stroke()
        ctx.set_line_width(outline_width)
        ctx.set_dash([1,0])  # set solid line pattern
        ctx.set_source_rgba(*black_rgba)
        return
    
    # Calculate figure height
    h = 2*border_width + 2*chain_height + interchain_space
    # Calculate chain width (in px)
    chain_width = w - (2*border_width)
    # Set chain A parameters
    a_x0 = border_width
    a_y0 = h/2 - 0.5*interchain_space - chain_height
    # Set chain B parameters
    b_x0 = a_x0
    b_y0 = h/2 + 0.5*interchain_space

    # Create a new surface (canvas) to draw on
    surface = cairo.SVGSurface(output_file, w, h)
    # Setting SVG unit
    surface.set_document_unit(document_unit)
    # Create a new context (pen) to draw with
    ctx = cairo.Context(surface)
    # Set pen width
    ctx.set_line_width(outline_width)
    # Set Font
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    ctx.set_font_size(font_size)
    
    # Find the most frequent interaction across crystal structures
    top_interaction = inter_contacts.groupby(['resi1', 'resi2','contact_type']).size().idxmax()
    top_interaction_count = inter_contacts.groupby(['resi1', 'resi2','contact_type']).size().max()

    # draw background
    if not transparent_bg:
        draw_background(ctx, white_rgba,w,h)

    # Draw Inter-molecular contacts
    for contact_type in list(contact_rgba.keys()):
        inter_contacts[(inter_contacts['contact_type'] == contact_type)].apply( draw_inter_contact, axis=1 )
    
    # Figure out which crystal structures were used used for this intermolecular interaction

    # Draw chain
    chain1_instances = list_chain_instances(inter_contacts, identifier='chain', mode='chain1')
    chain2_instances = list_chain_instances(inter_contacts, identifier='chain', mode='chain2')
    draw_chain(ctx, a_x0, a_y0, chain_width, chain_height, chain_instances=chain1_instances)
    draw_chain(ctx, b_x0, b_y0, chain_width, chain_height, chain_instances=chain2_instances)
    return



def main():
    parser = argparse.ArgumentParser(description='Draw PprA contact schematics.')
    parser.add_argument('-i', '--input', dest='input_csv', required=True, help='Path to the input CSV file')
    parser.add_argument('-o', '--output', dest='output_dir', default='./', help='Path to the output directory (default: ./)')
    args = parser.parse_args()

    input_csv = args.input_csv
    output_dir = os.path.abspath(args.output_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # load contacts that were output from pymol script
    contacts = load_interactions(input_csv)
    # convert all residue indices to match numbering of D.radiodurans PprA
    contacts = process_columns(contacts, homolog_alignment)

    # split contacts
    contacts_M = contacts[(contacts['interface'] == 'M')]# & (contacts['pdb'] == '6o5l')]
    contacts_S = contacts[(contacts['interface'] == 'S')]
    contacts_F = contacts[(contacts['interface'] == 'F')]
    contacts_P = contacts[(contacts['interface'] == 'P')]

    intra_contacts_filepath = os.path.join(output_dir, 'intra_contacts.svg')
    s_contacts_filepath = os.path.join(output_dir, 'S_contacts.svg')
    f_contacts_filepath = os.path.join(output_dir, 'F_contacts.svg')
    p_contacts_filepath = os.path.join(output_dir, 'P_contacts.svg')
    print(f"Generating {intra_contacts_filepath}")
    draw_intra_contacts(contacts_M, intra_contacts_filepath)
    print(f"Generating {s_contacts_filepath}")
    draw_inter_contacts(contacts_S, s_contacts_filepath)
    print(f"Generating {f_contacts_filepath}")
    draw_inter_contacts(contacts_F, f_contacts_filepath)
    print(f"Generating {p_contacts_filepath}")
    draw_inter_contacts(contacts_P, p_contacts_filepath)



main()

