#!/usr/bin/env python

import xml.etree.ElementTree as xml
import sys
import os
import re
import tkSimpleDialog
import tkMessageBox
import urllib2
import BeautifulSoup
import pymol
from pymol import cmd, setting

###Constants
##PDBTM
DB_ID_PDBTM="PDBTM"
NS_PDBTM="{http://pdbtm.enzim.hu}"
E_CHAIN="CHAIN"
E_SIDEDEF="SIDEDEFINITION"
E_REG="REGION"
A_CHAINID="CHAINID"
A_CHAINTYPE="TYPE"
##OPM
DB_ID_OPM="OPM"
OPM_SUBUNITS_WEBSITE="http://opm.phar.umich.edu/subunits.php"

##Color Coding
SEGMENT_HIGHLIGHTS = {
        'Helix': 'yellow', 'Unknown': 'green', 'Interfacial_helix': 'forest',
        'Side_1': 'red', 'Side_2' : 'blue', 'Membrane_loop': 'orange',
        'Outside': 'blue', 'Inside': 'red',
        'Beta-strand': 'yellow', 'Membrane_inside': 'deeppurple', 'Coil': 'yellow'
        }

##Nomenclature
SEGMENT_LABELS = {
        'H': 'Helix', 'U': 'Unknown', 'F': 'Interfacial_helix',
        '1': 'Side_1', '2' : 'Side_2', 'L': 'Membrane_loop',
        'I': 'Membrane_inside', 'C': 'Coil', 'B': 'Beta-strand'
        }
REVERSE_SIDES = { 'Outside': 'Inside', 'Inside': 'Outside' }

def __init__(self):
   self.menuBar.addmenuitem('Plugin', 'command',
                            'Membrane Protein Highlight',
                            label = 'Membrane Protein Highlight',
                            command = lambda s=self : highlight_membrane_dialog(s))

def highlight_membrane_dialog(app):
   pdbCode = tkSimpleDialog.askstring('Membrane Protein Highlight',
                                      'Please enter a 4-digit pdb code:',
                                      parent=app.root)
   if pdbCode is None:
       print 'No PDB code supplied'
   else:
       #For now, assume the protein is not loaded yet
       highlight_membrane(pdbCode, 0)
   
def highlight_membrane(pdbCode, loaded=0, db=DB_ID_PDBTM):
    loaded = int(loaded)
    if not loaded:
        cmd.fetch(pdbCode[0:4])
        
    if db==DB_ID_PDBTM:
        xml = get_pdbtm_xml(pdbCode[0:4])
        chains_dict = get_pdbtm_annotation(pdbCode[0:4], xml)
    else:
        if db==DB_ID_OPM:
            opm_file = get_opm_file(pdbCode[0:4])
            chains_dict = get_opm_annotation(opm_file)
        else:
            print 'Unkown database supplied: %s. Exiting' % (db)
            exit(1)
    highlight_molecule(chains_dict, pdbCode[0:4].lower(), loaded, pdbCode[4:5])
   
def main(sys_argv=sys.argv):
    #pdbCode = '1XFH'.lower()
    pdbCode = '2NWL'.lower()
    pymol.finish_launching()
    cmd.fetch(pdbCode)
    xml = get_pdbtm_xml(pdbCode)
    chains_dict = get_pdbtm_annotation(pdbCode, xml)
    highlight_molecule(chains_dict, pdbCode.lower())
 
def get_pdbtm_xml(arg_pdbid):
    """
    Make sure the XML file for the given pdbid is present and return its path
    Search order: 1:Locally in fetch_path, 2:PDBTM Remote
    """
    
    path = setting.get('fetch_path')    #Path also used by the 'fetch' command. If not set in pymolrc, it's probably ~
    localfn = os.path.join(path, arg_pdbid.upper() + '.xml')
    if not os.access(localfn, os.R_OK):
        #Retrieve file from remote
        try:
            xmlfile = urllib2.urlopen('http://pdbtm.enzim.hu/data/database/' + arg_pdbid[1:3] + '/' + arg_pdbid.lower() + '.xml')
            output = open(localfn,'wb')
            output.write(xmlfile.read())
            output.close()
        except:
            print 'Error during accession or retrieval of XML file' , sys.exc_info()[0]
            exit(1)
    return localfn
    
def get_opm_file(arg_pdbid):
    """
    Ensure the flatfile with annotations from OPM for the given PDB id is present
    and return the path to it.
    Search order: 1:Locally in fetch_path, 2:OPM Remote
    """
    
    path = setting.get('fetch_path')    #Path also used by the 'fetch' command. If not set in pymolrc, it's probably ~
    localfn = os.path.join(path, arg_pdbid.upper() + '.opm')
    if not os.access(localfn, os.R_OK):
        #Retrieve the file from OPM website
        try:
            #Retrieve all segments at once
            opm_segments_urlobj = urllib2.urlopen(OPM_SUBUNITS_WEBSITE)
            opm_segments_html = opm_segments_urlobj.read()
            #Get the relevant part from that site
            soup = BeautifulSoup.BeautifulSoup(opm_segments_html)
            opm_segments_div = str(soup.find("div", {"id": "body"}))
            
            if opm_segments_div is not None:
                print "Parsing OPM annotations"
                opm_segments = opm_segments_div.split('<br />')
                pdb_id = ''
                arg_pdbid_lines = []
                for segment in opm_segments:
                    #print 'Parsing %s' % (segment)
                    match = re.search('(\w{4}) <b>(\w)</b>', segment)
                    if match:
                        pdb_id = match.group(1)
                        if pdb_id == arg_pdbid.lower():
                            arg_pdbid_lines.append(segment)
                #Write all relevant lines to output file
                if len(arg_pdbid_lines) > 0:
                    output = open(localfn,'wb')
                    output.write('\n'.join(arg_pdbid_lines) + '\n')
                    output.close()
                else:
                    print 'Could not find any annotation for pdbid %s in OPM' % (arg_pdbid)
                    exit(1)
            else:
                print 'Error during accession or retrieval of OPM file. Cannot find segment annotations in page.'
                exit(1)
        except:
            print 'Error during accession or retrieval of OPM file'
            raise
    else:
        print 'Found local file'
    return localfn
    
def get_opm_annotation(arg_file):
    """
    Parse a OPM file to 3D structure annotation
    """
    chains_dict = {}
    
    with open(arg_file) as f:
        for line in f:
            match = re.match('(\w{4})\s+<b>(\w)</b>.+Segments:(.+)', line)
            if match:
                region_dict = {}
                pdbid = match.group(1)
                chain_id = match.group(2)
                segments = match.group(3)
                #Iterate all segment annotations
                for segment in re.finditer("(\d+)\s*\(\s*(\d+\s*-\s*\d+)\s*\)",segments):
                    #Too much whitespace will trip up PyMOL, so get rid of it
                    region_interval = re.sub('\s', '', segment.group(2))
                    #Use helix numbering from OPM (provided on group 1) and interval (at group 2)
                    region_dict[chain_id + '_' + 'Helix' + segment.group(1) + '_' + pdbid + '_' + DB_ID_OPM] = [SEGMENT_HIGHLIGHTS['Helix'], region_interval]
                chains_dict[chain_id] = region_dict
    
    return chains_dict

def get_pdbtm_annotation(arg_pdbid, arg_xml):
    """
    Parse a PDBTM XML-File to 3D-structure annotation for given PDB-ID
    """
    
    #Read xml file from PDBTM
    tree = xml.parse(arg_xml)
    root = tree.getroot()

    #Check if entry describes a membrane protein (Actually this should never happen since for non IMPs there is no xml to begin with)
    is_imp = root.attrib.get("TMP")
    if not is_imp == "yes":
        print "Unexpected error:", sys.exc_info()[0]
        tkMessageBox.showerror('Invalid Code',
                              'The PDB codes does not describe a membrane protein in PDBTM:' + arg_pdbid)
        
    #xml.dump(rootElem)
    #print list(rootElem)

    #Check whether a Side definition cross-reference to TOPDB is present
    #If so, use it to rename Side1/2 to Inside/Outside
    side_def = root.find(NS_PDBTM+E_SIDEDEF)
    if side_def is not None:
         SEGMENT_LABELS['1'] = side_def.attrib['Side1']
         SEGMENT_LABELS['2'] = REVERSE_SIDES[SEGMENT_LABELS['1']]

    #Get a list of all chain elements
    chains_list = root.findall(NS_PDBTM+E_CHAIN)
    if chains_list != None:
        
        chains_dict = {}
        
        for chain in chains_list:
            chain_id = chain.attrib.get(A_CHAINID)
            chain_type = chain.attrib.get(A_CHAINTYPE)
            
            if not (chain_type == "alpha" or chain_type == "beta"): 
                continue
            
            region_dict = {}
            helix_count = 0
            
            region_list = chain.findall(NS_PDBTM+E_REG)
            for region in region_list:
                region_atom_begin = region.attrib.get("pdb_beg")
                region_atom_end = region.attrib.get("pdb_end")
                region_type = region.attrib.get("type")
                region_label = SEGMENT_LABELS[parse_pdbtm_annotation(region_type)]
                
                #If we are processing a helix, add an entry with Helix\d, one for every helix
                if region_label == 'Helix':
                    region_dict[chain_id + '_' + region_label + str(helix_count) + '_' + arg_pdbid + '_' + DB_ID_PDBTM] = [SEGMENT_HIGHLIGHTS[region_label] , region_atom_begin + '-' + region_atom_end]
                    helix_count += 1
                #In all other casesm e.g. Side1 pool all the regions and make one big selection for all residues in Side1
                else:
                    if chain_id + '_' + region_label + '_' + arg_pdbid + '_' + DB_ID_PDBTM in region_dict:
                        region_dict[chain_id + '_' + region_label + '_' + arg_pdbid + '_' + DB_ID_PDBTM] = [SEGMENT_HIGHLIGHTS[region_label] , region_dict[chain_id + '_' + region_label + '_' + arg_pdbid + '_' + DB_ID_PDBTM][1] + '+' + region_atom_begin + '-' + region_atom_end]
                    else:
                        region_dict[chain_id + '_' + region_label + '_' + arg_pdbid + '_' + DB_ID_PDBTM] = [SEGMENT_HIGHLIGHTS[region_label] , '+' + region_atom_begin + '-' + region_atom_end]
                        
            chains_dict[chain_id] = region_dict
    else:
        print 'Exiting. No chains found for PDB-ID: ' + arg_pdbid
        exit(0)
        
    return chains_dict
    
def highlight_molecule(chains_dict, pdbCode, loaded=0, desiredChain=''):
    if not loaded:
        cmd.color('white', pdbCode)
    for chain in chains_dict:
        #If no desired chains is given, handle all, otherwise only the specified one
        if desiredChain == '' or desiredChain == chain: 
            for label, (col, region) in sorted(chains_dict[chain].items()):
                #Below statement uses Pymol's atom selection macro syntax: http://pymol.sourceforge.net/newman/user/S0220commands.html
                if cmd.count_atoms('/%s//%s/%s' % ( pdbCode, chain, region)) > 0:
                    cmd.select(label, '/%s//%s/%s' % ( pdbCode, chain, region))
                    cmd.color(col, label)
    
def parse_pdbtm_annotation(arg_type):
    """
    Parse the PDBTM type annotation to the common standard used for all databases.
    In this case there is nothing to do because we parse all other sources to the way PDBTM annotates
    """
    return arg_type   
    
if __name__=='__main__':
    main()
            
#Use hmem 1XFH to highlight structure 1XFH from pymol commandline
cmd.extend("hmem", highlight_membrane)
