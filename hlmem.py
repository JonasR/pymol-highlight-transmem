#!/usr/bin/env python

import xml.etree.ElementTree as xml
import sys
import os
import tkSimpleDialog
import tkMessageBox
import urllib2
from pymol import cmd, setting

###Constants
##PDBTM
NS_PDBTM="{http://pdbtm.enzim.hu}"
E_CHAIN="CHAIN"
E_REG="REGION"
A_CHAINID="CHAINID"
A_CHAINTYPE="TYPE"

##Color Coding
SEGMENT_HIGHLIGHTS = {'H': 'yellow', 'U': 'green', 'F': 'forest', '1': 'red', '2' : 'blue', 'L': 'orange' }
SEGMENT_LABELS = {'H': 'Helix', 'U': 'Unknown', 'F': 'Interfacial_helix', '1': 'Side_1', '2' : 'Side_2', 'L': 'Membrane_loop' }

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
       exit(0)
   
   #For now, assume the protein is not loaded yet
   highlight_membrane(pdbCode, 0)
   
def highlight_membrane(pdbCode, loaded=0):
    if not loaded:
        cmd.fetch(pdbCode)
    xml = get_pdbtm_xml(pdbCode)
    chains_dict = get_pdbtm_annotation(pdbCode, xml)
    highlight_molecule(chains_dict, pdbCode.lower(), loaded)
   
def main(sys_argv=sys.argv):
    #This will crash at some point, since pymol isn't running
    pdbCode = '1XFH'
    cmd.fetch(pdbCode)
    xml = get_pdbtm_xml(pdbCode)
    chains_dict = get_pdbtm_annotation(pdbCode, xml)
    highlight_molecule(chains_dict, pdbCode.lower(),0)
 
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
            xmlfile = urllib2.urlopen('http://pdbtm_data.enzim.hu/database/' + arg_pdbid[1:3] + '/' + arg_pdbid.lower() + '.xml')
            output = open(localfn,'wb')
            output.write(xmlfile.read())
            output.close()
        except:
            print 'Error during accession or retrieval of XML file' , sys.exc_info()[0]
            exit(1)
    return localfn
    
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

    #Get a list of all chain elements
    chains_list = root.findall(NS_PDBTM+E_CHAIN)
    if chains_list != None:
        
        chains_dict = {}
        
        for chain in chains_list:
            chain_id = chain.attrib.get(A_CHAINID)
            chain_type = chain.attrib.get(A_CHAINTYPE)
            
            if not chain_type == "alpha":   #At some later point add beta (maybe)
                continue
            
            region_dict = {}
            helix_count = 0
            
            region_list = chain.findall(NS_PDBTM+E_REG)
            for region in region_list:
                region_atom_begin = region.attrib.get("pdb_beg")
                region_atom_end = region.attrib.get("pdb_end")
                region_type = region.attrib.get("type")
                region_label = parse_pdbtm_annotation(region_type)
                
                #If we are processing a helix, add an entry with Helix\d, one for every helix
                if region_label == 'H':
                    region_dict[chain_id + '_' + SEGMENT_LABELS[region_label] + str(helix_count) + '_' + arg_pdbid] = [SEGMENT_HIGHLIGHTS[region_label] , region_atom_begin + '-' + region_atom_end]
                    helix_count += 1
                #In all other casesm e.g. Side1 pool all the regions and make one big selection for all residues in Side1
                else:
                    if chain_id + '_' + SEGMENT_LABELS[region_label] in region_dict:
                        region_dict[chain_id + '_' + SEGMENT_LABELS[region_label] + '_' + arg_pdbid] = [SEGMENT_HIGHLIGHTS[region_label] , region_dict[chain_id + '_' + SEGMENT_LABELS[region_label] + '_' + arg_pdbid][1] + '+' + region_atom_begin + '-' + region_atom_end]
                    else:
                        region_dict[chain_id + '_' + SEGMENT_LABELS[region_label] + '_' + arg_pdbid] = [SEGMENT_HIGHLIGHTS[region_label] , '+' + region_atom_begin + '-' + region_atom_end]
                        
            chains_dict[chain_id] = region_dict
    else:
        print 'Exiting. No chains found for PDB-ID: ' + arg_pdbid
        exit(0)
        
    return chains_dict
    
def highlight_molecule(chains_dict, pdbCode, loaded):
    if not loaded:
        cmd.do('color white, %s' % (pdbCode))
    for chain in chains_dict:
        for label, (col, region) in sorted(chains_dict[chain].items()):
            #Below statement uses Pymol's atom selection macro syntax: http://pymol.sourceforge.net/newman/user/S0220commands.html
            cmd.do('select %s, /%s//%s/%s' % (label, pdbCode, chain, region))
            cmd.do('color %s, %s' %(col, label))

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