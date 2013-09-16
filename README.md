pymol-highlight-transmem
========================
A pymol script to highlight transmembrane segments and related structural elements on a given PDB structure.
Currently, annotations are retrieved only from [PDBTM][1] and this script essentially reproduces the behaviour of the viewer on their website.

Install
-------
* For current session only, execute in Pymol CLI:

        import hlmem.py
* For permanent availability:
    * Start Pymol and go to Plugin -> Manage Plugins -> Install, then choose hlmem.py
    * Or dedicate a local folder as additional script ressource and include it in your pymolrc, as seen [here][4]

Usage
-----
* Pymol GUI: Plugin -> Membrane Protein Highlight
* Pymol CLI: hmem [pdbCode]

Planned improvements:
* Alternative annotation from [OPM][2]
* Use [TOPDB][3] cross-references in PDBTM (if available) to rename Side1/Side2 to inside/outside
* Allow preloaded structures (overwriting any conflicting present highlights) 

[1]: http://pdbtm.enzim.hu/
[2]: http://opm.phar.umich.edu/
[3]: http://topdb.enzim.hu/
[4]: http://www.pymolwiki.org/index.php/Git_install_scripts#Adding_Pymol-script-repo_to_PyMOL_search_path