pymol-highlight-transmem
========================
A [PyMOL][3] script to highlight transmembrane segments and related structural elements on a given PDB structure.
Annotations are retrieved from either [PDBTM][1] or [OPM][2]. Choosing PDBTM essentially reproduces the behaviour of the viewer on their website within PyMOL.

Note that OPM does not have explicit definitions for membrane sides. In addition there is no (accessible) distinction between interfacial, reentrant and transmembrane helices, so all are colored the same (yellow)
OPM only gives a non-redundant list of subunit annotations in one file. Therefore queuing for PDB code '1xfh' will not give a hit, since the representative structure for this cluster is '2nwl' and only this one appears in the subunits file.

Install
-------
* For current session only, execute in Pymol CLI:

        import hlmem
* For permanent availability:
    * Start Pymol and go to Plugin -> Manage Plugins -> Install, then choose hlmem.py
    * Or dedicate a local folder as additional script ressource and include it in your pymolrc, as seen [here][4]

Usage
-----
* Pymol GUI: Plugin -> Membrane Protein Highlight
* Pymol CLI: hmem [pdbCode [,loaded,(OPM|PDBTM)]]

If loaded=0 (default), the structure will be fetched before highlighting, otherwise only selections will be created and colouring performed.
Default highlights are from PDBTM, but OPM can be specified as well.

        hmem 2nwl,0,OPM     //Fetch 2nwl and highlight according to OPM
        hmem 2nwl           //Use already loaded 2nwl and highlight according to PDBTM

Planned improvements:
-----
* Offer to copy OPM annotations for redundant hits (annotations are only given for a non-redundant set of structures)
* Update GUI to new features available through 'hmem' on the commandline

References
-----
* Kozma,D. et al. (2013) _PDBTM: Protein Data Bank of transmembrane proteins after 8 years._ NAR [DOI][5]
* Lomize,M.A. et al. (2012) _OPM database and PPM web server: resources for positioning of proteins in membranes._ NAR [DOI][6]

[1]: http://pdbtm.enzim.hu/
[2]: http://opm.phar.umich.edu/
[3]: http://pymol.org/
[4]: http://www.pymolwiki.org/index.php/Git_install_scripts#Adding_Pymol-script-repo_to_PyMOL_search_path
[5]: http://dx.doi.org/10.1093/nar/gks1169
[6]: http://dx.doi.org/10.1093/nar/gkr703
