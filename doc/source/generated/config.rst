Section ``DEFAULT``
^^^^^^^^^^^^^^^^^^^

``file.input.iprscan``
    Raw output file from IPRScan that contains domain assignments for
    all the sequences we are interested in.
    
``file.input.sequences``
    A FASTA file containing all the sequences being analysed.
    
``file.mapping.gene_ontology``
    A file containing the Gene Ontology in OBO format.
    
    Default value: ``data/gene_ontology.obo``
    
``file.mapping.interpro2go``
    File containing the mapping of GO terms to InterPro entries, as
    downloaded from geneontology.org.
    
    Default value: ``data/interpro2go``
    
``file.mapping.interpro2name``
    File containing a tab-separated list of InterPro IDs and their
    corresponding human-readable descriptions. This can be constructed
    by the following command::
    
        $ bin/download_names.py | gzip -9 >data/names.dat.gz.
    
    Default value: ``data/names.dat.gz``
    
``file.mapping.interpro_parent_child``
    File containing the parent-child relationships of InterPro terms.
    
    Default value: ``data/ParentChildTreeFile.txt``
    
``folder.work``
    The working folder in which to put intermediary files.
    
    Default value: ``work``
    
``folder.output``
    The output folder in which to put the final results.
    
    Default value: ``work``
    
``sequence_id_regexp``
    Python regular expression that matches gene IDs from the sequence
    file. This is necessary if the IPRScan output uses a sequence ID
    that is only a part of the sequence ID in the input FASTA file. If
    it is empty, no ID transformation will be done, the sequence IDs
    in the input FASTA file will be matched intact to the IPRScan
    output. If it is not empty, it must be a valid Python regular
    expression with a *named* group "id" that matches the gene ID that
    is used in the IPRScan output. If you don't know what named groups
    are, check the documentation of the Python `re` module.
    
``untrusted_sources``
    Which assignment sources NOT to trust from InterPro? (space
    separated)
    
    Default value: ``HAMAP PatternScan FPrintScan Seg Coil``
    
``max_overlap``
    Maximum overlap allowed between assignments originating from the
    same data source.
    
    Default value: ``20``
    
Section ``analysis:blast_filter``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``min_seq_identity``
    Minimum sequence identity between two fragments to consider them
    as being in the same domain.
    
    Default value: ``45``
    
``min_alignment_length``
    Minimum alignment length between two fragments to consider them as
    being in the same domain. If normalization_method is not off, this
    must be the normalized alignment length threshold according to the
    chosen normalization method.
    
    Default value: ``0.7``
    
``max_e_value``
    Maximum E-value between two fragments in order to consider them as
    being in the same domain.
    
    Default value: ``1e-3``
    
``normalization_method``
    Normalization method to use for calculating normalized alignment
    length Must be one of: off, smaller, larger, query, hit.
    
    Default value: ``query``
    
Section ``analysis:find_domain_arch``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``min_novel_domain_size``
    A novel domain occur in at least this number of sequences.
    
    Default value: ``4``
    
Section ``analysis:find_unassigned``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``min_seq_length``
    Minimum number of amino acids in a sequence in order to consider
    it further (i.e. to calculate its unassigned fragments)
    
    Default value: ``30``
    
``min_fragment_length``
    Minimum number of amino acids in a sequence fragment in order to
    consider that fragment as a novel domain candidate.
    
    Default value: ``75``
    
``sequences_file``
    Input FASTA file containing all the sequences of the
    representative gene model being analysed.
    
    Default value: same as ``file.input.sequences``
    
Section ``analysis:iprscan_filter``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``e_value_thresholds``
    E-value thresholds to use when processing the initial InterPro
    file. This entry is a semicolon-separated list of source=threshold
    pairs, the given threshold will be used for the given data source.
    If an entry contains only a threshold, this will be considered as
    a default threshold for sources not explicitly mentioned here.
    
    Default value: ``1e-3;superfamily=inf;HMMPanther=inf;Gene3D=inf;HMMPIR=inf``
    
``interpro_parent_child_mapping``
    File containing the parent-child relationships of InterPro terms.
    
    Default value: same as ``file.mapping.interpro_parent_child``
    
``stages.1``
    These configuration keys specify which assignment sources are to
    be taken into account at each stage of the analysis. For more
    information about what these stages are, please refer to the
    documentation, especially the description of Step 2 in section
    "Steps of the GFam pipeline"
    
    Default value: ``ALL-HMMPanther-Gene3D``
    
Section ``analysis:jaccard``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``min_similarity``
    Minimum Jaccard similarity between the neighbour sets of two
    fragments in order to consider them as being in the same domain.
    
    Default value: ``0.66``
    
``assume_loops``
    Whether to assume that a protein is connected to itself or not.
    
    Default value: ``1``
    
``only_linked``
    Whether to consider only those protein pairs which are linked in
    the input file. If this is 1, protein pairs not in the input file
    will not be returned even if their Jaccard similarity is larger
    than the given threshold.
    
    Default value: ``1``
    
Section ``analysis:overrep``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``confidence``
    The p-value threshold in the hypergeometric test used in the
    overrepresentation analysis process.
    
    Default value: ``0.05``
    
``correction``
    The method used to account for multiple hypothesis testing. Valid
    choices: bonferroni (Bonferroni correction), Sidak (Sidak
    correction), fdr (Benjamini-Hochberg method), none (off). The
    Bonferroni and Sidak correction methods control the family-wise
    error rate (FWER), while the Benjamini-Hochberg method control the
    false discovery rate.
    
    Default value: ``fdr``
    
``min_term_size``
    The minimum number of annotated domains a GO term must have in
    order to be considered in the overrepresentation analysis.
    
    Default value: ``5``
    
Section ``generated``
^^^^^^^^^^^^^^^^^^^^^

``file.unassigned_fragments``
    File in which the unassigned sequence fragments are stored.
    
    Default value: ``%(folder.work)s/unassigned_fragments.ffa``
    
``file.valid_gene_ids``
    File containing a list of valid gene IDs (extracted from the input
    file)
    
    Default value: ``%(folder.work)s/gene_ids.txt``
    
``file.domain_architecture_details``
    File containing the detailed final domain architecture for each
    sequence.
    
    Default value: ``%(folder.output)s/domain_architecture_details.txt``
    
Section ``utilities``
^^^^^^^^^^^^^^^^^^^^^

``folder.blast``
    The folder containing the BLAST executables. It does not matter
    whether you have the old C-based or the newer C++-based tools,
    GFam can use both if you also have the ``legacy_blast.pl`` script
    that adapts the new tools to the command line syntax used by the
    older ones.
    
``util.formatdb``
    The path to formatdb. You may use the name of the folder
    containing the BLAST executables or the full path (including the
    name of the tool). If you have the newer, C++-based BLAST tools
    (which do not have formatdb), pass the name of the folder
    containing the BLAST executables here, and if you have the
    ``legacy_blast.pl`` script in the same folder (plus a working Perl
    setup), GFam will detect the situation and run ``legacy_blast.pl``
    accordingly.
    
    Default value: same as ``folder.blast``
    
``util.blastall``
    The path to blastall. You may use the name of the folder
    containing the BLAST executables or the full path (including the
    name of the tool). If you have the newer, C++-based BLAST tools
    (which do not have blastall), pass the name of the folder
    containing the BLAST executables here, and if you have the
    ``legacy_blast.pl`` script in the same folder (plus a working Perl
    setup), GFam will detect the situation and run ``legacy_blast.pl``
    accordingly.
    
    Default value: same as ``folder.blast``
    
