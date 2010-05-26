#####################################################################
## Files and directories used in this script

GO_TREE_FILE=data/gene_ontology.obo

TAIR_DOMAINS_FILE=data/tair9_iprscan1_0323_no_seg_no_coil.out.gz
TAIR_SEQUENCES_FILE=data/TAIR9_pep_20090619_representative_gene_model.txt

INTERPRO_NAMES_FILE=data/names.dat
INTERPRO_GO_MAPPING_FILE=data/interpro2go

BLAST_DIR=/home/local/tamas/src/blast/current/bin
SCPS_DIR=/home/local/tamas/bzr/clusterix/build/src/ui/cmdline

#####################################################################
## Analysis parameters

# Which assignment sources NOT to trust from InterPro? (space separated)
IGNORED_SOURCES=HAMAP PatternScan FPrintScan Seg Coil

# E-value threshold to use when processing the initial InterPro file
INTERPRO_ASSIGNMENT_E_VALUE_THRESHOLDS=1e-3;superfamily=inf;HMMPanther=inf;Gene3D=inf;HMMPIR=inf

# Unclassified sequence fragments

# Minimum number of amino acids in a sequence in order to consider
# it further (i.e. calculate its unassigned fragments)
MINIMUM_SEQ_LENGTH=30
# Minimum number of amino acids in a sequence fragment in order to
# consider that fragment as a domain candidate
MINIMUM_DOMAIN_CANDIDATE_LENGTH=75

# BLAST filtering step

# Minimum sequence identity between two fragments to consider them
# as being in the same domain
MINIMUM_SEQ_IDENTITY=45
# Minimum normalized alignment length between two fragments to consider them
# as being in the same domain
MINIMUM_NORMALIZED_ALIGNMENT_LENGTH=0.7
# Maximum E-value between two fragments in order to consider them
# as being in the same domain
MAXIMUM_E_VALUE=1e-3
# Minimum Jaccard similarity between the neighbour sets of two
# fragments in order to consider them as being in the same domain
MINIMUM_JACCARD_SIMILARITY=0.66

# Novel domain assignment step

# A novel domain occur in at least this number of sequences
MINIMUM_NOVEL_DOMAIN_SIZE=4

#####################################################################

PYTHON=python
FORMATDB=$(BLAST_DIR)/formatdb
BLASTALL=$(BLAST_DIR)/blastall
SCPS=$(SCPS_DIR)/clusterx

.PHONY: clean

all: domain_architectures.zip

coverage: work/coverage_comparison.txt

overrep: work/overrepresentation_analysis.txt

work/gene_ids.txt: $(TAIR_SEQUENCES_FILE)
	grep "^>" $< | cut -d ' ' -f 1 | tr -d '>' >$@

work/coverage_comparison.txt: $(TAIR_DOMAINS_FILE) work/filtered_assignments.txt work/gene_ids.txt gfam/scripts/coverage.py
	$(PYTHON) -m gfam.scripts.coverage -x "$(IGNORED_SOURCES)" --gene-ids work/gene_ids.txt $< | sort >work/coverage_all.txt
	$(PYTHON) -m gfam.scripts.coverage -x "$(IGNORED_SOURCES)" --gene-ids work/gene_ids.txt work/filtered_assignments.txt | sort >work/coverage_selected.txt
	join -j 1 -t ' ' work/coverage_all.txt work/coverage_selected.txt | sed -e 's/ /\t/g' >$@

work/filtered_assignments.txt: $(TAIR_DOMAINS_FILE) work/gene_ids.txt bin/assignment_source_filter.py data/ParentChildTreeFile.txt
	zcat $< | bin/assignment_source_filter.py -x "$(IGNORED_SOURCES)" -e "$(INTERPRO_ASSIGNMENT_E_VALUE_THRESHOLDS)" -i data/ParentChildTreeFile.txt --gene-ids work/gene_ids.txt >$@

work/unassigned_fragments.tab: work/filtered_assignments.txt $(TAIR_SEQUENCES_FILE) bin/find_unassigned.py
	bin/find_unassigned.py -l $(MINIMUM_SEQ_LENGTH) -f $(MINIMUM_DOMAIN_CANDIDATE_LENGTH) -S $(TAIR_SEQUENCES_FILE) <$< >$@

work/unassigned_fragments.ffa: work/unassigned_fragments.tab $(TAIR_SEQUENCES_FILE)
	$(PYTHON) -m gfam.scripts.seqslicer -i $< $(TAIR_SEQUENCES_FILE) >$@

work/unassigned_fragments.psq: work/unassigned_fragments.ffa
	cd work && $(FORMATDB) -n unassigned_fragments -i unassigned_fragments.ffa -o F

work/unassigned_fragment_matches.blast: work/unassigned_fragments.ffa work/unassigned_fragments.psq
	cd work && $(BLASTALL) -p blastp -d unassigned_fragments -i unassigned_fragments.ffa -m 8 -o unassigned_fragment_matches.blast

work/relevant_matches.blast: work/unassigned_fragment_matches.blast work/unassigned_fragments.ffa
	$(PYTHON) -m gfam.scripts.blast_filter -s $(MINIMUM_SEQ_IDENTITY) -l $(MINIMUM_NORMALIZED_ALIGNMENT_LENGTH) -n query -e $(MAXIMUM_E_VALUE) -S work/unassigned_fragments.ffa $< >$@

work/relevant_matches_jaccard.tab: work/relevant_matches.blast
	$(PYTHON) -m gfam.scripts.jaccard -l $< >$@

work/relevant_matches_cca.tab: work/relevant_matches_jaccard.tab
	$(SCPS) -m cca --param epsilon=$(MINIMUM_JACCARD_SIMILARITY) -f plain_compressed $< >$@

work/domain_architectures.tab: work/filtered_assignments.txt bin/find_domain_architecture.py work/relevant_matches_cca.tab $(TAIR_SEQUENCES_FILE) $(INTERPRO_NAMES_FILE)
	bin/find_domain_architecture.py -n $(INTERPRO_NAMES_FILE) -i data/ParentChildTreeFile.txt -d work/domain_architecture_details.txt -S $(TAIR_SEQUENCES_FILE) -s $(MINIMUM_NOVEL_DOMAIN_SIZE) - work/relevant_matches_cca.tab <$< >$@

work/overrepresentation_analysis.txt: work/domain_architectures.tab $(GO_TREE_FILE) $(INTERPRO_GO_MAPPING_FILE)
	$(PYTHON) -m gfam.scripts.overrep $(GO_TREE_FILE) $(INTERPRO_GO_MAPPING_FILE) $< >$@

domain_architectures.zip: work/domain_architectures.tab work/overrepresentation_analysis.txt
	zip $@ work/domain_architecture* work/overrepresentation_analysis.txt

clean:
	rm -rf work/
	mkdir work/