# LargeScaleMutagenesis

Data generation
    - The generation of the data requires a .txt file consisting of a list of PDB identifiers formatted as a single 'pdbid_chain' string per line (to be located in the '/List' directory).
    The second requirement are the PoPMuSiC output files corresponding to the entries of the afore-mentionned list (located in the '/Data/pop' directory).
    The whole pipeline for generating the data can then be achieved by first editing the 'config.py' with the relevant paths and then executing 'exec.py' (both are located in the '/Code' directory).
    See the 3 steps below for a more detailed explanation of the pipeline.
    - The lists for the datasets used in our paper ('Mpop' and 'Mexp') are located in the '/List' directory. 
    - Data examples are given for PDB id '1A9N'.

    1) Stability change computation:
    The first step in generating the data is to apply the PoPMuSiC algorithm (http://dezyme.com/) on the PDB structures (https://www.rcsb.org/) of a dataset of interest.
    The output contains, among other information, the predicted stability change for the 19 variants of each residue of the protein structures.
    Access to the version PoPMuSiCsym (https://doi.org/10.1016/j.ifacol.2015.05.068) will be available soon.

    2) Codon mapping:
    The second step consists in mapping the codons coding for the residues of each protein in the dataset of interest.
    First, 'get_seq_aa.py' generates a fasta file of the amino acid sequence of each protein.
    Then, 'get_seq_na.py' retreives the available nucleotide sequences of transcripts for each protein.
    Then, 'get_ali.py' aligns a translation of the retrieved nucleotide sequences with the amino acid sequence for each protein.
    Finally, 'get_ali_map.py' evaluates the best alignment (based on identity and similarity thresholds) and generates a file that maps the  ensembl identifier of the best alignment to the pdb identifier of its protein.

    3) Mutational ensembles definition:
    The last step is to define the mutational ensembles (i.e. single or multiple nucleotide substitution, at which position, ...) corresponding to each variants of the residues of each protein based on the mapped codons, and to integrate that information to the output files of PoPMuSiC.
    First, 'get_mut.py' generates the files that link the variants to their respective mutational ensembles.
    Then, 'get_xfiles.py' generates the final data files to be used in the analysis.

Data analysis
    - The data analysis simply requires to extract the relevant data from the previously generated .popx files.
    - The extra data required for some analyses (such as the codon usage bias data or the fitness data) is located in '/Extra'.
