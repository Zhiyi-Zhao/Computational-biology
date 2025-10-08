# Computational-biology
Computational biology note based on Wageningen University&Research's class. Course code:SSB34306
This is a course note writed by Zhiyi. While retaining the framework of the lectures and assignment content, some of the content of interest has been organized here.

# Skills that need to be mastered in advance
This course uses Jupyter for assignments, which uses two common cell types:
Code Cells: Used to write Python code, which generates output when run.
Markdown Cells: Used to write captions, formulas, images, tables, and more. It supports Markdown syntax and can directly render LaTeX mathematical formulas.

⭐⭐⭐**Course**:
- Biochemistry
- Molecular biology
- Evolutionaty biology

## Python
This lesson will only use simple Python code to perform formula calculations. If you don't want to learn, you can also do the calculations manually.

[Python--Basic Operators](https://www.learnpython.org/en/Basic_Operators)

[If you want to know more about python](https://wiki.python.org/moin/BeginnersGuide/NonProgrammers)

## Markdown
Markdown is a lightweight markup language that allows you to write documents in a plain text format that's easy to read and write. It also supports inserting rich web links and highlighting code blocks and their associated language.

[Why Markdown](https://www.markdownguide.org/getting-started/)\
[How Markdown](https://commonmark.org/help/)

## LaTeX
You can use LaTeX to write formulas in your daily homework, or you can use LaTeX to complete beautiful course papers.
[Learn Latex in 30 minute](https://www.overleaf.com/learn/latex/Learn_LaTeX_in_30_minutes)


# 1A Molecular biology primer: from gene to protein
## Lecture
- Central dogma\
[Central dogma](https://www.genome.gov/genetics-glossary/Central-Dogma). The central dogma of molecular biology is a theory stating that genetic information flows only in one direction, from DNA, to RNA, to protein, or RNA directly to protein.

Protein is more abundant in creatures, people also first sequenced protein by Edmann degradation method. But now, DNA sequence is becoming cheaper and cheaper, the data for DNA sequence is also increasing in a fast way.

The [DNA sequencing method](https://www.geneious.com/guides/introduction-to-dna-sequencing), we will learn how to analyze sequencing data in the coming course, but its also recommended to learn where the data come from.

- Genome\
In classical genetics, the genome of a diploid organism refers to a full set of chromosomes or genes in a gamete. A regular somatic cell contains two full sets of genomes(and a mitochondrial genome)\
In haploid organisms, including bacteria, archaea and viruses, a cell contains only a single set of the genome, usually in a single circular or contiguous linear DNA.\
Mitochondrial and Chloroplast genome are also present in single copy.\
Size of the genomes is variable and does not correlate with complexity of living beings. And small and simple does not mean harmless.

- Gene\
A gene is a programmable unit that can give rise to a multitude of products: protein and RNA products through (alternative) splicing and trans-splicing

## Assignment
- Reading Frame(RF), Open Reading Frame(ORF) and Coding Sequence(CDS)
  - A **Reading Frame** is how a genetic sequence is read during a translation process.
  - An **Open Reading Frame** is the part of a reading frame that has the potential to code for a protein or peptide and consists of a continuous stretch of codons that do not contain a stop codon (usually UAA, UAG or UGA). Note that in this case, the ORF will be interrupted with intervening sequences or introns.
  - The **Coding Sequence**, is the portion of a gene's DNA or RNA, composed **only of exons**, that codes for protein. The region is bounded at the 5' end with a start codon and at the 3' end with a stopcodon.
 
- Intron\
An intron is a nucleotide sequence within a gene in eukaryote. It is a noncoding sequence. During the final maturation of the RNA product, the introns are removed by splicing.

- Frameshift\
Frameshift is a shifting of the reading frame caused by insert nuclides not devided by 3. Frameshift mutation is a genetic mutation caused by indels (insertions or deletions) of a number of nucleotides in a DNA sequence that is not divisible by three. 

Comment of the assignment:
Use the [SMS website](https://www.bioinformatics.org/sms2/#:~:text=The%20Sequence%20Manipulation%20Suite%20is%20a%20collection%20of,for%20teaching%2C%20and%20for%20program%20and%20algorithm%20testing.) to find the ORF corresponding to the known peptide on the known DNA sequence. Noting:
- Continuous peptides do not equal continuous DNA sequences; they may be located in different ORFs and separated by introns.
- Knowing which DNA chain the ORF is located on means you don't have to worry about the ORF on the other chain.

# 1B Amino Acids

Only understand amino acids, can we understand protein. What happend in enzymes' action site, how secondary structure.

Remember 20 amino acids!
![Table of 20 amino acids](/AAKT.PNG)


# 2A Genome Coverage
## Lectures
Coverage: average number of times any given base in the genome is sequenced.
<center>

$$Coverage(a) = \frac{N \times L}{G}$$

</center>

- Number of reads -- N
- Read Length -- L
- Genome length -- G

Given coverage **a**, then the probability that base is sequenced n times or **a base not being sequenced(Gap)**:
<center>

$$P(n) = \frac{a^n e^{-a}}{n!}$$

</center>

By this fomula: P(0) = the probability of the gaps = $e^{-a}$, so for a genome of size G, the number of nucleotides in the gaps is: $P(0)G = e^{-a}G$

- Why Possion distribute and Why lambda=coverage:\
  Why does sequencing follow a Poisson distribution?

In next-generation sequencing, we have a genome of length L (e.g., a 3 Gb human genome).The sequencer randomly selects a DNA fragment from the library and generates a read.The position of each read on the genome is approximately random and independent.Therefore, for a given site on the genome:

The probability of each read covering that site is a constant: $p =\frac{ReadLength}{GenomeLength}$, if there are a total of 𝑁, then the coverage of this site **X ~ Binomial distribution**, When 𝑁 is large, p is small, and 𝑁𝑝=𝜆 (read length) is fixed, the binomial distribution approximates the Poisson Distribution.

H = number of contigs = NP(0) = $Ne^{-a}$
Compares to the calculation fomula of Gaps, the only difference is N(Number of reads) and G(Genome size)
But this is just calculate number, don't have any biology meaning, in real sequencing, some region is hard to sequence or assembly, like the highly repeat sequence. And the H can less than 1, but actually, it's not possible because creatures at least have 1 contig.

## Assignment

**Important Concepts**:
- Shotgun sequencing: The technology used to sequence complete genomes. Shotgun sequencing involves randomly breaking up or shredding a genomic or large DNA sequence into many small fragments (hence shotgun), sequence them and then use a computer to reassemble the original sequence by finding regions of overlap.
- Sequence read: The output of the sequencing process. A read refers to a string of nucleotides. This string is presented in FASTA or FASTQ format.
- Contig: When two or more sequence reads show a region of overlap, they can be merged into one single larger string. This merging by using overlap is the assembly process and the resulting string is called a contig and usually presented in fasta format.
- Scaffold: When the order of the contigs is known they can be connected in one long sequence with gaps in between the contigs. The length of the gaps can be estimated or are sometimes set with a fixed length. These ordered contigs are called scaffolds or supercontigs.
Paired-end or long_read data information or a reference genome is often used for scaffolding.
- N50: N50 is a statistic that defines the quality of an assembly. An assembly is a set of contigs that together represent a genome. Given a set of contigs, each with its own sequence length, N50 is defined as the shortest sequence length at 50% of the genome. A simple example of the calculation is found here.
- NG50：N50 is calculated in the context of the assembly size rather than the genome size. Therefore, comparisons of N50 values derived from assemblies of significantly different lengths are usually not informative, even if for the same genome.
To address this, the authors of the Assemblathon competition came up with a new measure called NG50. The NG50 statistic is the same as N50 except that it is 50% of the known or estimated genome size that must be of the NG50 length or longer. This allows for meaningful comparisons between different assemblies. In the typical case that the assembly size is not more than the genome size, the NG50 statistic will not be more than the N50 statistic.
- Coverage: This is a container term for a number of similar terms. Given a file with many sequence reads the average coverage can be calculated from the length of the genome (G), the number of reads in the file (N), and the average read length (L) of the reads, following the Lander-Waterman equation (N × L / G). In shotgun sequencing a high coverage is essential in order to find overlap between reads and repress sequence read errors as a result of the sequencing process which relates to the actual coverage at a particular position in a sequence.


**FASTA File**\
The FASTA format is a text-based format for representing either nucleotide sequences or amino acid (protein) sequences, in which nucleotides or amino acids are represented using single-letter codes.
A sequence begins with a greater-than character (">") followed by a description of the sequence (all in a single line). The lines immediately following the description line are the sequence representation, with one letter per amino acid or nucleic acid, and are typically no more than 80 characters in length.

For example:

\>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*

**FASTQ File**\
FASTQ File has four line-separated fields per sequence:

Field 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line). Field 2 is the raw sequence letters. Field 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again. Field 4 encodes the quality values for the sequence in Field 2, and must contain the same number of symbols as letters in the sequence. A FASTQ file containing a single sequence might look like this:

@SEQ_ID\
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\
+\
!''((((+))%%%++)(%%%%).1-+''))**55CCF>>>>>>CCCCCCC65

They have some difference on format: FASTA begins with ">", FASTQ begins with "@". And FASTQ file has "+" to link the coming line and has the information of quality value. FASTA files are more flexible, allowing to have the sequence in multiple lines while in FASTQ formats the sequence is restricted to only one line. And FASTQ file can give more information about the sequence, they have the record of quality. Therefore they also bigger in the size of files when they have same length of sequence.

More information about [FASTQ and Illumina](https://knowledge.illumina.com/software/general/software-general-reference_material-list/000002211)\
More information about [FASTQ and Q-score](https://en.wikipedia.org/wiki/FASTQ_format)

Transfer FASTQ into FASTA or combine FASTA and QUAL into FASTQ -- [GALAXY Platform](https://usegalaxy.org/) : Tools--FASTA/FASTQ\

**Pair-end sequence**\
The header of each sequence contains pertinent information about the read.\
'/1' at the end of the header stands for the forward read '/2' for the reverse read:\
Forward > @NC_000913.3:320000-380000-2388/1\
Reverse > @NC_000913.3:320000-380000-2388/2

**GALAXY Platform**\
[Galaxy](https://usegalaxy.org/) is a web-based, open-source bioinformatics analysis platform jointly developed by Pennsylvania State University and Johns Hopkins University. As an open platform, Galaxy integrates a large number of bioinformatics analysis tools, providing researchers with an easy-to-use bioinformatics analysis interface that **requires no programming knowledge**.

Galaxy Functions:
- Genomics: Gene assembly, variant detection, gene annotation, etc.
- Transcriptomics: RNA-seq analysis, differential expression analysis, etc.
- Proteomics: Mass spectrometry data analysis, protein identification and quantification, etc.
- Metabolomics: Metabolic pathway analysis, metabolite identification, etc.
- Systems biology: network construction and analysis, functional enrichment analysis, etc.

**Length of reads**:\
The sequencer model, the version of the sequencing chemistry and the number of cycles.\
The specific model of the instrument used (for example, HIseq, NovaSeq or MiSeq) will determine the length of reads produced.\
Illumina is a second generation technology, with a sequencing chemistry referred to as Sequencing by synthesis. Here, clusters of identical molecules are sequenced simultaneously. Each chemistry cycle adds a tagged base to the end of a sequence in a cluster followed by an imaging step where the newly added base is read. Sometimes not one base binds and thus gets out of sync with the other molecules in the cluster. More cycles therefore results in lower resolution, limiting the length of the reads that can be done.

Comment of the assignment:
- Use [GALAXY](https://usegalaxy.org/) to count the sequence number and the length of sequence([Line/Word/Character count](https://usegalaxy.org/?tool_id=wc_gnu&version=latest))
- Use [SPAdes](https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fnml%2Fspades%2Fspades%2F4.2.0%2Bgalaxy0&version=latest) for genome assembly
- Calculate H,N50,NG50

# 2B Genome Assembly & Annotation
## Lectures
**Course Aim**
- The differences between mapping and assembly:\
  Mapping = compare to a known genome\
  Assembly = build the genome from scratch
- The general strategy for genome assembly
  1. Find overlaps among reads
  2. Build a graph to visualize the connections
  3. Make the graph simpler
  4. Walk through the graph
- The specific approaches(Greedy overlapping) for genome assembly
- Compare the structure of contigs and scaffolds\
  Scaffold is made by the contigs placed in correct orientation and correct order, with the approximately correct distance sometimes including repeats.
- Define the scope of structural and functional annotations
  - Structural annotation: is there a gene? Three common methods:
    1. use staristical models to search for specific sequences that indicate the presence of a gene nearby, or statistial properties of the protein-coding sequence itself.
    2. Evidence-based approach: sequence mRNA
    3. Homology-based approach: align genome with known genes or proteins to find orthologs
  - Functional annotation: what's the function of the sequence?
    1. From protein to Gene:Experimental function of known proteins → mapped to corresponding genes
    2. Predict function from homology protein:If a gene has a high similarity to a known protein sequence, it is inferred that it may have the same or similar function.

[**SCS solution and Greedy overlapping**](https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_scs.pdf)

**Greedy overlapping based assembly**: Use heuristic approach to solve problems: find the best local choice at each step
1. calculate pairwise alignment of all reads vs all reads
2. Find the best matching read pair
3. Merge pair
4. Repeat until cannot continue

OLC(Overlap-layout-consensus):
1. Use dynamic programming to make a string comparison between the reads and obtaining the optimal alignment
2. Look for significant overlap allowing for a few mutations and deletions
3. Create an overlapping graph where every read is a vertex connected by an edge in case of significant overlap.
4. Sort edges by weight
5. Clean the redundancy

[de Bruijn graph](https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf)

**de Bruijn graph**:The de Bruijn graph transforms the overlap problem into a path problem:\
A graph is constructed using the (k-1) base overlaps between k-mers.\
Then, an Eulerian path is found within the graph to reconstruct the genome.

## Assignment
Comment of the assignment:\
Use Bowtie2 to mapping the reads(with ref_genome). Then use Jbrowse to visualize the results of alignment.\
We found Single Nucleotide Polymorphism (SNP), gaps(may means deletion or mutation on this site compare to ref_genome), 

# 3A Transcript analysis

**Course Aim:**
- types of RNA and main characteristics:\
  messenger RNA(main), code for proteins\
  transfer RNA(main), central to protein synthesis as adaptors between mRNA and amino acids\
  ribosomal RNA(main), form the basic structure of the ribosome and catalyze protein synthesis \
  microRNA\
  short interfering RNA\
  small nucleolar RNA\
  long non-coding RNA\
  piwi-interacting RNA etc.
- various processing steps involved between transcription and protein expression:\
  1. Pre-mRNA is spliced to form mature mRNA, capped and polyadenylated
  2. In eukaryotes, one gene can result in different mRNA transcripts; in bacteria, one mRNA can translate into multiple proteins
  3. Translation
- The use and limitation of measuring mRNA expression to understand biological process, including the influence of technical and biological noise:
  1. Usage: mRNA and protein often shows same trends, so we can use it for conditions comparison and time series comparison
  2. Many reason can cause difference: different genes,isoforms,tissues,developmental stages,cell cycle,circadian rhythm, individual cells,environment, and the results of mRNA synthesis and mRNA decay
- Explain the main technologies to measure transcript concentrations and list their specific advantages and disadvantages and main use:
  1. RT-qPCR: Low throughput
  2. Microarrays: pros: highly standardized, relatively cheap, samll data size--easy to handle; cons: gene sequence should be known, no position-specific information, can't detect new isoforms, not very quantitative: low dynamic range
  3. RNA-seq: untargeted, works for species without a sequenced genome, can identify alternatively spliced transcripts, can identify SNPAs in transcripts, high dynamic range(quantitative), no strand specificity
- Identify the challenges in mapping RNA-seq reads to reference genomes and describe the use of transcriptome assembly\
  challenges: RNA reads will span an intron on the genome; one exon can be part of multiple isoforms

**Pseudo-alignment**
It compares the **k-mers** in the reads with the **k-mers** set of the reference transcriptome to determine which transcripts each read is compatible with, without calculating its precise alignment position on the transcript.\
Why pseudo-alignment is faster:\
It eliminates base-level alignment and searching across splice junctions, instead using *k-mer hash matching* to search directly on the spliced ​​transcriptome, significantly reducing computational effort.

## Assignment
Comment of the assignment:\
Use Ensembl website download transcript data of yeast, analyze the mean length, from forward strand or reverse strand;\
Use IGV website to visualize transcript results, we can search a specific gene, from the transcript data, we can know that how many exons and introns we have. And we can have a brife idea of the level of the gene expression. In some region, there's no gene, but still have some expression, the reason can be: an unannotated gene, read-through of the upstream gene (transcription stop signal ignored), other forms of expression (enhancers can cause local transcription), contamination with genomic DNA. And we can find the alternative splicing. And we can use this tool to copy interest sequence.\
Then use SMART to analyse our interested sequence, to predict the protein structure.
 
# 3B Transcriptomics Quantification and differential expression(DE) analysis
## Lecture
- Normalizing for sequenceing depth:\
  **CPM**\
  **C**ount **P**er **M**ilion of mapped reads: counts scaled for total number of reads, relative measure for reads count where the total amount of reads is set to 1 milion to avoid samll numbers:
  <center>
 
  $$CPM =10^{6} \frac{ReadsPerGene}{TotalReads}$$

  </center> 

  **RPKM/FPKM**\
  **R**eads **P**er **K**ilo base of transcript per **M**illion mapped reads
    <center>

    $$RPKM =10^{6} \times \frac{ReadsPerGene}{TotalReads} \times 10^{3} \times \frac{1}{GeneLength} = 10^{9} \times\frac{ReadsPerGene}{GeneLength \times TotalReads} $$

    </center> 
  As to paired-end sequence, RPKM = 2 * FPKM

  **TPM**(Recommended, we can use this to compare both inside and between samples)\
  <center>

    $$TPM =10^{6} \times \frac{ReadsPerGene}{GeneLength} \times  \frac{1}{TotalLengthCorrectedCounts} = 10^{6} \times\frac{\frac{ReadsPerGene}{GeneLength}}{\sum\frac{ReadsPerGene}{GeneLength}} $$

    </center> 
    
 Different Expression:   t-test
 Volcano plot and cluster analysis, Go analysis(Gene Ontology) provides a system for hierarchically classifying genes or gene products into terms organized in a graph structure

## Assignment
Comment of the assignment:\
Analyze the transcript results, brife analyse which is expressed which not, and see is there difference between different expriment group.\
Use Goprofiler to see the different expression gene and its function.\
Use microarray to see common expression, by using SPELL website we can see the gene different expression under different condition in different published article.

# 4A MS-based proteomics  
## Lecture
- Why do we study proteins?
  1. It can tell us what happens now, what enzymes are currently active, which signals are being transduced.
  2. Transcriptomics sometimes have poor correlation with proteomics.
  3. Genomics can only tell what's the potential
- What is proteomics?
  The study of all proteins present in a sample at a given time.
- What are the challenges in analysing proteins?
- How do we do it?
- How can we manually give meaning to the spectra?
  

## Assignment
Comment of the assignment:\
By using MAmine3, we can find information like the elution time, the different isotopic peaks for the same peptide.\
We use MS-Viewer to compare the peptide with high score and low score.\
Calculate the peptide by hand

# 4B Identification and quantification of proteins in complex samples
**Course Aim:**
- Understanding how mass spectrometry data are analysed to produce protein indentifications
- Discover how proteins are quantified by mass spectrometry-based proteomics methods\
  Label free quantitation\
  Metabolic labeling\
  Isotopic or isobaric tags
- Reflect on how the quantity of the entity (proteins) is inferred by quantifying its components (peptides)
 
**FDR**： measures the error rate associated with a collection of PSMs:

$$
FDR(S_T) = \frac{N_d(S_T)}{N_t(S_T)}
$$

How to assemble identifed peptides into a list of proteins -- 1. grouping peptides; 2. grouping the protein side; 3.find the parsimony protein groouping solution

## Assignment
Comment of the assignment:\
Use OMSSA to see the result of MS, calculate the FDR, check the result of peptide hit

MS-MS have some limitation: The reason why these peptides were not discovered is basically that ionisation was too difficult. Peptides which are easier to measure are the peptides which are easily ionised. The ionization is a parameter highly influenced by the chemistry of the peptide (e.g. composition, length). Protein modification can also influent the result, so, we need to expand the search scope to take protein modifications into account. Hence, the false discovery rate goes up because the search space increased substantially.

# 5A Databases
## Lecture
- International databases:
  DNA/RNA: GenBnak for genetic sequences, SRA for next NGS data;\
  Protein: Uniprot for sequence information, PDB for structures;
  Metabolism: KEGG, BioCyc

## Assignment
**Eenzyme Code** consists of the letters "EC" followed by four numbers separated by periods. Those numbers represent a progressively finer classification of the enzyme. Preliminary EC numbers exist and have an 'n' as part of the fourth (serial) digit (e.g. EC 3.5.1.n3).

At NCBI you can start text searches which allows you to retrieve molecular biology data and bibliographic citations from the NCBI's integrated databases.

These include:

DNA sequences from GenBank, EMBL, and DDBJ
Protein sequences from SwissProt, PIR, PRF, PDB, and translated protein sequences from the DNA sequence databases.
Genome and chromosome mapping data.
Three-dimensional protein structures derived from PDB, and incorporated into NCBI's Molecular Modelling Database (MMDB)
PubMed bibliographic database containing citations for nearly 9 million biomedical articles from the National Library of Medicine's MEDLINE and pre-MEDLINE databases.
It has some useful features. For example:

Most of its records are linked to other records, both within a given database (such as PubMed) and between databases (related information shown at the right side of the screen).
Suggested Filters to refine your search
The sequence databases are linked to the Medline database such that it is possible to move from paper to sequence or vice versa seamlessly.

Comment of the assignment:\
Use Uniprot to search for a protein and see its linked desease in Disease & Variants, and download sequence data;\
Use NCBI to search a protein, and check the linked article in PubMed, we can also check the size of the protein and the function, and the site on the genome, and the intron and exon of a gene;\
Use KEGG to see how many pathway is related to the interested protein.\
Use PDB database to see the structure of the protein, like how many chains are contained in a specify protein.

# 5B Public resourses for genetic data & the need for data FAIRification
## Lecture
Nucleotide sequence database: DDBJ, EMBL, Genbank\
Uniprot: is made by two part: 1. Swiss-prot: the manually curated section; 2. TeEMBL: the automatically annotated section. this databse is focused on sequence data and functional annotaions\
PDB: focused on 3D structural data of proteins, nucleic acid and complexes

FAIR data: scientific data linked to machine readable metadata: Findable, Accessible, Interoperable, Reusable

## Assignment
[checklist](https://www.ebi.ac.uk/ena/browser/checklists) for how to make FAIR data
Comment of the assignment:\
Taking the data of diabetic patients as an example, the content of FAIRdata is shown.

# 6A Substitution Patterns & Matrices
## Lecture
**Course Aim:**
1. List reasons for why sequence alignments are done;
2. Compare the Jukes-Cantor and Kimura models for nucleotide substitution;\
   Jukes-Cantor:    $$d = \frac{-3}{4} \times ln(1-\frac{4}{3} \times D)$$, every mutation pair have same score \
   Kimura: take substitutions and transversions into consideration
3. Calculate log-adds ratios for amino acid substitutions;\
   $$log_2(\frac{Pr(x,y|R)}{Pr(x,y|U)}) = \sum log_2(\frac{q_{x,y}}{p_{x}p_{y}})$$
4. Describe why and how amino acid substitution matrices have been constructed\
**BLOSUM vs. PAM Matrices**

| Feature                        | BLOSUM62 (Default) | Other BLOSUM Matrices | PAM Matrices |
|--------------------------------|--------------------|----------------------|--------------|
| **Basis of Calculation**       | Derived from conserved blocks of protein sequences | BLOSUMX is derived from alignments with at least X% sequence identity | Based on evolutionary models of accepted mutations |
| **Best for**                   | Medium sequence divergence (~62% identity) | Higher X (e.g., BLOSUM80) for closely related sequences; Lower X (e.g., BLOSUM45) for distant sequences | Evolutionary distant sequences |
| **High Number (BLOSUM80 / PAM10)** | Detects close homologs | Sensitive for closely related sequences | Better for short evolutionary distances |
| **Low Number (BLOSUM45 / PAM250)** | Detects distant homologs | Good for highly divergent sequences | Suitable for long evolutionary distances |
| **Evolutionary Assumptions**   | Observed substitutions (not inferred from models) | More aggressive or conservative depending on X | Models mutations over time |   

**Key Differences**

1. **BLOSUM is Data-Driven, While PAM is Evolutionary**  
   - BLOSUM matrices are built from actual conserved regions, making them empirical.  
   - PAM matrices are based on a Markov model of accepted mutations over time.  

2. **BLOSUM is Better for Detecting Functional Similarities**  
   - Since BLOSUM is derived from conserved regions, it’s better at detecting proteins with similar functions.  

3. **PAM is Better for Studying Evolutionary Relationships**  
   - PAM matrices predict how sequences evolve over time, making them useful for phylogenetics.  

** Sequence alignments: comparisons of how each nucleotide or amino acid in two or more strings of DNA, RNA, or protein sequences with another 


## Assignment
Comment of the assignment:\
Link to a detailed calculation example

# 6B BLAST(Basic Local Alignment Search Tool)
## Lecture
**Learning goals**\
1. List the different tyoes of BLAST;
   - blastn n vs n
   - blastp p vs p
   - blastx n(trans) vs p
   - tblastn p vs n(trans)
   - tblastx n(trans) vs n(trans)  
2. Summarize the four steps of the BLAST algorithm;
   - List: 1. split sequence into k-mers; 2. calculate alignment score and keep only k-mers above certain threshold
   - Scan: 1. scan database for exact matches to each word from list step; 2. extend each hit in both directions calculate raw score until score drops a certain amount below highest score; 3. allows gaps to be filled as long as gap penalyies do not drop total score below a certain score
   - Report Raw scores (S) are converted to bit scores to allow comparisons between different searches, $S' = \frac{\lambda \times S - \ln(K)}{\ln(2)}$ 
3. Explain how changing BLAST parameters influence the algorithm and output;
   We can change the databse, targeted organism and protein matrix
4. Evaluate BLAST.


## Assignment
Comment of the assignment:\
Use BLAST in NCBI and learned how to evaluate the results and how to calculate these score.

The raw score (*S*) of an alignment depends on the number of matches/mismatches, the number of gaps and the length of the gaps, and can be calculated with the following formula:

$$
S = \sum (M_{i,j}) - cO - dG
$$

where $M_{i,j}$ is the sum of the scores for each match or mismatch (from the corresponding nucleotide or amino acid scoring matrices), *c* is the number of gaps, *O* is the penalty for the existence of a gap, *d* is the total length of gaps and *G* is the per-residue penalty for extending the gap.

Determine the **E-value** of the top hit that you found in **Question 4** with the following formula and show your calculations

$$
E = K \cdot m \cdot n \cdot e^{-\lambda S}
$$


Where *m* is the effective length of query and *n* is the effective length of database (i.e. sum of all sequences in the database). 


Like &lambda; and *K*, *m* and *n* are shown in the *Search summary* table.

Note the size of *n*. This is why a heuristic algorithm like BLAST is necessary!

## Lecture

## Assignment
Comment of the assignment:\

# 7 Protein Domains
## Lecture
**Domains** are distinct structural units of proteins that fold independently and have a hydrophobic core, they're the most basic functional units of proteins that are required for their activity.\
**Motifs** are found as patterns(specific or degenerate), which can define a domain, but also can be characteristic of a subset of a domain and found in DNA

Three general methods for domain identification:
- Patterns: highly conserved sequence motifs that are domain-unique
- Profiles: domain-specific position specific scoring matrices
- Hidden Markov Models(HMMs): probability-based models

Profiles, position-specific scoring matrices (PSSMs):\
Table of position-specific scores and gap penalties based on a multiple sequence aligment. Specific matrix for each position in multiple sequence alignment. PSSM is a dynamic, domain-specific matrix

PSI-BLAST(Position-Specific Iterative BLAST):\
While standard BLAST uses a fixed substitution matrix (such as BLOSUM62), PSI-BLAST automatically creates a position-specific scoring matrix (PSSM) for your query sequence and continuously updates it over multiple iterations. This allows starting with closely related homologs, it gradually captures more distantly related proteins (those with lower sequence similarity but related functions).

Profile Hidden Markov Models(HMMs):
The Profile HMM is a probabilistic model specifically designed to identify patterns in protein or DNA sequences and simulate evolutionary changes (substitutions, insertions, and deletions) that occur in a set of homologous sequences. In other words, while a PSSM only tells us "which amino acid is common at each position," an HMM can also tell us "what is likely to occur (insertions, deletions, or matches) from one position to the next in the sequence."

## Assignment
Comment of the assignment:\
Annotate two sequences using: SMART, PROSITE, Pfam/Interpro\
PSI-BLAST with RNAse 3 domain and check multiple sequence alignment for 'pattern'




## Lecture

## Assignment
Comment of the assignment:\

## Lecture

## Assignment
Comment of the assignment:\

## Lecture

## Assignment
Comment of the assignment:\

## Lecture

## Assignment
Comment of the assignment:\

## Lecture

## Assignment
Comment of the assignment:\

## Lecture

## Assignment
Comment of the assignment:\



## Lecture

## Assignment
Comment of the assignment:\




## Lecture

## Assignment
Comment of the assignment:\


## Lecture

## Assignment
Comment of the assignment:\


# List of all the bioinformatics tools
| Theme | Exercise | Tool | Introduction | Used for (in the course) | Type of data | Notes |
|---|---|---|---|---|---|---|
| Building Blocks of Life | 1A | [Sequence Manipulation Suite](https://www.bioinformatics.org/sms2/index.html) | SMS is for generating, formatting, and analyzing short DNA and protein sequences | Find the ORF | DNA/RNA/Protein sequence | Format conversion, sequence analysis, sequence figures |
| The OMICS | 2A | GALAXY | Galaxy is an open source, web-based platform for data intensive biomedical research | Genome Assembly (SPAdes), line/word/character count | FASTQ | — |
| The OMICS | 2B | [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example) | Bowtie2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences | Map reads against reference genome | FASTQ and ref_genome | in GALAXY environment |
| The OMICS | 2B | [Jbrowse](https://jbrowse.org/jb2/) |A genome browser that can run on the web,desktop, or embedded in application | display various genetic information on website and easy to share with others |  |  |
| The OMICS | 3A | [Ensembl](https://www.ensembl.org/index.html) |A public and open project providing access to genomes, annotations, tools and methods | here we can download transcipts data of yeast etc. |  |  |
| The OMICS | 3A | [yeastgenome](https://www.yeastgenome.org/genomesnapshot) |Saccharomyces cerevisiae Genome Snapshot | here we can check the information of yeast genome |  |  |
| The OMICS | 3A | [IGV](https://igv.org/app/) | A browser for visualize| use the IGV browser to visualize transcript result | .bam |  |
| The OMICS | 3A | [SMART](https://smart.embl.de/) |Protein prediction| predict protein structure from DNA sequence | .bam |  |
| The OMICS | 3B | [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) |performs functional enrichment analysis, also known as over-representation analysis (ORA) or gene set enrichment analysis| see the gene which have different expression and its function | |  |
| The OMICS | 3B | [SPELL](https://spell.yeastgenome.org/) |SPELL (Serial Pattern of Expression Levels Locator) is a query-driven search engine for large gene expression microarray compendia| see the gene different expression under different condition in different published article | |  |
| The OMICS | 4A | [MZmine3](https://github.com/mmattano/mzmine3) |mzmine is an open-source and platform-independent software for mass spectrometry (MS) data processing and visualization| | | in nuvolos environment |
| The OMICS | 4A | [MSViewer](https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msviewer) |Proteomics Data Visualization and Comparison| | | |
| The OMICS | 4B | [OMMSA](https://ssb-omssa.containers.wur.nl/) | The Open Mass Spectrometry Search Algorithm [OMSSA] is an efficient search engine for identifying MS/MS peptide spectra by searching libraries| | | |
| Sequence-defined properties | 7 | [SMART](https://smart.embl.de/) | An online resource for the identification and annotation of protein domains and the analysis of protein domain architectures| | | |
