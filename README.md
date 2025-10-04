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

- 







# List of all the bioinformatics tools
| Theme | Exercise | Tool | Introduction | Used for (in the course) | Type of data | Notes |
|---|---|---|---|---|---|---|
| Building Blocks of Life | 1A | [Sequence Manipulation Suite](https://www.bioinformatics.org/sms2/index.html) | SMS is for generating, formatting, and analyzing short DNA and protein sequences | Find the ORF | DNA/RNA/Protein sequence | Format conversion, sequence analysis, sequence figures |
| The OMICS | 2A | GALAXY | Galaxy is an open source, web-based platform for data intensive biomedical research | Genome Assembly (SPAdes), line/word/character count | FASTQ | — |
| The OMICS | 2B | [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example) | Bowtie2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences | Map reads against reference genome | FASTQ and ref_genome | in GALAXY environment |
| The OMICS | 2B | [Jbrowse]((https://jbrowse.org/jb2/) |A genome browser that can run on the web,desktop, or embedded in application | display various genetic information on website and easy to share with others |  |  |
| The OMICS | 3A | [Ensembl](https://www.ensembl.org/index.html) |A public and open project providing access to genomes, annotations, tools and methods | here we can download transcipts data of yeast etc. |  |  |
| The OMICS | 3A | [yeastgenome](https://www.yeastgenome.org/genomesnapshot) |Saccharomyces cerevisiae Genome Snapshot | here we can check the information of yeast genome |  |  |
| The OMICS | 3A | [IGV](https://igv.org/app/) | A browser for visualize| use the IGV browser to visualize transcript result | .bam |  |
| The OMICS | 3A | [SMART](https://smart.embl.de/) |Protein prediction| predict protein structure from DNA sequence | .bam |  |
