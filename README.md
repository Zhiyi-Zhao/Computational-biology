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
[Central dogma](https://www.genome.gov/genetics-glossary/Central-Dogma). The central dogma of molecular biology is a theory stating that genetic information flows only in one direction, from DNA, to RNA, to protein, or RNA directly to protein. Now its already expand like the picture below.
![Central Dogma](https://upload.wikimedia.org/wikipedia/commons/thumb/3/3a/Crick%27s_1958_central_dogma.svg/1024px-Crick%27s_1958_central_dogma.svg.png)

Protein is more abundant in creatures, people also first sequenced protein by Edmann degradation method. But now, DNA sequence is becoming cheaper and cheaper, the data for DNA sequence is also increasing in a fast way.

The [DNA sequencing method](https://www.geneious.com/guides/introduction-to-dna-sequencing), we will learn how to analyze sequencing data in the coming course, but its also recommended to learn where the data come from.\

- Genome\
In classical genetics, the genome of a diploid organism refers to a full set of chromosomes or genes in a gamete. A regular somatic cell contains two full sets of genomes(and a mitochondrial genome)\
In haploid organisms, including bacteria, archaea and viruses, a cell contains only a single set of the genome, usually in a single circular or contiguous linear DNA.\
Mitochondrial and Chloroplast genome are also present in single copy.\
Size of the genomes is variable and does not correlate with complexity of living beings. And small and simple does not mean harmless.

- Gene
A gene is a programmable unit that can give rise to a multitude of products: protein and RNA products through (alternative) splicing and trans-splicing

## Assignment
- Reading Frame(RF), Open Reading Frame(ORF) and Coding Sequence(CDS)
  - A **Reading Frame** is how a genetic sequence is read during a translation process.
  - An **Open Reading Frame** is the part of a reading frame that has the potential to code for a protein or peptide and consists of a continuous stretch of codons that do not contain a stop codon (usually UAA, UAG or UGA). Note that in this case, the ORF will be interrupted with intervening sequences or introns.
  - The **Coding Sequence**, is the portion of a gene's DNA or RNA, composed **only of exons**, that codes for protein. The region is bounded at the 5' end with a start codon and at the 3' end with a stopcodon.
 
- Intron
An intron is a nucleotide sequence within a gene in eukaryote. It is a noncoding sequence. During the final maturation of the RNA product, the introns are removed by splicing.

- Frameshift
Frameshift is a shifting of the reading frame caused by insert nuclides not devided by 3. Frameshift mutation is a genetic mutation caused by indels (insertions or deletions) of a number of nucleotides in a DNA sequence that is not divisible by three. 

Comment of the assignment:
Use the SMS website to find the ORF corresponding to the known peptide on the known DNA sequence. Noting:
- Continuous peptides do not equal continuous DNA sequences; they may be located in different ORFs and separated by introns.
- Knowing which DNA chain the ORF is located on means you don't have to worry about the ORF on the other chain.

| Theme | Exercise | Tool |Introduction |Used for | Type of data | Notes | Extension |
|-------|----------|------|-----|----------|--------------|-------|-----------|
|Buiding Blocks of Life|1A|[Sequence Manipulation Suite](https://www.bioinformatics.org/sms2/index.html)|SMS is for generating, formatting, and analyzing short DNA and protein sequences|Find the ORF|DNA/RNA/Protein sequence||Format conversion,sequence analysis,sequence figures|




