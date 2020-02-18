#Simplified procedure for NAD-tagSeq data analysis

1. Combine fastq files (pass & fail) to one fastq file.
In Windows OS CMD:  
    copy (file_name)_*.fastq mixed.fastq
In Linux OS: 
    cat (file_name)_*.fastq > mixed.fastq

2. Sort out the RNA with and without tag in the first 40 nt:
In Linux OS:
    Download main.py from our Git-Hub repository: https://github.com/rocketjishao/NAD-tagSeq/blob/master/main.py
    Install python (version 2.7.17): (http://ubuntuhandbook.org/index.php/2017/07/install-python-3-6-1-in-ubuntu-16-04-lts/) 
        $ sudo apt-get install python
        #$ python get-pip.py # pip install, optional
    Change directory to the file pathway for main.py; 
    Sort out the RNAs with and without tag RNA sequence, meanwhile cut off the first 40 nt of the tagged RNA by typing in:
        $ python main.py input_file.fastq tagged.fastq untagged.fastq
          # result file: tagged.fastq (as an example) and untagged.fastq
        
3. Minimap2 to analyze the RNA sequenced from Nanopore Direct RNA Sequencing:
In Linux OS:
    Install minimap2. (https://github.com/lh3/minimap2)
    Run Minimap2 by typing in:
        ./minimap2 -ax splice -uf -k14 reference.fa tagged.fastq > output.sam
          # result file: output.sam

4. Samtools to change the sam file to bam file and obtain its bam.bai file.
In Linux OS:
    Intall samtools. (https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856)
    Run Samtools by typing in (one by one):
        samtools view -bS output.sam > output.bam 
        samtools sort -O BAM -o output_sort.bam  output.bam
        samtools index output_sort.bam output_sort.bam.bai
        # result files: output.bam, output_sort.bam, output_sort.bam.bai

5. IGV to visualize the result
In Windows OS:
    Download IGV: (https://software.broadinstitute.org/software/igv/download)
    Download genome file from IGV for A. thaliana, human, mouse, or E.coli: Genome > Load Genome from Server > Select the genome file
    Import the bam and bam.bai to Windows OS, then: File > Load from File > Select the output.bam file
  
6. Use featureCounts to count each gene to the RNA reads of tagged and untagged counterparts, or to the tagged RNA in ADPRC- and ADPRC+ samples.
In Linux OS:
    And download gene annotation files in gtf format from Ensembl or GenBank (https://www.ncbi.nlm.nih.gov/genbank/), avoid UCSC
    Install featureCounts (http://subread.sourceforge.net/): sudo apt-get install subread 
    Run the command below:
        featureCounts -L -a gencode.vM23.annotation.gtf -o both tagged.sam untagged.sam
        # result files: both and both.summary
