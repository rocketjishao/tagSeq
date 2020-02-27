# Table of content

- [Demo files](#demo-files)
- [Softwares and code](#softwares-and-code)
- [Software installation and initiation](#software-installation-and-initiation)
- [NAD-tagSeq data analysis procedure](#nad-tagseq-data-analysis-procedure)


# Demo files			
       
|sofware|input_files|output_files|
|---|---|---|
| pycoQC | [sequencing_summary.txt](https://github.com/rocketjishao/NAD-tagSeq/blob/master/Rapiflex-PC_20191225_182450_FAL15529_minion_sequencing_run_1_sequencing_summary.tar.xz) | [pycoQC.html](https://github.com/rocketjishao/NAD-tagSeq/blob/master/pycoQC.html) ([(web browser)](https://rawcdn.githack.com/rocketjishao/NAD-tagSeq/37433efcfd6add36e27a77e0124571326b6ec05d/pycoQC.html)) |    
| Windows OS CMS | file_\*.fastq | mixed.fastq|  
| main.py | mixed.fastq | tagged.fastq; untagged.fastq |   
| minimap2 | tagged.fastq; reference_file [(A. thaliana TAIR10.fas](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas) or [mouse mm10.fa](https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz))  | tagged.sam |   
| samtools | tagged.sam | tagged.bam; tagged_sort.bam; tagged_sort.bam.bai |   
| IGV | genome files (TAIR10.genome or mm10.genome); tagged_sort.bam; tagged_sort.bam.bai | no |  
| featureCounts | annotation file ([TAIR10](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff) or [mouse](https://www.gencodegenes.org/mouse/)); tagged.sam; untagged.sam | both; both.summary |  



# Softwares and code
 (1) MinKNOW 19.6.8, with base-caller of Guppy embedded, from Oxford Nanopore Technology  
 (2) Ubuntu 18.04 LTS, Linux-based operating system (https://ubuntu.com/download)  
***The following code packages should be installed on Ubuntu18.04, except Integrative Genomics Viewer:***  
(3) Python 2.7 and 3.7 (http://www.python.org/downloads/)  
(4) Minicode3 (https://dos.conda.io/projects/conda/en/latest/user-guide/install/linux.html)             
(5) PycoQC (https://github.com/a-slide/pycoQC)  
(6) Homemade python script to sort out tagged and untagged RNA  (https://github.com/rocketjishao/NAD-tagSeq/blob/master/main.py)  
(7) Minimap2 (https://github.com/lh3/minimap2) to align the sequenced RNA to genome or transcriptome databases for interpretation of the RNA identities;  
(8) Samtools (http://samtools.sourceforge.net/) to translate the sam file to bam file and obtain its bam.bai file;  
(9) Integrative Genomics Viewer (https://software.broadinstitute.org/software/igv/) to visualize the RNA structures;  
(10) featureCounts (http://bioinf.wehi.edu.au/featureCounts/) to count the reads of tagged RNA in different samples.  

# Software installation and initiation
(1) python2.7 and python3.6 (Installed on Ubuntu:)
    
    $ sudo apt-get install python2.
          # Then type in password
    $ sudo apt-get install python-pip 
          # or try $ python get-pip.py

    $ sudo add-apt-repository ppa:jonathonf/python-3.6
          # Then type in password, or try $ sudo apt-get install python3

(2) Minicode3 (https://conda.io/projects/conda/en/latest/user-guide/install/linux.html), installed on Ubuntu18.04:  
    Download the installer: Miniconda installer for Linux.(https://docs.conda.io/en/latest/miniconda.html#linux-installers)  
    Verify your installer hashes, in a terminal window enter:  
        
    $ sha256sum Downloads/Miniconda3-latest-Linux-x86_64.sh
   In your terminal window, run Miniconda:  
        
    $ bash Downloads/Miniconda3-latest-Linux-x86_64.sh
   Follow the prompts on the installer screens.  
   If you are unsure about any setting, accept the defaults. You can change them later.  
   To make the changes take effect, type in the command below, or close and then re-open your terminal window. 
    
    $ source ~/.bashrc
       # Then you can see "(base)" in the front of the terminal command line. 
   Test your installation. In your terminal window, run the command:
   
    $ conda list
   A list of installed packages appears if it has been installed correctly.  
   To change the automatic conda activation because *auto_activate_base* is set to True. You can check this using the following command
             
    $ conda config --show | grep auto_activate_base
          
   To set it false

    $ conda config --set auto_activate_base False
    $ source ~/.bashrc
   To reactivate set it to True

    $ conda config --set auto_activate_base True
    $ source ~/.bashrc
      


(3) pycoQC: (https://a-slide.github.io/pycoQC/installation/)  
a. Create a clean virtual environment:  

    $ conda create -n pycoQC python=3.6

b. Install pycoQC with miniconda3:  

    $ conda install -c aleg pycoqc

c. Run pycoQC by the command:  

    $ pycoQC -f sequencing_summary.txt -o pycoQC_output.html
d. Quit conda 

    $ conda deactivate

(4) Minimap2 installation (https://github.com/lh3/minimap2) :    
    
    $ git clone https://github.com/lh3/minimap2
    $ cd minimap2 && make
    
(5) Samtools installation (https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856)
      
    $ sudo apt-get install -y samtools
        # then type in password
(6) Install featureCounts (http://subread.sourceforge.net/):  
    
    $ sudo apt-get install subread 
        # then type in password 

# NAD-tagSeq data analysis procedure

1. In MiniConda3 active virtual environment, 
   Run pycoQC to visualize the summary file generated from the sequencing and do the quality control analysis of the basecalling results:  
   Type in the command below. Open the html file with web browser to visualize the results.   
    
       $ pycoQC –f sequencing_summary.txt –o pycoQC.html

2. Combine fastq files (pass & fail) to one fastq file.  
    In Windows OS CMD:  
    
       $ copy file_*.fastq mixed.fastq
    In Linux OS: 
    
       $ cat file_*.fastq > mixed.fastq

3. Sort out the RNA with and without tag in the first 40 nt:
   Download main.py from our Git-Hub repository: https://github.com/rocketjishao/NAD-tagSeq/blob/master/main.py  
   Change directory to the file pathway for main.py; 
   Sort out the RNAs with and without tag RNA sequence by typing in:
        
       $ python main.py mixed.fastq tagged.fastq untagged.fastq
          # result files: tagged.fastq (as an example) and untagged.fastq
        
4. Minimap2 to analyze the RNA sequenced from Nanopore Direct RNA Sequencing: (In Linux OS):  
   Run Minimap2 by typing in:
        
       $ ./minimap2 -ax splice -uf -k14 reference.fa tagged.fastq > tagged.sam
          # reference file like TAIR10_chr_all.fa, result file is tagged.sam

5. Samtools to translate the sam file to bam file and obtain its bam.bai file. (In Linux OS):  
   Run Samtools by typing in (one by one):
    
       $ samtools view -bS tagged.sam > tagged.bam 
       $ samtools sort -O BAM -o tagged_sort.bam tagged.bam
       $ samtools index tagged_sort.bam
          # result files: tagged.bam, tagged_sort.bam, tagged_sort.bam.bai

6. IGV to visualize the result, in Windows OS:
    Download IGV: (https://software.broadinstitute.org/software/igv/download)
    Download genome file from IGV for A. thaliana, human, mouse, or E.coli:   
          Genome > Load Genome from Server > Select the genome file
    Import the bam and bam.bai to Windows OS, then:   
          File > Load from File > Select the tagged.bam file
  
7. Use featureCounts to count each gene to the RNA reads of tagged and untagged counterparts, or map each gene to the tagged RNA in ADPRC- and ADPRC+ samples. (In Linux OS:)  
   And download gene annotation files in gtf format from Ensembl or GenBank (https://www.ncbi.nlm.nih.gov/genbank/), avoid UCSC  
   Run the command below:  
        
       $ featureCounts -L -a annotation_file -o both tagged.sam untagged.sam
          # annotation file like TAIR10_GFF3_genes.gff, result files are both and both.summary
