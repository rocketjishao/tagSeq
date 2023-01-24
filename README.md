[![DOI](https://zenodo.org/badge/240644522.svg)](https://zenodo.org/badge/latestdoi/240644522)
# Note  
### This is used for processing CoA tagSeq and NAD tagSeq data.
### The demo files are mainly for Arabidopsis and data analysis of NAD tagSeq protocol   

# Table of content
- [Softwares and code](#softwares-and-code)
- [Software installation and configuration](#software-installation-and-configuration)
  * [python2.7 and python3.6](#python27-and-python36)
  * [Miniconda3](#miniconda3)
  * [pycoQC](#pycoqc)
  * [Minimap2](#minimap2)
  * [featureCounts](#featurecounts)
  * [Samtools](#samtools)
  * [IGV for Linux OS](#igv-for-linux-os)
  * [ont_fast5_api:](#ont-fast5-api-)
  * [ont-tombo:](#ont-tombo-)
- [Demo files](#demo-files)
- [tagSeq data analysis procedure](#tagseq-data-analysis-procedure)
  * [1. Run pycoQC in MiniConda3 active virtual environment](#1-run-pycoqc-in-miniconda3-active-virtual-environment)
  * [2. Combine fastq files to one fastq file](#2-combine-fastq-files-to-one-fastq-file)
  * [3. Sort out the RNA with and without tag in the first 50 nt](#3-sort-out-the-rna-with-and-without-tag-in-the-first-50-nt)
  * [4. Minimap2 to align the reads to reference sequence](#4-minimap2-to-align-the-reads-to-reference-sequence)
  * [5. Use featureCounts to count the aligned reads to genes](#5-use-featurecounts-to-count-the-aligned-reads-to-genes)
  * [6. Samtools to translate the sam file to bam file and obtain its bam.bai file](#6-samtools-to-translate-the-sam-file-to-bam-file-and-obtain-its-bambai-file)
  * [7. IGV to visualize the RNA structure](#7-igv-to-visualize-the-rna-structure)
  * [8. Electronic signal analysis for tagged model CoA-RNA](#8-electronic-signal-analysis-for-tagged-model-coa-rna)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


# Softwares and code
 (1) MinKNOW 19.6.8, with base-caller of Guppy embedded, from Oxford Nanopore Technology  
 (2) Ubuntu 18.04.3 LTS, Linux-based operating system (https://ubuntu.com/download)  
***The following code packages should be installed on Ubuntu**  
(3) Python 2.7 and 3.7 (http://www.python.org/downloads/)  
(4) Miniconda3 (Miniconda3-latest-Linux-x86_64.sh)(https://dos.conda.io/projects/conda/en/latest/user-guide/install/linux.html) for pycoQC uses;             
(5) PycoQC (https://github.com/a-slide/pycoQC) to analyze the basecalling results and data quality;  
(6) Homemade python script to sort out tagged and untagged RNA  (https://github.com/rocketjishao/tagSeq/blob/master/main.py)  
(7) Minimap2 (https://github.com/lh3/minimap2) to align the sequenced RNA to genome or transcriptome databases for interpretation of the RNA identities;  
(8) featureCounts 1.6.0 (http://bioinf.wehi.edu.au/featureCounts/) to map and count the reads of tagged RNA to genes in different samples.  
(9) Samtools 1.7 (http://samtools.sourceforge.net/) to translate the sam file to bam file and obtain its bam.bai file;  
(10) Integrative Genomics Viewer 2.7.2 (https://software.broadinstitute.org/software/igv/) to visualize the RNA structures;  
  
# Software installation and configuration
## python2.7 and python3.6
    
    $ sudo apt-get install python2
          # Then type in password
    $ sudo apt-get install python-pip 
          # or try $ python get-pip.py

    $ sudo add-apt-repository ppa:jonathonf/python-3.6
          # Then type in password, or try $ sudo apt-get install python3

## Miniconda3 
(https://conda.io/projects/conda/en/latest/user-guide/install/linux.html):  
   a. Download the installer: Miniconda installer for Linux.(https://docs.conda.io/en/latest/miniconda.html#linux-installers)  
   b. Verify your installer hashes, in a terminal window enter:  
        
    $ sha256sum Downloads/Miniconda3-latest-Linux-x86_64.sh
   c.In your terminal window, run Miniconda:  
        
    $ bash Downloads/Miniconda3-latest-Linux-x86_64.sh
   d.Follow the prompts on the installer screens.  
   e.If you are unsure about any setting, accept the defaults. You can change them later.  
   f.To make the changes take effect, type in the command below, or close and then re-open your terminal window. 
    
    $ source ~/.bashrc
       # Then you can see "(base)" in the front of the terminal command line. 
   g.Test your installation. In your terminal window, run the command:
   
    $ conda list
   h. A list of installed packages appears if it has been installed correctly.  
   i. To change the automatic conda activation because *auto_activate_base* is set to True. You can check this using the following command
             
    $ conda config --show | grep auto_activate_base
     
 Note: for activating miniconda3, type in
    $ source ~/miniconda3/bin/activate
     for inactivation, 
    $ source ~/.bashrc


## pycoQC   
(https://github.com/tleonardi/pycoQC)  
   a. Create a clean virtual environment (only needed for the 1st run):  

    $ conda create -n pycoQC python=3.6
    # Note: python 2 is not supported by pycoQC

   b. Install pycoQC with miniconda3:(only needed for the 1st run)  

    $ conda install -c aleg pycoqc

   c. Run pycoQC by the command:  

    $ pycoQC -f sequencing_summary.txt -o pycoQC_output.html
   d. Quit conda 

    $ conda deactivate
   e. To enter and exit conda for 2nd, 3rd,... run  
   To change the automatic conda activation because *auto_activate_base* is set to True. You can check this using the following command
             
    $ conda config --show | grep auto_activate_base
          
   To set it false

    $ conda config --set auto_activate_base False
    $ source ~/.bashrc
   To reactivate set it to True

    $ conda config --set auto_activate_base True
    $ source ~/.bashrc

## Minimap2 
(https://github.com/lh3/minimap2) :    
    
    $ git clone https://github.com/lh3/minimap2
    $ cd minimap2 && make
  
## featureCounts 
(http://subread.sourceforge.net/):  
    
    $ sudo apt-get install subread 
        # then type in password 

## Samtools 
(https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856)  
(http://www.htslib.org/download/)
      
    $ sudo apt-get install -y samtools
        # then type in password

## IGV for Linux OS  
   Download the IGV file: https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_Linux_2.8.0.zip;  
   Unzip the package;  
   In the terminal window, start IGV by the command line:

    $ java --module-path=lib -Xmx4g @igv.args --module=org.igv/org.broad.igv.ui.Main
   Download genome file from IGV for A. thaliana, human, mouse, or E.coli:   
          Genome > Load Genome from Server > Select the genome file  

## ont_fast5_api:
(https://github.com/nanoporetech/ont_fast5_api):  
   
    $ pip install ont-fast5-api


## ont-tombo:
(https://github.com/nanoporetech/tombo)
   Note: To install ont-tombo, conda is required and a virtual environment is recommened. 
   a. In your terminal window, run anaconda3:

    $ source ~/anaconda3/etc/profile.d/conda.sh
               $ conda activate 
               ## deactivate conda: $ conda deactivate
   b. Create a virtual environment called Python2 where python version is 2
    
    $ conda create --name python2 python=2
    $ source activate python2

   c. Tombo install:

    $ conda install -c bioconda ont-tombo

# Demo files			
       
|Step|sofware|input_files|output_files| demo files |
|---|---|---|---| ---|
|1| pycoQC | [sequencing_summary.txt]() | [pycoQC.html](https://rawcdn.githack.com/rocketjishao/NAD-tagSeq/37433efcfd6add36e27a77e0124571326b6ec05d/pycoQC.html) ([raw data](https://github.com/rocketjishao/NAD-tagSeq/blob/master/pycoQC.html)) | no |   
|2| Windows OS CMS | fastq files (ADPRC+_1.fastq,ADPRC+_2.fastq,ADPRC+_3.fastq; ADPRC-.fastq| ADPRC+.fastq, ADPRC-.fastq| [demo](https://github.com/rocketjishao/NAD-tagSeq/blob/Arabidopsis-only/Demo-total%20RNA_raw%20data.tar.xz)|  
|3| main.py | ADPRC+.fastq; ADPRC-.fastq | ADPRC+_tagged.fastq; ADPRC+_untagged.fastq; ADPRC-_tagged.fastq; ADPRC-_untagged.fastq |[deomo](https://github.com/rocketjishao/NAD-tagSeq/blob/Arabidopsis-only/Demo-tag_sorted.tar.xz)|   
|4| minimap2 | ADPRC+_tagged.fastq, ADPRC+_untagged.fastq; ADPRC-_tagged.fastq; ADPRC-_untagged.fastq; reference_file ([A. thaliana TAIR10.fas](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas)) | ADPRC+_tagged.sam, ADPRC+_untagged.sam; ADPRC-_tagged.sam; ADPRC-_untagged.sam |[demo](https://github.com/rocketjishao/NAD-tagSeq/blob/Arabidopsis-only/Demo-minimap2.tar.xz)|   
|5| featureCounts | ADPRC+_tagged.sam; ADPRC+_untagged.sam; ADPRC-_tagged.sam; ADPRC-_untagged.sam; annotation file ([TAIR10.gtf](ftp://ftp.ensemblgenomes.org/pub/release-46/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gtf.gz)) | all; all.summary |[demo](https://github.com/rocketjishao/NAD-tagSeq/blob/Arabidopsis-only/Demo-featureCounts.tar.xz)|  
|6| samtools | ADPRC+_tagged.sam | ADPRC+_tagged.bam; ADPRC+_tagged_sort.bam;[ADPRC+_tagged_sort.bam.bai |[demo](https://github.com/rocketjishao/NAD-tagSeq/blob/Arabidopsis-only/Demo-samtools.tar.xz)|   
|7| IGV | ADPRC+_tagged_sort.bam; ADPRC+_tagged_sort.bam.bai; ADPRC+_untagged_sort.bam; ADPRC+_untagged_sort.bam.bai; genome files (mm10.genome) | IGV figure |no|  




# tagSeq data analysis procedure

## 1. Run pycoQC in MiniConda3 active virtual environment
   To visualize the summary file generated from the sequencing and do the quality control analysis of the basecalling results:  
   Type in the command below. Open the html file with web browser to visualize the results.   
     
     $ pycoQC –f sequencing_summary.txt –o pycoQC.html

## 2. Combine fastq files to one fastq file  
   In Windows OS CMD:  
       
       $ copy ADPRC+_*.fastq ADPRC+.fastq
   In Linux OS: 
    
       $ cat ADPRC+_*.fastq > ADPRC+.fastq

## 3. Sort out the RNA with and without tag in the first 50 nt
   Download main.py from our Git-Hub repository: https://github.com/rocketjishao/NAD-tagSeq/blob/master/main.py  
   Change directory to the file pathway of main.py; 
   Sort out the RNAs with and without tag RNA sequence by typing in:
        
       $ python main.py ADPRC+.fastq ADPRC+_tagged.fastq ADPRC+_untagged.fastq
          # result files: ADPRC+_tagged.fastq and ADPRC+_untagged.fastq
        
## 4. Minimap2 to align the reads to reference sequence   
   Run Minimap2 for analyzing the Nanopore direct RNA sequencing data by typing in the command:
        Note: GNU make (C compiler) and zlib are needed.
        ## zlib: sudo apt-get install zlib1g-dev
        ## GNU make: $ sudo apt install build-essential 
              check gcc version $ gcc --version
       $ ./minimap2 -ax splice -uf -k14 reference.fa ADPRC+_tagged.fastq > ADPRC+_tagged.sam
          # reference file like TAIR10.fa, result file is ADPRC+_tagged.sam

## 5. Use featureCounts to count the aligned reads to genes
   Use simultaneously the tagged and untagged counterparts (or map each gene to the tagged RNA in ADPRC- and ADPRC+ samples.)  
   And download gene annotation files in gtf format from Ensembl or GenBank (https://www.ncbi.nlm.nih.gov/genbank/), avoid UCSC  
   Run the command below:  
        
       $ featureCounts -L -a annotation -o all ADPRC+_tagged.sam ADPRC+_untagged.sam ADPRC-_tagged.sam ADPRC-_untagged.sam
          # annotation file like TAIR10.gff, result files are all and all.summary

## 6. Samtools to translate the sam file to bam file and obtain its bam.bai file  
   Run Samtools by typing in (one by one):
    
       $ samtools view -bS ADPRC+_tagged.sam > ADPRC+_tagged.bam 
       $ samtools sort -O BAM -o ADPRC+_tagged_sort.bam ADPRC+_tagged.bam
       $ samtools index ADPRC+_tagged_sort.bam
          # result files: ADPRC+_tagged.bam, ADPRC+_tagged_sort.bam, ADPRC+_tagged_sort.bam.bai
       $ samtools stats ADPRC+_tagged.bam | grep '^SN' | cut -f 2-  
          # use this to visualize the # mismatches / bases mapped (cigar), which should be smaller than 0.25, indicating dismatched bases account for <20% and matched bases >80%
## 7. IGV to visualize the RNA structure  
   Import the bam and bam.bai to IGV by:   
          File > Load from File > Select the ADPRC+_tagged_sort.bam file
  
## 8. Electronic signal analysis for tagged model CoA-RNA
   a. To analyze the electronic signal of tagged CoA-RNA, firstly you need to get the fast5 file and figure out in which fast5 file the tagged CoA-RNA was located. In IGV, get the ID of the tagged model CoA-RNA (here it is "18cf1707...."), then sort out the tagged reads and copy their information (first two lines) to a text file by using grep, therefore you can see which fast5 file the tagged model CoA-RNA was in:  
        
    $ grep "18cf1707" ~/C44pos/44pos.tag.fastq -A 2 > tag1.txt
   b. Copy the fast5 file to the folder of "ont_fast5_api". Because a fast5 file contains 4000 reads, firstly separate multiple fast5 reads into single reads:
  
    $ multi_to_single_fast5 --input_path ~/ont_fast5_api/multi_reads --save_path ~/ont_fast5_api/single_readsc

   c. Transfer the fast5 files to another folder, like "'~/C44/fast5". repeat steps a-c to find out more fast5 files for tagged model CoA-RNA.
   
   d. Then run the command below, to requiggle electronic signal to according the sequence "~/minimap2-/model.fa":
     
     # tombo resquiggle ~/C44/fast5 ~/minimap2-/model.fa --processes 4 --num-most-common-errors 5

   e. Use the following commands to analyze the resultant mapping of electronic signals to tagged model CoA-RNA sequence:
    
     # tombo detect_modifications de_novo --fast5-basedirs /C44/fast5 --statistics-file-basename de_novo_testing --processes 4

     # tombo plot genome_locations --fast5-basedirs ~/C44/fast5 --genome-locations model:50 --num-bases 100

