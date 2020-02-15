# NAD-tagSeq

1. Combine fastq files (pass & fail) to one fastq file.
  In Windows OS CMD: (after changed to the file folder path) (file_folder_path>)copy (file_front)_*.fastq (mixed_file_name).fastq
  In Linux OS: cat (file_identical).fastq > (mixed_file_name).fastq

2. Find out the RNA with tag in the first 40 nt:
In Linux OS:
  Download main.py from my Git-Hub repository: https://github.com/rocketjishao/NAD-tagSeq/blob/master/main.py
  Install python (version 2.7.17): (http://ubuntuhandbook.org/index.php/2017/07/install-python-3-6-1-in-ubuntu-16-04-lts/)
      $ sudo apt-get update
      $ sudo apt-get install python
  Change directory to the file pathway for main.py; 
  Sort out the RNAs with tag RNA sequence, meanwhile cut off the first 40 nt by typing in:
      $ python main.py input_file.fastq output_file.fastq
        (Result file: tagged.fastq (as an example))
        
3. Minimap2 to analyze the RNA sequenced from Nanopore Direct RNA Sequencing:
In Linux OS:
  Download minimap2. (https://github.com/lh3/minimap2)
  Run Minimap2 by typing in:
    ./minimap2 -ax splice -uf -k14 reference.fa tagged.fastq > output.sam

4. Samtools to change the sam file to bam file and obtain its bam.bai file.
In Linux OS:
  Download samtools. (https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856)
  Run Samtools by typing in (one by one):
    samtools view -hb output.sam > output.bam 
    samtools sort -O BAM -T output.bam.temp -o output_sort.bam  output.bam
    samtools index output_sort.bam output_sort.bam.bai

5. IGV to visualize the result
In Windows OS:
  Download IGV: (https://software.broadinstitute.org/software/igv/download)
  Download genome file from IGV for A. thaliana, human, mouse, or E.coli: Genome > Load Genome from Server > Select the genome file
  Import the bam and bam.bai to Windows OS, then: File > Load from File > Select the output.bam file
  
6. Use htseq-count to do statistical analysis of each RNA of tagged and untagged counterparts in ADPRC- and ADPRC+ groups:

