# NAD-tagSeq
1.
Combine fastq files (pass & fail) to one fastq file.
In Windows OS CMD: (after changed to the file folder path) (file_folder_path>)copy (file_front)_*.fastq (mixed_file_name).fastq
In Linux OS: cat (file_identical).fastq > file_name.fastq

2. 
Use main.py to find out the RNA with tag in the first 40 nt
Result file: tagged.fa (as an example)

3.
Download minimap2 in Linux OS. (https://github.com/lh3/minimap2)
Download samtools in Linux OS. (https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856)

4. 
Code for minimap2 running:


