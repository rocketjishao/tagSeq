import os
import sys

# notify user
print("Please install minimap2 and python2, and download main.py code, along with genome files (model.fa & mm10.fa)")
print("Please include genome annotation files (lnc.gtf, mouse2.gtf, model.gtf")

# check if the correct number of arguments is passed
if len(sys.argv) != 2:
    print("Usage: python tagSeq.py fileA.fastq")
    sys.exit(1)

# get the file name from command line argument
fileA = sys.argv[1]

# check if the file exists
if not os.path.isfile(fileA):
    print("File not found!")
    sys.exit(1)

# run the commands
print("\r\n....1.start sorting out tagged and untagged reads\r\n")
os.system("python2 main.py " + fileA + " "+ fileA + ".tag.fastq " + fileA + ".untag.fastq")
print("\r\n....2.minimap2 to map tagged reads to tagged model CoA-RNA\r\n")
os.system("minimap2 -ax splice -uf -k14 model.fa " + fileA + ".tag.fastq > " + fileA + ".tag.model.sam")
print("\r\n....3.minimap2 to map tagged reads to mouse genome\r\n")
os.system("minimap2 -ax splice -uf -k14 mm10.fa " + fileA + ".tag.fastq > " + fileA + ".tag.mouse.sam")
print("\r\n....4.minimap2 to map untagged reads to tagged model CoA-RNA\r\n")
os.system("minimap2 -ax splice -uf -k14 model.fa " + fileA + ".untag.fastq > " + fileA + ".untag.model.sam")
print("\r\n....5.minimap2 to map untagged reads to mouse genome\r\n")
os.system("minimap2 -ax splice -uf -k14 mm10.fa " + fileA + ".untag.fastq > " + fileA + ".untag.mouse.sam")
print("\r\n....6.samtools to convert .sam files to generate .bam and .sort.bam files\r\n")
os.system("samtools view -bS " + fileA + ".tag.model.sam > " + fileA + ".tag.model.bam")
os.system("samtools sort -O BAM -o " + fileA + ".tag.model.sort.bam " + fileA + ".tag.model.bam")
os.system("samtools index " + fileA + ".tag.model.sort.bam")
os.system("samtools view -bS " + fileA + ".untag.model.sam > " + fileA + ".untag.model.bam")
os.system("samtools sort -O BAM -o " + fileA + ".untag.model.sort.bam " + fileA + ".untag.model.bam")
os.system("samtools index " + fileA + ".untag.model.sort.bam")
os.system("samtools view -bS " + fileA + ".tag.mouse.sam > " + fileA + ".tag.mouse.bam")
os.system("samtools sort -O BAM -o " + fileA + ".tag.mouse.sort.bam " + fileA + ".tag.mouse.bam")
os.system("samtools index " + fileA + ".tag.mouse.sort.bam")
os.system("samtools view -bS " + fileA + ".untag.mouse.sam > " + fileA + ".untag.mouse.bam")
os.system("samtools sort -O BAM -o " + fileA + ".untag.mouse.sort.bam " + fileA + ".untag.mouse.bam")
os.system("samtools index " + fileA + ".untag.mouse.sort.bam")
print("\r\n....7.featureCounts to map and count mapped reads to annotated gtf files (mouse, lncRNA, model.gtf)\r\n")
os.system("featureCounts -L -a model.gtf -o model " + fileA + ".tag.model.sam " + fileA + ".untag.model.sam")
os.system("featureCounts -L -a lnc.gtf -o lnc " + fileA + ".tag.mouse.sam " + fileA + ".untag.mouse.sam")
os.system("featureCounts -L -a mouse.gtf -o mouse " + fileA + ".tag.mouse.sam " + fileA + ".untag.mouse.sam")
