# **9577Q Final Assignment**  
# Comparing gene expresion in murine WT and Blimp-1 knockout plasma cells

For this final assignment, everything that was done to obtain and analyze the sample RNAseq data is listed below for reference. Each script that was written assumed that the user was using python 3. 

This assignment was designed to analyze data of nitrophenol (NP) specific high and low affinity IgG1 light zone germinal center B cells. The sample groups and n number per each group is listed below. Briefly, high, and low affinity IgG1 LZ GC B cells were FACS-sorted from 10 mice and pooled for each genotype. Illumina HiSeq 1500 deep sequencing was then completed on these B cells. For this proposal, I will only be analyzing two of the sample groups one from the high affinity, and the other being a low affinity B cell. The aim of this assignment is to determine how to compare the gene expression of two different populations of germinal center (GC) B cells and identify the read counts of transcripts of interest that are expressed.

## Downloading RNASeq FASTQ data from NCBI SRA
In order to download the data I first had to download the sratoolkit. After unzipping the file I was then able to fastq-dump the sequence data of High affinity1 and Low affinity1 nitrophenol (NP) specific IgG1 light zone germinal center B cells. 

```
$ wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/$ sdk/current/sratoolkit.current-ubuntu64.tar.gz
$ gunzip sratoolkit.2.8.2-1-ubuntu64.tar.gz 
$ tar -xf sratoolkit.2.8.2-1-ubuntu64.tar 
$ cd sratoolkit.2.8.2-1-ubuntu64/
$ cd bin
$ ./fastq-dump
$ ./fastq-dump -v SRR2558581
$ ./fastq-dump -v SRR2558584
```

In the end I will have obtained my two samples that I plan on analyzing which are listed below. High affinity1 (SRR2558581) was renamed to HA1.fastq and Low affinity1 (SRR2558584) was renamed to LA1.fastq.

Samples  | Accession ID | Read length | File Size | Number of reads
------------- | ---------| ----------- | --------- | ---------------
High affinity1 | SRR2558581 | 49 | 295.3M | 8505150
Low affinity1 | SRR2558584 | 49 | 343.6M | 9976173

## Determine the number of reads and whether it is paired or unpaired data
First I wanted to see how many reads each sample I was working with contained. To do this I used one simple Linux Commands:
```shell
$ cat HA1.fastq | wc -l
  34020600 # divide by four to get the number of reads
$ cat LA1.fastq | wc -l
  39904692 # divide by four to get the number of reads
```

For the sake of convenience I have added this information to the above table. The next thing I did was check to make sure that my records contained unpaired reads. Although you can assume it is unpaired since there was only one fastq file per sample, I double checked just in case. To do this I simply looked for duplicate read ids using the grep function. 
```
To check my HA1.fastq file:
grep @SRR2558581.1 HWI-1KL133:205:4:2302:14166:79223 HA1.fastq 

To check my LA1.fastq file:
grep @SRR2558584.1 HWI-1KL133:205:4:1101:1178:2072 LA1.fastq
```

## FastQC Quality Control check for FASTQ files
FastQC is an application that can take FASTQ files and run a series of tests to generate a comprehensive QC report. Since this program is written in java I first needed to check to see if java was downloaded on my computer before proceeding. Once I ensured java was installed I then downloaded and unzipped fastqc. After that I ran fastqc on my two FASTQ samples. 
```shell 
$ java -version
  # java was not installed so I downloaded it 
$ tar zxvf jre-8u73-linux-x64.tar.gz # unpack tarball and install java
$ unzip fastqc_v0.11.5.zip
$ cd FastQC/
$ chmod +x fastqc
$ cd ..
$ ./FastQC/fastqc HA1.fastq 
$ ./FastQC/fastqc LA1.fastq 
```

This will create two html files (titled HA1_fastqc.html, and LA1_fastqc.html) that will contain the fastqc report. Based on what was outputted we can see that our quality control looks good on both our samples and that there are no adapter sequences we need to worry about trimming off. 

## Creating a reference index to input into bowtie
Bowtie is a short read aligner that is capable of aligning large sets of short sequence reads in fasta and fastq format to a reference genome. They have an abundant source of built in reference indexes you can download, but for the purpose of this assignment I created my own that would reference transcript of [Mus musculus strain C57BL/6J](https://www.ncbi.nlm.nih.gov/genome/52?genome_assembly_id=279711). 


To create my reference index, first I had to download bowtie. I used the command line to download all my files:
```shell
$ pip install bowtie
$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.25_GRCm38.p5/GCF_000001635.25_GRCm38.p5_rna.fna.gz
  # fasta transcript sequence
$ gunzip GCF_000001635.25_GRCm38.p5/GCF_000001635.25_GRCm38.p5_rna.fna.gz
$ bowtie-build GCF_000001635.25_GRCm38.p5_rna.fna GRCm38.p5
  # builds reference index
```

The reference index created will have the leading string 'GRCm38.p5' with '.1.ebwt' added to the end followed by '2.ebwt' and so on and so forth. When refering back to this reference genome all you have to imput is 'GRCm38.p5' and the bowtie function will pull all files that lead with that string. 

## Explaning readcount_workflow.py - Determining reads of interest genes
I then created my reads of intrest readcount_workflow. The goal was to map my two different samples of high and low affinity murine B cells to the mouse genome, and determine which genes related to germinal center (GC) B cell development were being expressed. I had a total of 7 steps that I applied to my sample fastq files one at a time, and my reference genome. I implemented the sys module in order to call those two files by defining a sample as the argument 1 and the reference genome to index as argument 2 repectively. 

```
   sample_file = sys.argv[1]
   ref_index = sys.argv[2]

   NOTE: ref_index is where you are defining your argument to the reference genome created. AGAIN remember, you need to only specify the name you outputted as your bowtie build (in this case 'GRCm38.p5'). 
```

**Step 1:** This function was created to align my two fastq files to the reference genome. To do that I called the bowtie function by importing subprocess. I set up the script so that when this first function was called, the screen would print out "Aligning fastq file to reference genome..." If the script crashed this would allow you to see it was during this *align_sam* function. 
```
NOTE: 
  - This is the first function used that will require that you import subprocess as we are calling out to bowtie. 
  - Bowtie will allow you to multithread, meaning you can use more then one processor. I only have two processors available for my Ubuntu so I therefore set my *-p* to **2**. If more are available it would be beneficial to change this so that the output could be produced faster. 
```
  
**Step 2:** Next I had to convert my file from sam to bam format in order to utilize the sort command. Samtools is able to do this so again I called out to samtools through subprocess. I made sure to pipline my aligned sam output from my *align_sam* function created from the first step as the input for this function. The script was created to print out "Converting sam to bam file..." while this function is running. If the script crashed this again would allow you to see it was during this sam2bam function. 

**Step 3 :** After a bam file was generated I was then able to utilize the samtools sort function in order to sort my bam file according to the query template name (QNAME, in Col 1). So I piped the *sam2bam* output and inputted it into my *sort_bam* function. The script was created to print out "Sorting bam file..." so again you would know if there was a problem with this step. 
```
This was a tricky one to figure out since the bowtie sort function adds in the .bam extension on its own. To by pass this you had to make sure that when you called the input from this function in step 4 you called it (+ ".bam"). 
This was the benefit of having the print comments appear, since I was able to see that my file was looking for a *sorted_bam_file* but my file was actually called *sorted_bam_file.bam*
```

**Step 4 :** Next I wanted to filter out all the unmapped reads, secondary and supplementary alignments. This would leave you with only mapped, primary aligned reads. So i piped my output from *sort_bam* and inputted it into my *remove_unmapped* function. The screen would print, "Removing unmapped reads, secondary alignments and supplemental alignments from the bam files..." so you know this was the function that was currently running. 

**Step 5 :** To be able to understand and look at your mapped and aligned sequences you have to convert the bam file back to sam. Again I utilized samtools and piped the output from my mapped sorted bam file that was created from the *remove_unmapped* function and inputted that into this *bam2sam* function. The screen was set to print "Converting back to sam..." while this occured. 

**Step 6 :** Finally I wanted to determine the read counts for particular genes of interest. To do this first I had to look up the Accession ID of my genes of interest. Next I created a dictionary  which contained the Accession ID of my genes of interest as the key, and then set the value of each key to be 0. 

Next I had to pipe my mapped and sorted sam file that was outputted from my *bam2sam* file. Using this file I then scanned through looking for the start of a sam line (idicated by the start of the QNAME). In my case all my QNAME's started with SRR. I then looked in position 2 to find out my RNAME which is where the Accession ID's are stored. Therefore if my genes of interest were present in the mapped sorted sam file I was set to add one to the value of those specific keys. While this function ran the screen was set to print "Determining the read counts for genes of interest". 

For reference the gene name of my genes of interest are listed in the table below. These were chosen as they are specifically expressed on GC B cells. 

**Gene name** | **Accession ID** 
--------------| ----------- 
CD86 | NM_019388.3 
CD80 | NM_009855.2
CD40 | NM_011611.2 
PDL2 | NM_021396.2 
PDL1 | NM_021893.3 
ICOS-L | NM_015790.3 
Blimp1 | NM_007548.4
Bach2 | NM_001109661.1 
Bcl6 | NM_009744.4

**Step 7:** The final step that I need to do is write out the dictionary so that I had an output with my produced table. To do this I had to first create a file that would locate me in my current directory. Next I had to set my output to write to that fie, and output a header that would help readers understand what was in each column, and then the key and values produced from my dictionary. While this was running my screen was set to print "Writing out the dictionary". 
