import subprocess
import sys


def sam2bam(path):
    """
    Converts a sam file to a bam file
    @arg path:  Inputted path is in sam format
    @returns:   Outputted file is in bam format
    """
    bam_file = path.replace('.sam', '.bam')
    handle = open(bam_file, 'w')
    subprocess.check_call(['samtools', 'view', '-bS', path], stdout=handle)
    handle.close()
    return bam_file


def bam2sam(path):
    """
    Converts a bam file to a sam file
    @arg path:  Inputted path is in bam format
    @returns:   Outputted file in sam format
    """
    sam_file = path.replace('.bam', '.sam')
    handle = open(sam_file, 'w')
    subprocess.check_call(['samtools', 'view', '-h', '-o', sam_file, path],
                          stdout=handle)
    handle.close()
    return sam_file


def align_sam(index, path):
    """
    Aligns a fastq file to a reference genome
    @arg index: Reference index file with .ebwt1 extension
    @arg path:  Path is in fastq format
    @returns:   Aligned file in sam format
    """
    sam_file_aligned = path.replace('.fastq', '.sam')
    handle = open(sam_file_aligned, 'w')
    # Put the bowtie 1 command together...
    # -t                    Print out timestamps
    # -S                    SAM output
    # -p 2                  Number of cores to use for alignment
    p = subprocess.Popen(['bowtie', '-t', '-S', '-p 2', index, path],
                         stdout=subprocess.PIPE)
    # remove binary b produced by subprocess.Popen
    for line in p.stdout:
        ln = line.decode('ascii')
        handle.write(ln)

    handle.close()
    return sam_file_aligned


def sort_bam(bam_file_path):
    """
    Sorts aligned bam file by leftmost coordinates
    @arg bam_file_path: Calls on aligned bam file
    @returns:  Sorted bam file outputted
    """
    sorted_bam = bam_file_path.replace(".bam", "_sorted")
    p = subprocess.Popen(['samtools', 'sort', bam_file_path, sorted_bam],
                         stdout=subprocess.PIPE)
    # remove binary b produced by subprocess.Popen
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)

    return sorted_bam


def remove_unmapped(sorted_bam_path):
    """
    Filters out unmapped reads, secondary and
    supplemental alignments from the bam file
    @arg sorted_bam_path: Calls on sorted aligned bam file
    @returns:  Mapped sorted bam file outputted
    """
    sorted_bam = sorted_bam_path + ".bam"
    mapped_sorted_bam = sorted_bam.replace(".bam", "_mapped.bam")
    # Put the samtools view command together...
    # -b                 Output as bam
    # -F 904             Filters out unmapped, secondary and supplemental alignments
    # -o                 Designates output name
    p = subprocess.Popen(['samtools', 'view', '-b', '-F 0x904', '-o',
                          mapped_sorted_bam, sorted_bam_path], stdout=subprocess.PIPE)
    # remove binary b produced by subprocess.Popen    
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)
    return mapped_sorted_bam


sample_file = sys.argv[1]
ref_index = sys.argv[2]


# Step 1 :
print("Aligning fastq file to reference genome...")
sam_file = align_sam(ref_index, sample_file)


# Step 2 :
print("Converting sam to bam file...")
bam_file = sam2bam(sam_file)


# Step 3 :
print("Sorting bam file...")
sorted_bam_file = sort_bam(bam_file)


# Step 4 :
print("Removing unmapped reads, secondary alignments and "
      "supplemental alignments from the bam files...")
mapped_sorted_bam_file = remove_unmapped(sorted_bam_file + ".bam")
  # samtools sort function adds
  # the .bam extension so you MUST add bam here


# Step 5 :
print("Converting back to sam...")
mapped_sorted_sam_file = bam2sam(mapped_sorted_bam_file)


# Step 6 :
print("Determining the read counts for genes of interest...")

# build a dictionary with the Accession ID of genes of interest
GeneCounts = {
    'NM_019388.3': 0,
    'NM_009855.2': 0,
    'NM_011611.2': 0,
    'NM_021396.2': 0,
    'NM_021893.3': 0,
    'NM_015790.3': 0,
    'NM_007548.4': 0,
    'NM_001109661.1': 0,
    'NM_009744.4': 0
}

handle = open(mapped_sorted_sam_file, 'rU')

for line in handle:
    # find the Accession ID within the sam file
    if line.startswith('SRR'):
        genome_key = line.split('\t')[2]
        if genome_key in GeneCounts.keys():
            GeneCounts[genome_key] = GeneCounts[genome_key] + 1
        else:
            continue

handle.close()


# Step 7:
print("Writing out the dictionary")
# Ensures the output is written to your current working directory
gene_file = mapped_sorted_sam_file.split("/")
# Need to change the csv output file name to match sample name
gene_file[-1] = 'gene_counts.csv'
gene_file = "/".join(gene_file)

# Create an output file
output = open(gene_file, 'w')
header = 'RefID' + ',' + 'ReadCount' # Header for table created
output.write(header)
output.write('\n')
for ID, count in GeneCounts.items(): #iterate over the dictionary
    output.write(ID + ',' + str(count))
    output.write('\n')
output.close()

print("Done!! - You did it!!")
