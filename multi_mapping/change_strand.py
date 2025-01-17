import pysam 
import sys

def change_strand(input_bam, output_bam):
    with pysam.AlignmentFile(input_bam, "rb") as input_bamfile:
        header = input_bamfile.header
        with pysam.AlignmentFile(output_bam, "wb", header=header) as output_bamfile:
            for read in input_bamfile:
                if read.is_reverse:
                    read.is_reverse = False
                else:
                    read.is_reverse = True
                output_bamfile.write(read)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python change_strand.py <input_bam> <output_bam>")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_bam = sys.argv[2]
    change_strand(input_bam, output_bam)