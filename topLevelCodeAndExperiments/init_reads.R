library("Biostrings")
#R.utils::sourceDirectory("/Volumes/Personal/dabrams/thesis/polyester-master/R") #path for colby macs
#R.utils::sourceDirectory("/dabrams/thesis/polyester-master/R") #path for personal computer to Personal on colby server
R.utils::sourceDirectory("./polyester-master/R")#path for local copy

fold_changes = matrix(c(4,4,rep(1,18),1,1,4,4,rep(1,16)), nrow=20)

# FASTA annotation
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)

# subset the FASTA file to first 20 transcripts
small_fasta = fasta[1:20]


writeXStringSet(small_fasta, 'chr22_small.fa')

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(small_fasta) / 100)

# simulation call:

simulate_experiment('./data/input/chr22_small.fa', reads_per_transcript=readspertx, 
                    num_reps=c(10,10), fold_changes=fold_changes, outdir='simulated_reads') 

