###################################################################################################################
# R script to simulate rna-seq reads from hg38 transcriptome reference using polyester simulator                  #
###################################################################################################################
# Usage: R < polyester.r --no-save
# For hg38 transcriptome, needs ~ 16G memory
# Output fa files ~20G (unzipped)


library(Biostrings)
library(polyester)

#Set number of repeats to simulate
num_repeats = 5

#Set a seed
the_seed = 1

# Path to cdna file. 
# Used gffread and gtf reference found at /sudlab1/.../annotations/hg38_noalt_ensembl85/ensembl.dir/ to produce "transcripts.fa"
fastapath <-  '../../refiles/transcripts.fa'


# Read transcripts reference file into a variable, transcripts
transcripts = readDNAStringSet(fastapath)


# Split the fasta into chunks (to help with memory issues). 
chunks = 50
transcript_count = length(transcripts)
assignment = rep(1:chunks, each=ceiling(transcript_count/chunks))[1:transcript_count]


# Create a loop that loops from 1 -> Chunks
for(i in 1:chunks){
	
      # Creating a "sub" fasta file containing the transcripts needed for the chunk
      fasta_sub = transcripts[assignment==i]
      writeXStringSet(fasta_sub, 'fasta_sub.fasta')

      # set number of reads as a function of the length of the transcript * a fold a coverage
      readspertx = round(20 * width(fasta_sub) / 100)

      # adding gc_bias to sample
      # Count the number of transcripts in the sub_fasta chunk
      numtx = count_transcripts('fasta_sub.fasta')
	
      # Create a "readmatrix" that will hold the read counts for all the transcripts and for each sample
      # Assign initial values as readspertx (set above)
      readmat = matrix(readspertx, ncol=num_repeats, nrow=numtx)
      
      # Set the seed as the number given before the loop. This ensures each sample recieves the same GC-bias model 	
      set.seed(the_seed)
	
      # Assign a number between 0-7 which corresponds to a built-in GC-bias model	
      biases = sample(0:7, num_repeats, replace=TRUE)
	
      #Assign GC-bias models to readcounts using add_gc_bias function to the readmatrix created above
      readmat_biased = add_gc_bias(readmat, as.list(biases), fasta_sub)  
      
      # Remove any values which were very small (approaching 0 reads) or were "N/A" in the readmatrix, as this causes 
      # errors later on. 	
      readmat_biased[readmat_biased<0.0000000001] <- 0.000000001
      readmat_biased[is.numeric(readmat_biased) & is.na(readmat_biased)] <- 0.00001

      # Output the read counts for the chunk to a tab seperated file
      write.table(readmat_biased, file=paste0('simulated_reads/chunk_',i,'_readmatrix.txt'), row.names=FALSE, col.names=FALSE)

      # Run the read_simulator to output reads for each sample. Reads are saved in seperate "chunk" folders to be cat'ed 
      # together after the simulations
      simulate_experiment_countmat(fasta='fasta_sub.fasta',
			  paired=TRUE,
			  readmat=readmat_biased,
			  num_reps=c(num_repeats),
			  outdir=paste0('simulated_reads/chunk_', i))
}
