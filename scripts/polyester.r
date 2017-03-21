library(Biostrings)
library(polyester)

#Set num of repeats
num_repeats = 4
the_seed = 1

# Path to cdna file
fastapath <-  '../../refiles/transcripts.fa'

# Read DNA file into var transcripts
transcripts = readDNAStringSet(fastapath)

# splitting the fasta into chunks based on the number of trnascripts
chunks = 50
transcript_count = length(transcripts)
assignment = rep(1:chunks, each=ceiling(transcript_count/chunks))[1:transcript_count]

# Loop to run simulation on each section of the fasta file
for(i in 1:chunks){
      # Creating a subsection of the fasta file 
      fasta_sub = transcripts[assignment==i]
      writeXStringSet(fasta_sub, 'fasta_sub.fasta')

      # set Coverage
      readspertx = round(30 * width(fasta_sub) / 100)

      # adding gc_bias to sample (see ?add_gc_bias)
      numtx = count_transcripts('fasta_sub.fasta')
      readmat = matrix(readspertx, ncol=num_repeats, nrow=numtx)
      set.seed(the_seed)
      biases = sample(0:7, num_repeats, replace=TRUE)
      readmat_biased = add_gc_bias(readmat, as.list(biases), fasta_sub)  
      readmat_biased[readmat_biased<0.0000000001] <- 0.000000001
      readmat_biased[is.numeric(readmat_biased) & is.na(readmat_biased)] <- 0.00001

      #filename <- c('simulated_reads/chunk_',i,'_readmatrix.txt')
      write.table(readmat_biased, file=paste0('simulated_reads/chunk_',i,'_readmatrix.txt'), row.names=FALSE, col.names=FALSE)

      simulate_experiment_countmat(fasta='fasta_sub.fasta',
			  paired=TRUE,
			  readmat=readmat_biased,
			  num_reps=c(num_repeats),
			  outdir=paste0('simulated_reads/chunk_', i))
}
