##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################

from ruffus import *
from ruffus.combinatorics import product
import sys
import os
import sqlite3
import glob
import CGAT.Experiment as E
import CGAT.Sra as Sra
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineRnaseq as RnaSeq

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))

PARAMS["project_src"]=os.path.dirname(__file__)

# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.
RnaSeq.PARAMS = PARAMS

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

STRINGTIE_QUANT_FILES=["i_data.ctab", "e_data.ctab", "t_data.ctab",
                       "i2t.ctab", "e2t.ctab"]


# ---------------------------------------------------
@follows(mkdir("assembled_transcripts.dir"), mkdir("portcullis"))
@transform(["input_assemble.dir/*.bam",
            "input_assemble.dir/*.remote"],
           formatter(),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"])),
           "assembled_transcripts.dir/{basename[0]}.gtf.gz")
def assembleWithStringTie(infiles, outfile):

    infile, reference = infiles
    basefile = os.path.basename(infile)
    job_threads = PARAMS["stringtie_threads"]
    job_memory = PARAMS["stringtie_memory"]
    tmpfile = P.getTempFilename()

    statement =  '''
                    portcullis full 
                            -t 1
                            -o portcullis/%(basefile)s/
                            -r %(portcullis_bedref)s
                            -b 
                            %(portcullis_fastaref)s
                            %(infile)s;
                    checkpoint;
                    mv portcullis/%(basefile)s/portcullis.filtered.bam %(tmpfile)s;
		    rm -r portcullis/%(basefile)s/;
                    stringtie %(tmpfile)s
                           -p %(stringtie_threads)s
                           -G <(zcat %(reference)s)
                           %(stringtie_options)s
                           2> %(outfile)s.log
                   | gzip > %(outfile)s;
                   checkpoint;
                   rm %(tmpfile)s'''

    if infile.endswith(".remote"):
        token = glob.glob("gdc-user-token*")
        tmpfilename = P.getTempFilename()
        if len(token) > 0:
            token = token[0]
        else:
            token = None

        s, infile = Sra.process_remote_BAM(
            infile, token, tmpfilename,
            filter_bed=os.path.join(
                PARAMS["annotations_dir"],
                PARAMS["annotations_interface_contigs_bed"]))

        infile = " ".join(infile)
        statement = "; checkpoint ;".join(
            ["mkdir %(tmpfilename)s",
             s,
             statement,
             "rm -r %(tmpfilename)s"])

    P.run()


# ---------------------------------------------------
@follows(mkdir("final_genesets.dir"))
@merge([assembleWithStringTie,
       os.path.join(
           PARAMS["annotations_dir"],
           PARAMS["annotations_interface_geneset_all_gtf"])],
       "final_genesets.dir/agg-agg-agg.gtf.gz")
def mergeAllAssemblies(infiles, outfile):

    infiles = ["<(zcat %s)" % infile for infile in infiles]
    infiles, reference = infiles[:-1], infiles[-1]

    job_threads = PARAMS["stringtie_merge_threads"]
    job_memory = "12G"

    infiles = " ".join(infiles)

    statement = '''stringtie --merge
                             -G %(reference)s
                             -p %(stringtie_merge_threads)s
                             %(stringtie_merge_options)s
                             %(infiles)s
                            2> %(outfile)s.log
                   | python %(scriptsdir)s/gtf2gtf.py --method=sort
                           --sort-order=gene+transcript
                            -S %(outfile)s -L %(outfile)s.log'''

    P.run() 
 
@collate(assembleWithStringTie,
         regex("(.+)/(.+)-(.+)-(.+).gtf.gz"),
         add_inputs(os.path.join(
            PARAMS["annotations_dir"],
            PARAMS["annotations_interface_geneset_all_gtf"])),
         r"\1/\2-agg-agg.gtf.gz")
def merge_by_tissue(infiles, outfile):

    reference = "<(zcat %s)" % infiles[0][0]
    infiles = ["<(zcat %s)" % infile[0] for infile in infiles]

    job_threads = PARAMS["stringtie_merge_threads"]

    infiles = " ".join(infiles)

    statement = '''stringtie --merge
                             -G %(reference)s
                             -p %(stringtie_merge_threads)s
                             %(stringtie_merge_options)s
                             %(infiles)s
                            2> %(outfile)s.log
                   | python %(scriptsdir)s/gtf2gtf.py --method=sort
                           --sort-order=gene+transcript
                            -S %(outfile)s -L %(outfile)s.log'''

    P.run() 

@follows(mergeAllAssemblies, merge_by_tissue)
def Assembly():
    pass


# ---------------------------------------------------
@transform([assembleWithStringTie, mergeAllAssemblies, merge_by_tissue],
           suffix(".gtf.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"])),
           ".class.gz")
def classifyTranscripts(infiles, outfile):
    '''classify transcripts.
    '''
    to_cluster = True

    infile, reference = infiles

    counter = PARAMS['gtf2table_classifier']

    job_memory = "16G"

    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/gtf2table.py
           --counter=%(counter)s
           --reporter=transcripts
           --gff-file=%(reference)s
           --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()


# ---------------------------------------------------
@merge(classifyTranscripts,
       "transcript_class.load")
def loadTranscriptClassification(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).class.gz",
                         options="-i transcript_id -i gene_id"
                         "-i match_gene_id -i match_transcript_id"
                         "-i source --quick")


# ---------------------------------------------------
@follows(mkdir("utron_beds.dir"),classifyTranscripts)
@subdivide([assembleWithStringTie, mergeAllAssemblies, merge_by_tissue],
           regex("(.+)/(.+).gtf.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"]),
                      r"\1/\2.class.gz"),
           [r"utron_beds.dir/\2.all_utrons.bed.gz",
            r"utron_beds.dir/\2.partnered_utrons.bed.gz",
            r"utron_beds.dir/\2.novel_utrons.bed.gz"])
def find_utrons(infiles, outfiles):

    infile, reference, classfile = infiles

    job_memory="6G"

    all_out, part_out, novel_out = outfiles

    track = P.snip(all_out, ".all_utrons.bed.gz")

    statement = '''python %(scriptsdir)s/gtf2gtf.py -I %(infile)s
                             --method=sort
                             --sort-order=gene+transcript
                              -L %(track)s.log
                 | python %(project_src)s/utrons/find_utrons.py
                             --reffile=%(reference)s
                             --class-file=%(classfile)s
                             --outfile %(all_out)s
                             --partfile=%(part_out)s
                             --novel-file=%(novel_out)s
                              -L %(track)s.log'''

    P.run()


# ---------------------------------------------------
@transform(find_utrons,
           suffix(".bed.gz"),
           ".ids.gz")
def getUtronIds(infile, outfile):

    statement = '''zcat %(infile)s 
                 | cut -f 4
                 | sed 's/:/\\t/g'
                 | sort -u
                 | gzip > %(outfile)s'''

    P.run()


# ---------------------------------------------------
@collate(getUtronIds,
         regex(".+/(.+)\.(.+).ids.gz"),
         r"\2_ids.load")
def loadUtronIDs(infiles, outfile):

    header = "track,transcript_id"
    options = "-i track -i transcript_id"

    if not outfile == "all_utrons_ids.load":
        header += ",match_transcript_id"
        options += "-i match_transcript_id"

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)\..+\.ids.gz",
                         has_titles=False,
                         cat="track",
                         header=header,
                         options=options)


@follows(loadUtronIDs, loadTranscriptClassification)
def AnnotateAssemblies():
    pass


@follows(mkdir("export/indexed_gtfs.dir"))
@transform([mergeAllAssemblies,
            assembleWithStringTie],
           regex(".+/(.+).gtf.gz"),
           r"export/indexed_gtfs.dir/\1.gtf.gz")
def exportIndexedGTFs(infile, outfile):

    statement = '''zcat %(infile)s
                 | sort -k1,1 -k4,4n
                 | bgzip > %(outfile)s;

                 checkpoint;

                 tabix -p gff %(outfile)s'''

    P.run()


@follows(exportIndexedGTFs)
def export():
    pass


# ---------------------------------------------------
@follows(mkdir("salmon_index"),exportIndexedGTFs)
@transform("final_genesets.dir/*.gtf.gz",
           formatter(),
           "salmon_index/{basename[0]}.salmon.index")
def makeSalmonIndex(infile,outfile):
    # Long transcripts cause indexing to use lots of memory?
    job_memory="64G"
    job_threads=1

    gtf_basename = P.snip(os.path.basename(infile), ".gtf.gz")
    transcript_fasta = "salmon_index/" + gtf_basename + "transcripts.fa"
    fastaref =PARAMS["portcullis_fastaref"]
    index_options=PARAMS["salmon_indexoptions"]
    tmpfile = P.getTempFilename()
    


    statement = '''
    gunzip -c %(infile)s > %(tmpfile)s;
    gffread %(tmpfile)s -g %(fastaref)s -w %(transcript_fasta)s;
    salmon index
      -p %(job_threads)s
      %(index_options)s
      -t %(transcript_fasta)s
      -i %(outfile)s
      --perfectHash;
    rm %(tmpfile)s
    '''
    P.run()


#---------------------------------------------------------



@follows(mkdir("quantification.dir"),
         mergeAllAssemblies, makeSalmonIndex)
@product(["input_quantify.dir/*.bam",
          "input_quantify.dir/*.remote"],
         formatter(".+/(?P<TRACK>.+).(?:bam|remote)"),
         "final_genesets.dir/*.gtf.gz",
         formatter(".+/(?P<GENESET>.+).gtf.gz"),
         "quantification.dir/{TRACK[0][0]}_{GENESET[1][0]}")

def quantifyWithSalmon(infiles, outfile):
    '''Quantify existing samples against genesets'''
    job_threads=2
    job_memory="8G"

    infile, gtffile = infiles
    basefile = os.path.basename(infile)
    gtfbase = P.snip(os.path.basename(gtffile), ".gz")
    salmonIndex = "salmon_index/" + gtfbase + ".salmon.index"
    fastq1 = P.snip(outfile, "_agg-agg-agg")+".1.fastq"
    fastq2 = P.snip(outfile, "_agg-agg-agg")+".2.fastq"
    salmon_options=PARAMS["salmon_quantoptions"]

    statement = '''
    samtools fastq
         -1 %(fastq1)s
         -2 %(fastq2)s
         %(infile)s; 
    salmon quant -i %(salmonIndex)s
        --libType A
        -1 %(fastq1)s
        -2 %(fastq2)s
        -o %(outfile)s
        --threads %(job_threads)s
        %(salmon_options)s; 
    checkpoint; 
    mv %(outfile)s/quant.sf %(outfile)s.sf; 
    rm %(fastq1)s; rm %(fastq2)s 
    '''

    if infile.endswith(".remote"):
        token = glob.glob("gdc-user-token*")
        filename = "temp_bams/%s" % basefile
        tmpfilename = P.getTempFilename()
        if len(token) > 0:
            token = token[0]
        else:
            token = None

        s, infile = Sra.process_remote_BAM(
            infile, token, filename,
            filter_bed=os.path.join(
                PARAMS["annotations_dir"],
                PARAMS["annotations_interface_contigs_bed"]))

        infile = " ".join(infile)
        statement = "; checkpoint; ".join(
            ["mkdir %(filename)s",
             s,
             statement,
             "rm -r %(filename)s"])



    P.run()

#----------------------------------------------------
# ***********
# NOTE - last time I used this, the output database wasn't indexed correctly
# ***********
@follows(quantifyWithSalmon)
@merge("quantification.dir/*.sf", "salmon_quant.load")
def mergeAllQuants(infiles, outfile):
    job_memory="6G"
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="quantification.dir/(.*-.*-.*)_agg-agg-agg.sf",
                         options="-i Name -i Length -i EffectiveLength"
                         "-i TPM -i NumReads -i track"
                         "-i source --quick")



# ---------------------------------------------------
# Generic pipeline tasks
@follows(mergeAllQuants,
         AnnotateAssemblies,
         export,
         mergeAllQuants)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
