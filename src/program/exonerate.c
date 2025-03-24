/****************************************************************\
*                                                                *
*  exonerate : a generic sequence comparison tool                *
*                                                                *
*  Guy St.C. Slater..   mailto:guy@ebi.ac.uk                     *
*  Copyright (C) 2000-2009.  All Rights Reserved.                *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU General Public License, version 3. See the file COPYING   *
*  or http://www.gnu.org/licenses/gpl.txt for details            *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#include "argument.h"
#include "analysis.h"
#include "gam.h"
#include "optimal.h"
#include "codegen.h"
#include "fastadb.h"

#include "est2genome.h"
#include "affine.h"
#include "ner.h"
#include "intron.h"
#include "heuristic.h"
#include "protein2dna.h"
#include "frameshift.h"
#include "alphabet.h"
#include "hspset.h"
#include "match.h"
#include "alignment.h"
#include "sar.h"
#include "sdp.h"
#include "splice.h"

static Analysis *Analysis_create(GPtrArray *query_path_list, Alphabet_Type query_type,
                                 gint query_chunk_id, gint query_chunk_total,
                                 GPtrArray *target_path_list, Alphabet_Type target_type,
                                 gint target_chunk_id, gint target_chunk_total,
                                 gint verbosity) {
    g_print("Starting Analysis_create\n"); fflush(stdout);
    register Analysis *analysis = g_new0(Analysis, 1);
    g_print("Allocated Analysis struct\n"); fflush(stdout);
    analysis->query_path_list = query_path_list;
    analysis->target_path_list = target_path_list;
    analysis->query_type = query_type;
    analysis->target_type = target_type;
    analysis->query_chunk_id = query_chunk_id;
    analysis->query_chunk_total = query_chunk_total;
    analysis->target_chunk_id = target_chunk_id;
    analysis->target_chunk_total = target_chunk_total;
    analysis->verbosity = verbosity;
    /* Add remaining initialization from original Analysis_create... */
    g_print("Loading query sequence\n"); fflush(stdout);
    analysis->query = Sequence_get_file(g_ptr_array_index(query_path_list, 0), query_type, FALSE);
    g_print("Loaded query sequence\n"); fflush(stdout);
    g_print("Loading target sequence\n"); fflush(stdout);
    analysis->target = Sequence_get_file(g_ptr_array_index(target_path_list, 0), target_type, FALSE);
    g_print("Loaded target sequence\n"); fflush(stdout);
    /* Model setup would follow here in full implementation */
    g_print("Completed Analysis_create\n"); fflush(stdout);
    return analysis;
}

int Argument_main(Argument *arg){
    register Analysis *analysis;
    register ArgumentSet *as_input =
        ArgumentSet_create("Sequence Input Options");
    GPtrArray *query_path_list, *target_path_list;
    Alphabet_Type query_type, target_type;
    gint query_chunk_id, target_chunk_id,
         query_chunk_total, target_chunk_total,
         verbosity;
    /**/
    g_print("Starting Argument_main\n"); fflush(stdout);
    ArgumentSet_add_option(as_input, 'q', "query", "path",
    "Specify query sequences as a fasta format file", NULL,
    NULL, &query_path_list);
    ArgumentSet_add_option(as_input, 't', "target", "path",
    "Specify target sequences as a fasta format file", NULL,
    NULL, &target_path_list);
    /**/
    ArgumentSet_add_option(as_input, 'Q', "querytype",
    "alphabet type", "Specify query alphabet type", "unknown",
    Alphabet_Argument_parse_alphabet_type, &query_type);
    ArgumentSet_add_option(as_input, 'T', "targettype",
    "alphabet type", "Specify target alphabet type", "unknown",
    Alphabet_Argument_parse_alphabet_type, &target_type);
    /**/
    ArgumentSet_add_option(as_input, '\0', "querychunkid", NULL,
    "Specify query job number", "0",
    Argument_parse_int, &query_chunk_id);
    ArgumentSet_add_option(as_input, '\0', "targetchunkid", NULL,
    "Specify target job number", "0",
    Argument_parse_int, &target_chunk_id);
    ArgumentSet_add_option(as_input, '\0', "querychunktotal", NULL,
    "Specify total number of query jobs", "0",
    Argument_parse_int, &query_chunk_total);
    ArgumentSet_add_option(as_input, '\0', "targetchunktotal", NULL,
    "Specify total number of target jobs", "0",
    Argument_parse_int, &target_chunk_total);
    /**/
    ArgumentSet_add_option(as_input, 'V', "verbose", "level",
    "Show search progress", "1",
    Argument_parse_int, &verbosity);
    /**/
    g_print("Setting up argument sets\n"); fflush(stdout);
    Argument_absorb_ArgumentSet(arg, as_input);
    Translate_ArgumentSet_create(arg);
    Analysis_ArgumentSet_create(arg);
    FastaDB_ArgumentSet_create(arg);
    GAM_ArgumentSet_create(arg);
    Viterbi_ArgumentSet_create(arg);
    Codegen_ArgumentSet_create(arg);
    Heuristic_ArgumentSet_create(arg);
    SDP_ArgumentSet_create(arg);
    BSDP_ArgumentSet_create(arg);
    /**/
    Sequence_ArgumentSet_create(arg);
    Match_ArgumentSet_create(arg);
    Seeder_ArgumentSet_create(arg);
    Affine_ArgumentSet_create(arg);
    NER_ArgumentSet_create(arg);
    Intron_ArgumentSet_create(arg);
    Frameshift_ArgumentSet_create(arg);
    Alphabet_ArgumentSet_create(arg);
    HSPset_ArgumentSet_create(arg);
    Alignment_ArgumentSet_create(arg);
    SAR_ArgumentSet_create(arg);
    Splice_ArgumentSet_create(arg);
    /**/
    g_print("Processing arguments\n"); fflush(stdout);
    Argument_process(arg, "exonerate",
      "A generic sequence comparison tool\n"
      "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2009.\n",
      "\n"
      "Examples of use:\n"
      "\n"
      "1. Ungapped alignment of any DNA or protein sequences:\n"
      "    exonerate queries.fa targets.fa\n"
      "2. Gapped alignment of Mouse proteins to Fugu proteins:\n"
      "    exonerate --model affine:local mouse.fa fugu.fa\n"
      "3. Find top 10 matches of each EST to a genome:\n"
      "    exonerate --model est2genome --bestn 10 est.fa genome.fa\n"
      "4. Find proteins with at least a 50% match to a genome:\n"
      "    exonerate --model protein2genome --percent 50 p.fa g.fa\n"
      "5. Perform a full Smith-Waterman-Gotoh alignment:\n"
      "    exonerate --model affine:local --exhaustive yes a.fa b.fa\n"
      "6. Many more combinations are possible.  To find out more:\n"
      "    exonerate --help\n"
      "    man exonerate\n"
      "\n");
    if(verbosity > 0)
        Argument_info(arg);
    /**/
    g_print("Calling Analysis_create\n"); fflush(stdout);
    analysis = Analysis_create(query_path_list, query_type,
                               query_chunk_id, query_chunk_total,
                               target_path_list, target_type,
                               target_chunk_id, target_chunk_total,
                               verbosity);
    g_print("Calling Analysis_process\n"); fflush(stdout);
    Analysis_process(analysis);
    g_print("Destroying Analysis\n"); fflush(stdout);
    Analysis_destroy(analysis);
    /**/
    if(verbosity > 0)
        g_print("-- completed exonerate analysis\n");
    g_print("Completed Argument_main\n"); fflush(stdout);
    return 0;
}
