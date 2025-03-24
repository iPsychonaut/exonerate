/****************************************************************\
*                                                                *
*  Protein <-> Genome comparison model                           *
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

#include <string.h> /* For strlen() */

#include "protein2genome.h"
#include "phase.h"

Protein2Genome_Data *Protein2Genome_Data_create(
                         Sequence *query, Sequence *target){
    register Protein2Genome_Data *p2gd = g_new0(Protein2Genome_Data, 1);
    g_assert(query->alphabet->type == Alphabet_Type_PROTEIN);
    g_assert(target->alphabet->type == Alphabet_Type_DNA);
    Protein2DNA_Data_init(&p2gd->p2dd, query, target);
    if(!Protein2Genome_Data_get_Intron_Data(p2gd))
        Protein2Genome_Data_get_Intron_Data(p2gd)
                              = Intron_Data_create();
    return p2gd;
    }

void Protein2Genome_Data_destroy(Protein2Genome_Data *p2gd){
    if(Protein2Genome_Data_get_Intron_Data(p2gd)){
        Intron_Data_destroy(Protein2Genome_Data_get_Intron_Data(p2gd));
        Protein2Genome_Data_get_Intron_Data(p2gd) = NULL;
        }
    Protein2DNA_Data_clear(&p2gd->p2dd);
    g_free(p2gd);
    return;
    }

/**/

C4_Model *Protein2Genome_create(Affine_Model_Type type){
    g_message("Starting Protein2Genome_create with type %d", type);
    register C4_Model *model = Protein2DNA_create(type);
    register C4_Transition *match_transition;
    register C4_Model *phase_model;
    register Match *match;
    register gchar *name = g_strdup_printf("protein2genome:%s",
                                           Affine_Model_Type_get_name(type));
    g_assert(model);
    g_message("Renaming model to %s", name);
    C4_Model_rename(model, name);
    g_free(name);
    g_message("Opening model for modification");
    C4_Model_open(model);
    g_message("Selecting MATCH transition");
    match_transition = C4_Model_select_single_transition(model,
                                                 C4_Label_MATCH);
    g_assert(match_transition);
    /* Add phased intron model */
    match = match_transition->label_data;
    g_message("Creating phase model");
    phase_model = Phase_create(NULL, match, FALSE, TRUE);
    g_message("Inserting phase model");
    C4_Model_insert(model, phase_model, match_transition->input,
                                        match_transition->output);
    g_message("Destroying phase model");
    C4_Model_destroy(phase_model);
    /**/
    g_message("Closing model");
    C4_Model_close(model);
    g_message("Completed Protein2Genome_create");
    return model;
}

