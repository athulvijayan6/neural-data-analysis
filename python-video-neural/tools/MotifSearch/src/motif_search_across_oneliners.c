
/*
Matches all the refs in the reffile ionfo with all the queries in the queyfile info. This is different from what we used to have as an input we a <query-ref filelist>. This is done for time saving whihc is wasted in loading the file.

%%%%%%%%%%%%%%%%
Perticularly It is assumed here that query list and ref list are the same files because it is meant to test across ragas.
Similarlly this code is tuned for specific needs. All this tuning is commented by *****. Uncommenr them for general use.
%%%%%%%%%%%%%%%%

Written by: Shrey Dutta
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "rlcs.h"
#include "types.h"
#include "maxmin.h"

void motif_search_in_reference(const SegmentPTR motif_ptr,const SegmentPTR reference_ptr, const char *  ctrl_file, SpotVectorPTR spotted_potential_motif_regions_ptr);

int main(int argc, char* argv[]) {
  
  int i,j; // simple iterators.
  
  if (argc != 4+1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: %s <ctrl_file> <query_fileinfo> <ref-fileinfo> <out_filename>\n", argv[0]);
    exit(-1);
  }

  char *  ctrl_file = argv[1];
  //  char *  query_reference_list = argv[2];
  char *  query_list = argv[2];
  char *  reference_list = argv[3];
  char *  out_filename = argv[4];


  FILEPTR fquery = fopen(query_list, "r");
  assert(fquery != NULL);
  FILEPTR fref = fopen(reference_list, "r");
  assert(fref != NULL);
  

  int match_num_query,
      match_num_reference;
  
  fscanf(fquery, "%d", &match_num_query);
  fscanf(fref, "%d", &match_num_reference);
 

  assert(match_num_query == match_num_reference); // **** comment this when query and ref list are different ***
  int match_num = (match_num_query *  match_num_reference);
  
  match_num = match_num -  match_num_query; // **** comment  when query and ref list are different ***
  
  char query_path[1000][1000];
  char reference_path[1000][1000];
  
  for (i = 0; i < match_num_query; i++){
    fscanf(fquery, "%s", query_path[i]);
    printf( "%s\n", query_path[i]);
  }
  fclose(fquery);

  for (i = 0; i < match_num_reference; i++){
    fscanf(fref, "%s", reference_path[i]);
    printf( "%s\n", reference_path[i]);
  }
  fclose(fref);

  
  SpotVectorPPTR spotted_potential_motif_regions_all_pptr =  (SpotVectorPPTR) malloc(sizeof(SpotVectorPTR)*match_num);
   
  // creating and allocating spot vectors
  for(i = 0; i < match_num; i++){
    spotted_potential_motif_regions_all_pptr[i] = create_SpotVectorPTR();
    // simply allocating 10. It will be rellocated inside if necessary.
    alloc_members_SpotVectorPTR(spotted_potential_motif_regions_all_pptr[i], 10);
  }

  // loading the data.
  printf("Loading the data....\n");
  SegmentPTR query_ptr[match_num_query];
  SegmentPTR reference_ptr[match_num_reference];

  printf("Loading the queries...\n");
  for (i =0; i< match_num_query; i++) {
    printf("%s...loading...\n", query_path[i]);
    query_ptr[i] = create_SegmentPTR();
    fill_SegmentPTR_from_file(query_path[i], query_ptr[i]);
  }
  printf("Loading the reference...\n");
  for (i =0; i < match_num_reference; i++) {
    printf("%s...loading...\n", reference_path[i]);
    reference_ptr[i] = create_SegmentPTR();
    fill_SegmentPTR_from_file(reference_path[i], reference_ptr[i]);
  }


  // saving the indexes
  int mi = 0;
  int qrinds[match_num][2];
  for (i = 0; i < match_num_query; i++) {
    for (j = 0; j < match_num_reference; j++) {

      if (i == j) continue; // *** remove this line if query list and ref list are different ***  

      qrinds[mi][0] = i; // indexes
      qrinds[mi][1] = j; // used for matching
      ++mi;
    }
  }
  assert(mi ==  match_num);

  //TODO: This loop can be made parallel
  for (i = 0; i < match_num; i++) {
    
    fprintf(stdout, "Match Iter:\t%d\n", i);
    
    motif_search_in_reference(query_ptr[qrinds[i][0]], reference_ptr[qrinds[i][1]], ctrl_file, spotted_potential_motif_regions_all_pptr[i]);
  }
  

  // freeing the segments
  for (i =0; i< match_num_query; i++) {
    empty_SegmentPTR(query_ptr[i]);
    destroy_SegmentPTR(query_ptr[i]);
  }
  printf("Loading the queries...\n");
  for (i =0; i < match_num_reference; i++) {
    empty_SegmentPTR(reference_ptr[i]);
    destroy_SegmentPTR(reference_ptr[i]);
  }

  printf("\n!!FInished!!!!\n");
  printf("\nSaving all the spots...\n");

  //TODO: No need to print here you can call it later.

  char buff[100];
  strcpy(buff, out_filename);
  FILEPTR fout_regions = fopen(strcat(buff, ".pmrmod"), "wb");

  printf_SpotVectorPPTR(fout_regions, spotted_potential_motif_regions_all_pptr, match_num);

  printf("\n!!!Saving Done!!!\n");

  fclose(fout_regions);
  //TODO: this freeing gives some error. Look into it.

  // freeing and destroyign spot vectors
  /*for(i = 0; i < match_num; i++){
    free_members_SpotVectorPTR(spotted_potential_motif_regions_all_pptr[i]);
    destroy_SpotVectorPTR(spotted_potential_motif_regions_all_pptr[i]);

    free_members_SpotVectorPTR(spotted_motifs_all_pptr[i]);
    destroy_SpotVectorPTR(spotted_motifs_all_pptr[i]);
  }
  */


  return 0;
}




void motif_search_in_reference(const SegmentPTR query_ptr,const SegmentPTR reference_ptr, const char *  ctrl_file, SpotVectorPTR spotted_potential_motif_regions_ptr) {
 
  /*******Approach*********
   ******First-Pass********
   * - for each voiced part
   **** - ignore if that voice part is too small compared to the
   ****   queried motif (in terms of saddle_points).
   **** - if it is little small then padd it with zeros
   ****   till it becomes equal or greater than the motif.
   **** - Do windowed RLCS using saddle points for that part and save the results.
   * - End For
   ***********************
   *****Second-Pass*******
   * - for each group from the first pass
   **** - Do windowed RLCS using entire pitch for that part and save the results.
   ***********************/
  
  
  int vi;   // voice part iterator
  int qvi; // voice part iterator of query
  int ri;  // potential motif region iterator
  
  printf("Setting RLCS parameters...\n");
  ParametersRLCSPTR parameters_rlcs_ptr = create_ParametersRLCSPTR();
  set_ParametersRLCSPTR(ctrl_file, parameters_rlcs_ptr); //setting rlcs parameters from the control file

  printf_ParametersRLCSPTR(stdout, parameters_rlcs_ptr);
  printf("RLCS parameters setting finished...\n");

  // first pass to find potential motif regions

  // created just to pass during the first pass of rlcs.
  /* DoubleMatrixPTR dummy_shift_points_ptr = create_DoubleMatrixPTR(); */
  /* dummy_shift_points_ptr->size = 0; */
  /* dummy_shift_points_ptr->dim = 0; */



  /* for (vi = 0; vi < reference_ptr->roi_ptr->size;  vi++) { */
  /*   printf("Voiced part: %d\n", vi); */

  /*   // Find the potential motif reagion in each voiced part */
  /*   int start = arg_find_smallest_num_after_D(reference_ptr->saddle_points_ptr->x_ptr, reference_ptr->roi_ptr->vec[vi].start); */
  /*   int end = arg_find_largest_num_before_D(reference_ptr->saddle_points_ptr->x_ptr, reference_ptr->roi_ptr->vec[vi].end); */

  /*   if (vi+1 < reference_ptr->roi_ptr->size && start >= reference_ptr->roi_ptr->vec[vi+1].start) */
  /*     continue; */

  /*   RLCS_window(query_ptr->saddle_points_ptr->y_ptr, reference_ptr->saddle_points_ptr->y_ptr, parameters_rlcs_ptr, 0, query_ptr->saddle_points_ptr->y_ptr->size-1, start, end, dummy_shift_points_ptr, parameters_rlcs_ptr->hop_size_fp, reference_ptr->path, query_ptr->path, spotted_potential_motif_regions_ptr);  */
  /* } */

  /* destroy_DoubleMatrixPTR(dummy_shift_points_ptr); */



  for (qvi =0; qvi < query_ptr->roi_ptr->size; qvi++) {
    int start_q = query_ptr->roi_ptr->vec[qvi].start;
    int end_q = query_ptr->roi_ptr->vec[qvi].end;
    
    // rejecting every query  i.e. less than 1.2 sec long
    if (((end_q*196 - start_q*196)/44100.0f) < 1.2){
      printf("******************* rejecting query ***************\n");

      continue;  // *** very much tuned for prticular needs. change this for general perpose
    }
    
    for (vi = 0; vi < reference_ptr->roi_ptr->size; vi++) {
      printf("Query Voiced part: %d\tVoiced part: %d\n", qvi, vi);

      int start = reference_ptr->roi_ptr->vec[vi].start;
      int end =  reference_ptr->roi_ptr->vec[vi].end;
         
      RLCS_window(query_ptr->feat_mat_ptr, reference_ptr->feat_mat_ptr, parameters_rlcs_ptr, start_q, end_q, start, end, reference_ptr->saddle_points_ptr->x_ptr, parameters_rlcs_ptr->hop_size_sp, reference_ptr->path, query_ptr->path, spotted_potential_motif_regions_ptr, 0, query_ptr->saddle_points_ptr, reference_ptr->saddle_points_ptr); 
    }
  }


  /* // Correcting start and ends */
  /* for (ri = 0; ri < spotted_potential_motif_regions_ptr->size; ri++) { */
  /*   printf("Region num: %d, best member score: %f\n", ri, spotted_potential_motif_regions_ptr->vec[ri].best_member_ptr->score); */

  /*   // Changing start and in the entire spot to correspond with the original */
  /*   // input file instaed of the saddle file. */
    
   
  /*   // group start and group end */
  /*   spotted_potential_motif_regions_ptr->vec[ri].group_start = reference_ptr->saddle_points_ptr->x_ptr->mat[spotted_potential_motif_regions_ptr->vec[ri].group_start][0]; */
  /*   spotted_potential_motif_regions_ptr->vec[ri].group_end = reference_ptr->saddle_points_ptr->x_ptr->mat[spotted_potential_motif_regions_ptr->vec[ri].group_end][0]; */

  /*   // best member start and end */
  /*   /\* spotted_potential_motif_regions_ptr->vec[ri].best_member_ptr->start = reference_ptr->saddle_points_ptr->x_ptr->mat[spotted_potential_motif_regions_ptr->vec[ri].best_member_ptr->start][0]; *\/ */
  /*   /\* spotted_potential_motif_regions_ptr->vec[ri].best_member_ptr->end = reference_ptr->saddle_points_ptr->x_ptr->mat[spotted_potential_motif_regions_ptr->vec[ri].best_member_ptr->end][0]; *\/ */
    
  /*   // all members start and end */
  /*   int mi; */
  /*   for (mi = 0; mi < spotted_potential_motif_regions_ptr->vec[ri].member_vec_ptr->size; mi++) { */
  /*     spotted_potential_motif_regions_ptr->vec[ri].member_vec_ptr->vec[mi].start = reference_ptr->saddle_points_ptr->x_ptr->mat[spotted_potential_motif_regions_ptr->vec[ri].member_vec_ptr->vec[mi].start][0]; */
  /*     spotted_potential_motif_regions_ptr->vec[ri].member_vec_ptr->vec[mi].end = reference_ptr->saddle_points_ptr->x_ptr->mat[spotted_potential_motif_regions_ptr->vec[ri].member_vec_ptr->vec[mi].end][0]; */
  /*   }     */
  /* } */

  destroy_ParametersRLCSPTR(parameters_rlcs_ptr);
}
