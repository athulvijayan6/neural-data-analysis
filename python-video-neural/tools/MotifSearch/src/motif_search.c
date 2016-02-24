
/*
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

void motif_search_in_reference(const SegmentPTR motif_ptr,const SegmentPTR reference_ptr, const String ctrl_file, SpotVectorPTR spotted_potential_motif_regions_ptr, SpotVectorPTR spotted_motifs_ptr);

int main(int argc, char* argv[]) {
  
  int i; // simple iterators.

  if (argc != 3+1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: motif_search <ctrl_file> <query_reference_list> <out_filename>\n");
    exit(-1);
  }

  String ctrl_file = argv[1];
  String query_reference_list = argv[2];
  String out_filename = argv[3];


  FILEPTR fin = fopen(query_reference_list, "r");
  assert(fin != NULL);

  int match_num;
  
  fscanf(fin, "%d", &match_num);

  char query_path[1000][1000];
  char reference_path[1000][1000];
  
  for (i = 0; i < match_num; i++){
    fscanf(fin, "%s%s", query_path[i], reference_path[i]);
    printf( "%s\t%s\n", query_path[i], reference_path[i]);
  }
  fclose(fin);

  
  SpotVectorPPTR spotted_potential_motif_regions_all_pptr =  (SpotVectorPPTR) malloc(sizeof(SpotVectorPTR)*match_num);
  SpotVectorPPTR spotted_motifs_all_pptr = (SpotVectorPPTR) malloc(sizeof(SpotVectorPTR)*match_num);


  
  // creating and allocating spot vectors
  for(i = 0; i < match_num; i++){
    spotted_potential_motif_regions_all_pptr[i] = create_SpotVectorPTR();
    // simply allocating 10. It will be rellocated inside if necessary.
    alloc_members_SpotVectorPTR(spotted_potential_motif_regions_all_pptr[i], 10);

    spotted_motifs_all_pptr[i] = create_SpotVectorPTR();
    // simply allocating 10. It will be rellocated inside if necessary.
    alloc_members_SpotVectorPTR(spotted_motifs_all_pptr[i], 10);
  }
  
  
  //TODO: This loop can be made parallel
  for (i = 0; i < match_num; i++) {
    
    fprintf(stdout, "Match Iter:\t%d\n", i);

    SegmentPTR query_ptr = create_SegmentPTR();
    SegmentPTR reference_ptr = create_SegmentPTR();

    fill_SegmentPTR_from_file(query_path[i], query_ptr);
    fill_SegmentPTR_from_file(reference_path[i], reference_ptr);

    //spotted_potential_motif_regions_ptr[i] = create_SpotVectorPTR();
    //spotted_motifs_ptr[i] = create_SpotVectorPTR();

    motif_search_in_reference(query_ptr, reference_ptr, ctrl_file, spotted_potential_motif_regions_all_pptr[i], spotted_motifs_all_pptr[i]);

    //  TODO: check whether empty_Segment is working here or not.
    destroy_SegmentPTR(query_ptr);
    destroy_SegmentPTR(reference_ptr);
  }

  printf("\n!!FInished!!!!\n");
  printf("\nSaving all the spots...\n");

  //TODO: No need to print here you can call it later.

  char buff[100];
  strcpy(buff, out_filename);
  FILEPTR fout_regions = fopen(strcat(buff, ".pmr"), "wb");
  strcpy(buff, out_filename);
  FILEPTR fout_motifs = fopen(strcat(buff, ".spm"), "wb");
 


  printf_SpotVectorPPTR(fout_regions, spotted_potential_motif_regions_all_pptr, match_num);
  printf_SpotVectorPPTR(fout_motifs, spotted_motifs_all_pptr, match_num);



  printf("\n!!!Saving Done!!!\n");

  fclose(fout_regions);
  fclose(fout_motifs);



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




void motif_search_in_reference(const SegmentPTR query_ptr,const SegmentPTR reference_ptr, const String ctrl_file, SpotVectorPTR spotted_potential_motif_regions_ptr, SpotVectorPTR spotted_motifs_ptr) {
 
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
  int ri;  // potential motif region iterator
  
  
  printf("Setting RLCS parameters...\n");
  ParametersRLCSPTR parameters_rlcs_ptr = create_ParametersRLCSPTR();
  set_ParametersRLCSPTR(ctrl_file, parameters_rlcs_ptr); //setting rlcs parameters from the control file

  printf_ParametersRLCSPTR(stdout, parameters_rlcs_ptr);
  printf("RLCS parameters setting finished...\n");

  // first pass to find potential motif regions

  // created just to pass during the first pass of rlcs.
  DoubleMatrixPTR dummy_shift_points_ptr = create_DoubleMatrixPTR();
  dummy_shift_points_ptr->size = 0;
  dummy_shift_points_ptr->dim = 0;
  for (vi = 0; vi < reference_ptr->roi_ptr->size;  vi++) {
    printf("Voiced part: %d\n", vi);

    // Find the potential motif reagion in each voiced part
    int start = arg_find_smallest_num_after_D(reference_ptr->saddle_points_ptr->x_ptr, reference_ptr->roi_ptr->vec[vi].start);
    int end = arg_find_largest_num_before_D(reference_ptr->saddle_points_ptr->x_ptr, reference_ptr->roi_ptr->vec[vi].end);

    if (vi+1 < reference_ptr->roi_ptr->size && start >= reference_ptr->roi_ptr->vec[vi+1].start)
      continue;
    
    RLCS_window(query_ptr->saddle_points_ptr->y_ptr, reference_ptr->saddle_points_ptr->y_ptr, parameters_rlcs_ptr, 0, query_ptr->saddle_points_ptr->y_ptr->size-1, start, end, dummy_shift_points_ptr, parameters_rlcs_ptr->hop_size_fp, reference_ptr->path, query_ptr->path, spotted_potential_motif_regions_ptr); 
  }

  destroy_DoubleMatrixPTR(dummy_shift_points_ptr);

  // second pass to find motifs in the potential motif regions.
  for (ri = 0; ri < spotted_potential_motif_regions_ptr->size; ri++) {
    printf("Region num: %d, best member score: %f\n", ri, spotted_potential_motif_regions_ptr->vec[ri].best_member_ptr->score);

    // ignoring groups with bad match
    if (spotted_potential_motif_regions_ptr->vec[ri].best_member_ptr->score < 0.50) {
      printf("score very less\n");
      continue;
    }
    // Find the motifs in each potential motif region
    int start = reference_ptr->saddle_points_ptr->x_ptr->mat[spotted_potential_motif_regions_ptr->vec[ri].group_start][0];
    int end = reference_ptr->saddle_points_ptr->x_ptr->mat[spotted_potential_motif_regions_ptr->vec[ri].group_end][0];

    RLCS_window(query_ptr->feat_mat_ptr, reference_ptr->feat_mat_ptr, parameters_rlcs_ptr, 0, query_ptr->feat_mat_ptr->size-1, start, end, reference_ptr->saddle_points_ptr->x_ptr, parameters_rlcs_ptr->hop_size_sp, reference_ptr->path, query_ptr->path, spotted_motifs_ptr);
  }
  destroy_ParametersRLCSPTR(parameters_rlcs_ptr);
}
