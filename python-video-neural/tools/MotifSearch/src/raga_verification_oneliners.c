
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
#include "yamlwrapper.h"


#define MAX_COHORTS 50
#define MAX_RECS 500

// TODO: run it with sample inputs and see.


void motif_search_in_reference(const SegmentPTR motif_ptr,const SegmentPTR reference_ptr, const char * ctrl_file, SpotVectorPTR spotted_potential_motif_regions_ptr);

int main(int argc, char* argv[]) {
  
  int i,j; // simple iterators.
  
  
  
  if (argc != 6+1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: %s <ctrl_file> <ragaId2CohortsId.yaml> <refId2rec_train> <refId2rec_test> <ragaId_list>  <out_filename>\n", argv[0]); 
    exit(-1);
  }
  
  
  // getting input arguments 
  char * ctrl_file = argv[1];
  char * ragaId2CohortsId_file = argv[2]; // yaml file
  char * ragaId2Rec_train_file = argv[3]; // yaml file
  char * ragaId2Rec_test_file = argv[4];  // yaml file
  char * ragaId_list_file = argv[5];        // ragas to be tested
                                          // their data info should be there in train 
                                          // and test files
  char * out_filename = argv[6];          
  
  /*
  // for debugging...
  char ctrl_file[100] = "ctrl_files/ctrl_file_across_oneline"; 
  char ragaId2CohortsId_file[100] = "RagaId2CohortsId.yaml";
  char ragaId2Rec_train_file[100] =  "RagaId2RecPaths_train.yaml";
  char ragaId2Rec_test_file[100] = "RagaId2RecPaths_test.yaml";
  char ragaId_list_file[100] = "ragaId_list_shrey";
  char out_filename[100] = "hist_scores_out.yaml";
  */
  

  // no need to read the ctrl_file
  // reading all the yaml files
  FILE * ragaId2CohortsId_fp = fopen(ragaId2CohortsId_file, "r");
  assert(ragaId2CohortsId_fp != NULL);

  FILE * ragaId2Rec_train_fp = fopen(ragaId2Rec_train_file, "r");
  assert(ragaId2Rec_train_fp != NULL);

  FILE * ragaId2Rec_test_fp = fopen(ragaId2Rec_test_file, "r");
  assert(ragaId2Rec_test_fp != NULL);

  // reading raga list file
  FILE * ragaId_list_fp = fopen(ragaId_list_file, "r");
  assert(ragaId_list_fp != NULL);


  // getting yamlfiles to a data structure
  IntVector * ragaId2CohortsId = yaml_read_file_type1(ragaId2CohortsId_fp); 
  StringVector * ragaId2Rec_train = yaml_read_file_type2(ragaId2Rec_train_fp);
  StringVector * ragaId2Rec_test = yaml_read_file_type2(ragaId2Rec_test_fp);
  
  // getting raga_list to array ds
  IntVectorPTR ragaIds_ptr = create_IntVectorPTR();
  alloc_members_IntVectorPTR(ragaIds_ptr, 500);
  {
    char * buff =  (char *) malloc(500*sizeof(char));    
    while(fgets(buff, 500, ragaId_list_fp) != NULL) {
      assert(ragaIds_ptr -> size < 500);
      ragaIds_ptr -> vec[ragaIds_ptr -> size] = atoi(buff);
      ragaIds_ptr -> size += 1;
    }
    free(buff);
  }

  

  /* Writing mind:
    - for every raga in ragaId
    - get the recordings from the test file one by one
    - test each recording against train of that raga and all its cohorts
      (test is to be defined)
    - define test
       - first do oneliners match and get all the scores of all the members
       - do histogram of all the scores
       - first verify manually scores histogram of raga against its cohorts.
       - if that makes sense then come up with a measure.               
    - output should be yaml file in this format:
       (ragaId cannot be more than 500. So we exploit this and use modulo thing 
       to encode for diff recording)
       -$ <ragaId 1 + (500 * 0)>: 
       -$ - <-1 test rec path>
       -$ - <0 followed by scores with raga itself as a string with spaces>
       -$ - <cohortId 1 followed by scores with cohort 1 ... >
       -$ - <cohortId 2 followed by scores with cohort 2 ... >
       -$ <ragaId 1 + (500 * 1)>:
       -$ - <-1 test rec path>
       -$ - <0 followed by scores with raga itself as a string with spaces>
       -$ - <cohortId 1 followed by scores with cohort 1 ... >
       -$ - <cohortId 2 followed by scores with cohort 2 ... >
       -$ <ragaId 1 + (500 * k1)>:
       -$ - <-1 test rec path>
       -$ - <0 followed by scores with raga itself as a string with spaces>
       -$ - <cohortId 1 followed by scores with cohort 1 ... >
       -$ - <cohortId 2 followed by scores with cohort 2 ... >
       -$ <ragaId N + (500 * 0)>:
       -$ - <-1 test rec path>
       -$ - <0 scores with raga itself as a string with spaces>
       -$ - <cohortId 1 followed by scores with cohort 1 ... >
       -$ - <cohortId M followed by scores with cohort M ... >
       -$ <ragaId N + (500 * kN)>:
       -$ - <-1 test rec path>
       -$ - <0 scores with raga itself as a string with spaces>
       -$ - <cohortId 1 followed by scores with cohort 1 ... >
       -$ - <cohortId M followed by scores with cohort M ... >
   */

  /* Things to do:
    - for every raga get the cohorts.
    - get the train and test recordings for raga and cohorts
    - do motif search of a raga test rec with raga itself. Save all the scores.
    - then motif search of raga test rec with all the cohorts. Save all the scores.
  */



  DoubleVectorPTR raga_scores[ragaIds_ptr -> size][MAX_RECS][MAX_COHORTS]; // 2nd dim for recordings
                                                              // 3rd dim for raga.. cohorts 
 
  int raga_i;
  for (raga_i = 0; raga_i < ragaIds_ptr -> size; raga_i++) {
    int ragaId = ragaIds_ptr -> vec[raga_i];
    
    // for every raga get the cohorts.
    // get the train and test recordings for raga and cohorts
    assert(ragaId2CohortsId[ragaId].size != 0); // cohorts are there
    assert(ragaId2Rec_train[ragaId].size != 0); // train recs are there
    assert(ragaId2Rec_test[ragaId].size != 0);  // test recs are there 

    // create allocate raga_scores[][] some space
    {
      int i;
      for (i = 0; i < MAX_RECS; i++) {
        int j;
	for (j = 0; j < MAX_COHORTS; j++) {
	  raga_scores[raga_i][i][j] = create_DoubleVectorPTR();
	  alloc_members_DoubleVectorPTR(raga_scores[raga_i][i][j], 100);
	}
      }
    }

    
    
    // do motif search of a raga test rec with raga itself. Save all the scores.
    // then motif search of raga test rec with all the cohorts. Save all the scores.
    int test_reci; 
    for (test_reci = 0; test_reci < ragaId2Rec_test[ragaId].size; test_reci++) {
      //   test_reci = 1; //remove this
      SegmentPTR test_segment_ptr = create_SegmentPTR();
      fill_SegmentPTR_from_file(ragaId2Rec_test[ragaId].vec[test_reci].data, test_segment_ptr);    
            
      {
	// test this segment file against all trian files of raga itself
	int train_reci;
	for (train_reci = 0; train_reci < ragaId2Rec_train[ragaId].size; train_reci++) {
	  //	    train_reci = 11;  // remove this
	  SegmentPTR train_raga_segment_ptr = create_SegmentPTR();
	  fill_SegmentPTR_from_file(ragaId2Rec_train[ragaId].vec[train_reci].data, train_raga_segment_ptr);

	  SpotVectorPTR spot_vec_ptr =  create_SpotVectorPTR();
	  // simply allocating 10. It will be rellocated inside if necessary.
	  alloc_members_SpotVectorPTR(spot_vec_ptr, 100);   


	  motif_search_in_reference(test_segment_ptr, train_raga_segment_ptr, ctrl_file, spot_vec_ptr);


	  // if spot vec is empty continue
	  if (spot_vec_ptr -> size == 0) {
	    free_members_SpotVectorPTR(spot_vec_ptr);
	    destroy_SpotVectorPTR(spot_vec_ptr);	      
	    continue;
	  }

	  // save against raga scores
	  // get all the number of members to get number of scores
	  int scores_num = 0;
	  {
	    int i;
	    for (i = 0; i < spot_vec_ptr -> size; i++)
	      scores_num += spot_vec_ptr -> vec[i].member_vec_ptr -> size;
	  }

	  assert(scores_num != 0); // scores_num cannot be zero

	  int cid = 0;
	  if (raga_scores[raga_i][test_reci][cid] -> size + scores_num > raga_scores[raga_i][test_reci][cid] -> max_size)
	    realloc_members_DoubleVectorPTR(raga_scores[raga_i][test_reci][cid], raga_scores[raga_i][test_reci][cid] -> max_size * 2 + scores_num);
	    

	  // finally save/append the scores
	  {
	    int i;
	    int s = raga_scores[raga_i][test_reci][cid] -> size;
	    for (i = 0; i < spot_vec_ptr -> size; i++) {
	      int j;
	      for (j = 0; j < spot_vec_ptr -> vec[i].member_vec_ptr -> size; j++) {
		raga_scores[raga_i][test_reci][cid] -> vec[s] = spot_vec_ptr -> vec[i].member_vec_ptr -> vec[j].score;
		s += 1;
	      }
	      raga_scores[raga_i][test_reci][cid] -> size = s;
	    }
	  }

	  // freeing ...
	  free_members_SpotVectorPTR(spot_vec_ptr);
	  destroy_SpotVectorPTR(spot_vec_ptr);

	  empty_SegmentPTR(train_raga_segment_ptr);
	  destroy_SegmentPTR(train_raga_segment_ptr);
	}
      }
          
      
      {
	// test this segment file against all trian files of all cohorts ragas
	int cohort_i;
	for (cohort_i = 0; cohort_i < ragaId2CohortsId[ragaId].size; cohort_i++) {
	  //  cohort_i=1; // remove this
	  int cohortId = ragaId2CohortsId[ragaId].vec[cohort_i]; 
	  int train_reci;
	  for (train_reci = 0; train_reci < ragaId2Rec_train[cohortId].size; train_reci++) {
            //train_reci =8; // remove this
	    SegmentPTR train_raga_segment_ptr = create_SegmentPTR();
	    fill_SegmentPTR_from_file(ragaId2Rec_train[cohortId].vec[train_reci].data, train_raga_segment_ptr);
	    
	    SpotVectorPTR spot_vec_ptr =  create_SpotVectorPTR();
	    // simply allocating 10. It will be rellocated inside if necessary.
	    alloc_members_SpotVectorPTR(spot_vec_ptr, 100);   
	    

	    motif_search_in_reference(test_segment_ptr, train_raga_segment_ptr, ctrl_file, spot_vec_ptr);


	    // if spot vec is empty continue
	    if (spot_vec_ptr -> size == 0) {
	      free_members_SpotVectorPTR(spot_vec_ptr);
	      destroy_SpotVectorPTR(spot_vec_ptr);	      
	      continue;
	    }

	    // save against raga scores
	    // get all the number of members to get number of scores
	    int scores_num = 0;
	    {
	      int i;
	      for (i = 0; i < spot_vec_ptr -> size; i++)
		scores_num += spot_vec_ptr -> vec[i].member_vec_ptr -> size;
	    }

	    assert(scores_num != 0); // scores_num cannot be zero

	    // reallocate more space if less
            int cid = cohort_i + 1;
	    if (raga_scores[raga_i][test_reci][cid] -> size + scores_num > raga_scores[raga_i][test_reci][cid] -> max_size)
	      realloc_members_DoubleVectorPTR(raga_scores[raga_i][test_reci][cid], raga_scores[raga_i][test_reci][cid] -> max_size * 2 + scores_num);
	    

	    // finally save/append the scores
	    {
	      int i;
	      int s = raga_scores[raga_i][test_reci][cid] -> size;
	      for (i = 0; i < spot_vec_ptr -> size; i++) {
		int j;
		for (j = 0; j < spot_vec_ptr -> vec[i].member_vec_ptr -> size; j++) {
		  raga_scores[raga_i][test_reci][cid] -> vec[s] = spot_vec_ptr -> vec[i].member_vec_ptr -> vec[j].score;
		  s += 1;
		}
		raga_scores[raga_i][test_reci][cid] -> size = s;
	      }
	    }

	    // freeing ...
	    free_members_SpotVectorPTR(spot_vec_ptr);
	    destroy_SpotVectorPTR(spot_vec_ptr);

	    empty_SegmentPTR(train_raga_segment_ptr);
	    destroy_SegmentPTR(train_raga_segment_ptr);
            //break;
	  }
	  //  break;
	}
      }
      

      // freeing ...
      empty_SegmentPTR(test_segment_ptr);
      destroy_SegmentPTR(test_segment_ptr);
      // break;
    }
    //break; 
  }

  // print raga_scores in a yaml file as described above
  FILE * out_fp = fopen(out_filename, "w");

  {
    int i;
    for (i = 0; i < ragaIds_ptr -> size; i++) {
      int j;
      for (j = 0; j < ragaId2Rec_test[ragaIds_ptr -> vec[i]].size; j++) {
	// check if this recording is valid
	//	if (raga_scores[i][j][0] -> size == 0)
	//  break;

	int key = ragaIds_ptr -> vec[i] + (500 * j);
        fprintf(out_fp, "%d:\n", key);
	fprintf(out_fp, "- -1 %s\n", ragaId2Rec_test[ragaIds_ptr ->vec[i]].vec[j].data);
	int k;
	for (k = 0; k < ragaId2CohortsId[ragaIds_ptr -> vec[i]].size + 1; k++) {
	  
	  if (k == 0) { 
	    fprintf(out_fp, "- 0 ");
	  } else {
	    fprintf(out_fp, "- %d ", ragaId2CohortsId[ragaIds_ptr->vec[i]].vec[k-1]);
	  }	  

	  int l;
	  for (l = 0; l <  raga_scores[i][j][k] -> size; l++)
	    fprintf(out_fp, "%lf ", raga_scores[i][j][k] -> vec[l]);
	  fprintf(out_fp, "\n");
	}
      }
    }
  }

  



  // freeing ...
  fclose(out_fp);

  {
    int i;
    for (i = 0; i < ragaIds_ptr -> size; i++) {
      int j;
      for (j = 0; j < MAX_RECS; j++) {
	int k;
	for (k = 0; k < MAX_COHORTS; k++) {
	  free_members_DoubleVectorPTR(raga_scores[i][j][k]);
	  destroy_DoubleVectorPTR(raga_scores[i][j][k]);
	}
      }
    }
  }
  

  free_members_IntVectorPTR(ragaIds_ptr);
  destroy_IntVectorPTR(ragaIds_ptr);
 
  yaml_delete_type1_ds(ragaId2CohortsId);
  yaml_delete_type2_ds(ragaId2Rec_train);
  yaml_delete_type2_ds(ragaId2Rec_test);

  fclose(ragaId2CohortsId_fp);
  fclose(ragaId2Rec_train_fp);
  fclose(ragaId2Rec_test_fp);
  fclose(ragaId_list_fp);

  return 0;
}




void motif_search_in_reference(const SegmentPTR query_ptr,const SegmentPTR reference_ptr, const char * ctrl_file, SpotVectorPTR spotted_potential_motif_regions_ptr) {
 
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

  // directly second pass to find potential motif regions

  for (qvi =0; qvi < query_ptr->roi_ptr->size; qvi++) {
    //qvi = 15;
    int start_q = query_ptr->roi_ptr->vec[qvi].start;
    int end_q = query_ptr->roi_ptr->vec[qvi].end;
    
    // rejecting every query  i.e. less than 1.2 sec long
    if (((end_q*196 - start_q*196)/44100.0f) < 1.2){
      printf("******************* rejecting query ***************\n");

      continue;  // *** very much tuned for prticular needs. change this for general perpose
    }
    
    for (vi = 0; vi < reference_ptr->roi_ptr->size; vi++) {

      /*
      //      vi = 10;
      if (qvi ==15 && vi == 10) {

	printf("I am in\n");
      }
      */

      printf("Query Voiced part: %d\tVoiced part: %d\n", qvi, vi);
      int start = reference_ptr->roi_ptr->vec[vi].start;
      int end =  reference_ptr->roi_ptr->vec[vi].end;
         
      RLCS_window(query_ptr->feat_mat_ptr, reference_ptr->feat_mat_ptr, parameters_rlcs_ptr, start_q, end_q, start, end, reference_ptr->saddle_points_ptr->x_ptr, parameters_rlcs_ptr->hop_size_sp, reference_ptr->path, query_ptr->path, spotted_potential_motif_regions_ptr, 1, query_ptr->saddle_points_ptr, reference_ptr->saddle_points_ptr); 
      //break;
    }
    //  break;
  }

  destroy_ParametersRLCSPTR(parameters_rlcs_ptr);
}
