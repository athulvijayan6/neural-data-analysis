
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
#include "yamlwrapper.h"


#define MAX_COHORTS 10
#define MAX_RECS 100

// TODO: run it with sample inputs and see.



void motif_search_in_reference(const SegmentPTR query_ptr,const SegmentPTR reference_ptr, const char * ctrl_file, MemberVectorPTR member_vec_ptr);

void motif_search_in_reference_2(const SegmentPTR query_ptr,const SegmentPTR reference_ptr, const char * ctrl_file, MemberVectorPTR member_vec_ptr);

int main(int argc, char* argv[]) {



  
  int i,j; // simple iterators.
 
   
  if (argc != 5+1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: %s <ctrl_file> <ragaId2CohortsId.yaml> <refId2rec_train> <ragaId_list>  <out_filename>\n", argv[0]); 
    exit(-1);
  }

  // getting input arguments 
  char * ctrl_file = argv[1];
  char * ragaId2CohortsId_file = argv[2]; // yaml file
  char * ragaId2Rec_train_file = argv[3]; // yaml file
  //  char * ragaId2Rec_test_file = argv[4];  // yaml file
  char * ragaId_list_file = argv[4];        // ragas to be tested
                                          // their data info should be there in train 
                                          // and test files
  char * out_filename = argv[5];          



  /*
  // getting input arguments 
  char * ctrl_file = "ctrl_files/ctrl_file_slcs";
  char * ragaId2CohortsId_file =  "RagaId2CohortsId.yaml";
  char * ragaId2Rec_train_file =  "RagaId2RecPaths_train.yaml";
//  char * ragaId2Rec_test_file =  "RagaId2RecPaths_test.yaml";
  char * ragaId_list_file =  "ragaId_list_shrey" ;
  char * out_filename = "hist_scores_out_slcs.yaml";          
  */
    

  

  // no need to read the ctrl_file
  // reading all the yaml files
  FILE * ragaId2CohortsId_fp = fopen(ragaId2CohortsId_file, "r");
  assert(ragaId2CohortsId_fp != NULL);

  FILE * ragaId2Rec_train_fp = fopen(ragaId2Rec_train_file, "r");
  assert(ragaId2Rec_train_fp != NULL);

  //  FILE * ragaId2Rec_test_fp = fopen(ragaId2Rec_test_file, "r");
  //assert(ragaId2Rec_test_fp != NULL);

  // reading raga list file
  FILE * ragaId_list_fp = fopen(ragaId_list_file, "r");
  assert(ragaId_list_fp != NULL);


  // getting yamlfiles to a data structure
  IntVector * ragaId2CohortsId = yaml_read_file_type1(ragaId2CohortsId_fp); 
  StringVector * ragaId2Rec_train = yaml_read_file_type2(ragaId2Rec_train_fp);
  //StringVector * ragaId2Rec_test = yaml_read_file_type2(ragaId2Rec_test_fp);
  
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

  // allocating space
  DoubleVectorPTR *raga_scores = (DoubleVectorPTR*) malloc(ragaIds_ptr->size*sizeof(DoubleVectorPTR));
  
  // 3rd dim for raga.. cohorts 

  
  int raga_i;
  for (raga_i = 0; raga_i < ragaIds_ptr -> size; raga_i++) {
    int ragaId = ragaIds_ptr -> vec[raga_i];

    // for every raga get the cohorts.
    // get the train and test recordings for raga and cohorts
    assert(ragaId2CohortsId[ragaId].size != 0); // cohorts are there
    assert(ragaId2Rec_train[ragaId].size != 0); // train recs are there
    //   assert(ragaId2Rec_test[ragaId].size != 0);  // test recs are there 

    // create allocate raga_scores[] some space
    raga_scores[raga_i] = create_DoubleVectorPTR();
    alloc_members_DoubleVectorPTR(raga_scores[raga_i], 100);

      

    // compare train of a raga with all train of its cohorts
    int q_train_reci; 
    for (q_train_reci = 0; q_train_reci < ragaId2Rec_train[ragaId].size; q_train_reci++) { 
      {

	SegmentPTR q_train_segment_ptr = create_SegmentPTR();
	fill_SegmentPTR_from_file(ragaId2Rec_train[ragaId].vec[q_train_reci].data, q_train_segment_ptr);


	// test this segment file against all trian files of all cohorts ragas
	int cohort_i;
	for (cohort_i = 0; cohort_i < ragaId2CohortsId[ragaId].size; cohort_i++) {
	  //  cohort_i=1; // remove this
	  int cohortId = ragaId2CohortsId[ragaId].vec[cohort_i]; 
	  int r_train_reci;
          // make this parallel
          #pragma omp parallel for
	  for (r_train_reci = 0; r_train_reci < ragaId2Rec_train[cohortId].size; r_train_reci++) {

	    SegmentPTR r_train_segment_ptr = create_SegmentPTR();
	    fill_SegmentPTR_from_file(ragaId2Rec_train[cohortId].vec[r_train_reci].data, r_train_segment_ptr);
	    
	    MemberVectorPTR member_vec_ptr =  create_MemberVectorPTR();
	    // simply allocating 10. It will be rellocated inside if necessary.
	    alloc_members_MemberVectorPTR(member_vec_ptr, 100);   

	    motif_search_in_reference_2(q_train_segment_ptr, r_train_segment_ptr, ctrl_file, member_vec_ptr);
	    
	    
	    // if spot vec is empty continue
	    if (member_vec_ptr -> size == 0) {
	      free_members_MemberVectorPTR(member_vec_ptr);
	      destroy_MemberVectorPTR(member_vec_ptr);	      
	      continue;
	    }
	    
	    // save against raga scores
	    // get all the number of members to get number of scores
	    int scores_num = member_vec_ptr->size;
	    
	    assert(scores_num != 0); // scores_num cannot be zero
	    
            int cid = cohort_i + 1;

	    if (raga_scores[raga_i] -> size + scores_num > raga_scores[raga_i] -> max_size)
	      realloc_members_DoubleVectorPTR(raga_scores[raga_i], raga_scores[raga_i] -> max_size * 2 + scores_num);
	    
	    
	 
	    // finally save/append the scores
	    {
	      int j;
	      for (j = 0; j < member_vec_ptr -> size; j++) {
		if (member_vec_ptr -> vec[j].alignment_vec_ptr->size != 0) {
		  raga_scores[raga_i] -> vec[raga_scores[raga_i]->size] = member_vec_ptr -> vec[j].score;
		  raga_scores[raga_i] -> size += 1;
		}
	      }  
	    }
	 


   
	    // freeing ...
	    free_members_MemberVectorPTR(member_vec_ptr);
	    destroy_MemberVectorPTR(member_vec_ptr);
	    
	    empty_SegmentPTR(r_train_segment_ptr);
	    destroy_SegmentPTR(r_train_segment_ptr);
            //break;
	  }
	  //  break;
	}

	// freeing ...
	empty_SegmentPTR(q_train_segment_ptr);
	destroy_SegmentPTR(q_train_segment_ptr);
	// break;
      }
      
      
    }
    //break; 
  }
  
  
  // print raga_scores in a yaml file as described above
  FILE * out_fp = fopen(out_filename, "w");
  
  {
    int i;
    for (i = 0; i < ragaIds_ptr -> size; i++) {
      int key = ragaIds_ptr -> vec[i];

      fprintf(out_fp, "%d:\n", key);
      fprintf(out_fp, "- ");
      int l;
      for (l = 0; l <  raga_scores[i] -> size; l++)
	fprintf(out_fp, "%lf ", raga_scores[i] -> vec[l]);
      fprintf(out_fp, "\n");
    }
  }


  // freeing ...
  fclose(out_fp);

  {
    int i;
    for (i = 0; i < ragaIds_ptr -> size; i++) {
      free_members_DoubleVectorPTR(raga_scores[i]);
      destroy_DoubleVectorPTR(raga_scores[i]);
    }
  }

  free(raga_scores);
  
  free_members_IntVectorPTR(ragaIds_ptr);
  destroy_IntVectorPTR(ragaIds_ptr);
 
  yaml_delete_type1_ds(ragaId2CohortsId);
  yaml_delete_type2_ds(ragaId2Rec_train);
  //  yaml_delete_type2_ds(ragaId2Rec_test);

  fclose(ragaId2CohortsId_fp);
  fclose(ragaId2Rec_train_fp);
  //fclose(ragaId2Rec_test_fp);
  fclose(ragaId_list_fp);

  return 0;
}




void motif_search_in_reference(const SegmentPTR query_ptr,const SegmentPTR reference_ptr, const char * ctrl_file, MemberVectorPTR member_vec_ptr) {
 
  int vi;   // voice part iterator
  int qvi; // voice part iterator of query
  int ri;  // potential motif region iterator
  
  printf("Setting SLCS parameters...\n");
  ParametersRLCSPTR parameters_slcs_ptr = create_ParametersRLCSPTR();
  set_ParametersRLCSPTR(ctrl_file, parameters_slcs_ptr); //setting rlcs parameters from the control file

  printf_ParametersRLCSPTR(stdout, parameters_slcs_ptr);
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
         


      //  checking for space in member_vec_ptr
      if (member_vec_ptr->size == member_vec_ptr->max_size)
	realloc_members_MemberVectorPTR(member_vec_ptr, 2*member_vec_ptr->max_size);
      int si = member_vec_ptr->size;
      SLCS_perform(query_ptr->feat_mat_ptr, reference_ptr->feat_mat_ptr, parameters_slcs_ptr, start_q, end_q, start, end, &(member_vec_ptr -> vec[si]));
      member_vec_ptr -> size += 1;
      //break;
    }
    //  break;
  }
  destroy_ParametersRLCSPTR(parameters_slcs_ptr);
}



void motif_search_in_reference_2(const SegmentPTR query_ptr,const SegmentPTR reference_ptr, const char * ctrl_file, MemberVectorPTR member_vec_ptr) {
 
  
  printf("Setting SLCS parameters...\n");
  ParametersRLCSPTR parameters_slcs_ptr = create_ParametersRLCSPTR();
  set_ParametersRLCSPTR(ctrl_file, parameters_slcs_ptr); //setting rlcs parameters from the control file

  printf_ParametersRLCSPTR(stdout, parameters_slcs_ptr);
  printf("RLCS parameters setting finished...\n");
    
  // rejecting every query  i.e. less than 1.2 sec long
  if (((query_ptr->feat_mat_ptr->size*196)/44100.0f) < 1.2){
    printf("******************* rejecting query ***************\n");

    return;  // *** very much tuned for prticular needs. change this for general perpose
  }
    


  //  checking for space in member_vec_ptr
  if (member_vec_ptr->size == member_vec_ptr->max_size)
    realloc_members_MemberVectorPTR(member_vec_ptr, 2*member_vec_ptr->max_size);

  int si = member_vec_ptr->size;
  //printf("hi1");
  SLCS_hard_perform(query_ptr->feat_mat_ptr, reference_ptr->feat_mat_ptr, parameters_slcs_ptr, 0, query_ptr->feat_mat_ptr->size-1, 0, reference_ptr->feat_mat_ptr->size-1, &(member_vec_ptr -> vec[si]));
  //  printf("hi2");
  member_vec_ptr -> size += 1;


  destroy_ParametersRLCSPTR(parameters_slcs_ptr);
}
