
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

int main(int argc, char* argv[]) {



  
  int i,j; // simple iterators.


  
  if (argc != 6+1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: %s <ctrl_file> <ragaId2ioi_test.yaml> <refId2rec_train> <refId2rec_test> <ragaID_list> <out_filename>\n", argv[0]); 
    exit(-1);
  }

  char *  ctrl_file = argv[1];
  char *  ragaId2ioi_test_file = argv[2]; // yaml file
  char *  ragaId2Rec_train_file = argv[3]; // yaml file
  char *  ragaId2Rec_test_file = argv[4];  // yaml file
  char *  ragaId_list_file = argv[5];        // ragas to be tested
  char *  outfilename = argv[6];

 


  /* 
  // getting input arguments 
  char * ctrl_file = "ctrl_files/ctrl_file_slcs";
  char * ragaId2ioi_test_file =  "RagaId2ioi_test.yaml";
  char * ragaId2Rec_train_file =  "RagaId2RecPaths_train.yaml";
  char * ragaId2Rec_test_file =  "RagaId2RecPaths_test.yaml";
  char * ragaId_list_file =  "ragaId_list_shrey" ;
  char * outfilename = "outfile_temp.yaml";          
  */
    


  // no need to read the ctrl_file
  // reading all the yaml files

  FILE * ragaId2ioi_test_fp = fopen(ragaId2ioi_test_file, "r");
  assert(ragaId2ioi_test_fp != NULL);
  
  FILE * ragaId2Rec_train_fp = fopen(ragaId2Rec_train_file, "r");
  assert(ragaId2Rec_train_fp != NULL);

  FILE * ragaId2Rec_test_fp = fopen(ragaId2Rec_test_file, "r");
  assert(ragaId2Rec_test_fp != NULL);


  // reading raga list file
  FILE * ragaId_list_fp = fopen(ragaId_list_file, "r");
  assert(ragaId_list_fp != NULL);


  // getting yamlfiles to a data structure
  IntVector * ragaId2ioi_test = yaml_read_file_type1(ragaId2ioi_test_fp); 
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

  short temp_indexes[30000][3];
  int z_i = 0;
  int ragai;

  MemberVectorPTR member_vec_ptr = create_MemberVectorPTR();
  alloc_members_MemberVectorPTR(member_vec_ptr, 1000);

  for (ragai = 0; ragai < ragaIds_ptr->size; ragai++){
    int ragaId = ragaIds_ptr->vec[ragai];
    int t_i;

    for (t_i = 0; t_i < ragaId2ioi_test[ragaId].size; t_i++ ) {
      int test_i = ragaId2ioi_test[ragaId].vec[t_i];      
      SegmentPTR test_segment_ptr = create_SegmentPTR();
      fill_SegmentPTR_from_file(ragaId2Rec_test[ragaId].vec[test_i].data, test_segment_ptr);    

      MemberVector member_vec[ragaId2Rec_train[ragaId].size];
      int train_i;
      #pragma omp parallel for
      for (train_i = 0; train_i < ragaId2Rec_train[ragaId].size; train_i++) {      

	SegmentPTR train_raga_segment_ptr = create_SegmentPTR();
	fill_SegmentPTR_from_file(ragaId2Rec_train[ragaId].vec[train_i].data, train_raga_segment_ptr);

	alloc_members_MemberVectorPTR(&member_vec[train_i], 100);	
        
	motif_search_in_reference(test_segment_ptr, train_raga_segment_ptr, ctrl_file, &member_vec[train_i]);
        int j;  
	for(j = 0; j < member_vec[train_i].size; j++) {
	  strcpy(member_vec[train_i].vec[j].querypath, ragaId2Rec_test[ragaId].vec[test_i].data);
	  strcpy(member_vec[train_i].vec[j].refpath, ragaId2Rec_train[ragaId].vec[train_i].data);
	}

	empty_SegmentPTR(train_raga_segment_ptr);
	destroy_SegmentPTR(train_raga_segment_ptr);
      }
      
      // get the member with highest score
      int h_i, h_j;
      double max_score = -100.0;
      int m_i, m_j;
      for (h_i = 0; h_i < ragaId2Rec_train[ragaId].size; h_i++) {
	for(h_j = 0; h_j < member_vec[h_i].size; h_j++){
	  // to check validy
          if (member_vec[h_i].vec[h_j].alignment_vec_ptr->size!=0)
	    if (member_vec[h_i].vec[h_j].score > max_score) {
	      max_score = member_vec[h_i].vec[h_j].score;
	      m_i = h_i;
	      m_j = h_j;
	    }
	}
      }

      int si = member_vec_ptr->size;
      copy_MemberPTR(&member_vec_ptr->vec[si], &member_vec[m_i].vec[m_j]);
      member_vec_ptr->size += 1;

      for (h_i = 0; h_i < ragaId2Rec_train[ragaId].size; h_i++) {
	free_members_MemberVectorPTR(&member_vec[h_i]);	
      }

      empty_SegmentPTR(test_segment_ptr);
      destroy_SegmentPTR(test_segment_ptr);
    }


  }
  
  FILE * fout = fopen(outfilename, "w");
  printf_MemberVectorPTR(fout, member_vec_ptr, 0, 0);
  fclose(fout);


  free_members_MemberVectorPTR(member_vec_ptr);
  destroy_MemberVectorPTR(member_vec_ptr);

  free_members_IntVectorPTR(ragaIds_ptr);
  destroy_IntVectorPTR(ragaIds_ptr);
 
  yaml_delete_type1_ds(ragaId2ioi_test);
  yaml_delete_type2_ds(ragaId2Rec_train);
  yaml_delete_type2_ds(ragaId2Rec_test);

  fclose(ragaId2ioi_test_fp);
  fclose(ragaId2Rec_train_fp);
  fclose(ragaId2Rec_test_fp);
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
      // printf("###### %f \n", member_vec_ptr -> vec[si].score);
      member_vec_ptr -> size += 1;
      //break;
    }
    //  break;
  }
  destroy_ParametersRLCSPTR(parameters_slcs_ptr);
}
