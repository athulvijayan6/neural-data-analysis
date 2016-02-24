/* Description: 
 * Creation Date: 9 March 2014
 * Author: Shrey Dutta
 */


#include "rlcs.h"
//#define min3(a,b,c) a<b?(a<c?a:c):(b<c?b:c)


void RLCS_window(const DoubleMatrixPTR  query_feat_ptr, const DoubleMatrixPTR reference_feat_ptr, const ParametersRLCSPTR parameters_rlcs_ptr, const int start_q, const int end_q, const int start, const int end, const DoubleMatrixPTR shift_points_ptr, const int hop_size, const char *  reference_path, const char *  query_path, SpotVectorPTR spot_vec_ptr, const int isFilter, const SaddlePointsPTR query_sp_ptr, const SaddlePointsPTR reference_sp_ptr) {

  /************* Approach *****************
   * - first check whether the reference region given by start and end
   *   is too small or not in comparison with the query.
   * - if the reference region is too small, do not do anything, just return
   * - if it is not that small, then padd the region by zeros.
   * - Call the helper RLCS_window0 function to do windowed rlcs.
   *****************************************
   */ 


  // how large query is from the reference region.
  //  int padd_length = query_feat_ptr->size - (end - start + 1)>0?query_feat_ptr->size - (end - start + 1) :0;
  int   q_size  = end_q - start_q + 1;
  int padd_length = q_size - (end - start + 1)>0?q_size - (end - start + 1) :0;

  // if reference region is too small the return
  if (padd_length >  parameters_rlcs_ptr->safe_region_factor * (q_size)) {
    printf("Reference is too small to compare.\n");
    return;
  }
  RLCS_window0(query_feat_ptr, reference_feat_ptr, parameters_rlcs_ptr, start_q, end_q, start, end, padd_length, shift_points_ptr, hop_size, reference_path, query_path, spot_vec_ptr, isFilter, query_sp_ptr, reference_sp_ptr);
}

void RLCS_window0(const DoubleMatrixPTR query_feat_ptr, const DoubleMatrixPTR reference_feat_ptr, const ParametersRLCSPTR parameters_rlcs_ptr, const int start_q, const int end_q, const int start, const int end, const int padd_length, const DoubleMatrixPTR shift_points_ptr, const int hop_size, const char *  reference_path, const char *  query_path, SpotVectorPTR spot_vec_ptr, const int isFilter, const SaddlePointsPTR query_sp_ptr, const SaddlePointsPTR reference_sp_ptr) {
  
  // HELPER FUNCTION FOR RLCS_window

  /*************************************
   * - if padd_length is zero, extract the reference region given by start and end from the reference.
   * - if padd_length is non-zero, extract the reference region (make a copy) and padd it with those many zero vectors of the same dimension.
   * - for each windowed region call RLCS
   *************************************/
  int q_size = end_q - start_q + 1;
  printf("query length: %d, roi region: [%d %d]\n", q_size, start, end);
  printf("roi region length: %d, padd length: %d\n", end-start+1, padd_length);


  DoubleMatrixPTR region_reference_ptr = create_DoubleMatrixPTR(); // feature array to store the reference region.


  //  if padd_length is zero, extrat the reference region given by start and end from the reference.
  int region_alloc_flag  =  0;
  if (!padd_length) {
    // just passing the reference instead of copying because it is fast.
    region_reference_ptr->mat = &(reference_feat_ptr->mat[start]);
    region_reference_ptr->dim = reference_feat_ptr->dim;
    region_reference_ptr->size = end-start+1;
    region_reference_ptr->max_size = end-start+1;
   
    region_alloc_flag = 0;
    // roi is not important for this purpose
  } else {
    // extract the region and make a copy of it then padd that many zeros.
    alloc_members_DoubleMatrixPTR(region_reference_ptr, end-start+1+padd_length, reference_feat_ptr->dim);
    
    // copy the region of reference  from start to end.
    copy_DoubleMatrixPTR_region(region_reference_ptr, reference_feat_ptr, start, end);
    
    // padd zeros to make it equal to motif length
    int i, j;
    for (i = end-start+1; i < region_reference_ptr->max_size; i++) {
      // add a zero vector of the same dimension
      for (j = 0; j < region_reference_ptr->dim; j++)
	region_reference_ptr->mat[i][j] = 0;
      
      region_reference_ptr->size += 1;
    }
    region_alloc_flag = 1;
  }

  DoubleMatrixPTR region_query_ptr = create_DoubleMatrixPTR();
  region_query_ptr->mat = &(query_feat_ptr->mat[start_q]);
  region_query_ptr->dim = query_feat_ptr->dim;
  region_query_ptr->size = end_q-start_q+1;
  region_query_ptr->max_size = end_q-start_q+1;

  SaddlePointsPTR region_reference_sp_ptr = create_SaddlePointsPTR();
  SaddlePointsPTR region_query_sp_ptr = create_SaddlePointsPTR();

  int sp_start = arg_find_smallest_num_after_D(reference_sp_ptr->x_ptr, start);
  int sp_end = arg_find_largest_num_before_D(reference_sp_ptr->x_ptr, end);

  alloc_members_SaddlePointsPTR(region_reference_sp_ptr, sp_end -  sp_start + 1);
  
  copy_SaddlePointsPTR_region(region_reference_sp_ptr, reference_sp_ptr, sp_start, sp_end);
  
  // correcting indexes
  int ii;
  for (ii = 0; ii < region_reference_sp_ptr->size; ii++)
    region_reference_sp_ptr->x_ptr->mat[ii][0] = region_reference_sp_ptr->x_ptr->mat[ii][0] - start;
 
    // for query
  int sp_start_q = arg_find_smallest_num_after_D(query_sp_ptr->x_ptr, start_q);
  int sp_end_q = arg_find_largest_num_before_D(query_sp_ptr->x_ptr, end_q);


  alloc_members_SaddlePointsPTR(region_query_sp_ptr, sp_end_q - sp_start_q + 1);
  
  copy_SaddlePointsPTR_region(region_query_sp_ptr, query_sp_ptr, sp_start_q, sp_end_q);
  
  // correcting indexes
  
  for (ii = 0; ii < region_query_sp_ptr->size; ii++)
    region_query_sp_ptr->x_ptr->mat[ii][0] = region_query_sp_ptr->x_ptr->mat[ii][0] - start_q;


  
  int start_w; // start and end of the analysis window. 
  int end_w;

  // defining analysis window length
  double window_size = floor(q_size * parameters_rlcs_ptr->reference_window_size_ratio);


  // if shift_points.size = 0 that means that there are no shift points given and only worry about hop size.
  // if shift points are given that means analysis window should be placed at shift points and hop size referes to hopping the shift points.
  int spi;
  if (shift_points_ptr->size == 0)
    spi = -1;
  else 
    // gives the index of closest shift point from the start.
    spi = arg_find_smallest_num_after_D(shift_points_ptr, start);

  // setting the initial start end of the window.
  start_w = 0;
  end_w = start_w + window_size - 1 > region_reference_ptr->size ? region_reference_ptr->size - 1 : start_w + window_size - 1;  


  // run rlcs on each window.
  /*****************
   * - for each window. 
   ***** - get the score, cost, war, waq, diag matrices from rlcs.
   ***** - get the best sequence from this result and add it to member array
   * - from the member array group the member into some groups and append
   *   these groups to spots (output).
   ************************
   */


  if (parameters_rlcs_ptr->cpu == 1) {
    // code for doing rlcs for each window using  cpu   

    OutputRLCSPTR rlcs_out_ptr = create_OutputRLCSPTR();    // for storing the rlcs output matrices.
    MemberVectorPTR member_vec_ptr = create_MemberVectorPTR();  //to store the members  before grouping them;
    alloc_members_MemberVectorPTR(member_vec_ptr, 500);
    //int c__i = 0;
    while (end_w < region_reference_ptr->size) {
      // for debugging
      // if (c__i == 55) {
      //	printf("I am in c__i\n");
      //}
      //printf("w_i : %d\n", c__i); c__i++;
      
      if (!strcmp(parameters_rlcs_ptr->distname, "dtw")) {
        // calculate length in terms of saddle points

	int rlcontext = parameters_rlcs_ptr->rlcontext;

	int sp_startr = arg_find_smallest_num_after_D(region_reference_sp_ptr->x_ptr, start_w);
	int sp_endr = arg_find_largest_num_before_D(region_reference_sp_ptr->x_ptr, end_w);

	//	int sp_startq = arg_find_smallest_num_after_D(query_sp_ptr->x_ptr, start_q);
	//int sp_endq = arg_find_largest_num_before_D(region_reference_sp_ptr->x_ptr, end_w);
  

	alloc_members_OutputRLCSPTR(rlcs_out_ptr, (sp_endr - 2*rlcontext) - sp_startr + 1 + 1, (sp_end_q - 2*rlcontext) - sp_start_q + 1 + 1);
    
      } else {
	alloc_members_OutputRLCSPTR(rlcs_out_ptr, end_w - start_w + 1  + 1, region_query_ptr->size + 1); // +1 for additional zero vector. 
      }

      // find the output matrices from RLCS or RLCS_MOD
      RLCS(region_query_ptr, region_reference_ptr, start_w, end_w,  parameters_rlcs_ptr, rlcs_out_ptr, isFilter, region_query_sp_ptr, region_reference_sp_ptr);

      // find the best sequence and append it to member_vec_ptr.
      if (!strcmp(parameters_rlcs_ptr->distname, "dtw")){
	//	int start_w_sp = arg_find_smallest_num_after_D(region_reference_sp_ptr->x_ptr, start_w);
 	//find_best_sequence(rlcs_out_ptr, parameters_rlcs_ptr, start_q, start+ start_w_sp, padd_length, member_vec_ptr);
	// or this would do

	int start_final_r = arg_find_smallest_num_after_D(reference_sp_ptr->x_ptr, start+start_w);
        int start_final_q = arg_find_smallest_num_after_D(query_sp_ptr -> x_ptr, start_q);
	find_best_sequence(rlcs_out_ptr, parameters_rlcs_ptr, start_final_q, start_final_r, padd_length, member_vec_ptr);
      }else{

	find_best_sequence(rlcs_out_ptr, parameters_rlcs_ptr, start_q, start+ start_w, padd_length, member_vec_ptr);

      }

      start_w = spi == -1 ? start_w + hop_size: shift_points_ptr->mat[spi+hop_size][0]-start;
      end_w = start_w + window_size - 1;
      spi = spi == -1 ? -1 : spi + hop_size;

      free_members_OutputRLCSPTR(rlcs_out_ptr);
    }

    // group members


    group_members(member_vec_ptr, parameters_rlcs_ptr, reference_path, query_path, spot_vec_ptr);



    //now write this func(done) and dtw func(TODO) and correct where you make rlcs_out_matrix (done); aasta la vista baby!!!
    if (!strcmp(parameters_rlcs_ptr->distname, "dtw")) {
      // the start and end and all the refsindexes and queries are interms of saddle points we need to change it
      // to the indexes of the main feat file.
      // we have to take care of query indexes and ref indexes even after that because they wont be continous 
      // and will be only with respect to saddle points.
      correct_SpotVec_for_dtw_dist(spot_vec_ptr, reference_sp_ptr, query_sp_ptr);
    }


    // do not free members of MemberVectorPTR because those access by some members of spot_vec_ptr, after group_members were called.
    // free_members_MemberVectorPTR(member_vec_ptr);

    // code changed again now you can free

    free_members_MemberVectorPTR(member_vec_ptr);
    


    // but you can destroy the pointer.
    destroy_MemberVectorPTR(member_vec_ptr);
    destroy_OutputRLCSPTR(rlcs_out_ptr);

  } else if (parameters_rlcs_ptr->cpu > 1){
    // TODO: for more than one cpu
    
  } else {
    // TODO: for doing it on a gpu

  }

  if (region_alloc_flag)
    free_members_DoubleMatrixPTR(region_reference_ptr);




  destroy_DoubleMatrixPTR(region_reference_ptr);
  destroy_DoubleMatrixPTR(region_query_ptr);
  //printf("out of here\n");


  
  free_members_SaddlePointsPTR(region_reference_sp_ptr);
  free_members_SaddlePointsPTR(region_query_sp_ptr);
  destroy_SaddlePointsPTR(region_reference_sp_ptr);
  destroy_SaddlePointsPTR(region_query_sp_ptr);
}




void RLCS_perform(const DoubleMatrixPTR  query_feat_ptr, const DoubleMatrixPTR reference_feat_ptr, const ParametersRLCSPTR parameters_rlcs_ptr, const int start_q, const int end_q, const int start, const int end, MemberPTR member_ptr, const int isFilter, const SaddlePointsPTR query_sp_ptr, const SaddlePointsPTR reference_sp_ptr) {

  // how large query is from the reference region.
  //  int padd_length = query_feat_ptr->size - (end - start + 1)>0?query_feat_ptr->size - (end - start + 1) :0;
  int   q_size  = end_q - start_q + 1;
  int padd_length = q_size - (end - start + 1)>0?q_size - (end - start + 1) :0;

  // if reference region is too small the return
  if (padd_length >  parameters_rlcs_ptr->safe_region_factor * (q_size)) {
    printf("Reference is too small to compare.\n");
    return;
  }
  RLCS_perform0(query_feat_ptr, reference_feat_ptr, parameters_rlcs_ptr, start_q, end_q, start, end, padd_length, member_ptr, isFilter, query_sp_ptr, reference_sp_ptr);
}

void RLCS_perform0(const DoubleMatrixPTR query_feat_ptr, const DoubleMatrixPTR reference_feat_ptr, const ParametersRLCSPTR parameters_rlcs_ptr, const int start_q, const int end_q, const int start, const int end, const int padd_length, MemberPTR member_ptr, const int isFilter, const SaddlePointsPTR query_sp_ptr, const SaddlePointsPTR reference_sp_ptr) {
  
  // HELPER FUNCTION FOR 

  /*************************************
   * - if padd_length is zero, extract the reference region given by start and end from the reference.
   * - if padd_length is non-zero, extract the reference region (make a copy) and padd it with those many zero vectors of the same dimension.
   * - for each windowed region call RLCS
   *************************************/
  int q_size = end_q - start_q + 1;
  printf("query length: %d, roi region: [%d %d]\n", q_size, start, end);
  printf("roi region length: %d, padd length: %d\n", end-start+1, padd_length);


  DoubleMatrixPTR region_reference_ptr = create_DoubleMatrixPTR(); // feature array to store the reference region.


  //  if padd_length is zero, extrat the reference region given by start and end from the reference.
  int region_alloc_flag  =  0;
  if (!padd_length) {
    // just passing the reference instead of copying because it is fast.
    region_reference_ptr->mat = &(reference_feat_ptr->mat[start]);
    region_reference_ptr->dim = reference_feat_ptr->dim;
    region_reference_ptr->size = end-start+1;
    region_reference_ptr->max_size = end-start+1;
   
    region_alloc_flag = 0;
    // roi is not important for this purpose
  } else {
    // extract the region and make a copy of it then padd that many zeros.
    alloc_members_DoubleMatrixPTR(region_reference_ptr, end-start+1+padd_length, reference_feat_ptr->dim);
    
    // copy the region of reference  from start to end.
    copy_DoubleMatrixPTR_region(region_reference_ptr, reference_feat_ptr, start, end);
    
    // padd zeros to make it equal to motif length
    int i, j;
    for (i = end-start+1; i < region_reference_ptr->max_size; i++) {
      // add a zero vector of the same dimension
      for (j = 0; j < region_reference_ptr->dim; j++)
	region_reference_ptr->mat[i][j] = 0;
      
      region_reference_ptr->size += 1;
    }
    region_alloc_flag = 1;
  }

  DoubleMatrixPTR region_query_ptr = create_DoubleMatrixPTR();
  region_query_ptr->mat = &(query_feat_ptr->mat[start_q]);
  region_query_ptr->dim = query_feat_ptr->dim;
  region_query_ptr->size = end_q-start_q+1;
  region_query_ptr->max_size = end_q-start_q+1;


  OutputRLCSPTR rlcs_out_ptr = create_OutputRLCSPTR();    // for storing the rlcs output matrices.
  alloc_members_OutputRLCSPTR(rlcs_out_ptr, region_reference_ptr->size+1, region_query_ptr->size +1);
  
  //alloc_members_MemberPTR(member_ptr, 1);


  SaddlePointsPTR region_reference_sp_ptr = create_SaddlePointsPTR();
  SaddlePointsPTR region_query_sp_ptr = create_SaddlePointsPTR();


  if (isFilter){
    int sp_start = arg_find_smallest_num_after_D(reference_sp_ptr->x_ptr, start);
    int sp_end = arg_find_largest_num_before_D(reference_sp_ptr->x_ptr, end);
  
    alloc_members_SaddlePointsPTR(region_reference_sp_ptr, sp_end -  sp_start + 1);
  
    copy_SaddlePointsPTR_region(region_reference_sp_ptr, reference_sp_ptr, sp_start, sp_end);
  
    // correcting indexes
    int ii;
    for (ii = 0; ii < region_reference_sp_ptr->size; ii++)
      region_reference_sp_ptr->x_ptr->mat[ii][0] = region_reference_sp_ptr->x_ptr->mat[ii][0] - start;
    
    // for query


    int sp_start_q = arg_find_smallest_num_after_D(query_sp_ptr->x_ptr, start_q);
    int sp_end_q = arg_find_largest_num_before_D(query_sp_ptr->x_ptr, end_q);


    alloc_members_SaddlePointsPTR(region_query_sp_ptr, sp_end_q - sp_start_q + 1);
  
    copy_SaddlePointsPTR_region(region_query_sp_ptr, query_sp_ptr, sp_start_q, sp_end_q);
  
  // correcting indexes
  
    for (ii = 0; ii < region_query_sp_ptr->size; ii++)
      region_query_sp_ptr->x_ptr->mat[ii][0] = region_query_sp_ptr->x_ptr->mat[ii][0] - start_q;
  }

  // dummy vars
  //SaddlePointsPTR dummy1 = create_SaddlePointsPTR();
  //SaddlePointsPTR dummy2 = create_SaddlePointsPTR();

  RLCS(region_query_ptr, region_reference_ptr, 0, region_reference_ptr->size-1,  parameters_rlcs_ptr, rlcs_out_ptr, isFilter, region_query_sp_ptr, region_reference_sp_ptr);
  
  if (isFilter){
    free_members_SaddlePointsPTR(region_reference_sp_ptr);
    free_members_SaddlePointsPTR(region_query_sp_ptr);
  }


  destroy_SaddlePointsPTR(region_reference_sp_ptr);
  destroy_SaddlePointsPTR(region_query_sp_ptr);

  //destroying dummy vars
  //  destroy_SaddlePointsPTR(dummy1);
  //destroy_SaddlePointsPTR(dummy2);

  find_best_sequence_rlcs(rlcs_out_ptr, parameters_rlcs_ptr, start_q, start, padd_length, member_ptr);

  free_members_OutputRLCSPTR(rlcs_out_ptr);
  destroy_OutputRLCSPTR(rlcs_out_ptr);

  if (region_alloc_flag)
    free_members_DoubleMatrixPTR(region_reference_ptr);

  destroy_DoubleMatrixPTR(region_reference_ptr);
  destroy_DoubleMatrixPTR(region_query_ptr);
}





void RLCS(const DoubleMatrixPTR query_feat_ptr, const DoubleMatrixPTR region_reference_ptr, const int start, const int end, const ParametersRLCSPTR parameters_rlcs_ptr, OutputRLCSPTR rlcs_out_ptr, const int isFilter, const SaddlePointsPTR query_sp_ptr, const SaddlePointsPTR region_reference_sp_ptr) {
  
  assert(end < region_reference_ptr->size);

  int i,j; // iterators

  // initilizing 1st col to zeros
  for (i = 0; i < rlcs_out_ptr->row; i++) {
    rlcs_out_ptr->cost_ptr->mat[i][0] = 0;
    rlcs_out_ptr->cost_actual_ptr->mat[i][0] = 0;
    rlcs_out_ptr->war_ptr->mat[i][0] = 0;
    rlcs_out_ptr->waq_ptr->mat[i][0] = 0;
    rlcs_out_ptr->score_ptr->mat[i][0] = 0;
  }
  // initizing 1st row to zeros
  for (j = 0; j < rlcs_out_ptr->col; j++) {
    rlcs_out_ptr->cost_ptr->mat[0][j] = 0;
    rlcs_out_ptr->cost_actual_ptr->mat[0][j] = 0;
    rlcs_out_ptr->war_ptr->mat[0][j] = 0;
    rlcs_out_ptr->waq_ptr->mat[0][j] = 0;
    rlcs_out_ptr->score_ptr->mat[0][j] = 0;
  }

  
  for (i = 1; i < rlcs_out_ptr->row; i++) {

    for (j = 1; j < rlcs_out_ptr->col; j++) {
      
      double dist = compute_distance(region_reference_ptr, query_feat_ptr, start, i-1, j-1, parameters_rlcs_ptr, query_sp_ptr, region_reference_sp_ptr);
     
      if (dist < parameters_rlcs_ptr->Td) {
	rlcs_out_ptr->diag_ptr->mat[i][j] = 1; // for cross

        rlcs_out_ptr->cost_ptr->mat[i][j] = rlcs_out_ptr->cost_ptr->mat[i-1][j-1] + (1 - dist/parameters_rlcs_ptr->Td);
	rlcs_out_ptr->cost_actual_ptr->mat[i][j] = rlcs_out_ptr->cost_actual_ptr->mat[i-1][j-1] + 1;

	if (!strcmp(parameters_rlcs_ptr->algo, "rlcs")) {
	  rlcs_out_ptr->war_ptr->mat[i][j] = rlcs_out_ptr->war_ptr->mat[i-1][j-1] + 1;
	  rlcs_out_ptr->waq_ptr->mat[i][j] = rlcs_out_ptr->waq_ptr->mat[i-1][j-1] + 1;
	} else if (!strcmp(parameters_rlcs_ptr->algo, "rlcs_mod")) {
	  rlcs_out_ptr->war_ptr->mat[i][j] = rlcs_out_ptr->war_ptr->mat[i-1][j-1] + (1 - dist/parameters_rlcs_ptr->Td);
	  rlcs_out_ptr->waq_ptr->mat[i][j] = rlcs_out_ptr->waq_ptr->mat[i-1][j-1] + (1 - dist/parameters_rlcs_ptr->Td);
	}
      } else if (rlcs_out_ptr->cost_ptr->mat[i-1][j] > rlcs_out_ptr->cost_ptr->mat[i][j-1]) {
	
	rlcs_out_ptr->diag_ptr->mat[i][j] = 2;  // for up

	rlcs_out_ptr->cost_ptr->mat[i][j] = rlcs_out_ptr->cost_ptr->mat[i-1][j];
	rlcs_out_ptr->cost_actual_ptr->mat[i][j] = rlcs_out_ptr->cost_actual_ptr->mat[i-1][j];
	

	rlcs_out_ptr->waq_ptr->mat[i][j] = rlcs_out_ptr->waq_ptr->mat[i-1][j];
	if (rlcs_out_ptr->war_ptr->mat[i-1][j] > 0)
	  rlcs_out_ptr->war_ptr->mat[i][j] = rlcs_out_ptr->war_ptr->mat[i-1][j] + 1;
	else
	  rlcs_out_ptr->war_ptr->mat[i][j] = 0;

      } else {
	rlcs_out_ptr->diag_ptr->mat[i][j] = 3; // for left

	rlcs_out_ptr->cost_ptr->mat[i][j] = rlcs_out_ptr->cost_ptr->mat[i][j-1];
	rlcs_out_ptr->cost_actual_ptr->mat[i][j] = rlcs_out_ptr->cost_actual_ptr->mat[i][j-1];

	rlcs_out_ptr->war_ptr->mat[i][j] = rlcs_out_ptr->war_ptr->mat[i][j-1];
	if (rlcs_out_ptr->waq_ptr->mat[i][j-1] > 0)
	  rlcs_out_ptr->waq_ptr->mat[i][j] = rlcs_out_ptr->waq_ptr->mat[i][j-1] + 1;
	else
	  rlcs_out_ptr->waq_ptr->mat[i][j] = 0;
      }

      if (!strcmp(parameters_rlcs_ptr->algo, "rlcs")) {
	
	if (rlcs_out_ptr->war_ptr->mat[i][j] && rlcs_out_ptr->waq_ptr->mat[i][j] && rlcs_out_ptr->cost_ptr->mat[i][j] > parameters_rlcs_ptr->rho*query_feat_ptr->size) {

	  rlcs_out_ptr->score_ptr->mat[i][j] = ( (parameters_rlcs_ptr->beta * rlcs_out_ptr->cost_ptr->mat[i][j]/rlcs_out_ptr->war_ptr->mat[i][j]) + ((1-parameters_rlcs_ptr->beta) * rlcs_out_ptr->cost_ptr->mat[i][j]/rlcs_out_ptr->waq_ptr->mat[i][j]) ) * (rlcs_out_ptr->cost_ptr->mat[i][j]/query_feat_ptr->size);
	
	} else
	  rlcs_out_ptr->score_ptr->mat[i][j] = 0;  

      } else if (!strcmp(parameters_rlcs_ptr->algo, "rlcs_mod")) {
	
	if (rlcs_out_ptr->war_ptr->mat[i][j] && rlcs_out_ptr->waq_ptr->mat[i][j] && rlcs_out_ptr->cost_actual_ptr->mat[i][j] > parameters_rlcs_ptr->rho*query_feat_ptr->size) {
  
	  rlcs_out_ptr->score_ptr->mat[i][j] = ( (parameters_rlcs_ptr->beta * rlcs_out_ptr->cost_ptr->mat[i][j]/rlcs_out_ptr->war_ptr->mat[i][j]) + ((1-parameters_rlcs_ptr->beta) * rlcs_out_ptr->cost_ptr->mat[i][j]/rlcs_out_ptr->waq_ptr->mat[i][j]) ) * ((rlcs_out_ptr->cost_ptr->mat[i][j] + rlcs_out_ptr->cost_ptr->mat[i][j]) / (rlcs_out_ptr->cost_actual_ptr->mat[i][j] + query_feat_ptr->size));
	
	} else
	  rlcs_out_ptr->score_ptr->mat[i][j] = 0;  
      }


      if (isFilter && rlcs_out_ptr->score_ptr->mat[i][j] != 0) {
	int x = i, y = j;
	while (x != 0 && y != 0 && rlcs_out_ptr->cost_ptr->mat[x][y] != 0) {
	  if (rlcs_out_ptr->diag_ptr->mat[x][y] == 1)
	    x -= 1, y -= 1;
	  else if (rlcs_out_ptr->diag_ptr->mat[x][y] == 2)
	    x -= 1;
	  else if (rlcs_out_ptr->diag_ptr->mat[x][y] == 3)
	    y -= 1;
	}

	int f_i = x==0?1:x;
	int f_j = y==0?1:y;

	// (not really) f_i and f_j represent the first index where score value is non-zero
	int s_r = arg_find_smallest_num_after_D(region_reference_sp_ptr->x_ptr, start + f_i-1);
	int e_r = arg_find_largest_num_before_D(region_reference_sp_ptr->x_ptr, start + i-1);
	int s_q = arg_find_smallest_num_after_D(query_sp_ptr->x_ptr, 0 + f_j-1);
	int e_q = arg_find_largest_num_before_D(query_sp_ptr->x_ptr, 0 + j-1);
      
	double rho =0.0;
	if (s_r != e_r && s_q != e_q) {
          // Changed 9 Feb 2015
	  /* double muQ = mean_of_first_diff(query_sp_ptr->y_ptr, s_q, e_q); */
	  /* double sigmaQ = std_of_first_diff(query_sp_ptr->y_ptr, s_q, e_q); */

	  /* double muR = mean_of_first_diff(region_reference_sp_ptr->y_ptr, s_r, e_r); */
	  /* double sigmaR = std_of_first_diff(region_reference_sp_ptr->y_ptr, s_r, e_r); */

	  double muQ = slope_of_best_linear_fit(query_sp_ptr->y_ptr, s_q, e_q); // used as slope of linear trend.
	  double sigmaQ = std_of_slope(query_sp_ptr->y_ptr, s_q, e_q); // used as std of slope of linear trend.

	  double muR = slope_of_best_linear_fit(region_reference_sp_ptr->y_ptr, s_r, e_r);
	  double sigmaR = std_of_slope(region_reference_sp_ptr->y_ptr, s_r, e_r);

	
	  double Z1 = muQ * sigmaR;
	  double Z2 = muR * sigmaQ;
      
	  if (Z1 < 0 && Z2 < 0)
	    rho  = max(Z1,Z2)/min(Z1,Z2);
	  else if (Z1 == 0 || Z2 == 0)
	    rho  = 0;
	  else
	    rho = min(Z1,Z2) / max(Z1,Z2);
	}
	
	rlcs_out_ptr->score_ptr->mat[i][j] = 0.6*rlcs_out_ptr->score_ptr->mat[i][j] + 0.4*rho;	
      }

    }  
    
  }
}




void find_best_sequence_rlcs(const OutputRLCSPTR rlcs_out_ptr, const ParametersRLCSPTR  parameters_rlcs_ptr, const int query_index_correction,  const int reference_index_correction, const int padd_length, MemberPTR member_ptr) {

  IntVectorPTR max_score_index_ptr = create_IntVectorPTR();
  alloc_members_IntVectorPTR(max_score_index_ptr, 2);


  //  argmax_matrix_D_gen(rlcs_out_ptr->score_ptr, min_row, min_col, max_row, max_col, max_score_index_ptr);
  argmax_matrix_D_gen(rlcs_out_ptr->score_ptr, 0, 0, rlcs_out_ptr->row - padd_length -1, rlcs_out_ptr->col-1, max_score_index_ptr);
  
  assert(max_score_index_ptr->vec[0] < rlcs_out_ptr->row - padd_length);

  
  if (rlcs_out_ptr->score_ptr->mat[max_score_index_ptr->vec[0]][max_score_index_ptr->vec[1]] <= parameters_rlcs_ptr->seqFilterTd) {
    free_members_IntVectorPTR(max_score_index_ptr);
    destroy_IntVectorPTR(max_score_index_ptr);
    return;
  }

  int i;

  //  int x = max_score_index_ptr->vec[0] >= rlcs_out_ptr->row - padd_length ? rlcs_out_ptr->row - padd_length -1 : max_score_index_ptr->vec[0]; 
  int x = max_score_index_ptr->vec[0];
  int y = max_score_index_ptr->vec[1];
  
  free_members_IntVectorPTR(max_score_index_ptr);
  destroy_IntVectorPTR(max_score_index_ptr);
  

  member_ptr->score = rlcs_out_ptr->score_ptr->mat[x][y];
  
  member_ptr->cost = rlcs_out_ptr->cost_ptr->mat[x][y];
  member_ptr->cost_actual  = rlcs_out_ptr->cost_actual_ptr->mat[x][y];
  
  member_ptr->war = rlcs_out_ptr->war_ptr->mat[x][y];
  
  member_ptr->waq = rlcs_out_ptr->waq_ptr->mat[x][y];


  if (member_ptr->alignment_vec_ptr->max_size < rlcs_out_ptr->row + rlcs_out_ptr->col)
    realloc_members_AlignmentVectorPTR(member_ptr->alignment_vec_ptr, rlcs_out_ptr->row + rlcs_out_ptr->col);
  

  while (x != 0 && y != 0 && rlcs_out_ptr->cost_ptr->mat[x][y] != 0) {
    i = member_ptr->alignment_vec_ptr->size;

    member_ptr->alignment_vec_ptr->vec[i].reference = x - 1 + reference_index_correction;

    member_ptr->alignment_vec_ptr->vec[i].query = y - 1 + query_index_correction;

    if (rlcs_out_ptr->diag_ptr->mat[x][y] == 1)
      x -= 1, y -= 1;
    else if (rlcs_out_ptr->diag_ptr->mat[x][y] == 2)
      x -= 1;
    else if (rlcs_out_ptr->diag_ptr->mat[x][y] == 3)
      y -= 1;

    // remove this for debugging
    
    //printf("%d %d\n", x-1+reference_index_correction, y-1+query_index_correction);
    
    /// till hgere

    member_ptr->alignment_vec_ptr->size = i + 1;
  }
  assert(member_ptr->alignment_vec_ptr->size <= member_ptr->alignment_vec_ptr->max_size);

  // reverse the alignment to make it correct;
  int align_size = member_ptr->alignment_vec_ptr->size;
  int align_middle_index = align_size % 2 ? floor((align_size-1)/2) - 1 : floor((align_size-1)/2);
  for (i = 0; i <= align_middle_index; i++) {
    // swapping two numbers
    swap_int(&(member_ptr->alignment_vec_ptr->vec[i].query), &(member_ptr->alignment_vec_ptr->vec[align_size-i-1].query));
    swap_int(&(member_ptr->alignment_vec_ptr->vec[i].reference), &(member_ptr->alignment_vec_ptr->vec[align_size-i-1].reference));
  }

  // Store as sample number
  member_ptr->start = member_ptr->alignment_vec_ptr->vec[0].reference;
  member_ptr->end = member_ptr->alignment_vec_ptr->vec[align_size-1].reference;
  member_ptr->start_q = member_ptr->alignment_vec_ptr -> vec[0].query;
  member_ptr->end_q = member_ptr->alignment_vec_ptr -> vec[align_size - 1].query;
}


void find_best_sequence(const OutputRLCSPTR rlcs_out_ptr, const ParametersRLCSPTR  parameters_rlcs_ptr, const int query_index_correction,  const int reference_index_correction, const int padd_length,  MemberVectorPTR member_vec_ptr) {



  IntVectorPTR max_score_index_ptr = create_IntVectorPTR();
  alloc_members_IntVectorPTR(max_score_index_ptr, 2);


  //  argmax_matrix_D_gen(rlcs_out_ptr->score_ptr, min_row, min_col, max_row, max_col, max_score_index_ptr);
  argmax_matrix_D_gen(rlcs_out_ptr->score_ptr, 0, 0, rlcs_out_ptr->row - padd_length -1, rlcs_out_ptr->col-1, max_score_index_ptr);
  
  assert(max_score_index_ptr->vec[0] < rlcs_out_ptr->row - padd_length);


  if (rlcs_out_ptr->score_ptr->mat[max_score_index_ptr->vec[0]][max_score_index_ptr->vec[1]] <= parameters_rlcs_ptr->seqFilterTd) {
    free_members_IntVectorPTR(max_score_index_ptr);
    destroy_IntVectorPTR(max_score_index_ptr);
    return;
  }

  int i;




  //  int x = max_score_index_ptr->vec[0] >= rlcs_out_ptr->row - padd_length ? rlcs_out_ptr->row - padd_length -1 : max_score_index_ptr->vec[0]; 
  int x = max_score_index_ptr->vec[0];
  int y = max_score_index_ptr->vec[1];

  free_members_IntVectorPTR(max_score_index_ptr);
  destroy_IntVectorPTR(max_score_index_ptr);

  assert(member_vec_ptr->size <= member_vec_ptr->max_size);
  // reallocating the member if max_size less.
  if (member_vec_ptr->size == member_vec_ptr->max_size)
    realloc_members_MemberVectorPTR(member_vec_ptr, 2*(member_vec_ptr->max_size));
  

  int mi = member_vec_ptr->size;

  // incrementing the mem_vec_ptr size by because one entry is goign to be added in it.
  member_vec_ptr->size = member_vec_ptr->size + 1;


  member_vec_ptr->vec[mi].score = rlcs_out_ptr->score_ptr->mat[x][y];

  member_vec_ptr->vec[mi].cost = rlcs_out_ptr->cost_ptr->mat[x][y];
  member_vec_ptr->vec[mi].cost_actual  = rlcs_out_ptr->cost_actual_ptr->mat[x][y];

  member_vec_ptr->vec[mi].war = rlcs_out_ptr->war_ptr->mat[x][y];

  member_vec_ptr->vec[mi].waq = rlcs_out_ptr->waq_ptr->mat[x][y];

 

  if (member_vec_ptr->vec[mi].alignment_vec_ptr->max_size < rlcs_out_ptr->row + rlcs_out_ptr->col)
    realloc_members_AlignmentVectorPTR(member_vec_ptr->vec[mi].alignment_vec_ptr, rlcs_out_ptr->row + rlcs_out_ptr->col);
  

  while (x != 0 && y != 0 && rlcs_out_ptr->cost_ptr->mat[x][y] != 0) {
    i = member_vec_ptr->vec[mi].alignment_vec_ptr->size;

    member_vec_ptr->vec[mi].alignment_vec_ptr->vec[i].reference = x - 1 + reference_index_correction;

    member_vec_ptr->vec[mi].alignment_vec_ptr->vec[i].query = y - 1 + query_index_correction;

    if (rlcs_out_ptr->diag_ptr->mat[x][y] == 1)
      x -= 1, y -= 1;
    else if (rlcs_out_ptr->diag_ptr->mat[x][y] == 2)
      x -= 1;
    else if (rlcs_out_ptr->diag_ptr->mat[x][y] == 3)
      y -= 1;

    // remove this for debugging
    
    //printf("%d %d\n", x-1+reference_index_correction, y-1+query_index_correction);
    
    /// till hgere

    member_vec_ptr->vec[mi].alignment_vec_ptr->size = i + 1;
  }
  assert(member_vec_ptr->vec[mi].alignment_vec_ptr->size <= member_vec_ptr->vec[mi].alignment_vec_ptr->max_size);

  // reverse the alignment to make it correct;
  int align_size = member_vec_ptr->vec[mi].alignment_vec_ptr->size;
  int align_middle_index = align_size % 2 ? floor((align_size-1)/2) - 1 : floor((align_size-1)/2);
  for (i = 0; i <= align_middle_index; i++) {
    // swapping two numbers
    swap_int(&(member_vec_ptr->vec[mi].alignment_vec_ptr->vec[i].query), &(member_vec_ptr->vec[mi].alignment_vec_ptr->vec[align_size-i-1].query));
    swap_int(&(member_vec_ptr->vec[mi].alignment_vec_ptr->vec[i].reference), &(member_vec_ptr->vec[mi].alignment_vec_ptr->vec[align_size-i-1].reference));
  }

  // Store as sample number
  member_vec_ptr->vec[mi].start = member_vec_ptr->vec[mi].alignment_vec_ptr->vec[0].reference;
  member_vec_ptr->vec[mi].end = member_vec_ptr->vec[mi].alignment_vec_ptr->vec[align_size-1].reference;
  member_vec_ptr->vec[mi].start_q = member_vec_ptr->vec[mi].alignment_vec_ptr -> vec[0].query;
  member_vec_ptr->vec[mi].end_q = member_vec_ptr->vec[mi].alignment_vec_ptr -> vec[align_size - 1].query;
}

void group_members(const MemberVectorPTR member_vec_ptr, const ParametersRLCSPTR paramters_rlcs_ptr, const char *  reference_path, const char *  query_path,  SpotVectorPTR spot_vec_ptr) {
  
  if (!member_vec_ptr -> size) return;

  // sort members on the basis of their start timings.
  qsort(member_vec_ptr->vec, member_vec_ptr->size, sizeof(Member), members_comp);


  int i; // iterator over members;

  // Group Criteria: if the start of any member is greter than maximum end 
  // seen thus far. Then entire chunk goes to one group and the process is 
  // repeated.
  double max_end = member_vec_ptr->vec[0].end;  
  double max_end_q = member_vec_ptr->vec[0].end_q;
  double best_score = member_vec_ptr->vec[0].score;
  int best_score_index = 0;
  int from = 0; // index od the begging member of the group;

  // starting i from 1 is very important because i = 0 is already taken care of in max_end and best_score before
  for (i = 1; i < member_vec_ptr->size; i++) {
    if (member_vec_ptr->vec[i].start > max_end) {
    
      // check if there is space in spots vector to store this spot
      // and allocate more space if there is not much space.
      if (spot_vec_ptr->size == spot_vec_ptr->max_size)
	realloc_members_SpotVectorPTR(spot_vec_ptr, 2*(spot_vec_ptr->max_size)); // a<<1 means a*2


      int si = spot_vec_ptr->size;

      // clean up the existing members
      //      free_members_MemberVectorPTR(spot_vec_ptr->vec[si].member_vec_ptr);

      {
	int j;
	if (spot_vec_ptr -> vec[si].member_vec_ptr -> size + ((i-1) - from + 1) > spot_vec_ptr -> vec[si].member_vec_ptr -> max_size)
	  realloc_members_MemberVectorPTR(spot_vec_ptr -> vec[si].member_vec_ptr, 2*(spot_vec_ptr -> vec[si].member_vec_ptr -> max_size) + ((i-1) - from + 1));
	for (j = 0; j< (i-1) - from + 1; j++) {
	  copy_MemberPTR(&(spot_vec_ptr->vec[si].member_vec_ptr->vec[j]), &(member_vec_ptr->vec[from + j]));
          spot_vec_ptr->vec[si].member_vec_ptr->size += 1;
	} 
      }
      // spot_vec_ptr->vec[si].member_vec_ptr->vec = &(member_vec_ptr->vec[from]);
      assert(spot_vec_ptr->vec[si].member_vec_ptr->size == (i-1) - from + 1);
      //spot_vec_ptr->vec[si].member_vec_ptr->max_size = (i-1) - from + 1;
      strcpy(spot_vec_ptr->vec[si].reference_path, reference_path);
      strcpy(spot_vec_ptr->vec[si].query_path, query_path);
      spot_vec_ptr->vec[si].group_start = member_vec_ptr->vec[from].start;
      spot_vec_ptr->vec[si].group_end = max_end;
      spot_vec_ptr->vec[si].group_start_q = member_vec_ptr->vec[from].start_q;
      spot_vec_ptr->vec[si].group_end_q = max_end_q;
      
      copy_MemberPTR(spot_vec_ptr->vec[si].best_member_ptr, &(member_vec_ptr->vec[best_score_index]));

      spot_vec_ptr->size = si + 1;
      
      //for the next group
      from = i;
      max_end = member_vec_ptr->vec[from].end;
      max_end_q = member_vec_ptr->vec[from].end_q;
      best_score = member_vec_ptr->vec[from].score;
      best_score_index = from;
    } else {
      // update best_score, best_member_index, max_end if needed
      if (member_vec_ptr->vec[i].end > max_end) 
	max_end = member_vec_ptr->vec[i].end;

      if (member_vec_ptr->vec[i].end_q > max_end_q) 
	max_end_q = member_vec_ptr->vec[i].end_q;

      
      if (member_vec_ptr->vec[i].score > best_score){
	best_score = member_vec_ptr->vec[i].score;
	best_score_index = i;
      }
    }
  }
  
 
  // making the last group...
  if (spot_vec_ptr->size == spot_vec_ptr->max_size)
    realloc_members_SpotVectorPTR(spot_vec_ptr, 2*(spot_vec_ptr->max_size)); // a<<1 means a*2
   
  int si = spot_vec_ptr->size;

  // clean up the existing members
  //free_members_MemberVectorPTR(spot_vec_ptr->vec[si].member_vec_ptr);


  {
    int j;
    if (spot_vec_ptr -> vec[si].member_vec_ptr -> size + ((i-1) - from + 1) > spot_vec_ptr -> vec[si].member_vec_ptr -> max_size)
      realloc_members_MemberVectorPTR(spot_vec_ptr -> vec[si].member_vec_ptr, 2*(spot_vec_ptr -> vec[si].member_vec_ptr -> max_size) + ((i-1) - from + 1));
    for (j = 0; j< (i-1) - from + 1; j++) {
      copy_MemberPTR(&(spot_vec_ptr->vec[si].member_vec_ptr->vec[j]), &(member_vec_ptr->vec[from + j]));
      spot_vec_ptr->vec[si].member_vec_ptr->size += 1;
    } 
  }
  // spot_vec_ptr->vec[si].member_vec_ptr->vec = &(member_vec_ptr->vec[from]);
  assert(spot_vec_ptr->vec[si].member_vec_ptr->size == (i-1) - from + 1);
  //spot_vec_ptr->vec[si].member_vec_ptr->max_size = (i-1) - from + 1;

  strcpy(spot_vec_ptr->vec[si].reference_path, reference_path);
  strcpy(spot_vec_ptr->vec[si].query_path, query_path);
  spot_vec_ptr->vec[si].group_start = member_vec_ptr->vec[from].start;
  spot_vec_ptr->vec[si].group_end = max_end;
  spot_vec_ptr->vec[si].group_start_q = member_vec_ptr->vec[from].start_q;
  spot_vec_ptr->vec[si].group_end_q = max_end_q;
  
  copy_MemberPTR(spot_vec_ptr->vec[si].best_member_ptr, &(member_vec_ptr->vec[best_score_index]));  

  spot_vec_ptr->size = si + 1;
}


void SLCS_perform(const DoubleMatrixPTR  query_feat_ptr, const DoubleMatrixPTR reference_feat_ptr, const ParametersRLCSPTR parameters_slcs_ptr, const int start_q, const int end_q, const int start, const int end, MemberPTR member_ptr) {

  // how large query is from the reference region.
  //  int padd_length = query_feat_ptr->size - (end - start + 1)>0?query_feat_ptr->size - (end - start + 1) :0;
  int   q_size  = end_q - start_q + 1;
  int padd_length = q_size - (end - start + 1)>0?q_size - (end - start + 1) :0;

  // if reference region is too small the return
  if (padd_length >  parameters_slcs_ptr->safe_region_factor * (q_size)) {
    printf("Reference is too small to compare.\n");
    return;
  }
  SLCS_perform0(query_feat_ptr, reference_feat_ptr, parameters_slcs_ptr, start_q, end_q, start, end, padd_length, member_ptr);
}

void SLCS_perform0(const DoubleMatrixPTR query_feat_ptr, const DoubleMatrixPTR reference_feat_ptr, const ParametersRLCSPTR parameters_slcs_ptr, const int start_q, const int end_q, const int start, const int end, const int padd_length, MemberPTR member_ptr) {
  
  // HELPER FUNCTION FOR 

  /*************************************
   * - if padd_length is zero, extract the reference region given by start and end from the reference.
   * - if padd_length is non-zero, extract the reference region (make a copy) and padd it with those many zero vectors of the same dimension.
   * - for each windowed region call RLCS
   *************************************/
  int q_size = end_q - start_q + 1;
  printf("query length: %d, roi region: [%d %d]\n", q_size, start, end);
  printf("roi region length: %d, padd length: %d\n", end-start+1, padd_length);


  DoubleMatrixPTR region_reference_ptr = create_DoubleMatrixPTR(); // feature array to store the reference region.


  //  if padd_length is zero, extrat the reference region given by start and end from the reference.
  int region_alloc_flag  =  0;
  if (!padd_length) {
    // just passing the reference instead of copying because it is fast.
    region_reference_ptr->mat = &(reference_feat_ptr->mat[start]);
    region_reference_ptr->dim = reference_feat_ptr->dim;
    region_reference_ptr->size = end-start+1;
    region_reference_ptr->max_size = end-start+1;
   
    region_alloc_flag = 0;
    // roi is not important for this purpose
  } else {
    // extract the region and make a copy of it then padd that many zeros.
    alloc_members_DoubleMatrixPTR(region_reference_ptr, end-start+1+padd_length, reference_feat_ptr->dim);
    
    // copy the region of reference  from start to end.
    copy_DoubleMatrixPTR_region(region_reference_ptr, reference_feat_ptr, start, end);
    
    // padd zeros to make it equal to motif length
    int i, j;
    for (i = end-start+1; i < region_reference_ptr->max_size; i++) {
      // add a zero vector of the same dimension
      for (j = 0; j < region_reference_ptr->dim; j++)
	region_reference_ptr->mat[i][j] = 0;
      
      region_reference_ptr->size += 1;
    }
    region_alloc_flag = 1;
  }

  DoubleMatrixPTR region_query_ptr = create_DoubleMatrixPTR();
  region_query_ptr->mat = &(query_feat_ptr->mat[start_q]);
  region_query_ptr->dim = query_feat_ptr->dim;
  region_query_ptr->size = end_q-start_q+1;
  region_query_ptr->max_size = end_q-start_q+1;


  OutputSLCSPTR slcs_out_ptr = create_OutputSLCSPTR();    // for storing the rlcs output matrices.
  alloc_members_OutputSLCSPTR(slcs_out_ptr, region_reference_ptr->size+1, region_query_ptr->size +1);
  
  //  alloc_members_MemberPTR(member_ptr, 1);
  
  Segment_LCS_v2_2(region_query_ptr, region_reference_ptr, 0, region_reference_ptr->size-1,  parameters_slcs_ptr, slcs_out_ptr);

  find_best_sequence_slcs(slcs_out_ptr, parameters_slcs_ptr, start_q, start, padd_length, member_ptr);

  free_members_OutputSLCSPTR(slcs_out_ptr);
  destroy_OutputSLCSPTR(slcs_out_ptr);

  if (region_alloc_flag)
    free_members_DoubleMatrixPTR(region_reference_ptr);

  destroy_DoubleMatrixPTR(region_reference_ptr);
  destroy_DoubleMatrixPTR(region_query_ptr);
}




void SLCS_hard_perform(const DoubleMatrixPTR  query_feat_ptr, const DoubleMatrixPTR reference_feat_ptr, const ParametersRLCSPTR parameters_slcs_ptr, const int start_q, const int end_q, const int start, const int end, MemberPTR member_ptr) {

  // how large query is from the reference region.
  //  int padd_length = query_feat_ptr->size - (end - start + 1)>0?query_feat_ptr->size - (end - start + 1) :0;
  int   q_size  = end_q - start_q + 1;
  int padd_length = q_size - (end - start + 1)>0?q_size - (end - start + 1) :0;

  // if reference region is too small the return
  if (padd_length >  parameters_slcs_ptr->safe_region_factor * (q_size)) {
    printf("Reference is too small to compare.\n");
    return;
  }




  SLCS_hard_perform0(query_feat_ptr, reference_feat_ptr, parameters_slcs_ptr, start_q, end_q, start, end, padd_length, member_ptr);
}

void SLCS_hard_perform0(const DoubleMatrixPTR query_feat_ptr, const DoubleMatrixPTR reference_feat_ptr, const ParametersRLCSPTR parameters_slcs_ptr, const int start_q, const int end_q, const int start, const int end, const int padd_length, MemberPTR member_ptr) {
  
  // HELPER FUNCTION FOR 

  /*************************************
   * - if padd_length is zero, extract the reference region given by start and end from the reference.
   * - if padd_length is non-zero, extract the reference region (make a copy) and padd it with those many zero vectors of the same dimension.
   * - for each windowed region call RLCS
   *************************************/
  int q_size = end_q - start_q + 1;
  printf("query length: %d, roi region: [%d %d]\n", q_size, start, end);
  printf("roi region length: %d, padd length: %d\n", end-start+1, padd_length);


  DoubleMatrixPTR region_reference_ptr = create_DoubleMatrixPTR(); // feature array to store the reference region.

  

  //  if padd_length is zero, extrat the reference region given by start and end from the reference.
  int region_alloc_flag  =  0;
  if (!padd_length) {
    // just passing the reference instead of copying because it is fast.

   
    region_reference_ptr->mat = &(reference_feat_ptr->mat[start]);
    region_reference_ptr->dim = reference_feat_ptr->dim;
    region_reference_ptr->size = end-start+1;
    region_reference_ptr->max_size = end-start+1;
   
    region_alloc_flag = 0;
    // roi is not important for this purpose
  } else {
    // extract the region and make a copy of it then padd that many zeros.
    alloc_members_DoubleMatrixPTR(region_reference_ptr, end-start+1+padd_length, reference_feat_ptr->dim);
    
    // copy the region of reference  from start to end.
    copy_DoubleMatrixPTR_region(region_reference_ptr, reference_feat_ptr, start, end);
    
    // padd zeros to make it equal to motif length
    int i, j;
    for (i = end-start+1; i < region_reference_ptr->max_size; i++) {
      // add a zero vector of the same dimension
      for (j = 0; j < region_reference_ptr->dim; j++)
	region_reference_ptr->mat[i][j] = 0;
      
      region_reference_ptr->size += 1;
    }
    region_alloc_flag = 1;
  }

  DoubleMatrixPTR region_query_ptr = create_DoubleMatrixPTR();
  region_query_ptr->mat = &(query_feat_ptr->mat[start_q]);
  region_query_ptr->dim = query_feat_ptr->dim;
  region_query_ptr->size = end_q-start_q+1;
  region_query_ptr->max_size = end_q-start_q+1;


  OutputSLCSPTR slcs_out_ptr = create_OutputSLCSPTR();    // for storing the rlcs output matrices.
  alloc_members_OutputSLCSPTR(slcs_out_ptr, region_reference_ptr->size+1, region_query_ptr->size +1);
  
  //  alloc_members_MemberPTR(member_ptr, 1);

  Segment_LCS(region_query_ptr, region_reference_ptr, 0, region_reference_ptr->size-1,  parameters_slcs_ptr, slcs_out_ptr);

  find_best_sequence_slcs(slcs_out_ptr, parameters_slcs_ptr, start_q, start, padd_length, member_ptr);

  free_members_OutputSLCSPTR(slcs_out_ptr);
  destroy_OutputSLCSPTR(slcs_out_ptr);

  if (region_alloc_flag)
    free_members_DoubleMatrixPTR(region_reference_ptr);

  destroy_DoubleMatrixPTR(region_reference_ptr);
  destroy_DoubleMatrixPTR(region_query_ptr);
}


void Segment_LCS(const DoubleMatrixPTR query_feat_ptr, const DoubleMatrixPTR region_reference_ptr, const int start, const int end, const ParametersRLCSPTR parameters_slcs_ptr, OutputSLCSPTR slcs_out_ptr) {
  
  assert(end < region_reference_ptr->size);

  int i,j; // iterators
  int qsize = query_feat_ptr->size;

  // initilizing 1st col to zeros
  for (i = 0; i < slcs_out_ptr->row; i++) {
    slcs_out_ptr->cost_ptr->mat[i][0] = 0;
    slcs_out_ptr->cost_segment_ptr->mat[i][0] = 0;
    slcs_out_ptr->score_ptr->mat[i][0] = 0;
  }
  // initizing 1st row to zeros
  for (j = 0; j < slcs_out_ptr->col; j++) {
    slcs_out_ptr->cost_ptr->mat[0][j] = 0;
    slcs_out_ptr->cost_segment_ptr->mat[0][j] = 0;
    slcs_out_ptr->score_ptr->mat[0][j] = 0;
  }

  



  
  for (i = 1; i < slcs_out_ptr->row; i++) {

    for (j = 1; j < slcs_out_ptr->col; j++) {

      double dist = compute_distance_slcs(region_reference_ptr, query_feat_ptr, start, i-1, j-1, parameters_slcs_ptr);
     
      if (dist < parameters_slcs_ptr->Td) {

	slcs_out_ptr->diag_ptr->mat[i][j] = 1; // for cross

        slcs_out_ptr->cost_ptr->mat[i][j] = slcs_out_ptr->cost_ptr->mat[i-1][j-1] + (1 - dist/parameters_slcs_ptr->Td);

	slcs_out_ptr->cost_segment_ptr->mat[i][j] = slcs_out_ptr->cost_segment_ptr->mat[i-1][j-1] + (1 - dist/parameters_slcs_ptr->Td);
       
	slcs_out_ptr->score_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i-1][j-1];

      } else if (slcs_out_ptr->cost_ptr->mat[i-1][j] > slcs_out_ptr->cost_ptr->mat[i][j-1]) {
	
	slcs_out_ptr->diag_ptr->mat[i][j] = 2;  // for up

	slcs_out_ptr->cost_ptr->mat[i][j] = slcs_out_ptr->cost_ptr->mat[i-1][j];
	
	slcs_out_ptr->cost_segment_ptr->mat[i][j] = 0;

	if (slcs_out_ptr->diag_ptr->mat[i-1][j] == 1) 
	  slcs_out_ptr->score_ptr->mat[i][j] = (slcs_out_ptr->score_ptr->mat[i-1][j] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[i-1][j]))/SQR(qsize);
	else
	  slcs_out_ptr->score_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i-1][j];

      } else {
	slcs_out_ptr->diag_ptr->mat[i][j] = 3; // for left

	slcs_out_ptr->cost_ptr->mat[i][j] = slcs_out_ptr->cost_ptr->mat[i][j-1];
	slcs_out_ptr->cost_segment_ptr->mat[i][j] = 0;

	if (slcs_out_ptr->diag_ptr->mat[i][j-1] == 1) 
	  slcs_out_ptr->score_ptr->mat[i][j] = (slcs_out_ptr->score_ptr->mat[i][j-1] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[i][j-1]))/SQR(qsize);
	else
	  slcs_out_ptr->score_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i][j-1];
      }
    }  // end for    
  } // end for
  
  // compute score for last column and row
  int lri = slcs_out_ptr->row; // last row index of score matrix
  int lci = slcs_out_ptr->col;

  for (i = 0; i < lri; i++) {   
    // last column
    if (slcs_out_ptr->diag_ptr->mat[i][lci-1] == 1) 
      slcs_out_ptr->score_ptr->mat[i][lci] = (slcs_out_ptr->score_ptr->mat[i][lci-1] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[i][lci-1]))/SQR(qsize);
    else
      slcs_out_ptr->score_ptr->mat[i][lci] = slcs_out_ptr->score_ptr->mat[i][lci-1];
    
    slcs_out_ptr->diag_ptr->mat[i][lci] = 3; // left
  }

  for (i = 0; i < lci; i++) {
    // last row
    if (slcs_out_ptr->diag_ptr->mat[lri-1][i] == 1) 
      slcs_out_ptr->score_ptr->mat[lri][i] = (slcs_out_ptr->score_ptr->mat[lri-1][i] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[lri-1][i]))/SQR(qsize);
    else
      slcs_out_ptr->score_ptr->mat[lri][i] = slcs_out_ptr->score_ptr->mat[lri-1][i];

    slcs_out_ptr->diag_ptr->mat[lri][i] = 2; // up

  }

  // last cell // it will be computed properly in the upper loop
  slcs_out_ptr->score_ptr->mat[lri][lci] = slcs_out_ptr->score_ptr->mat[lri-1][lci];
  slcs_out_ptr->diag_ptr->mat[lri][lci] = 2; // up

} //end function



void Segment_LCS_v2(const DoubleMatrixPTR query_feat_ptr, const DoubleMatrixPTR region_reference_ptr, const int start, const int end, const ParametersRLCSPTR parameters_slcs_ptr, OutputSLCSPTR slcs_out_ptr) {
  double delta = 0.5;
  assert(end < region_reference_ptr->size);

  int i,j; // iterators
  int qsize = query_feat_ptr->size;

  // initilizing 1st col to zeros
  for (i = 0; i < slcs_out_ptr->row; i++) {
    slcs_out_ptr->cost_ptr->mat[i][0] = 0;
    slcs_out_ptr->cost_segment_ptr->mat[i][0] = 0;
    slcs_out_ptr->score_ptr->mat[i][0] = 0;
    slcs_out_ptr->adder_ptr->mat[i][0] = 0;
  }
  // initizing 1st row to zeros
  for (j = 0; j < slcs_out_ptr->col; j++) {
    slcs_out_ptr->cost_ptr->mat[0][j] = 0;
    slcs_out_ptr->cost_segment_ptr->mat[0][j] = 0;
    slcs_out_ptr->score_ptr->mat[0][j] = 0;
    slcs_out_ptr->adder_ptr->mat[0][j] = 0;
  }


  for (i = 1; i < slcs_out_ptr->row; i++) {

    for (j = 1; j < slcs_out_ptr->col; j++) {

      double dist = compute_distance_slcs(region_reference_ptr, query_feat_ptr, start, i-1, j-1, parameters_slcs_ptr);
     
      if (dist < parameters_slcs_ptr->Td) {

	slcs_out_ptr->diag_ptr->mat[i][j] = 1; // for cross

        slcs_out_ptr->cost_ptr->mat[i][j] = slcs_out_ptr->cost_ptr->mat[i-1][j-1] + (1 - dist/parameters_slcs_ptr->Td);

	slcs_out_ptr->cost_segment_ptr->mat[i][j] = slcs_out_ptr->cost_segment_ptr->mat[i-1][j-1] + (1 - dist/parameters_slcs_ptr->Td);
       

	// filling adder 
	double aa = slcs_out_ptr->adder_ptr->mat[i-1][j-1];
	double bb = slcs_out_ptr->adder_ptr->mat[i-1][j];
	double cc = slcs_out_ptr->adder_ptr->mat[i][j-1];

        double mabc = MYMAX(aa,bb,cc);
        if (mabc==-1)
	  slcs_out_ptr->adder_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i-1][j-1];
	else
	  slcs_out_ptr->adder_ptr->mat[i][j] = mabc;

	slcs_out_ptr->score_ptr->mat[i][j] = 0; // diag scores to zero

      } else if (slcs_out_ptr->cost_ptr->mat[i-1][j] > slcs_out_ptr->cost_ptr->mat[i][j-1]) {
	
	slcs_out_ptr->diag_ptr->mat[i][j] = 2;  // for up

	slcs_out_ptr->cost_ptr->mat[i][j] = slcs_out_ptr->cost_ptr->mat[i-1][j];
	
	slcs_out_ptr->cost_segment_ptr->mat[i][j] = slcs_out_ptr->cost_segment_ptr->mat[i-1][j] - delta;
	
	if (slcs_out_ptr->cost_segment_ptr->mat[i][j] < 0)
	  slcs_out_ptr->cost_segment_ptr->mat[i][j] = 0;

	// filling adder ptr
        if (slcs_out_ptr->cost_segment_ptr->mat[i][j] <= 0.5)
	  slcs_out_ptr->adder_ptr->mat[i][j] = -1;
	else {
	  double aa = slcs_out_ptr->adder_ptr->mat[i-1][j-1];
	  double bb = slcs_out_ptr->adder_ptr->mat[i-1][j];
	  double cc = slcs_out_ptr->adder_ptr->mat[i][j-1];
	  
	  double mabc = MYMAX(aa,bb,cc);
	  slcs_out_ptr->adder_ptr->mat[i][j] = mabc;
	}
	
	// filling score_ptr
	if (slcs_out_ptr->diag_ptr->mat[i-1][j] == 1) 
	  slcs_out_ptr->score_ptr->mat[i][j] = (slcs_out_ptr->adder_ptr->mat[i-1][j] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[i-1][j]))/SQR(qsize);
	else
	  slcs_out_ptr->score_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i-1][j];

      } else {
	
	slcs_out_ptr->diag_ptr->mat[i][j] = 3;  // for left

	slcs_out_ptr->cost_ptr->mat[i][j] = slcs_out_ptr->cost_ptr->mat[i][j-1];
	
	slcs_out_ptr->cost_segment_ptr->mat[i][j] = slcs_out_ptr->cost_segment_ptr->mat[i][j-1] - delta;
	
	if (slcs_out_ptr->cost_segment_ptr->mat[i][j] < 0)
	  slcs_out_ptr->cost_segment_ptr->mat[i][j] = 0;

	// filling adder ptr
        if (slcs_out_ptr->cost_segment_ptr->mat[i][j] <= 0.5)
	  slcs_out_ptr->adder_ptr->mat[i][j] = -1;
	else {
	  double aa = slcs_out_ptr->adder_ptr->mat[i-1][j-1];
	  double bb = slcs_out_ptr->adder_ptr->mat[i-1][j];
	  double cc = slcs_out_ptr->adder_ptr->mat[i][j-1];
	  
	  double mabc = MYMAX(aa,bb,cc);
	  slcs_out_ptr->adder_ptr->mat[i][j] = mabc;
	}
	
	// filling score_ptr
	if (slcs_out_ptr->diag_ptr->mat[i][j-1] == 1) 
	  slcs_out_ptr->score_ptr->mat[i][j] = (slcs_out_ptr->adder_ptr->mat[i][j-1] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[i][j-1]))/SQR(qsize);
	else
	  slcs_out_ptr->score_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i][j-1];
      }
    }  // end for    
  } // end for
  
  // compute score for last column and row
  int lri = slcs_out_ptr->row; // last row index of score matrix
  int lci = slcs_out_ptr->col;

  // fill the adder last row col with -1 /// no need to fill
  // also compute score at the same time 
  for (i = 0; i < lri; i++) {   
    // last column
    if (slcs_out_ptr->diag_ptr->mat[i][lci-1] == 1) 
      slcs_out_ptr->score_ptr->mat[i][lci] = (slcs_out_ptr->adder_ptr->mat[i][lci-1] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[i][lci-1]))/SQR(qsize);
    else
      slcs_out_ptr->score_ptr->mat[i][lci] = slcs_out_ptr->score_ptr->mat[i][lci-1];
    
    slcs_out_ptr->diag_ptr->mat[i][lci] = 3; // left
  }

  for (i = 0; i < lci; i++) {
    // last row
    if (slcs_out_ptr->diag_ptr->mat[lri-1][i] == 1) 
      slcs_out_ptr->score_ptr->mat[lri][i] = (slcs_out_ptr->adder_ptr->mat[lri-1][i] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[lri-1][i]))/SQR(qsize);
    else
      slcs_out_ptr->score_ptr->mat[lri][i] = slcs_out_ptr->score_ptr->mat[lri-1][i];

    slcs_out_ptr->diag_ptr->mat[lri][i] = 2; // up
  }

  // last cell // it will be computed properly in the upper loop
  slcs_out_ptr->score_ptr->mat[lri][lci] = slcs_out_ptr->score_ptr->mat[lri-1][lci];
  slcs_out_ptr->diag_ptr->mat[lri][lci] = 2; // up
} //end function



void Segment_LCS_v2_2(const DoubleMatrixPTR query_feat_ptr, const DoubleMatrixPTR region_reference_ptr, const int start, const int end, const ParametersRLCSPTR parameters_slcs_ptr, OutputSLCSPTR slcs_out_ptr) {
  
  printf("LCSS v2.2\n");

  double delta = 0.5;
  assert(end < region_reference_ptr->size);

  int i,j; // iterators
  int qsize = query_feat_ptr->size;

  // initilizing 1st col to zeros
  for (i = 0; i < slcs_out_ptr->row; i++) {
    //    slcs_out_ptr->cost_ptr->mat[i][0] = 0;
    slcs_out_ptr->cost_segment_ptr->mat[i][0] = 0;
    slcs_out_ptr->score_ptr->mat[i][0] = 0;
    slcs_out_ptr->adder_ptr->mat[i][0] = 0;
  }
  // initizing 1st row to zeros
  for (j = 0; j < slcs_out_ptr->col; j++) {
    //slcs_out_ptr->cost_ptr->mat[0][j] = 0;
    slcs_out_ptr->cost_segment_ptr->mat[0][j] = 0;
    slcs_out_ptr->score_ptr->mat[0][j] = 0;
    slcs_out_ptr->adder_ptr->mat[0][j] = 0;
  }
  
  //double madd =-1;

  for (i = 1; i < slcs_out_ptr->row; i++) {

    for (j = 1; j < slcs_out_ptr->col; j++) {

      double dist = compute_distance_slcs(region_reference_ptr, query_feat_ptr, start, i-1, j-1, parameters_slcs_ptr);
     
      if (dist < parameters_slcs_ptr->Td) {

	slcs_out_ptr->diag_ptr->mat[i][j] = 1; // for cross

	//  slcs_out_ptr->cost_ptr->mat[i][j] = slcs_out_ptr->cost_ptr->mat[i-1][j-1] + (1 - dist/parameters_slcs_ptr->Td);

	slcs_out_ptr->cost_segment_ptr->mat[i][j] = slcs_out_ptr->cost_segment_ptr->mat[i-1][j-1] + (1 - dist/parameters_slcs_ptr->Td);
       

	slcs_out_ptr->score_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i-1][j-1]; // diag scores to zero

      } else if (slcs_out_ptr->cost_segment_ptr->mat[i-1][j] > slcs_out_ptr->cost_segment_ptr->mat[i][j-1]) {
	
	slcs_out_ptr->diag_ptr->mat[i][j] = 2;  // for up

	//slcs_out_ptr->cost_ptr->mat[i][j] = slcs_out_ptr->cost_ptr->mat[i-1][j];
	
	slcs_out_ptr->cost_segment_ptr->mat[i][j] = slcs_out_ptr->cost_segment_ptr->mat[i-1][j] - delta;
	
	if (slcs_out_ptr->cost_segment_ptr->mat[i][j] < 0)
	  slcs_out_ptr->cost_segment_ptr->mat[i][j] = 0;

	
	// filling score_ptr
	if (slcs_out_ptr->diag_ptr->mat[i-1][j] == 1) 
	  slcs_out_ptr->score_ptr->mat[i][j] = (slcs_out_ptr->adder_ptr->mat[i-1][j] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[i-1][j]))/SQR(qsize);
	else
	  slcs_out_ptr->score_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i-1][j];

      } else {
	
	slcs_out_ptr->diag_ptr->mat[i][j] = 3;  // for left

	//slcs_out_ptr->cost_ptr->mat[i][j] = slcs_out_ptr->cost_ptr->mat[i][j-1];
	
	slcs_out_ptr->cost_segment_ptr->mat[i][j] = slcs_out_ptr->cost_segment_ptr->mat[i][j-1] - delta;
	
	if (slcs_out_ptr->cost_segment_ptr->mat[i][j] < 0)
	  slcs_out_ptr->cost_segment_ptr->mat[i][j] = 0;

	
	// filling score_ptr
	if (slcs_out_ptr->diag_ptr->mat[i][j-1] == 1) 
	  slcs_out_ptr->score_ptr->mat[i][j] = (slcs_out_ptr->adder_ptr->mat[i][j-1] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[i][j-1]))/SQR(qsize);
	else
	  slcs_out_ptr->score_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i][j-1];
      }
      

      double aa = slcs_out_ptr->adder_ptr->mat[i-1][j-1];
      double bb = slcs_out_ptr->adder_ptr->mat[i-1][j];
      double cc = slcs_out_ptr->adder_ptr->mat[i][j-1];
      
      double mabc = MYMAX(aa,bb,cc);
      // filling adder ptr
      if (mabc == -1 && slcs_out_ptr->diag_ptr->mat[i][j] == 1)
	slcs_out_ptr->adder_ptr->mat[i][j] = slcs_out_ptr->score_ptr->mat[i-1][j-1];
      else if (slcs_out_ptr->cost_segment_ptr->mat[i][j] <= 0.5 && slcs_out_ptr->diag_ptr->mat[i][j] != 1)
	slcs_out_ptr->adder_ptr->mat[i][j] = -1;
      else {
	slcs_out_ptr->adder_ptr->mat[i][j] = mabc;
      }
    }  // end for    
  } // end for
  
  // compute score for last column and row
  int lri = slcs_out_ptr->row; // last row index of score matrix
  int lci = slcs_out_ptr->col;

  // fill the adder last row col with -1 /// no need to fill
  // also compute score at the same time 
  for (i = 0; i < lri; i++) {   
    // last column
    if (slcs_out_ptr->diag_ptr->mat[i][lci-1] == 1) 
      slcs_out_ptr->score_ptr->mat[i][lci] = (slcs_out_ptr->adder_ptr->mat[i][lci-1] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[i][lci-1]))/SQR(qsize);
    else
      slcs_out_ptr->score_ptr->mat[i][lci] = slcs_out_ptr->score_ptr->mat[i][lci-1];
    
    slcs_out_ptr->diag_ptr->mat[i][lci] = 3; // left
  }

  for (i = 0; i < lci; i++) {
    // last row
    if (slcs_out_ptr->diag_ptr->mat[lri-1][i] == 1) 
      slcs_out_ptr->score_ptr->mat[lri][i] = (slcs_out_ptr->adder_ptr->mat[lri-1][i] * SQR(qsize) + SQR(slcs_out_ptr->cost_segment_ptr->mat[lri-1][i]))/SQR(qsize);
    else
      slcs_out_ptr->score_ptr->mat[lri][i] = slcs_out_ptr->score_ptr->mat[lri-1][i];

    slcs_out_ptr->diag_ptr->mat[lri][i] = 2; // up
  }

  // last cell // it will be computed properly in the upper loop
  slcs_out_ptr->score_ptr->mat[lri][lci] = slcs_out_ptr->score_ptr->mat[lri-1][lci];
  slcs_out_ptr->diag_ptr->mat[lri][lci] = 2; // up
} //end function




void find_best_sequence_slcs(const OutputSLCSPTR slcs_out_ptr, const ParametersRLCSPTR  parameters_slcs_ptr, const int query_index_correction,  const int reference_index_correction, const int padd_length, MemberPTR member_ptr) {


  IntVectorPTR max_score_index_ptr = create_IntVectorPTR();
  alloc_members_IntVectorPTR(max_score_index_ptr, 2);


  //  argmax_matrix_D_gen(rlcs_out_ptr->score_ptr, min_row, min_col, max_row, max_col, max_score_index_ptr);
  argmax_matrix_D_gen(slcs_out_ptr->score_ptr, 0, 0, slcs_out_ptr->row+1 - padd_length -1, slcs_out_ptr->col + 1 -1, max_score_index_ptr);
  
  assert(max_score_index_ptr->vec[0] < slcs_out_ptr->row+1 - padd_length);


  if (slcs_out_ptr->score_ptr->mat[max_score_index_ptr->vec[0]][max_score_index_ptr->vec[1]] <= parameters_slcs_ptr->seqFilterTd) {
    free_members_IntVectorPTR(max_score_index_ptr);
    destroy_IntVectorPTR(max_score_index_ptr);
    return;
  }

  int i;

  //  int x = max_score_index_ptr->vec[0] >= rlcs_out_ptr->row - padd_length ? rlcs_out_ptr->row - padd_length -1 : max_score_index_ptr->vec[0]; 
  int x = max_score_index_ptr->vec[0];
  int y = max_score_index_ptr->vec[1];

  free_members_IntVectorPTR(max_score_index_ptr);
  destroy_IntVectorPTR(max_score_index_ptr);


  // fill member

  member_ptr->score = slcs_out_ptr->score_ptr->mat[x][y];

  //  printf("%f\n", member_ptr->score); // comment it

  // check from here
  // also reduce x, y accordingly for proper color info using cost segment.
  if (slcs_out_ptr->diag_ptr->mat[x][y] ==1) {
    member_ptr->cost = slcs_out_ptr->cost_segment_ptr->mat[x-1][y-1];
    x -= 1; y-=1;
  }
  else if (slcs_out_ptr->diag_ptr->mat[x][y] ==2) {
    member_ptr->cost = slcs_out_ptr->cost_segment_ptr->mat[x-1][y];
    x-=1;
  }
  else { 
    member_ptr->cost = slcs_out_ptr->cost_segment_ptr->mat[x][y-1];
    y-=1;
  }
  member_ptr->cost_actual  = -1;

  member_ptr->war = -1;

  member_ptr->waq = -1;

  if (member_ptr->alignment_vec_ptr->max_size < slcs_out_ptr->row + slcs_out_ptr->col)
    realloc_members_AlignmentVectorPTR(member_ptr->alignment_vec_ptr, slcs_out_ptr->row + slcs_out_ptr->col);


  //  int seg_start = 0;
  while (x != 0 && y != 0) {
    i = member_ptr->alignment_vec_ptr->size;
    
    member_ptr->alignment_vec_ptr->vec[i].reference = x - 1 + reference_index_correction;
    
    member_ptr->alignment_vec_ptr->vec[i].query = y - 1 + query_index_correction;
    member_ptr->alignment_vec_ptr->vec[i].color = slcs_out_ptr->cost_segment_ptr->mat[x][y];

    if (slcs_out_ptr->diag_ptr->mat[x][y] == 1)
      x -= 1, y -= 1;
    else if (slcs_out_ptr->diag_ptr->mat[x][y] == 2)
      x -= 1;
    else if (slcs_out_ptr->diag_ptr->mat[x][y] == 3)
      y -= 1;

    member_ptr->alignment_vec_ptr->size += 1;
  }

  assert(member_ptr->alignment_vec_ptr->size <= member_ptr->alignment_vec_ptr->max_size);

  // reverse the alignment to make it correct;
  int align_size = member_ptr->alignment_vec_ptr->size;
  int align_middle_index = align_size % 2 ? floor((align_size-1)/2) - 1 : floor((align_size-1)/2);
  for (i = 0; i <= align_middle_index; i++) {
    // swapping two numbers
    swap_int(&(member_ptr->alignment_vec_ptr->vec[i].query), &(member_ptr->alignment_vec_ptr->vec[align_size-i-1].query));
    swap_int(&(member_ptr->alignment_vec_ptr->vec[i].reference), &(member_ptr->alignment_vec_ptr->vec[align_size-i-1].reference));
    swap_float(&(member_ptr->alignment_vec_ptr->vec[i].color), &(member_ptr->alignment_vec_ptr->vec[align_size-i-1].color));
  }

  // Store as sample number
  member_ptr->start = member_ptr->alignment_vec_ptr->vec[0].reference;
  member_ptr->end = member_ptr->alignment_vec_ptr->vec[align_size-1].reference;
  member_ptr->start_q = member_ptr->alignment_vec_ptr->vec[0].query;
  member_ptr->end_q = member_ptr->alignment_vec_ptr->vec[align_size - 1].query;
}

double DTWDistance(const DoubleMatrixPTR query_feat_ptr, const DoubleMatrixPTR region_reference_ptr, const ParametersRLCSPTR parameters_slcs_ptr) {

  
  int n = region_reference_ptr->size;
  int m = query_feat_ptr->size;
  
  double **DTW = (double**)malloc((n+1)*sizeof(double*));
  {
    int i;
    for (i=0; i<n+1; i++)
      DTW[i] = (double*) malloc((m+1)*sizeof(double));
  }


  int i,j;
  for (i = 1; i < n+1; i++)
    DTW[i][0] = 100000;
  for (j = 1; j < m+1; j++)
    DTW[0][j] = 100000;
  
  DTW[0][0] = 0;

  for (i=1; i < n+1; i++) {
    for (j=1; j < m+1; j++) {
      double cost = compute_distance_slcs(region_reference_ptr, query_feat_ptr, 0, i-1, j-1, parameters_slcs_ptr);
   
      double mymin = 10000;
      if (DTW[i-1][j]<DTW[i][j-1] && DTW[i-1][j] < DTW[i-1][j-1])
	mymin = DTW[i-1][j];
      else if (DTW[i][j-1]<DTW[i-1][j] && DTW[i][j-1] < DTW[i-1][j-1])
	mymin = DTW[i][j-1];
      else 
	mymin = DTW[i-1][j-1];

      DTW[i][j] = cost + mymin;
      //      printf("%f ", DTW[i][j]);
    }
    // printf("\n");
  }

  double ret =  1.0 - DTW[n][m]/(n+m);


  {
    int i;
    for (i=0; i<n+1; i++)
      free(DTW[i]);
  }
  free(DTW);

  return ret;
}



int members_comp(const void* a, const void* b) {
    
  if ( ((Member*)a)->start ==  ((Member*)b)->start ) return 0;
  else if (((Member*)a)->start >  ((Member*)b)->start) return 1;
  else return -1;

}

double compute_distance(const DoubleMatrixPTR region_reference_ptr, const DoubleMatrixPTR query_feat_ptr, const int startr,const int i, const int j, const ParametersRLCSPTR parameters_rlcs_ptr,  const SaddlePointsPTR query_sp_ptr, const SaddlePointsPTR region_reference_sp_ptr){
  // dimension should be same;
  assert(region_reference_ptr->dim == query_feat_ptr->dim);

  double ret = 0;

  if (!strcmp(parameters_rlcs_ptr-> distname, "cubic")) {
    int new_i = startr + i;
    double diff;
    int k;

    for (k = 0; k < region_reference_ptr->dim; k++) {
      diff = abs(region_reference_ptr->mat[new_i][k] - query_feat_ptr->mat[j][k]);
      if (diff < 3*parameters_rlcs_ptr->semitone)
	ret += pow(diff,3)/(pow(3*parameters_rlcs_ptr->semitone,3));
      else
	ret += 1;
    }
    
    ret /= region_reference_ptr->dim;
    
    return ret;
  } else if (!strcmp(parameters_rlcs_ptr-> distname, "dtw")) { // TODO: NOT WORKING 
    // everything is set start filling this
    
    // Start from correcting this line.
    // int s_r = arg_find_smallest_num_after_D(region_reference_sp_ptr->x_ptr, start + f_i-1);

    int rlcontext = parameters_rlcs_ptr->rlcontext;

    int sp_startr = arg_find_smallest_num_after_D(region_reference_sp_ptr->x_ptr, startr);
    int dtw_r_st = region_reference_sp_ptr->x_ptr->mat[(sp_startr + i)][0];
    int dtw_r_ed = region_reference_sp_ptr->x_ptr->mat[(sp_startr + i) + 2*rlcontext][0];

    //    int sp_startq = arg_find_smallest_num_after_D(region_reference_sp_ptr->x_ptr, 0); // no need for this it will be zero.
    int dtw_q_st = query_sp_ptr->x_ptr->mat[(0 + j)][0];
    int dtw_q_ed = query_sp_ptr->x_ptr->mat[(0 + j) + 2*rlcontext][0];

    return dtw_dist(region_reference_ptr, query_feat_ptr, dtw_r_st, dtw_r_ed, dtw_q_st, dtw_q_ed, parameters_rlcs_ptr->distname2, 1); // last arg for flag for normalizing.

  }
  return -1;
}



double compute_distance_slcs(const DoubleMatrixPTR region_reference_ptr, const DoubleMatrixPTR query_feat_ptr, const int startr,const int i, const int j, const ParametersRLCSPTR parameters_rlcs_ptr){
  // dimension should be same;
  assert(region_reference_ptr->dim == query_feat_ptr->dim);

  double ret = 0;

  if (!strcmp(parameters_rlcs_ptr-> distname, "cubic")) {
    int new_i = startr + i;
    double diff;
    int k;

    for (k = 0; k < region_reference_ptr->dim; k++) {
      diff = abs(region_reference_ptr->mat[new_i][k] - query_feat_ptr->mat[j][k]);
      if (diff < 3*parameters_rlcs_ptr->semitone)
	ret += pow(diff,3)/(pow(3*parameters_rlcs_ptr->semitone,3));
      else
	ret += 1;
    }
    
    ret /= region_reference_ptr->dim;
    
    return ret;
  }
  return -1;
}




double dtw_dist(const DoubleMatrixPTR region_reference_ptr, const DoubleMatrixPTR query_feat_ptr, const int dtw_r_st, const int dtw_r_ed, const int dtw_q_st, const int dtw_q_ed, const char * distname, const int isNorm) {
  //TODO: complete the function using (intergrating )sankalp's code

  return -1;
}

void correct_SpotVec_for_dtw_dist(SpotVectorPTR spot_vec_ptr, SaddlePointsPTR reference_sp_ptr, SaddlePointsPTR query_sp_ptr) {
  // read for its use where this function is called in RLCS function.
  // This function needs to be reviewed. There could be mistakes.
  int i;

  for (i = 0; i < spot_vec_ptr->size; i++) {
    spot_vec_ptr -> vec[i].group_start = reference_sp_ptr -> x_ptr->mat[spot_vec_ptr ->vec[i].group_start][0];
    spot_vec_ptr ->vec[i].group_end = reference_sp_ptr -> x_ptr->mat[spot_vec_ptr ->vec[i].group_end][0];
    spot_vec_ptr ->vec[i].group_start_q = query_sp_ptr -> x_ptr->mat[spot_vec_ptr ->vec[i].group_start_q][0];
    spot_vec_ptr ->vec[i].group_end_q = query_sp_ptr -> x_ptr->mat[spot_vec_ptr ->vec[i].group_end_q][0];

    int j;
    for (j = 0; j < spot_vec_ptr ->vec[i].member_vec_ptr->size; j++) {
      spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].start = reference_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].start][0];
      spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].end = reference_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].end][0];
      spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].start_q = query_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].start_q][0];
      spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].end_q = query_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].end_q][0];
      
      int k;
      for (k = 0; k < spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].alignment_vec_ptr->size; k++) {
	spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].alignment_vec_ptr->vec[k].query = query_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].alignment_vec_ptr->vec[k].query][0];
	spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].alignment_vec_ptr->vec[k].reference = reference_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].member_vec_ptr->vec[j].alignment_vec_ptr->vec[k].reference][0];
      }
    }

    spot_vec_ptr -> vec[i].best_member_ptr->start = reference_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].best_member_ptr->start][0];
    spot_vec_ptr -> vec[i].best_member_ptr->end = reference_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].best_member_ptr->end][0];
    spot_vec_ptr -> vec[i].best_member_ptr->start_q = query_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].best_member_ptr->start_q][0];
    spot_vec_ptr -> vec[i].best_member_ptr->end_q = query_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].best_member_ptr->end_q][0];
      
      int k;
      for (k = 0; k < spot_vec_ptr -> vec[i].best_member_ptr->alignment_vec_ptr->size; k++) {
	spot_vec_ptr -> vec[i].best_member_ptr->alignment_vec_ptr->vec[k].query = query_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].best_member_ptr->alignment_vec_ptr->vec[k].query][0];
	spot_vec_ptr -> vec[i].best_member_ptr->alignment_vec_ptr->vec[k].reference = reference_sp_ptr -> x_ptr->mat[spot_vec_ptr -> vec[i].best_member_ptr->alignment_vec_ptr->vec[k].reference][0];
      }
  }

}


double mean_of_first_diff(DoubleMatrixPTR mat_ptr, int start, int end) {
  int i; 
  double mean = 0.0;
  for (i = start; i <= end-1; i++) {
    mean = mean + (mat_ptr->mat[i+1][0] - mat_ptr->mat[i][0]);
  }

  mean  =  mean/((end-1) - start +1);

  return mean;
}

double std_of_first_diff(DoubleMatrixPTR mat_ptr, int start, int end) {
  int i; 
  double mean = mean_of_first_diff(mat_ptr, start, end);
  double std = 0;
  for (i = start; i <= end-1; i++) {
    double temp = (mat_ptr->mat[i+1][0] - mat_ptr->mat[i][0]) - mean;
    std = std + (temp*temp);
  }

  std  =  std/((end-1) - start +1);
  
  return sqrt(std);
}


void LinReg (const DoubleMatrixPTR mat_ptr, const int start, const int numPts,double  *intercept, double *slope) {
  double               sumX = 0.0, sumY = 0.0, sqSumX = 0.0, sqSumY = 0.0, sqSumXY = 0.0;
  double              m, c;  // y = mx + c
  double              norm;
  int                 i, N;
  N = start+numPts < mat_ptr->size ? numPts : mat_ptr->size - start;
  for (i = start; i < (start+numPts) && i < mat_ptr->size; i++) {
    sumY += mat_ptr->mat[i][0];
    sqSumY += (mat_ptr->mat[i][0] * mat_ptr->mat[i][0]);
    sqSumXY += (mat_ptr->mat[i][0] * (i+1));
    sumX += (i+1);
    sqSumX += ((i+1)*(i+1));
  }
  norm = (N * sqSumX) - (sumX * sumX);
  if (norm != 0) {
    c = ((sqSumX * sumY) - (sumX * sqSumXY))/norm;
    m = ((N * sqSumXY) - (sumX * sumY))/norm;
  } else {
    m = 0;
    c = 0;
  }
  
  *intercept = c;
  *slope = m;
}

double slope_of_best_linear_fit(const DoubleMatrixPTR mat_ptr, const int start, const int end) { // used as slope of linear trend.
  double slope, intercept;
  LinReg(mat_ptr, start, end-start+1, &intercept, &slope);
  return slope;  
}

void avg_error(const DoubleMatrixPTR mat_ptr, const int start, const int numPts, const double intercept, const double slope, double *error) {
  double sqerror_sum=0.0;
  int i;
  for (i = start; i < start+numPts; i++) 
    sqerror_sum += ((slope*i + intercept) - mat_ptr->mat[i][0]) * ((slope*i + intercept) - mat_ptr->mat[i][0]);
  
  sqerror_sum /= numPts;
  *error = sqrt(sqerror_sum);
}

double std_of_slope(const DoubleMatrixPTR mat_ptr, const int start, const int end) {  // used as std of slope of linear trend.
  double slope, intercept, error;
  LinReg(mat_ptr, start, end-start+1, &intercept, &slope);
  avg_error(mat_ptr, start, end-start+1, intercept, slope, &error);
  return error;
}



double min(double x, double y)  {
  return (x<=y ? x : y);
}

double  max(double x, double y)  {
  return x>=y ? x : y;
}
