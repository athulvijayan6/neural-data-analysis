/* Description: 
 * Creation Date: 9 March 2014
 * Author: Shrey Dutta
 */


#include "types.h"


// Never use datastructer as it is. Always use the ____PTR version to make less errors.
// - whenever you create something always destroy it when use it over.
// - whenever you allocate something always free it when use it over.
// Double pointer is ___PPTR


//TODO: write printf functions for all the sructures.

RangePTR create_RangePTR() {return (RangePTR) malloc(sizeof(Range));}
PointPTR create_PointPTR() {return (PointPTR) malloc(sizeof(PointPTR));}
AlignmentPTR create_AlignmentPTR() {return (AlignmentPTR) malloc(sizeof(AlignmentPTR));}
OutputRLCSPTR create_OutputRLCSPTR() {return (OutputRLCSPTR) malloc(sizeof(OutputRLCS));}
OutputSLCSPTR create_OutputSLCSPTR() {return (OutputSLCSPTR) malloc(sizeof(OutputSLCS));}
ParametersRLCSPTR create_ParametersRLCSPTR() {return (ParametersRLCSPTR) malloc(sizeof(ParametersRLCS));}
IntVectorPTR create_IntVectorPTR() {return (IntVectorPTR) malloc(sizeof(IntVector));}
DoubleVectorPTR create_DoubleVectorPTR() {return (DoubleVectorPTR) malloc(sizeof(DoubleVector));}
StringVectorPTR create_StringVectorPTR() {return (StringVectorPTR) malloc(sizeof(StringVector));}
RangeVectorPTR create_RangeVectorPTR() {return (RangeVectorPTR) malloc(sizeof(RangeVector));}
PointVectorPTR create_PointVectorPTR() {return (PointVectorPTR) malloc(sizeof(PointVector));}
AlignmentVectorPTR create_AlignmentVectorPTR() {return (AlignmentVectorPTR) malloc(sizeof(AlignmentVector));}
IntMatrixPTR create_IntMatrixPTR() {return (IntMatrixPTR) malloc(sizeof(IntMatrix));}
DoubleMatrixPTR create_DoubleMatrixPTR() {return (DoubleMatrixPTR) malloc(sizeof(DoubleMatrix));}
CharMatrixPTR create_CharMatrixPTR() {return (CharMatrixPTR) malloc(sizeof(CharMatrix));}
StringMatrixPTR create_StringMatrixPTR() {return (StringMatrixPTR) malloc(sizeof(StringMatrix));}
SaddlePointsPTR create_SaddlePointsPTR() {return (SaddlePointsPTR) malloc(sizeof(SaddlePoints));}
SegmentPTR create_SegmentPTR() {return (SegmentPTR) malloc(sizeof(Segment));}
MemberPTR create_MemberPTR() {return (MemberPTR) malloc(sizeof(Member));} 
MemberVectorPTR create_MemberVectorPTR() {return (MemberVectorPTR) malloc(sizeof(MemberVector));}
SpotPTR create_SpotPTR() {return (SpotPTR) malloc(sizeof(Spot));}
SpotVectorPTR create_SpotVectorPTR() {return (SpotVectorPTR) malloc(sizeof(SpotVector));}


void destroy_RangePTR(RangePTR temp) {free(temp);}
void destroy_PointPTR(PointPTR temp) {free(temp);}
void destroy_AlignmentPTR(AlignmentPTR temp) {free(temp);}
void destroy_OutputRLCSPTR(OutputRLCSPTR temp) {free(temp);}
void destroy_OutputSLCSPTR(OutputSLCSPTR temp) {free(temp);}
void destroy_ParametersRLCSPTR(ParametersRLCSPTR temp) {free(temp);}
void destroy_IntVectorPTR(IntVectorPTR temp) {free(temp);}
void destroy_DoubleVectorPTR(DoubleVectorPTR temp) {free(temp);}
void destroy_StringVectorPTR(StringVectorPTR temp) {free(temp);}
void destroy_RangeVectorPTR(RangeVectorPTR temp) {free(temp);}
void destroy_PointVectorPTR(PointVectorPTR temp) {free(temp);}
void destroy_AlignmentVectorPTR(AlignmentVectorPTR temp) {free(temp);}
void destroy_IntMatrixPTR(IntMatrixPTR temp) {free(temp);}
void destroy_DoubleMatrixPTR(DoubleMatrixPTR temp) {free(temp);}
void destroy_CharMatrixPTR(CharMatrixPTR temp) {free(temp);}
void destroy_StringMatrixPTR(StringMatrixPTR temp) {free(temp);}
void destroy_SaddlePointsPTR(SaddlePointsPTR temp) {free(temp);}
void destroy_SegmentPTR(SegmentPTR temp) {free(temp);}
void destroy_MemberPTR(MemberPTR temp) {free(temp);}
void destroy_MemberVectorPTR(MemberVectorPTR temp) {free(temp);}
void destroy_SpotPTR(SpotPTR temp) {free(temp);}
void destroy_SpotVectorPTR(SpotVectorPTR temp) {free(temp);}




void alloc_members_OutputRLCSPTR(OutputRLCSPTR rlcs_out_ptr, const int row, const int col) {
  rlcs_out_ptr->row = row;
  rlcs_out_ptr->col = col;
  
  rlcs_out_ptr->score_ptr = create_DoubleMatrixPTR();
  rlcs_out_ptr->cost_ptr = create_DoubleMatrixPTR();
  rlcs_out_ptr->cost_actual_ptr = create_DoubleMatrixPTR();
  rlcs_out_ptr->war_ptr = create_DoubleMatrixPTR();
  rlcs_out_ptr->waq_ptr = create_DoubleMatrixPTR();
  rlcs_out_ptr->diag_ptr = create_DoubleMatrixPTR();

  alloc_members_DoubleMatrixPTR(rlcs_out_ptr->score_ptr, row, col);
  alloc_members_DoubleMatrixPTR(rlcs_out_ptr->cost_ptr, row, col);
  alloc_members_DoubleMatrixPTR(rlcs_out_ptr->cost_actual_ptr, row, col);
  alloc_members_DoubleMatrixPTR(rlcs_out_ptr->war_ptr, row, col);
  alloc_members_DoubleMatrixPTR(rlcs_out_ptr->waq_ptr, row, col);
  alloc_members_DoubleMatrixPTR(rlcs_out_ptr->diag_ptr, row, col);

  rlcs_out_ptr->score_ptr->size = row; 
  rlcs_out_ptr->cost_ptr->size = row; 
  rlcs_out_ptr->cost_actual_ptr->size = row; 
  rlcs_out_ptr->war_ptr->size = row; 
  rlcs_out_ptr->waq_ptr->size = row; 
  rlcs_out_ptr->diag_ptr->size = row; 
}



void alloc_members_OutputSLCSPTR(OutputSLCSPTR slcs_out_ptr, const int row, const int col) {
  slcs_out_ptr->row = row;
  slcs_out_ptr->col = col;
  
  slcs_out_ptr->score_ptr = create_DoubleMatrixPTR();
  slcs_out_ptr->cost_ptr = create_DoubleMatrixPTR();
  slcs_out_ptr->cost_segment_ptr = create_DoubleMatrixPTR();
  slcs_out_ptr->diag_ptr = create_DoubleMatrixPTR();
  slcs_out_ptr->adder_ptr = create_DoubleMatrixPTR();

  alloc_members_DoubleMatrixPTR(slcs_out_ptr->score_ptr, row+1, col+1);
  alloc_members_DoubleMatrixPTR(slcs_out_ptr->cost_ptr, row, col);
  alloc_members_DoubleMatrixPTR(slcs_out_ptr->cost_segment_ptr, row, col);
  alloc_members_DoubleMatrixPTR(slcs_out_ptr->diag_ptr, row+1, col+1);
  alloc_members_DoubleMatrixPTR(slcs_out_ptr->adder_ptr, row, col);

  slcs_out_ptr->score_ptr->size = row+1; 
  slcs_out_ptr->cost_ptr->size = row; 
  slcs_out_ptr->cost_segment_ptr->size = row; 
  slcs_out_ptr->diag_ptr->size = row+1; 
  slcs_out_ptr->adder_ptr->size = row; 
}



void alloc_members_MemberPTR(MemberPTR member_ptr, int align_size) {
  member_ptr->alignment_vec_ptr = create_AlignmentVectorPTR();
  alloc_members_AlignmentVectorPTR(member_ptr->alignment_vec_ptr, align_size);
}


void alloc_members_SpotPTR(SpotPTR spot_ptr, const int members_num) {
  spot_ptr->reference_path = (char*) malloc(sizeof(char) * 500);
  spot_ptr->query_path = (char*) malloc(sizeof(char) * 500);

  spot_ptr->member_vec_ptr = create_MemberVectorPTR();
  alloc_members_MemberVectorPTR(spot_ptr->member_vec_ptr, members_num);

  spot_ptr->best_member_ptr = create_MemberPTR();
  alloc_members_MemberPTR(spot_ptr->best_member_ptr, 1);
}


void alloc_members_IntVectorPTR(IntVectorPTR vec_ptr, const int max_size) {
  vec_ptr->vec = (int*) malloc(sizeof(int)*max_size);
  vec_ptr->size = 0;
  vec_ptr->max_size = max_size;
}

void alloc_members_DoubleVectorPTR(DoubleVectorPTR vec_ptr, const int max_size) {
  vec_ptr->vec = (double*) malloc(sizeof(double)*max_size);
  vec_ptr->size = 0;
  vec_ptr->max_size = max_size;
}

void alloc_members_StringVectorPTR(StringVectorPTR vec_ptr, const int max_size) {
  vec_ptr->vec = (String*) malloc(sizeof(String)*max_size);
  vec_ptr->size = 0;
  vec_ptr->max_size = max_size;
}


void alloc_members_RangeVectorPTR(RangeVectorPTR vec_ptr, const int max_size) {
  vec_ptr->vec = (Range*) malloc(sizeof(Range)*max_size);
  vec_ptr->size = 0;
  vec_ptr->max_size = max_size;
}

void alloc_members_PointVectorPTR(PointVectorPTR vec_ptr, const int max_size) {
  vec_ptr->vec = (Point*) malloc(sizeof(Point)*max_size);
  vec_ptr->size = 0;
  vec_ptr->max_size = max_size;
}

void  alloc_members_AlignmentVectorPTR(AlignmentVectorPTR vec_ptr, const int max_size) {
  vec_ptr->vec = (Alignment*) malloc(sizeof(Alignment)*max_size);
  vec_ptr->size = 0;
  vec_ptr->max_size = max_size;
}

// Matrices
void  alloc_members_IntMatrixPTR(IntMatrixPTR mat_ptr, const int max_size, const int dim) {
  mat_ptr->mat = (int**) malloc(sizeof(int*) * max_size);
  
  int i;
  for (i = 0; i< max_size; i++) 
    mat_ptr->mat[i] = (int*) malloc(sizeof(int) * dim);

  mat_ptr->dim = dim;
  mat_ptr->size = 0;
  mat_ptr->max_size = max_size;
}

void alloc_members_DoubleMatrixPTR(DoubleMatrixPTR mat_ptr, const int max_size, const int dim) {
  mat_ptr->mat = (double**) malloc(sizeof(double*) * max_size);

  int i;  

  for (i = 0; i< max_size; i++) 
    mat_ptr->mat[i] = (double*) malloc(sizeof(double) * dim);

  mat_ptr->dim = dim;
  mat_ptr->size = 0;
  mat_ptr->max_size = max_size;
}


void alloc_members_CharMatrixPTR(CharMatrixPTR mat_ptr, const int max_size, const int dim) {
  mat_ptr->mat = (char**) malloc(sizeof(char*) * max_size);

  int i;  

  for (i = 0; i< max_size; i++) 
    mat_ptr->mat[i] = (char*) malloc(sizeof(char) * dim);

  mat_ptr->dim = dim;
  mat_ptr->size = 0;
  mat_ptr->max_size = max_size;
}


void alloc_members_StringMatrixPTR(StringMatrixPTR mat_ptr, const int max_size, const int dim) {
  mat_ptr->mat = (String**) malloc(sizeof(String*) * max_size);

  int i;  

  for (i = 0; i< max_size; i++) 
    mat_ptr->mat[i] = (String*) malloc(sizeof(String) * dim);

  mat_ptr->dim = dim;
  mat_ptr->size = 0;
  mat_ptr->max_size = max_size;
}


// Special structs 
void  alloc_members_MemberVectorPTR(MemberVectorPTR vec_ptr, const int max_size) {
  vec_ptr->vec = (Member*) malloc(sizeof(Member) * max_size);
  int i;
  for (i = 0; i < max_size; i++) {
    alloc_members_MemberPTR(&(vec_ptr->vec[i]), 1);
  }
  
  vec_ptr->size = 0;
  vec_ptr->max_size = max_size;
}

void alloc_members_SpotVectorPTR(SpotVectorPTR  vec_ptr, const int max_size) {
  vec_ptr->vec = (Spot*) malloc(sizeof(Spot)*max_size);

  int i;
  for (i = 0; i < max_size; i++) {
    alloc_members_SpotPTR(&(vec_ptr->vec[i]), 1);
  }
  vec_ptr->size = 0;
  vec_ptr->max_size = max_size;
}

void alloc_members_SaddlePointsPTR(SaddlePointsPTR saddle_points_ptr, int max_size) {
  saddle_points_ptr->size = max_size;
  
  saddle_points_ptr->x_ptr = create_DoubleMatrixPTR();
  saddle_points_ptr->y_ptr = create_DoubleMatrixPTR();

  alloc_members_DoubleMatrixPTR(saddle_points_ptr->x_ptr, max_size, 1);
  alloc_members_DoubleMatrixPTR(saddle_points_ptr->y_ptr, max_size, 1);

  saddle_points_ptr->x_ptr->size = max_size;
  saddle_points_ptr->y_ptr->size = max_size;
}



// rellocating structs
void realloc_members_IntVectorPTR(IntVectorPTR  vec_ptr, const int max_size) {
  assert(max_size > vec_ptr->max_size); 

  vec_ptr->vec = (int*) realloc(vec_ptr->vec, sizeof(int)*max_size);
  vec_ptr->max_size = max_size;
  vec_ptr->size = vec_ptr->size;
}

void realloc_members_DoubleVectorPTR(DoubleVectorPTR vec_ptr, const int max_size) {

  assert(max_size > vec_ptr->max_size); 

  vec_ptr->vec = (double*) realloc(vec_ptr->vec, sizeof(double)*max_size);
  vec_ptr->max_size = max_size;
  vec_ptr->size = vec_ptr->size;
}


void realloc_members_StringVectorPTR(StringVectorPTR vec_ptr, const int max_size) {

  assert(max_size > vec_ptr->max_size); 

  vec_ptr->vec = (String*) realloc(vec_ptr->vec, sizeof(String)*max_size);
  vec_ptr->max_size = max_size;
  vec_ptr->size = vec_ptr->size;
}


void realloc_members_RangeVectorPTR(RangeVectorPTR  vec_ptr, const int max_size) {
 
  assert(max_size > vec_ptr->max_size); 

  vec_ptr->vec = (Range*) realloc(vec_ptr->vec, sizeof(Range)*max_size);
  vec_ptr->max_size = max_size;
  vec_ptr->size = vec_ptr->size;
}

void realloc_members_AlignmentVectorPTR(AlignmentVectorPTR  vec_ptr, const int max_size) {

  assert(max_size > vec_ptr->max_size); 

  vec_ptr->vec = (Alignment*) realloc(vec_ptr->vec, sizeof(Alignment)*max_size);
  vec_ptr->max_size = max_size;
  vec_ptr->size = vec_ptr->size;
  
}

void  realloc_members_IntMatrixPTR(IntMatrixPTR  mat_ptr, const int max_size) {

  assert(max_size > mat_ptr->max_size);

  mat_ptr->mat = (int**) realloc(mat_ptr->mat, sizeof(int*) * max_size);
  
  int i;
  for (i = mat_ptr->max_size; i < max_size; i++) 
    mat_ptr->mat[i] = (int*) malloc(sizeof(int) * mat_ptr->dim);

  mat_ptr->max_size = max_size;
  mat_ptr->size = mat_ptr->size;
  mat_ptr->dim = mat_ptr->dim;

}

void realloc_members_DoubleMatrixPTR(DoubleMatrixPTR mat_ptr, const int max_size) {

  assert(max_size > mat_ptr->max_size); 

  mat_ptr->mat = (double**) realloc(mat_ptr->mat, sizeof(double*) * max_size);
  
  int i;
  for (i = mat_ptr->max_size; i< max_size; i++) 
    mat_ptr->mat[i] = (double*) malloc(sizeof(double) * mat_ptr->dim);

  mat_ptr->max_size = max_size;
  mat_ptr->size = mat_ptr->size;
  mat_ptr->dim = mat_ptr->dim;
}

void realloc_members_CharMatrixPTR(CharMatrixPTR mat_ptr, const int max_size) {

  assert(max_size > mat_ptr->max_size); 

  mat_ptr->mat = (char**) realloc(mat_ptr->mat, sizeof(char*) * max_size);
  
  int i;
  for (i = mat_ptr->max_size; i< max_size; i++) 
    mat_ptr->mat[i] = (char*) malloc(sizeof(char) * mat_ptr->dim);

  mat_ptr->max_size = max_size;
  mat_ptr->size = mat_ptr->size;
  mat_ptr->dim = mat_ptr->dim;
}


void realloc_members_StringMatrixPTR(StringMatrixPTR mat_ptr, const int max_size) {

  assert(max_size > mat_ptr->max_size); 

  mat_ptr->mat = (String**) realloc(mat_ptr->mat, sizeof(String*) * max_size);
  
  int i;
  for (i = mat_ptr->max_size; i< max_size; i++) 
    mat_ptr->mat[i] = (String*) malloc(sizeof(String) * mat_ptr->dim);

  mat_ptr->max_size = max_size;
  mat_ptr->size = mat_ptr->size;
  mat_ptr->dim = mat_ptr->dim;
}



void realloc_members_MemberVectorPTR(MemberVectorPTR  vec_ptr, const int max_size) {

  assert(max_size > vec_ptr->max_size); 

  vec_ptr->vec = (Member*) realloc(vec_ptr->vec, sizeof(Member) * max_size);

  int i;
  for (i = vec_ptr->max_size; i <  max_size; i++) {
    alloc_members_MemberPTR(&(vec_ptr->vec[i]), 1);
  }
  
  vec_ptr->max_size = max_size;
  vec_ptr->size = vec_ptr->size;
}


void realloc_members_SpotVectorPTR(SpotVectorPTR  vec_ptr, const int max_size) {

  assert(max_size > vec_ptr->max_size); 

  vec_ptr->vec = (Spot*) realloc(vec_ptr->vec, sizeof(Spot)*max_size);
  int i;

  for(i = vec_ptr->max_size; i < max_size; i++) {
    alloc_members_SpotPTR(&(vec_ptr->vec[i]), 1);
  }

  vec_ptr->max_size = max_size;
  vec_ptr->size = vec_ptr->size;
}


//free structs

void free_members_MemberPTR(MemberPTR member_ptr) {
  free_members_AlignmentVectorPTR(member_ptr->alignment_vec_ptr);
  destroy_AlignmentVectorPTR(member_ptr->alignment_vec_ptr);
}

void free_members_SpotPTR(SpotPTR spot_ptr) {
  
  //  printf("1\n");
  free_members_MemberVectorPTR(spot_ptr->member_vec_ptr);
  //printf("2\n");
  destroy_MemberVectorPTR(spot_ptr->member_vec_ptr);
  //printf("3\n");
  free_members_MemberPTR(spot_ptr->best_member_ptr);
  //printf("4\n");
  destroy_MemberPTR(spot_ptr->best_member_ptr);
  //printf("5\n");
  free(spot_ptr->query_path);
  //printf("6\n");
  free(spot_ptr->reference_path);
  //printf("7\n");
}
 
void free_members_IntVectorPTR(IntVectorPTR  vec_ptr) {
  free(vec_ptr->vec);
  vec_ptr->size = 0;
  vec_ptr->max_size = 0;
}

void free_members_DoubleVectorPTR(DoubleVectorPTR  vec_ptr) {
  free(vec_ptr->vec);
  vec_ptr->size = 0;
  vec_ptr->max_size = 0;
}


void free_members_StringVectorPTR(StringVectorPTR  vec_ptr) {
  free(vec_ptr->vec);
  vec_ptr->size = 0;
  vec_ptr->max_size = 0;
}


void free_members_RangeVectorPTR(RangeVectorPTR  vec_ptr) {
  free(vec_ptr->vec);
  vec_ptr->size = 0;
  vec_ptr->max_size = 0;
}

void free_members_PointVectorPTR(PointVectorPTR  vec_ptr) {
  free(vec_ptr->vec);
  vec_ptr->size = 0;
  vec_ptr->max_size = 0;
}

void free_members_AlignmentVectorPTR(AlignmentVectorPTR  vec_ptr) {
  free(vec_ptr->vec);
  vec_ptr->size = 0;
  vec_ptr->max_size = 0;
}

// Matrices
void free_members_IntMatrixPTR(IntMatrixPTR  mat_ptr) {
  
  int i;
  for (i = 0; i< mat_ptr->max_size; i++) 
    free(mat_ptr->mat[i]);
  free(mat_ptr->mat);

  mat_ptr->dim = 0;
  mat_ptr->size = 0;
  mat_ptr->max_size = 0;
}

void free_members_DoubleMatrixPTR(DoubleMatrixPTR  mat_ptr) {
  
  int i;
  for (i = 0; i< mat_ptr->max_size; i++)
    free(mat_ptr->mat[i]);
   
  free(mat_ptr->mat);

  mat_ptr->dim = 0;
  mat_ptr->size = 0;
  mat_ptr->max_size = 0;
}


void free_members_CharMatrixPTR(CharMatrixPTR  mat_ptr) {
  
  int i;
  for (i = 0; i< mat_ptr->max_size; i++)
    free(mat_ptr->mat[i]);

  free(mat_ptr->mat);

  mat_ptr->dim = 0;
  mat_ptr->size = 0;
  mat_ptr->max_size = 0;
}


void free_members_StringMatrixPTR(StringMatrixPTR  mat_ptr) {
  
  int i;
  for (i = 0; i< mat_ptr->max_size; i++)
    free(mat_ptr->mat[i]);

  free(mat_ptr->mat);

  mat_ptr->dim = 0;
  mat_ptr->size = 0;
  mat_ptr->max_size = 0;
}



void free_members_MemberVectorPTR(MemberVectorPTR  vec_ptr) {
  int i = 0;
  
  for(i = 0; i < vec_ptr->max_size;  i++)
    free_members_MemberPTR(&(vec_ptr->vec[i]));

  free(vec_ptr->vec);
  vec_ptr->size = 0;
  vec_ptr->max_size = 0;
}

void free_members_SpotVectorPTR(SpotVectorPTR  vec_ptr) {
  int i = 0;

  for (i = 0; i < vec_ptr->max_size; i++)
    free_members_SpotPTR(&(vec_ptr->vec[i]));

  free(vec_ptr->vec);
  vec_ptr->size = 0;
  vec_ptr->max_size = 0;
}

void free_members_OutputRLCSPTR(OutputRLCSPTR rlcs_out_ptr) {

  rlcs_out_ptr->score_ptr->size = 0; 
  rlcs_out_ptr->cost_ptr->size = 0; 
  rlcs_out_ptr->cost_actual_ptr->size = 0; 
  rlcs_out_ptr->war_ptr->size = 0; 
  rlcs_out_ptr->waq_ptr->size = 0; 
  rlcs_out_ptr->diag_ptr->size = 0; 
 
  free_members_DoubleMatrixPTR(rlcs_out_ptr->diag_ptr);
  free_members_DoubleMatrixPTR(rlcs_out_ptr->score_ptr);
  free_members_DoubleMatrixPTR(rlcs_out_ptr->cost_ptr);
  free_members_DoubleMatrixPTR(rlcs_out_ptr->cost_actual_ptr);
  free_members_DoubleMatrixPTR(rlcs_out_ptr->war_ptr);
  free_members_DoubleMatrixPTR(rlcs_out_ptr->waq_ptr);

  destroy_DoubleMatrixPTR(rlcs_out_ptr->score_ptr);
  destroy_DoubleMatrixPTR(rlcs_out_ptr->cost_ptr);
  destroy_DoubleMatrixPTR(rlcs_out_ptr->cost_actual_ptr);
  destroy_DoubleMatrixPTR(rlcs_out_ptr->war_ptr);
  destroy_DoubleMatrixPTR(rlcs_out_ptr->waq_ptr);
  destroy_DoubleMatrixPTR(rlcs_out_ptr->diag_ptr);

  rlcs_out_ptr->row = 0;
  rlcs_out_ptr->col = 0;
}


void free_members_OutputSLCSPTR(OutputSLCSPTR slcs_out_ptr) {

  slcs_out_ptr->score_ptr->size = 0; 
  slcs_out_ptr->cost_ptr->size = 0; 
  slcs_out_ptr->cost_segment_ptr->size = 0; 
  slcs_out_ptr->diag_ptr->size = 0; 
  slcs_out_ptr->adder_ptr->size = 0; 


  free_members_DoubleMatrixPTR(slcs_out_ptr->cost_ptr); 
  free_members_DoubleMatrixPTR(slcs_out_ptr->cost_segment_ptr);
  free_members_DoubleMatrixPTR(slcs_out_ptr->diag_ptr);
  free_members_DoubleMatrixPTR(slcs_out_ptr->score_ptr);
  free_members_DoubleMatrixPTR(slcs_out_ptr->adder_ptr);
  
  destroy_DoubleMatrixPTR(slcs_out_ptr->score_ptr);
  destroy_DoubleMatrixPTR(slcs_out_ptr->cost_ptr);
  destroy_DoubleMatrixPTR(slcs_out_ptr->cost_segment_ptr);
  destroy_DoubleMatrixPTR(slcs_out_ptr->diag_ptr);
  destroy_DoubleMatrixPTR(slcs_out_ptr->adder_ptr);

  slcs_out_ptr->row = 0;
  slcs_out_ptr->col = 0;
}



void free_members_SaddlePointsPTR(SaddlePointsPTR saddle_points_ptr) {
  saddle_points_ptr->x_ptr->size = 0;
  saddle_points_ptr->y_ptr->size = 0;
 
  free_members_DoubleMatrixPTR(saddle_points_ptr->x_ptr);
  free_members_DoubleMatrixPTR(saddle_points_ptr->y_ptr);

  destroy_DoubleMatrixPTR(saddle_points_ptr->x_ptr);
  destroy_DoubleMatrixPTR(saddle_points_ptr->y_ptr);

  saddle_points_ptr->size = 0;
}

//Spots
void printf_SpotVectorPTR(FILEPTR fout, const SpotVectorPTR  spot_vec_ptr) {
  int i;
  fprintf(fout, "%d\n", spot_vec_ptr->size);
  for (i = 0; i < spot_vec_ptr->size; i++) {
    printf_SpotPTR(fout, &(spot_vec_ptr->vec[i]));  
  }
}

void printf_MemberVectorPTR(FILEPTR fout, const MemberVectorPTR member_vec_ptr, const int space_num, const int start_count) {
  int i; 
  for (i =0+start_count; i < member_vec_ptr->size+start_count; i++) {
    fprintf(fout, "%*s_%d: \n", space_num + (int) strlen("Member"), "Member", i+1);
    fprintf(fout, "%*s\n", space_num + (int) strlen("   -  "), "   -  ");
    
    printf_MemberPTR(fout, &(member_vec_ptr->vec[i]), space_num+strlen("   -  "));
  }
}


// TO print all at once.
void printf_SpotVectorPPTR(FILEPTR fout, const SpotVectorPPTR  spot_vec_pptr, const int match_num) {
  int i,j;
  int tot_size =0;

  for (i = 0; i < match_num; i++)
    tot_size += spot_vec_pptr[i]->size;

  fprintf(fout, "%d\n", tot_size);

  for (i = 0; i < match_num; i++) 
    for (j = 0; j < spot_vec_pptr[i]->size; j++) 
      printf_SpotPTR(fout, &(spot_vec_pptr[i]->vec[j]));
}

void scanf_SpotVectorPTR(FILEPTR  fin, SpotVectorPTR  spot_vec_ptr) {
  int max_size;
  fscanf(fin, "%d", &max_size);

  if(spot_vec_ptr->max_size < max_size)
    realloc_members_SpotVectorPTR(spot_vec_ptr, max_size);


  int i;
  for (i = 0; i < max_size; i++){
    scanf_SpotPTR(fin, &(spot_vec_ptr->vec[i]));
    spot_vec_ptr -> size += 1;
  }

}

void scanf_MemberVectorPTR(FILEPTR  fin, MemberVectorPTR  member_vec_ptr) {
  int max_size;
  fscanf(fin, "%d", &max_size);
  
  if(member_vec_ptr->max_size < max_size)
    realloc_members_MemberVectorPTR(member_vec_ptr, max_size);

  int i;
  for (i = 0; i < member_vec_ptr->max_size; i++)
    scanf_MemberPTR(fin, &(member_vec_ptr->vec[i]));
  member_vec_ptr -> size = member_vec_ptr -> max_size;
}


void printf_SpotPTR(FILEPTR  fout, const SpotPTR  spot_ptr) {

  int  j;
    //, k; 

  fprintf(fout, "Group-Begin\n"); // group starts

  fprintf(fout, "query_path:%s\n", spot_ptr->query_path);
  fprintf(fout, "reference_path:%s\n", spot_ptr->reference_path);
  fprintf(fout, "query_length:%d\n", spot_ptr->query_length);
  fprintf(fout, "group_start:%d\n",spot_ptr->group_start);
  fprintf(fout, "group_end:%d\n",spot_ptr->group_end);
  fprintf(fout, "group_start_q:%d\n",spot_ptr->group_start_q);
  fprintf(fout, "group_end_q:%d\n",spot_ptr->group_end_q);
  
  fprintf(fout, "members_size:%d\n", spot_ptr->member_vec_ptr->size);

  for (j = 0; j < spot_ptr->member_vec_ptr->size; j++) {
    fprintf(fout, "Member-%d\n", j+1);

    printf_MemberPTR(fout, &(spot_ptr->member_vec_ptr->vec[j]), 0);
    /* fprintf(fout, "start:%d\n",spot_ptr->member_vec_ptr->vec[j].start); */
    /* fprintf(fout, "end:%d\n",spot_ptr->member_vec_ptr->vec[j].end); */
    /* fprintf(fout, "start_q:%d\n",spot_ptr->member_vec_ptr->vec[j].start_q); */
    /* fprintf(fout, "end_q:%d\n",spot_ptr->member_vec_ptr->vec[j].end_q); */
    /* fprintf(fout, "score:%f\n",spot_ptr->member_vec_ptr->vec[j].score); */
    /* fprintf(fout, "cost:%f\n",spot_ptr->member_vec_ptr->vec[j].cost); */
    /* fprintf(fout, "cost_actual:%f\n",spot_ptr->member_vec_ptr->vec[j].cost_actual); */
    /* fprintf(fout, "war:%f\n",spot_ptr->member_vec_ptr->vec[j].war); */
    /* fprintf(fout, "waq:%f\n",spot_ptr->member_vec_ptr->vec[j].waq); */

    /* fprintf(fout, "alignment_size:%d\n", spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr->size); */
    /* fprintf(fout, "query reference\n"); */
    /* for (k = 0; k < spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr->size; k++) */
    /*   fprintf(fout, "%d %d\n", spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr->vec[k].query, spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr->vec[k].reference); */
  }// j loop ends

  fprintf(fout, "Best Member\n");
  printf_MemberPTR(fout, spot_ptr->best_member_ptr, 0);
  /* fprintf(fout, "start:%d\n",spot_ptr->best_member_ptr->start); */
  /* fprintf(fout, "end:%d\n",spot_ptr->best_member_ptr->end); */
  /* fprintf(fout, "start_q:%d\n",spot_ptr->best_member_ptr->start_q); */
  /* fprintf(fout, "end_q:%d\n",spot_ptr->best_member_ptr->end_q); */
  /* fprintf(fout, "score:%f\n",spot_ptr->best_member_ptr->score); */
  /* fprintf(fout, "cost:%f\n",spot_ptr->best_member_ptr->cost); */
  /* fprintf(fout, "cost_actual:%f\n",spot_ptr->best_member_ptr->cost_actual); */
  /* fprintf(fout, "war:%f\n",spot_ptr->best_member_ptr->war); */
  /* fprintf(fout, "waq:%f\n",spot_ptr->best_member_ptr->waq); */

  /* fprintf(fout, "alignment_size:%d\n", spot_ptr->best_member_ptr->alignment_vec_ptr->size); */
  /* fprintf(fout, "query reference\n"); */
  /* for (k = 0; k < spot_ptr->best_member_ptr->alignment_vec_ptr->size; k++) */
    /* fprintf(fout, "%d %d\n", spot_ptr->best_member_ptr->alignment_vec_ptr->vec[k].query, spot_ptr->best_member_ptr->alignment_vec_ptr->vec[k].reference); */
}


void printf_MemberPTR(FILEPTR  fout, const MemberPTR  member_ptr, const int space_num) {
  int k;
  fprintf(fout,"%*s: %s\n", space_num + (int) strlen("querypath"),  "querypath", member_ptr->querypath);
  fprintf(fout,"%*s: %s\n", space_num + (int) strlen("refpath"),  "refpath", member_ptr->refpath);
  fprintf(fout,"%*s: %d\n", space_num + (int) strlen("start"),  "start", member_ptr->start);
  fprintf(fout,"%*s: %d\n", space_num + (int) strlen("end"), "end", member_ptr->end);
  fprintf(fout,"%*s: %d\n", space_num + (int) strlen("start_q"),  "start_q", member_ptr->start_q);
  fprintf(fout,"%*s: %d\n", space_num + (int) strlen("end_q"), "end_q", member_ptr->end_q);
  fprintf(fout,"%*s: %lf\n", space_num + (int) strlen("score"), "score", member_ptr->score);
  fprintf(fout,"%*s: %lf\n", space_num + (int) strlen("cost"), "cost", member_ptr->cost);
  fprintf(fout,"%*s: %lf\n", space_num + (int) strlen("cost_actual"), "cost_actual", member_ptr->cost_actual);
  fprintf(fout,"%*s: %lf\n", space_num + (int) strlen("war"), "war", member_ptr->war);
  fprintf(fout,"%*s: %lf\n", space_num + (int) strlen("waq"), "waq", member_ptr->waq);

  fprintf(fout,"%*s: %d\n", space_num + (int) strlen("alignment_size"), "alignment_size", member_ptr->alignment_vec_ptr->size);

  fprintf(fout,"%*s: \n", space_num + (int) strlen("alignemnt"), "alignment");
  fprintf(fout,"%*s\n", space_num + (int) strlen("   -  "), "   -  ");
  fprintf(fout,"%*s: ", space_num + (int) strlen("      query"), "      query");
  for (k = 0; k < member_ptr->alignment_vec_ptr->size; k++)
    fprintf(fout, "%d ", member_ptr->alignment_vec_ptr->vec[k].query);
  fprintf(fout, "\n");

  fprintf(fout,"%*s: ", space_num + (int) strlen("      reference"), "      reference");
  for (k = 0; k < member_ptr->alignment_vec_ptr->size; k++)
    fprintf(fout, "%d ", member_ptr->alignment_vec_ptr->vec[k].reference);
  fprintf(fout, "\n");

  fprintf(fout,"%*s: ", space_num + (int) strlen("      color"), "      color");
  for (k = 0; k < member_ptr->alignment_vec_ptr->size; k++)
    fprintf(fout, "%f ", member_ptr->alignment_vec_ptr->vec[k].color);
  fprintf(fout, "\n");
}


void scanf_SpotPTR(FILEPTR fin, SpotPTR spot_ptr) {

  int  j;
  //, k;
  char* line = NULL;
  size_t len = 0;
  ssize_t read = 0;

  // these two are used for temp storage of mem size and alignment size.
  int members_size;
  //  int align_size;


  assert(fin != NULL); // input stream invalid


  advance_through_wspace(fin);
  // line by line reading of the file.
  read = getline(&line, &len, fin);

  assert(read != -1);

  assert(strcmp(line, "Group-Begin\n") == 0);

  strcpy(spot_ptr->query_path, read_spot_entry(fin, "query_path"));
  strcpy(spot_ptr->reference_path, read_spot_entry(fin, "reference_path"));
  spot_ptr->query_length = atoi(read_spot_entry(fin, "query_length"));
  spot_ptr->group_start = atoi(read_spot_entry(fin, "group_start"));
  spot_ptr->group_end = atoi(read_spot_entry(fin, "group_end"));
  spot_ptr->group_start_q = atoi(read_spot_entry(fin, "group_start_q"));
  spot_ptr->group_end_q = atoi(read_spot_entry(fin, "group_end_q"));

  members_size = atoi(read_spot_entry(fin, "members_size"));

  //spot_ptr->member_vec_ptr = create_MemberVectorPTR();
  if (spot_ptr->member_vec_ptr->max_size < members_size)
    realloc_members_MemberVectorPTR(spot_ptr->member_vec_ptr, members_size);
  
  spot_ptr->member_vec_ptr->size = members_size;
  for (j = 0; j < spot_ptr->member_vec_ptr->size; j++) {
    char buff[50]; // used to store buffer strings
    advance_through_wspace(fin);
    read = getline(&line, &len, fin);
    assert(read != -1);

    sprintf(buff, "Member-%d\n", j+1);
    assert(strcmp(line, buff) == 0);


    scanf_MemberPTR(fin, &(spot_ptr->member_vec_ptr->vec[j]));

    /* spot_ptr->member_vec_ptr->vec[j].start = atoi(read_spot_entry(fin, "start")); */
    /* spot_ptr->member_vec_ptr->vec[j].end = atoi(read_spot_entry(fin, "end")); */
    /* spot_ptr->member_vec_ptr->vec[j].score = strtod(read_spot_entry(fin, "score"), NULL); */
    /* spot_ptr->member_vec_ptr->vec[j].cost = strtod(read_spot_entry(fin, "cost"), NULL); */
    /* spot_ptr->member_vec_ptr->vec[j].cost_actual = strtod(read_spot_entry(fin, "cost_actual"), NULL); */
    /* spot_ptr->member_vec_ptr->vec[j].war = strtod(read_spot_entry(fin, "war"), NULL); */
    /* spot_ptr->member_vec_ptr->vec[j].waq = strtod(read_spot_entry(fin, "waq"), NULL); */
    /* align_size = atoi(read_spot_entry(fin, "alignment_size")); */

    /* //spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr = create_AlignmentVectorPTR(); */
    /* if (spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr->max_size < align_size) */
    /*   realloc_members_AlignmentVectorPTR(spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr, align_size); */


    /* spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr->size = align_size; */
    
    /* // to read the query reference heading */
    /* advance_through_wspace(fin); */
    /* read = getline(&line, &len, fin); */
    /* assert(read != -1); */

    /* assert(strcmp(line, "query reference\n") == 0);  */

    /* for (k = 0; k < spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr->size; k++) { */
    /*   advance_through_wspace(fin); */
    /*   read = getline(&line, &len, fin); */
    /*   assert(read = -1); */

    /*   char* tokens[2]; */
    /*   get_tokens(line,tokens, " "); */
    /*   spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr->vec[k].query = atoi(tokens[0]); */
    /*   spot_ptr->member_vec_ptr->vec[j].alignment_vec_ptr->vec[k].reference = atoi(tokens[1]);  */
    /* }// k loop ends */
  }// j loop ends


  advance_through_wspace(fin);
  read = getline(&line, &len, fin);
  assert(read != -1);
  assert(strcmp(line, "Best Member\n") == 0); 

  scanf_MemberPTR(fin, spot_ptr->best_member_ptr);
  /* spot_ptr->best_member_ptr->start = atoi(read_spot_entry(fin, "start")); */
  /* spot_ptr->best_member_ptr->end = atoi(read_spot_entry(fin, "end")); */
  /* spot_ptr->best_member_ptr->score = strtod(read_spot_entry(fin, "score"), NULL); */
  /* spot_ptr->best_member_ptr->cost = strtod(read_spot_entry(fin, "cost"), NULL); */
  /* spot_ptr->best_member_ptr->cost_actual = strtod(read_spot_entry(fin, "cost_actual"), NULL); */
  /* spot_ptr->best_member_ptr->war = strtod(read_spot_entry(fin, "war"), NULL); */
  /* spot_ptr->best_member_ptr->waq = strtod(read_spot_entry(fin, "waq"), NULL); */
  /* align_size = atoi(read_spot_entry(fin, "alignment_size")); */

  /* //spot_ptr->best_member_ptr->alignment_vec_ptr = create_AlignmentVectorPTR(); */
  /* if (spot_ptr->best_member_ptr->alignment_vec_ptr->max_size < align_size) */
  /*   realloc_members_AlignmentVectorPTR(spot_ptr->best_member_ptr->alignment_vec_ptr, align_size); */

  /* spot_ptr->best_member_ptr->alignment_vec_ptr->size = align_size; */
  /* // to read the query reference heading */
  /* advance_through_wspace(fin); */
  /* read = getline(&line, &len, fin); */
  /* assert(read != -1); */
  /* assert(strcmp(line, "query reference\n") == 0);  */

  /* for (k = 0; k < spot_ptr->best_member_ptr->alignment_vec_ptr->size; k++) { */
  /*   advance_through_wspace(fin); */
  /*   read = getline(&line, &len, fin); */
  /*   assert(read != 1); */
  /*   char* tokens[2]; */
  /*   get_tokens(line,tokens, " "); */
  /*   spot_ptr->best_member_ptr->alignment_vec_ptr->vec[k].query = atoi(tokens[0]); */
  /*   spot_ptr->best_member_ptr->alignment_vec_ptr->vec[k].reference = atoi(tokens[1]); */
  /* } */
}



void scanf_MemberPTR(FILEPTR fin, MemberPTR member_ptr) {
  ssize_t read = 0;
  size_t len = 0;
  char* line = NULL;
  int k;

  member_ptr->start = atoi(read_spot_entry(fin, "start"));
  member_ptr->end = atoi(read_spot_entry(fin, "end"));
  member_ptr->start_q = atoi(read_spot_entry(fin, "start_q"));
  member_ptr->end_q = atoi(read_spot_entry(fin, "end_q"));
  member_ptr->score = strtod(read_spot_entry(fin, "score"), NULL);
  member_ptr->cost = strtod(read_spot_entry(fin, "cost"), NULL);
  member_ptr->cost_actual = strtod(read_spot_entry(fin, "cost_actual"), NULL);
  member_ptr->war = strtod(read_spot_entry(fin, "war"), NULL);
  member_ptr->waq = strtod(read_spot_entry(fin, "waq"), NULL);
  int align_size = atoi(read_spot_entry(fin, "alignment_size"));

  //spot_ptr->best_member_ptr->alignment_vec_ptr = create_AlignmentVectorPTR();
  if (member_ptr->alignment_vec_ptr->max_size < align_size)
    realloc_members_AlignmentVectorPTR(member_ptr->alignment_vec_ptr, align_size);

  member_ptr->alignment_vec_ptr->size = align_size;
  // to read the query reference heading
  advance_through_wspace(fin);
  read = getline(&line, &len, fin);
  assert(read != -1);
  assert(strcmp(line, "query reference\n") == 0); 

  for (k = 0; k < member_ptr->alignment_vec_ptr->size; k++) {
    advance_through_wspace(fin);
    read = getline(&line, &len, fin);
    assert(read != 1);
    char* tokens[2];
    get_tokens(line,tokens, " ");
    member_ptr->alignment_vec_ptr->vec[k].query = atoi(tokens[0]);
    member_ptr->alignment_vec_ptr->vec[k].reference = atoi(tokens[1]);
  }
}


void printf_DoubleMatrixPTR(FILEPTR fout, DoubleMatrixPTR mat){
  fprintf(fout, "%d %d\n", mat->size, mat->dim);
  int i,j;
  for (i = 0; i < mat->size; i++){
    for (j = 0; j < mat->dim; j++)
      fprintf(fout, "%lf ", mat->mat[i][j]);
    fprintf(fout, "\n");
  }
}


void printf_CharMatrixPTR(FILEPTR fout, CharMatrixPTR mat){
  fprintf(fout, "%d %d\n", mat->size, mat->dim);
  int i,j;
  for (i = 0; i < mat->size; i++){
    for (j = 0; j < mat->dim; j++)
      fprintf(fout, "%c ", mat->mat[i][j]);
    fprintf(fout, "\n");
  }
}

void printf_StringMatrixPTR(FILEPTR fout, StringMatrixPTR mat){
  fprintf(fout, "%d %d\n", mat->size, mat->dim);
  int i,j;
  for (i = 0; i < mat->size; i++){
    for (j = 0; j < mat->dim; j++)
      fprintf(fout, "%s ", mat->mat[i][j].data);
    fprintf(fout, "\n");
  }
}



void get_tokens(char *  line, char** tokens, const char *  delimiters) {
  /* Parameters:
   * line- input string   
   * tokens- array of string where the tokens needs to be stored
   * delimeters- delimeters for tokenizations (default value is ":")
   * 
   * Description:
   * Tokenize the string "line" spereated by "delimeters" and store the 
   * tokens in "tokens".
   */

  int i = 0;
  tokens[i] = strtok(line, delimiters);
  while (tokens[i++] != NULL) {
    assert(i<3);
    tokens[i] = strtok(NULL, delimiters); // splitting;
  }
}

int advance_through_wspace(FILEPTR fin) {
  int c = 0;
  while(isspace(fpeek(fin))) {
    fseek(fin, 1, SEEK_CUR);
    ++c;
  }
  return c;
}

int fpeek(FILEPTR fin)
{
    int c;
    c = fgetc(fin);
    ungetc(c, fin);
    return c;
}


char *  read_spot_entry (FILEPTR fin, char *  entry_name) {
  /* Parameters:
   * fin- FILE pointer poiting to a Spot file.
   * entry_name- name of the entry.
   * 
   * Description:
   * Read an entry from the Spot file and return the value associated with
   * that entry as string.
   */
  charPTR line = NULL;
  size_t len = 0;
  ssize_t read = 0;
  char* tokens[5];// (char**) malloc(sizeof(charPTR)*3);  
  //tokens[0] = (charPTR) malloc(sizeof(char)*200);
  //tokens[1] = (charPTR) malloc(sizeof(char)*200);
  

  advance_through_wspace(fin);
  read = getline(&line, &len, fin);
  assert(read != 1);

  get_tokens(line, tokens, ":"); // to get the tokens by ':'

  //  printf("t:%s\te:%s\n", tokens[0], entry_name);

  assert(strcmp(tokens[0], entry_name) == 0);

  charPTR  ret =  tokens[1];
  
  //free(tokens[0]);
  //free(tokens[2]);
  //free(tokens);

  ret[strlen(ret)-1] = '\0';
  return ret;
}


//copy_functions
void copy_MemberPTR(MemberPTR  dest_ptr, const MemberPTR  src_ptr) {
  dest_ptr->start = src_ptr->start;
  dest_ptr->end = src_ptr->end;
  dest_ptr->start_q = src_ptr->start_q;
  dest_ptr->end_q = src_ptr->end_q;
  dest_ptr->score  = src_ptr->score;
  dest_ptr->cost   = src_ptr->cost;
  dest_ptr->cost_actual  = src_ptr->cost_actual;
  dest_ptr->war  = src_ptr->war;
  dest_ptr->waq  = src_ptr->waq;
  strcpy(dest_ptr->querypath, src_ptr->querypath);
  strcpy(dest_ptr->refpath, src_ptr->refpath);


  if (dest_ptr->alignment_vec_ptr->max_size < src_ptr->alignment_vec_ptr->max_size)
    realloc_members_AlignmentVectorPTR(dest_ptr->alignment_vec_ptr, src_ptr->alignment_vec_ptr->max_size);

  copy_AlignmentVectorPTR(dest_ptr->alignment_vec_ptr, src_ptr->alignment_vec_ptr);
}


void copy_AlignmentPTR(AlignmentPTR  dest_ptr, const AlignmentPTR  src_ptr) {
  dest_ptr->query = src_ptr->query;
  dest_ptr->reference = src_ptr->reference;
  dest_ptr->color = src_ptr->color;

}

void copy_AlignmentVectorPTR(AlignmentVectorPTR  dest_ptr, const AlignmentVectorPTR src_ptr) {
  int i;
  for (i = 0; i < src_ptr->size; i++)
    copy_AlignmentPTR(&(dest_ptr->vec[i]), &(src_ptr->vec[i]));
  dest_ptr->size = src_ptr->size;
}


void copy_DoubleMatrixPTR(DoubleMatrixPTR  dest_ptr, const DoubleMatrixPTR  src_ptr) {
  copy_DoubleMatrixPTR_region(dest_ptr, src_ptr, 0, src_ptr->size - 1);
}

// TODO: make copy region version for all the types.
void copy_DoubleMatrixPTR_region(DoubleMatrixPTR  dest_ptr, const DoubleMatrixPTR  src_ptr, const int start, const int end) {
  int i, j;
  for (i = 0; i < end-start+1; i++)
    for (j = 0; j < src_ptr->dim; j++) {     
      //      printf("%d, %d\n", i, j);
      dest_ptr -> mat[i][j] = src_ptr -> mat[start+i][j];
    }

  dest_ptr -> size = end-start+1;
}

void copy_SaddlePointsPTR_region(SaddlePointsPTR dest_ptr, const SaddlePointsPTR src_ptr, const int start, const int end) {
  copy_DoubleMatrixPTR_region(dest_ptr->x_ptr, src_ptr->x_ptr, start, end);
  copy_DoubleMatrixPTR_region(dest_ptr->y_ptr, src_ptr->y_ptr, start, end);

  dest_ptr -> size = end-start-1;
}


// swap types
void swap_int(int * a, int * b) {
  *a = *a  +  *b;
  *b = *a  -  *b;
  *a = *a  -  *b;
} 

// swap types
void swap_double(double * a, double * b) {
  *a = *a  +  *b;
  *b = *a  -  *b;
  *a = *a  -  *b;
} 


// swap types
void swap_float(float *a, float *b) {
  *a = *a  +  *b;
  *b = *a  -  *b;
  *a = *a  -  *b;
} 


// ParametersRLCS
void set_ParametersRLCSPTR(const char *  ctrl_file, ParametersRLCSPTR parameters_rlcs_ptr) {
  assign_default_values_ParametersRLCSPTR(parameters_rlcs_ptr);
  FILEPTR fin =  fopen(ctrl_file, "r");

  

  //assert(fin!=NULL);


  if (fin == NULL) {
    printf("This is tragedy\n");
    exit(-1);
  }

  char dummy[100];
  int ret =  fscanf(fin, "%s", dummy);

  while(ret != EOF) {
    if (strcmp(dummy, "safe_region_factor") == 0) {
      ret = fscanf(fin, "%lf", &(parameters_rlcs_ptr->safe_region_factor));
      assert(ret != EOF);
    } else if (strcmp(dummy, "reference_window_size_ratio") == 0) {
      ret = fscanf(fin, "%lf", &(parameters_rlcs_ptr->reference_window_size_ratio));
      assert(ret != EOF);
    } else if (strcmp(dummy, "Td") == 0) {
      ret = fscanf(fin, "%lf", &(parameters_rlcs_ptr->Td));
      assert(ret != EOF);
    } else if (strcmp(dummy, "rho") == 0) {
      ret = fscanf(fin, "%lf", &(parameters_rlcs_ptr->rho));
      assert(ret != EOF);
    } else if (strcmp(dummy, "beta") == 0) {
      ret = fscanf(fin, "%lf", &(parameters_rlcs_ptr->beta));
      assert(ret != EOF);
    } else if (strcmp(dummy, "seqFilterTd") == 0) {
      ret = fscanf(fin, "%lf", &(parameters_rlcs_ptr->seqFilterTd));
      assert(ret != EOF);
    } else if (strcmp(dummy, "Fs") == 0) {
      ret = fscanf(fin, "%lf", &(parameters_rlcs_ptr->Fs));
      assert(ret != EOF);
    } else if (strcmp(dummy, "frame_shift") == 0) {
      ret = fscanf(fin, "%lf", &(parameters_rlcs_ptr->frame_shift));
      assert(ret != EOF);
    } else if (strcmp(dummy, "algo") == 0) {
      ret = fscanf(fin, "%s", parameters_rlcs_ptr->algo);
      assert(ret != EOF);
    } else if (strcmp(dummy, "cpu") == 0) {
      ret = fscanf(fin, "%d", &(parameters_rlcs_ptr->cpu));
      assert(ret != EOF);
    } else if (strcmp(dummy, "distname") == 0) {
      ret = fscanf(fin, "%s", parameters_rlcs_ptr->distname);
      assert(ret != EOF);
    } else if (strcmp(dummy, "hop_size_fp") == 0) {
      ret = fscanf(fin, "%d", &(parameters_rlcs_ptr->hop_size_fp));
      assert(ret != EOF);
    } else if (strcmp(dummy, "hop_size_sp") == 0) {
      ret = fscanf(fin, "%d", &(parameters_rlcs_ptr->hop_size_sp));
      assert(ret != EOF);
    } else if (strcmp(dummy, "semitone") == 0) {
      ret = fscanf(fin, "%d", &(parameters_rlcs_ptr->semitone));
      assert(ret != EOF);
    } else if (strcmp(dummy, "rlcontext") == 0) {
      ret = fscanf(fin, "%d", &(parameters_rlcs_ptr->rlcontext));
      assert(ret != EOF);
    } else if (strcmp(dummy, "distname2") == 0) {
      ret = fscanf(fin, "%s", parameters_rlcs_ptr->distname2);
      assert(ret != EOF);
    } else if (strcmp(dummy, "END") == 0) {
      printf("All parameters read.\n");
    } else {
      printf("There is some problem with parameter %s.\n", dummy);
    }
    ret =  fscanf(fin, "%s", dummy);
  }

  fclose(fin); 
}

// write this func in .h also
void printf_ParametersRLCSPTR(FILEPTR fout, const ParametersRLCSPTR parameters_rlcs_ptr) {
  fprintf(fout, "safe_region_factor %lf\n", parameters_rlcs_ptr->safe_region_factor);
  fprintf(fout, "reference_window_size_ratio %lf\n", parameters_rlcs_ptr->reference_window_size_ratio);
  fprintf(fout, "Td %lf\n", parameters_rlcs_ptr->Td);
  fprintf(fout, "rho %lf\n", parameters_rlcs_ptr->rho);
  fprintf(fout, "beta %lf\n", parameters_rlcs_ptr->beta);
  fprintf(fout, "seqFilterTd %lf\n", parameters_rlcs_ptr->seqFilterTd);
  fprintf(fout, "Fs %lf\n", parameters_rlcs_ptr->Fs);
  fprintf(fout, "frame_shift %lf\n", parameters_rlcs_ptr->frame_shift);
  fprintf(fout, "algo %s\n", parameters_rlcs_ptr->algo);
  fprintf(fout, "cpu %d\n", parameters_rlcs_ptr->cpu);
  fprintf(fout, "distname %s\n", parameters_rlcs_ptr->distname);
  fprintf(fout, "hop_size_fp %d\n", parameters_rlcs_ptr->hop_size_fp);
  fprintf(fout, "hop_size_sp %d\n", parameters_rlcs_ptr->hop_size_sp);
  fprintf(fout, "semitone %d\n", parameters_rlcs_ptr->semitone);
  fprintf(fout, "rlcontext %d\n", parameters_rlcs_ptr->rlcontext);
  fprintf(fout, "distname2 %s\n", parameters_rlcs_ptr->distname2);
}

void assign_default_values_ParametersRLCSPTR(ParametersRLCSPTR parameters_rlcs_ptr) {
  parameters_rlcs_ptr -> safe_region_factor = 0.5;
  parameters_rlcs_ptr -> reference_window_size_ratio = 1.5;
  parameters_rlcs_ptr -> Td = 0.45;
  parameters_rlcs_ptr -> rho = 0;
  parameters_rlcs_ptr -> seqFilterTd = 0.0;
  parameters_rlcs_ptr -> beta = 0.5;
  parameters_rlcs_ptr -> Fs = 44100;
  parameters_rlcs_ptr -> frame_shift = 441;
  strcpy(parameters_rlcs_ptr -> algo, "rlcs_mod");
  parameters_rlcs_ptr -> cpu = 1;
  parameters_rlcs_ptr -> hop_size_fp = 1;
  parameters_rlcs_ptr -> hop_size_sp = 1;
  strcpy(parameters_rlcs_ptr -> distname,"cubic");
  parameters_rlcs_ptr -> semitone = 10;
  parameters_rlcs_ptr -> rlcontext = 0; // this is same as just distname
  strcpy(parameters_rlcs_ptr -> distname2, "cubic");
  
}
 
// Segment
void fill_SegmentPTR_from_file(char *  file_path, SegmentPTR  segment_ptr) {

  segment_ptr->path = file_path;

  segment_ptr->feat_mat_ptr = create_DoubleMatrixPTR();

  fill_DoubleMatrixPTR_from_file(file_path, segment_ptr->feat_mat_ptr);

  char dummy[500], dummy2[500];
  find_dir_in_path(file_path, dummy);
  find_filename_in_path_without_extn(file_path, dummy2);
  strcat(dummy, "/");
  strcat(dummy, dummy2);
  strcat(dummy, ".roi");
  
  segment_ptr->roi_ptr = create_RangeVectorPTR();
  
  fill_RangeVectorPTR_from_file(dummy, segment_ptr->roi_ptr);

  strcpy(dummy, file_path);

  //TODO: for a time-being make sure segment is made with sp's.
  // Also write a gen_sp function (convert vig's code) and call that fucn insead.
  //char dummy[500], dummy2[50];
  find_dir_in_path(file_path, dummy);
  find_filename_in_path_without_extn(file_path, dummy2);
  strcat(dummy, "/");
  strcat(dummy, dummy2);
  strcat(dummy, ".sp");
  
  segment_ptr->saddle_points_ptr = create_SaddlePointsPTR();
  fill_SaddlePointsPTR_from_file(dummy, segment_ptr->saddle_points_ptr);
}


void fill_DoubleMatrixPTR_from_file(const char *  file_path, DoubleMatrixPTR feat_mat_ptr) {

  int i, j;
  char * buffer;

  // do not forget to free buffer
  buffer = read_file_into_buffer(file_path);
  int bi = 0;

  int size, dim;
  int offset; // to store the offset.
  

  // getting teh row and col of feat_mat;
  sscanf(&buffer[bi], "%d%d%n", &size, &dim, &offset);
  bi += offset;
  //TODO: change it to free and reallocate
  alloc_members_DoubleMatrixPTR(feat_mat_ptr, size, dim);

  feat_mat_ptr->size = feat_mat_ptr->max_size;

  for (i = 0; i < size; i++)
    for (j = 0; j < dim; j++){
      sscanf(&buffer[bi], "%lf%n", &(feat_mat_ptr->mat[i][j]), &offset);
      bi += offset;
    }


  // freeing buffer.
  free(buffer);
}


void fill_SaddlePointsPTR_from_file(const char *  file_path, SaddlePointsPTR  saddle_points_ptr) {
  int i;
  char *  buffer;
  buffer = read_file_into_buffer(file_path);
  int bi = 0;

  saddle_points_ptr->x_ptr = create_DoubleMatrixPTR();
  saddle_points_ptr->y_ptr = create_DoubleMatrixPTR();

  // it is basically a vector.
  saddle_points_ptr->x_ptr->dim = 1;
  saddle_points_ptr->y_ptr->dim = 1;

  int max_size;
  int offset;
  sscanf(&buffer[bi], "%d%*d%n", &max_size, &offset);
  bi += offset;

  alloc_members_DoubleMatrixPTR(saddle_points_ptr->x_ptr, max_size, 1);
  alloc_members_DoubleMatrixPTR(saddle_points_ptr->y_ptr, max_size, 1);

  for (i =0; i < max_size; i++) {
    sscanf(&buffer[bi], "%lf%lf%n", &(saddle_points_ptr->x_ptr->mat[i][0]), &(saddle_points_ptr->y_ptr->mat[i][0]), &offset);
    bi += offset;
  }

  saddle_points_ptr->x_ptr->size = max_size;
  saddle_points_ptr->y_ptr->size = max_size;
  saddle_points_ptr->size = max_size;  
  free(buffer);
}

void fill_RangeVectorPTR_from_file(const char *  file_path, RangeVectorPTR range_vec_ptr) {

  int i;
  char *  buffer;
  buffer = read_file_into_buffer(file_path);
  int bi = 0;

  int rsize;
  int offset;
  sscanf(&buffer[bi], "%d%*d%n", &rsize, &offset);
  bi += offset;
  
  alloc_members_RangeVectorPTR(range_vec_ptr, rsize);

  range_vec_ptr->size = range_vec_ptr->max_size;

  double tempa, tempb;
  for (i = 0; i < rsize; i++) {
    
    sscanf(&buffer[bi], "%lf%lf%n", &tempa, &tempb, &offset);

    range_vec_ptr->vec[i].start = (int)tempa;
    range_vec_ptr->vec[i].end   = (int)tempb;
    
    bi += offset;
  }

  free(buffer);
}


char *  read_file_into_buffer(const char *  file_path) { 
  
  printf("\nReading file: %s\n", file_path);
  
  charPTR buffer; 

  FILEPTR fin = fopen(file_path, "rb");  
  assert(fin != NULL);

  // getting the size of the file in bytes
  fseek(fin, 0L, SEEK_END);
  long sz = ftell(fin);
  fseek(fin, 0L, SEEK_SET);
  
  
  printf("Size of the file in bytes is %ld\n", sz);

  // allocating  buffer with a size to hold the file contents.
  buffer = (char*) malloc(sizeof(char) * sz);

  // reading the file in the buffer
  fread(buffer, 1, sz, fin);
  fclose(fin);

    return buffer;
}


void empty_SegmentPTR(SegmentPTR segment_ptr) {
  segment_ptr->path = NULL;

  empty_DoubleMatrixPTR(segment_ptr->feat_mat_ptr);
  destroy_DoubleMatrixPTR(segment_ptr->feat_mat_ptr);

  empty_RangeVectorPTR(segment_ptr->roi_ptr);
  destroy_RangeVectorPTR(segment_ptr->roi_ptr);

  empty_SaddlePointsPTR(segment_ptr->saddle_points_ptr);
  destroy_SaddlePointsPTR(segment_ptr->saddle_points_ptr);
}

void empty_DoubleMatrixPTR(DoubleMatrixPTR mat_ptr) {
  free_members_DoubleMatrixPTR(mat_ptr);
}

void empty_RangeVectorPTR(RangeVectorPTR vec_ptr) {
  free_members_RangeVectorPTR(vec_ptr);
}

void empty_SaddlePointsPTR(SaddlePointsPTR  saddle_points_ptr) {
  free_members_DoubleMatrixPTR(saddle_points_ptr->x_ptr);
  destroy_DoubleMatrixPTR(saddle_points_ptr->x_ptr);

  free_members_DoubleMatrixPTR(saddle_points_ptr->y_ptr);
  destroy_DoubleMatrixPTR(saddle_points_ptr->y_ptr);
}


// char *  Manipulations 

void find_dir_in_path(const char *  file_path, char *  dir) {
  int len = strlen(file_path);
  int i = len-1;
  while(i >= 0) {
    if (file_path[i] == '/')
      break;
    --i;
  }
  int j;
  for (j = 0; j <= i-1; j++) {
    dir[j] = file_path[j];
  }
  dir[j] = '\0';
}

void find_filename_in_path(const char *  file_path, char *  filename) {
  int len = strlen(file_path);
  int st_i = len-1;
  while(st_i >= 0) {
    if (file_path[st_i] == '/')
      break;
    --st_i;
  }
  st_i = st_i < 0 ? 0: st_i;

  int j;
  for (j = 0; j < len - (st_i+1) + 1 ; j++) {
    filename[j] = file_path[j+ st_i+1];
  }
  
  filename[j] = '\0';
    
}

void find_filename_in_path_without_extn(const char *  file_path, char *  filename) {
  int len = strlen(file_path);
  int st_i = len-1;
  while(st_i >= 0) {
    if (file_path[st_i] == '/')
      break;
    --st_i;
  }

  int ed_i = len-1;
  while(ed_i >= st_i) {
    if (file_path[ed_i] == '.')
      break;
    --ed_i;
  }

  int j;
  for (j = 0; j < ed_i-1 - (st_i+1) +1; j++) {
    filename[j] = file_path[j+ st_i+1];
  }
  
  filename[j] = '\0';
 
}

