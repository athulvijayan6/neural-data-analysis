/* Description: Contains the structures defines for usage and the methods
 * based on them
 * Creation Date: 9 March 2014
 * Author: Shrey Dutta
 */


#ifndef TYPES_H
#define TYPES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>


// Creating aliases for the pointers and double pointers of the standard 
// types for clear and easy usage.
typedef float* floatPTR;
typedef float** floatPPTR;

typedef double* doublePTR;
typedef double** doublePPTR;

typedef int* intPTR;
typedef int** intPPTR;

typedef char* charPTR;
typedef char** charPPTR;

typedef FILE* FILEPTR;
typedef FILE** FILEPPTR;


// Simple Creation of Data Structures
typedef struct Range {
  int start;
  int end;
} Range;
typedef Range* RangePTR;
typedef Range** RangePPTR;


typedef struct Point{
  double x;
  double y;
} Point;
typedef Point* PointPTR;
typedef Point** PointPPTR;

typedef struct Alignment { // query reference allignment.
  int query;
  int reference;
  float color;
} Alignment;
typedef Alignment* AlignmentPTR;
typedef Alignment** AlignmentPPTR;


// This DS contains the variables needed for tuning RLCS. These are supplied
// in a ctrl file.
typedef struct ParametersRLCS {

  // This parameter gives the factor which determine the safe range for the
  // size of teh reference. Too small refernce should be compared with
  // query.
  // It is used like this.
  //  if (padd_length >  parameters_rlcs.safe_region_factor * (query_feat_ptr->size)) return;
  double safe_region_factor;

  // defining analysis window length
  // double window_size = floor(query_feat_ptr->size * parameters_rlcs_ptr->reference_window_size_ratio);
  double reference_window_size_ratio;

  //if (dist < parameters_rlcs_ptr->Td) {
  //	rlcs_out_ptr->direction[i][j] = 1; // for cross
  double Td;

  // for score update;
  // determines how much query should be matrched.
  double rho;

  // weight on WAR and WAQ
  double beta;
  
  // threshold on score to be considered for a sequence.
  // if (rlcs_out_ptr->score[max_score_index[0]][max_score_index[1]] < parameters_rlcs_ptr->seqFilterTd) return;
  double seqFilterTd;

  // Sampling frequency;
  double Fs;

  // the frame shift;
  double frame_shift;


  char algo[50];
  
  int cpu;

  // hop size in first pass
  int hop_size_fp;

  // hop size in second pass
  int hop_size_sp;

  char distname[50];

  // how many units makes a semiton
  int semitone;

  // (valid when dist is dtw) how much right-left-context for dtw distance.
  int rlcontext;

  // (valid when dist is dtw) tells what distance function should be used by dtw.
  char distname2[50];

} ParametersRLCS; 
typedef ParametersRLCS* ParametersRLCSPTR;
typedef ParametersRLCS** ParametersRLCSPPTR;


// Vectors
typedef struct IntVector {
  int * vec;
  int size;
  int max_size;
} IntVector;
typedef IntVector* IntVectorPTR;
typedef IntVector** IntVectorPPTR;


typedef struct DoubleVector {
  double * vec;
  int size;
  int max_size;
} DoubleVector;
typedef DoubleVector* DoubleVectorPTR;
typedef DoubleVector** DoubleVectorPPTR;


typedef struct String { 
  char data[500];
} String;


typedef struct StringVector {
  String * vec;
  int size;
  int max_size;
} StringVector;
typedef StringVector* StringVectorPTR;
typedef StringVector** StringVectorPPTR;


typedef struct RangeVector {
  Range* vec;
  int size;
  int max_size;
} RangeVector;
typedef RangeVector* RangeVectorPTR;
typedef RangeVector** RangeVectorPPTR;


typedef struct PointVector {
  Point * vec;
  int size;
  int max_size;
} PointVector;
typedef PointVector* PointVectorPTR;
typedef PointVector** PointVectorPPTR;


typedef struct AlignmentVector {
  Alignment * vec;
  int size;
  int max_size;
} AlignmentVector;
typedef AlignmentVector* AlignmentVectorPTR;
typedef AlignmentVector** AlignmentVectorPPTR;

// Matrices
typedef struct IntMatrix {
  int ** mat;
  int size;
  int dim;
  int max_size;
} IntMatrix;
typedef IntMatrix* IntMatrixPTR;
typedef IntMatrix** IntMatrixPPTR;


typedef struct DoubleMatrix {
  double ** mat;
  int dim;
  int size;
  int max_size;
} DoubleMatrix;
typedef DoubleMatrix* DoubleMatrixPTR;
typedef DoubleMatrix** DoubleMatrixPPTR;


typedef struct CharMatrix {
  char ** mat;
  int dim;
  int size;
  int max_size;
} CharMatrix;
typedef CharMatrix* CharMatrixPTR;
typedef CharMatrix** CharMatrixPPTR;



typedef struct StringMatrix {
  String ** mat;
  int dim;
  int size;
  int max_size;
} StringMatrix;
typedef StringMatrix* StringMatrixPTR;
typedef StringMatrix** StrinMatrixPPTR;



// Special structs
typedef struct SaddlePoints {
  // no need for this double matrix.
  // Int vector should have also sufficed.
  // But RLCS_window only accepts double matrix and we need to pass
  // saddle point in the first pass. That is why!!!
  DoubleMatrixPTR x_ptr;
  DoubleMatrixPTR y_ptr;
  int size;
} SaddlePoints;
typedef SaddlePoints* SaddlePointsPTR;
typedef SaddlePoints** SaddlePointsPPTR;


// structure representing music segment.
typedef struct Segment {
  char *  path;               // id of this piece.
  DoubleMatrixPTR feat_mat_ptr;        // features.
  SaddlePointsPTR saddle_points_ptr;           // saddle points
  RangeVectorPTR roi_ptr;
} Segment;
typedef Segment* SegmentPTR;
typedef Segment** SegmentPPTR;

// Member of a group which represents the Spot
typedef struct Member {
  char querypath[500]; 
  char refpath[500];

  int start;
  int end;
  int start_q;
  int end_q;

  double score;
  double cost;
  double cost_actual; // useless for slsc
  double war; // useless for slcs
  double waq; // useless for scls
  AlignmentVectorPTR alignment_vec_ptr; // for slcs
                                        // it will contain the segment 
                                        // region only
}  Member;
typedef Member* MemberPTR;
typedef Member** MemberPPTR;

typedef struct MemberVector {
  Member* vec;
  int size;
  int max_size;
} MemberVector;
typedef MemberVector* MemberVectorPTR;
typedef MemberVector** MemberVectorPPTR;


// data structure to store the spotted motif group.
typedef struct Spot {
  char *  reference_path;         
  char *  query_path;
  int query_length; // This does not work
  int  group_start;
  int  group_end;
  int group_start_q;
  int group_end_q;
  MemberVectorPTR member_vec_ptr; // for all the members of the group
  
  MemberPTR best_member_ptr; // the best member of all;
} Spot;
typedef Spot* SpotPTR;
typedef Spot** SpotPPTR;

typedef struct SpotVector {
  Spot * vec;
  int size;
  int max_size;
} SpotVector;
typedef SpotVector* SpotVectorPTR;
typedef SpotVector** SpotVectorPPTR;


// To store the output of RLCS
typedef struct OutputRLCS {
  DoubleMatrixPTR score_ptr;
  DoubleMatrixPTR  cost_ptr;
  DoubleMatrixPTR  cost_actual_ptr;
  DoubleMatrixPTR  war_ptr;
  DoubleMatrixPTR  waq_ptr;
  DoubleMatrixPTR  diag_ptr;

  int row;
  int col;
} OutputRLCS;
typedef OutputRLCS* OutputRLCSPTR;
typedef OutputRLCS** OutputRLCSPPTR;


// To store the output of SLCS
typedef struct OutputSLCS {
  DoubleMatrixPTR score_ptr; // have 1 extra dimesion to row, col
  DoubleMatrixPTR  cost_ptr;
  DoubleMatrixPTR  cost_segment_ptr;
  DoubleMatrixPTR  adder_ptr;
  DoubleMatrixPTR  diag_ptr; // have 1 extra dimension to row col

  int row;
  int col;
} OutputSLCS;
typedef OutputSLCS* OutputSLCSPTR;
typedef OutputSLCS** OutputSLCSPPTR;




// Creation of Pointers. We only work with pointers.
RangePTR create_RangePTR();
PointPTR create_PointPTR();
AlignmentPTR create_AlignmentPTR();
OutputRLCSPTR create_OutputRLCSPTR();
OutputSLCSPTR create_OutputSLCSPTR();
ParametersRLCSPTR create_ParametersRLCSPTR();
IntVectorPTR create_IntVectorPTR();
DoubleVectorPTR create_DoubleVectorPTR();
StringVectorPTR create_StringVectorPTR();
RangeVectorPTR create_RangeVectorPTR();
PointVectorPTR create_PointVectorPTR();
AlignmentVectorPTR create_AlignmentVectorPTR();
IntMatrixPTR create_IntMatrixPTR();
DoubleMatrixPTR create_DoubleMatrixPTR();
CharMatrixPTR create_CharMatrixPTR();
StringMatrixPTR create_StringMatrixPTR();
SaddlePointsPTR create_SaddlePointsPTR();
SegmentPTR create_SegmentPTR();
MemberPTR create_MemberPTR();
MemberVectorPTR create_MemberVectorPTR();
SpotPTR create_SpotPTR();
SpotVectorPTR create_SpotVectorPTR();

// Destroying them
void destroy_RangePTR(RangePTR temp);
void destroy_PointPTR(PointPTR temp);
void destroy_AlignmentPTR(AlignmentPTR temp);
void destroy_OutputRLCSPTR(OutputRLCSPTR temp);
void destroy_OutputSLCSPTR(OutputSLCSPTR temp);
void destroy_ParametersRLCSPTR(ParametersRLCSPTR temp);
void destroy_IntVectorPTR(IntVectorPTR temp);
void destroy_DoubleVectorPTR(DoubleVectorPTR temp);
void destroy_StringVectorPTR(StringVectorPTR temp);
void destroy_RangeVectorPTR(RangeVectorPTR temp);
void destroy_PointVectorPTR(PointVectorPTR temp);
void destroy_AlignmentVectorPTR(AlignmentVectorPTR temp);
void destroy_IntMatrixPTR(IntMatrixPTR temp);
void destroy_DoubleMatrixPTR(DoubleMatrixPTR temp);
void destroy_CharMatrixPTR(CharMatrixPTR temp);
void destroy_StringMatrixPTR(StringMatrixPTR temp);
void destroy_SaddlePointsPTR(SaddlePointsPTR temp);
void destroy_SegmentPTR(SegmentPTR temp);
void destroy_MemberPTR(MemberPTR temp);
void destroy_MemberVectorPTR(MemberVectorPTR temp);
void destroy_SpotPTR(SpotPTR temp);
void destroy_SpotVectorPTR(SpotVectorPTR temp);

// alloc members of different structures
void alloc_members_OutputRLCSPTR(OutputRLCSPTR rlcs_out_ptr, const int row, const int col);
void alloc_members_OutputSLCSPTR(OutputSLCSPTR slcs_out_ptr, const int row, const int col);
void alloc_members_SpotPTR(SpotPTR spot_ptr, const int members_num);
void alloc_members_IntVectorPTR(IntVectorPTR vec_ptr, const int max_size);
void alloc_members_DoubleVectorPTR(DoubleVectorPTR vec_ptr, const int max_size);
void alloc_members_StringVectorPTR(StringVectorPTR vec_ptr, const int max_size);
void alloc_members_RangeVectorPTR(RangeVectorPTR vec_ptr, const int max_size);
void alloc_members_AlignmentVectorPTR(AlignmentVectorPTR vec_ptr, const int max_size);
void alloc_members_IntMatrixPTR(IntMatrixPTR mat_ptr, const int max_size, const int dim);
void alloc_members_DoubleMatrixPTR(DoubleMatrixPTR mat_ptr, const int max_size, const int dim);
void alloc_members_CharMatrixPTR(CharMatrixPTR mat_ptr, const int max_size, const int dim);
void alloc_members_StringMatrixPTR(StringMatrixPTR mat_ptr, const int max_size, const int dim);
void alloc_members_MemberVectorPTR(MemberVectorPTR vec_ptr, const int max_size);
void alloc_members_SpotVectorPTR(SpotVectorPTR vec_ptr, const int max_size);
void alloc_members_SaddlePointsPTR(SaddlePointsPTR vec_ptr, const int max_size);
void alloc_members_MemberPTR(MemberPTR member_ptr, int align_size);

// realloction of memory of the members of some structers
void realloc_members_IntVectorPTR(IntVectorPTR vec_ptr, const int max_size);
void realloc_members_DoubleVectorPTR(DoubleVectorPTR vec_ptr, const int max_size);
void realloc_members_StringVectorPTR(StringVectorPTR vec_ptr, const int max_size);
void realloc_members_RangeVectorPTR(RangeVectorPTR vec_ptr, const int max_size);
void realloc_members_AlignmentVectorPTR(AlignmentVectorPTR vec_ptr, const int max_size);
void realloc_members_IntMatrixPTR(IntMatrixPTR mat_ptr, const int max_size);
void realloc_members_DoubleMatrixPTR(DoubleMatrixPTR mat_ptr, const int max_size);
void realloc_members_CharMatrixPTR(CharMatrixPTR mat_ptr, const int max_size);
void realloc_members_StringMatrixPTR(StringMatrixPTR mat_ptr, const int max_size);
void realloc_members_MemberVectorPTR(MemberVectorPTR vec_ptr, const int max_size);
void realloc_members_SpotVectorPTR(SpotVectorPTR vec_ptr, const int max_size);

// freeing the memory of the structures.
void free_members_OutputRLCSPTR(OutputRLCSPTR rlcs_out_ptr);
void free_members_OutputSLCSPTR(OutputSLCSPTR slcs_out_ptr);
void free_members_SpotPTR(SpotPTR spot_ptr);
void free_members_IntVectorPTR(IntVectorPTR vec_ptr);
void free_members_DoubleVectorPTR(DoubleVectorPTR vec_ptr);
void free_members_StringVectorPTR(StringVectorPTR vec_ptr);
void free_members_RangeVectorPTR(RangeVectorPTR vec_ptr);
void free_members_AlignmentVectorPTR(AlignmentVectorPTR vec_ptr);
void free_members_IntMatrixPTR(IntMatrixPTR mat_ptr);
void free_members_DoubleMatrixPTR(DoubleMatrixPTR mat_ptr);
void free_members_CharMatrixPTR(CharMatrixPTR mat_ptr);
void free_members_StringMatrixPTR(StringMatrixPTR mat_ptr);
void free_members_MemberVectorPTR(MemberVectorPTR vec_ptr);
void free_members_SpotVectorPTR(SpotVectorPTR vec_ptr);
void free_members_SaddlePointsPTR(SaddlePointsPTR vec_ptr);

// Routines to read and write Spots from/to a file.
// TODO: write print and scan function for all types
void printf_SpotVectorPTR(FILEPTR fout, const SpotVectorPTR spot_vec_ptr);

// TODO: verify whether this PPTR func is required or not.
void printf_SpotVectorPPTR(FILEPTR fout, const SpotVectorPPTR spot_vec_pptr, const int match_num);
void scanf_SpotVectorPTR(FILEPTR fin, SpotVectorPTR spot_vec);
void printf_SpotPTR(FILEPTR fout, const SpotPTR spot);
void scanf_SpotPTR(FILEPTR fin, SpotPTR spot);

void printf_MemberVectorPTR(FILEPTR fout, const MemberVectorPTR member_vec_ptr, const int space_num, const int start_count);
void scanf_MemberVectorPTR(FILEPTR fin, MemberVectorPTR member_vec_ptr);
void printf_MemberPTR(FILEPTR fout, const MemberPTR member_ptr, const int space_num);
void scanf_MemberPTR(FILEPTR fin, MemberPTR member_ptr);

void printf_DoubleMatrixPTR(FILEPTR fin, DoubleMatrixPTR mat);
void printf_CharMatrixPTR(FILEPTR fin, CharMatrixPTR mat);
void printf_StringMatrixPTR(FILEPTR fin, StringMatrixPTR mat);


void get_tokens(char *  line, char** tokens, const char *  delimiters);
int advance_through_wspace(FILEPTR fin);
int fpeek(FILEPTR fin);
char *  read_spot_entry(FILEPTR fin, char *  entry_name);

// copy functions
// TODO: make copy function for all types.
void copy_MemberPTR(MemberPTR dest_ptr, const MemberPTR src_ptr);
void copy_AlignmentPTR(AlignmentPTR dest_ptr, const AlignmentPTR src_ptr);
void copy_AlignmentVectorPTR(AlignmentVectorPTR dest_ptr, const AlignmentVectorPTR src_ptr);
void copy_DoubleMatrixPTR(DoubleMatrixPTR dest_ptr, const DoubleMatrixPTR src_ptr);
void copy_DoubleMatrixPTR_region(DoubleMatrixPTR dest_ptr, const DoubleMatrixPTR src_ptr, const int start, const int end);
void copy_SaddlePointsPTR_region(SaddlePointsPTR dest_ptr, const SaddlePointsPTR src_ptr, const int start, const int end);
// Swap types
// TODO: to write for all types.
void swap_int(int * a, int * b);
void swap_double(double * a, double * b);
void swap_float(float *a, float *b);

// ParametrsRLCS: to set this structure from ctrl file and
// assigning default values
void set_ParametersRLCSPTR(const char *  ctrl_file, ParametersRLCSPTR parameters_rlcs_ptr);
void assign_default_values_ParametersRLCSPTR(ParametersRLCSPTR parameters_rlcs_ptr);
void printf_ParametersRLCSPTR(FILEPTR fout, const ParametersRLCSPTR parameters_rlcs_ptr);

// Segments: routines to fill DS by reading from a file.
void fill_SegmentPTR_from_file(char *  file_path, SegmentPTR segment_ptr);
void fill_DoubleMatrixPTR_from_file(const char *  file_path, DoubleMatrixPTR feat_mat_ptr);
void fill_SaddlePointsPTR_from_file(const char *  file_path, SaddlePointsPTR  saddle_points_ptr);
void fill_RangeVectorPTR_from_file(const char *  file_path, RangeVectorPTR range_vec_ptr);

// empting the Structures
void empty_SegmentPTR(SegmentPTR segment_ptr);
void empty_DoubleMatrixPTR(DoubleMatrixPTR mat_ptr);
void empty_RangeVectorPTR(RangeVectorPTR vec_ptr);
void empty_SaddlePointsPTR(SaddlePointsPTR  saddle_points_ptr);
// reading a file and and saving to the buffer.
char *  read_file_into_buffer(const char *  file_path);

// char *  manipulations 
void find_dir_in_path(const char *  file_path, char *  dir);
void find_filename_in_path_without_extn(const char *  file_path, char *  filename);
void find_filename_in_path(const char *  file_path, char *  filename);


#endif
