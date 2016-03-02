#include "sp.h"


void find_voiced_parts(DoubleMatrixPTR pitch_ptr, DoubleMatrixPTR voiced_parts_ptr, int v_th, int s_th) {

  // pitch should be of one dimension.
  assert(pitch_ptr->dim == 1);

  int st = 0, ed = 0;
 
  while (ed < pitch_ptr->size && pitch_ptr->mat[ed][0] == 0) {
    ++st;
    ++ed;
  }  
  int flag_voice = 0;
  while (ed < pitch_ptr->size) {
    if (pitch_ptr->mat[ed][0]) {
      ++ed;
      flag_voice = 1;
    } else {
      flag_voice = 0;
      int sil_st = ed;
      ++ed;
      while(ed < pitch_ptr->size && pitch_ptr->mat[ed][0] == 0)
	++ed;
      
      int sil_len = ed-1 - sil_st + 1;
      
      if (ed>=pitch_ptr->size || sil_len > s_th){ 
	int v_len = sil_st-1 - st + 1;
	if(v_len > v_th) {  
	  if (voiced_parts_ptr->size==voiced_parts_ptr->max_size)
	    realloc_members_DoubleMatrixPTR(voiced_parts_ptr, 2*voiced_parts_ptr->max_size);

	  voiced_parts_ptr->mat[voiced_parts_ptr->size][0] = st;

	  voiced_parts_ptr->mat[voiced_parts_ptr->size][1] = sil_st-1; // sil_st is ed + 1;

	  voiced_parts_ptr->size += 1;
	}	
	st = ed;
      }
    }
  }

  // for the last voiced part
  if (flag_voice) {
    int v_len = ed - 1 - st + 1;
    if (v_len > v_th) {
      if (voiced_parts_ptr->size==voiced_parts_ptr->max_size)
	realloc_members_DoubleMatrixPTR(voiced_parts_ptr, 2*voiced_parts_ptr->max_size);

      voiced_parts_ptr->mat[voiced_parts_ptr->size][0] = st;

      voiced_parts_ptr->mat[voiced_parts_ptr->size][1] = ed-1; // sil_st is ed + 1;

      voiced_parts_ptr->size += 1;
    }
  }
  
 



}
