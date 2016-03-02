#include "types.h"

// finds the start and end indices of the voiced part in a DoubleMatrix of dim 2. The threshold of the silence and voiced parts is decided by v_th andf s_th
void find_voiced_parts(DoubleMatrixPTR pitch_ptr, DoubleMatrixPTR voiced_parts_ptr, int v_th, int s_th);
