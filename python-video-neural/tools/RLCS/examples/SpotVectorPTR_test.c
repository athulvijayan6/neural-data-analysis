#include "types.h"


int main(int argc, char* argv[]) {
  SpotVectorPTR spot_vec_ptr =  create_SpotVectorPTR();
  alloc_members_SpotVectorPTR(spot_vec_ptr, 1);
  FILEPTR fin  =  fopen(argv[1],"r");

  scanf_SpotVectorPTR(fin, spot_vec_ptr);
  printf_SpotVectorPTR(stdout, spot_vec_ptr);

  fclose(fin);
  free_members_SpotVectorPTR(spot_vec_ptr);
  destroy_SpotVectorPTR(spot_vec_ptr);

  return 0;
}
