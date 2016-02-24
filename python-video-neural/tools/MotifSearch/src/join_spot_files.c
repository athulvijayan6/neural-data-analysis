#include "types.h"
/* Join spot files
 */



// my main
int main(int argc, char* argv[]) {
  
  if (argc <= 3 +1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: %s <give the spotfiles to join seperated by space> <name of the joined files> ", argv[0]);
    exit(-1);
  }
  
  int num_spot_files = argc - 1 - 1;
  int i;



  SpotVectorPPTR spots_all_pptr = (SpotVectorPPTR) malloc(sizeof(SpotVectorPTR)*num_spot_files);


  
  // creating and allocating spot vectors
  for(i = 0; i < num_spot_files; i++){
    spots_all_pptr[i] = create_SpotVectorPTR();
    // simply allocating 10. It will be rellocated inside if necessary.
    alloc_members_SpotVectorPTR(spots_all_pptr[i], 10);

  }
  

  //TODO: This loop can be made parallel
  for (i = 0; i < num_spot_files; i++) {
    
    fprintf(stdout, "Spot_file_no:\t%d\n", i);
    FILEPTR fin  =  fopen(argv[i+1], "r");
    scanf_SpotVectorPTR(fin, spots_all_pptr[i]);
  }

  FILEPTR fout = fopen(argv[argc-1], "w");
  printf_SpotVectorPPTR(fout, spots_all_pptr, num_spot_files);

  return 0;
}
