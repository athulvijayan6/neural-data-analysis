#include "types.h"

int main(int argc, char* argv[]) {

  ParametersRLCSPTR parameters_rlcs_ptr = create_ParametersRLCSPTR();
  
  set_ParametersRLCSPTR(argv[1], parameters_rlcs_ptr);

  //  FILEPTR fcopy = fopen("ctrl_file_sample","w");

  // printf_ParametersRLCSPTR(fcopy, parameters_rlcs_ptr);
  printf_ParametersRLCSPTR(stdout, parameters_rlcs_ptr);
  
  // fclose(fcopy);

  destroy_ParametersRLCSPTR(parameters_rlcs_ptr);
  return 0;
}
