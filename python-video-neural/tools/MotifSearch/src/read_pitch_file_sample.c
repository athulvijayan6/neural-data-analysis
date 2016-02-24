#include "types.h"

int main(int argc, char* argv[]) {
  
  char * buffer = read_file_into_buffer(argv[1]);
  int buflen  =  strlen(buffer);
  int bi = 0;

  int offset;
  DoubleMatrixPTR time_pitch_ptr = create_DoubleMatrixPTR();
  alloc_members_DoubleMatrixPTR(time_pitch_ptr, 100, 2);
  int i= 0;

  while (bi < buflen) {
    
    if (time_pitch_ptr->max_size <= i)
      realloc_members_DoubleMatrixPTR(time_pitch_ptr, time_pitch_ptr->max_size*2);
    
    sscanf(&buffer[bi], "%lf%lf%n", &(time_pitch_ptr->mat[i][0]), &(time_pitch_ptr->mat[i][1]), &offset);
    printf("%lf\t%lf\n", time_pitch_ptr->mat[i][0], time_pitch_ptr->mat[i][1]);    

    bi += offset;
    
    // to remove white spaces
    while(isspace(buffer[bi]))
    bi++;

    ++i;
  }
  time_pitch_ptr->size = i;

  free_members_DoubleMatrixPTR(time_pitch_ptr);
  destroy_DoubleMatrixPTR(time_pitch_ptr);

  free(buffer);

  return  0;
}
