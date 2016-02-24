/* Description: 
 * Creation Date: 24 March 2014
 * Author: Shrey Dutta
 */

#include "types.h"
#include "sp.h"

int main (int argc, char* argv[]) {
  if (argc < 4){
    printf("Usage: get_roi pathOfFileContainingListOfFiles voiceThreshold silenceThreshold outPath\n");
    exit(-1);
  }
  printf("Inside\n");
  char *  list_of_files = argv[1];
  double v_th = strtod(argv[2], NULL);
  double s_th = strtod(argv[3], NULL);
  char *   out_path = argv[4];
  // Reading the filelist into an array.
  FILEPTR fin = fopen(list_of_files, "r");

  assert(fin != NULL);

  int num_of_files;
  
  fscanf(fin, "%d", &num_of_files);

  char files[5000][500];
  int i;

  /*
  for (i = 0; i < num_of_files; i++){
    fscanf(fin, "%s", files[i]);
    printf("%s\n", files[i]);
  }
 */ 



  // This way is better than above
  // This we can also have tonic in the list files.
  i=0;
  char str[500];
  while(fgets(str, 500, fin)) {
    if (strlen(str) < 2)
      continue;
    // gettin first word
    char ch;
    int ci = 0;
    //files[i]
    while (str[ci] != '\0') {
      if (str[ci] == '\n' || str[ci] == ' ') {
        files[i][ci] = '\0';
        ++i;
        break;
      }
      files[i][ci] = str[ci];
      ++ci;
    }
  }
  
  
  assert(i == num_of_files);

  fclose(fin);
  // Reading over;


  DoubleMatrixPTR roi_ptr = create_DoubleMatrixPTR();
  DoubleMatrixPTR mat_ptr =  create_DoubleMatrixPTR();
  for (i=0; i < num_of_files; i++) {
    alloc_members_DoubleMatrixPTR(roi_ptr, 50, 2);
    //alloc_members_DoubleMatrixPTR(mat_ptr, 50, 1);
    fill_DoubleMatrixPTR_from_file(files[i], mat_ptr);
   
    find_voiced_parts(mat_ptr, roi_ptr, v_th, s_th);


    //setting the filename 
    char dummy[500], dummy2[100];
    strcpy(dummy, out_path);
    find_filename_in_path_without_extn(files[i], dummy2);
    strcat(dummy, dummy2);
    strcat(dummy, ".roi");


    FILEPTR fout = fopen(dummy, "w");


    printf_DoubleMatrixPTR(fout, roi_ptr);
    
    fclose(fout);
    free_members_DoubleMatrixPTR(roi_ptr);
    free_members_DoubleMatrixPTR(mat_ptr);
  }
  
  destroy_DoubleMatrixPTR(roi_ptr);
  destroy_DoubleMatrixPTR(mat_ptr);
  return 0;
}

