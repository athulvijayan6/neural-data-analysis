#include "types.h"
/* this will parse the spot vector and prepare some numbers of list
 * equal to the number of unique motifs queried. The file will contain the
 * info of the matches in below given format. This file is to be processed 
 * further to get the desired result.
 * <query_file_name> <ref_file_name> <start_sample> <end_sample> <score>
 */

// my main
int main(int argc, char* argv[]) {
  
  if (argc != 2 +1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: %s <result_file> <out_dir_path>", argv[0]);
    exit(-1);
  }
  
  int top_n = atoi(argv[2]);

  SpotVectorPTR spot_vec_ptr = create_SpotVectorPTR();
  alloc_members_SpotVectorPTR(spot_vec_ptr, 500);
  
  // not sure whether the read binary mode works.
  FILEPTR fin;
  if ((fin = fopen(argv[1], "rb")) == NULL) {
    printf("ERROR: Cannot open file %s\n", argv[1]);
    exit(-1);
  }
  
  // reading the result into spot vector
  printf("Reading spot_vector...\n");
  scanf_SpotVectorPTR(fin, spot_vec_ptr);
  printf("Reading over, %d spots read\n", spot_vec_ptr->size);
  // closing the file
  fclose(fin);

   
  // sorting the spot_vector_ptr w.r.t. to best_member score in descending order
  //   printf("Sorting spots...\n"); 
  //  qsort(spot_vec_ptr->vec, spot_vec_ptr->size, sizeof(Spot), spot_comp_func); 
  //   printf("Sorting over\n"); 
  
  
  char queries[200][500];
  int query_ctr = 0;
  {
    // for query 
    int i;
    for (i = 0; i < spot_vec_ptr -> size; i++) { 
      int j = 0; 
      for (j = 0; j < query_ctr; j++) 
	if (strcmp(spot_vec_ptr->vec[i].query_path, queries[j]) == 0) 
	  break; 
      // if current query is not found in queries then add it to the list. */
      if (j == query_ctr) { 
	strcpy(queries[query_ctr], spot_vec_ptr->vec[i].query_path); 
	++query_ctr; 
      } 
    } 

 
    char dummy[100];
    strcpy(dummy, argv[1]);
    strcat(dummy, "_query.map");
    FILEPTR querymap = fopen(dummy, "w");

 
    for (i = 0; i <  query_ctr; i++){ 
      printf("%d)\t%s\n", i, queries[i]); 
      fprintf(querymap, "%d %s\n", i, queries[i]); 
    }

    fclose(querymap);

      
  }


  // loop over and do something
  int i;

  FILEPTR ch_wave_out = fopen("ch_wave_out_file", "w");

    
  for (i = 0; i < spot_vec_ptr -> size; i++) {

    int j;
    for (j = 0; j < query_ctr; j++) {
      if (strcmp(spot_vec_ptr->vec[i].query_path, queries[j]) == 0)
	break;
    }
    assert(j!=query_ctr);

    char queryraganame[500], dum[500], dum2[500];
    find_dir_in_path(spot_vec_ptr->vec[i].query_path, dum);
    find_filename_in_path(dum, dum2);
    find_dir_in_path(dum, dum);
    find_filename_in_path(dum, queryraganame);

    char refraganame[500], dump[500], dump2[500];
    find_dir_in_path(spot_vec_ptr->vec[i].reference_path, dump);
    //find_filename_in_path(dum, dum2);
    //find_dir_in_path(dum, dum);
    find_filename_in_path(dump, refraganame);
      
    char outfilename[500];
    strcpy(dump, argv[2]);
    strcat(dump, "/");
    strcat(dump, queryraganame);
    strcat(dump, "/");
    sprintf(dump, "%s%d_", dump, j+1);
    strcat(dump, refraganame);
    strcat(dump, "_");
    sprintf(dump2, "%lf", spot_vec_ptr->vec[i].best_member_ptr->score);
    strcat(dump, dump2);
    strcat(dump, ".wav");
    strcpy(outfilename, dump);

    //for ref_wave
    find_filename_in_path_without_extn(spot_vec_ptr->vec[i].reference_path, dump);
    find_dir_in_path(spot_vec_ptr->vec[i].reference_path, dump2);
    find_dir_in_path(dump2, dump2);
    find_dir_in_path(dump2, dump2);
    strcat(dump2, "/");
    strcat(dump2, "recordings/");
    strcat(dump2, refraganame);
    strcat(dump2, "/");
    strcat(dump2, dump);
    strcat(dump2, ".wav");
    
      
    int from =  spot_vec_ptr->vec[i].best_member_ptr->start * 196;
    int to =  spot_vec_ptr->vec[i].best_member_ptr->end * 196;
    fprintf(ch_wave_out, "ch_wave -from %d -to %d -o \"%s\" \"%s\"\n", from, to, outfilename, dump2);      
      
  }

  fclose(ch_wave_out);
  return 0;
}
