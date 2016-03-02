#include "types.h"
/* this will parse the spot vector and prepare some numbers of list
 * equal to the number of unique motifs queried. The file will contain the
 * info of the matches in below given format. This file is to be processed 
 * further to get the desired result.
 * <query_file_name> <ref_file_name> <start_sample> <end_sample> <score>
 */

//declaring functions
int spot_comp_func(const void* a, const void* b);



// my main
int main(int argc, char* argv[]) {
  
  if (argc != 1 +1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: %s <result_file>", argv[0]);
    exit(-1);
  }
  

  MemberVectorPTR member_vec_ptr = create_MemberVectorPTR();
  alloc_members_MemberVectorPTR(member_vec_ptr, 500);
  
  // not sure whether the read binary mode works.
  FILEPTR fin;
  if ((fin = fopen(argv[1], "rb")) == NULL) {
    printf("ERROR: Cannot open file %s\n", argv[1]);
    exit(-1);
  }

  // reading the result into spot vector
  printf("Reading spot_vector...\n");
  scanf_MemberVectorPTR(fin, member_vec_ptr);
  printf("Reading over, %d members read\n", member_vec_ptr->size);
  // closing the file
  fclose(fin);

  // loop over and do something
  int i;  
  printf("Finding details of hits in each ref...\n");

  char  outfilename[200];
  int r = 0;
  int last_start = 99999999;
  FILEPTR fout = fopen("zzzzzmm----dummy-----mmmzzzz", "r");
  for (i = 0; i < member_vec_ptr -> size; i++) {
    if (member_vec_ptr->vec[i].end < last_start) {
      r = r+1;
      if (fout!=NULL) 
	fclose(fout);
      strcpy(outfilename, argv[1]);
      strcat(outfilename, ".parsed.");
      char dum[10];
      sprintf(dum, "%d", r);
      strcat(outfilename, dum);
      fout = fopen(outfilename, "w");
      if (fout!=NULL) printf("%s\n",outfilename);
    }
    

    fprintf(fout, "%d %d %f %f %f %f %f\n",
	    member_vec_ptr -> vec[i].start,
	    member_vec_ptr -> vec[i].end,
	    member_vec_ptr -> vec[i].score,
            member_vec_ptr -> vec[i].cost,
	    member_vec_ptr -> vec[i].cost_actual,
	    member_vec_ptr -> vec[i].war,
	    member_vec_ptr -> vec[i].waq
   	    );
    last_start = member_vec_ptr->vec[i].start;
  }


  //destroing all the fps
  free_members_MemberVectorPTR(member_vec_ptr);
  destroy_MemberVectorPTR(member_vec_ptr);

  return 0;
}


// define functions
int spot_comp_func(const void* a, const void* b) {
  if (((Spot*) a) -> best_member_ptr -> score == ((Spot*) b) -> best_member_ptr -> score) 
    return 0;
  else if (((Spot*) a) -> best_member_ptr -> score > ((Spot*) b) -> best_member_ptr -> score) 
    return -1; // this is done for desceding order
  else 
    return 1; // this is done for descending order
}
