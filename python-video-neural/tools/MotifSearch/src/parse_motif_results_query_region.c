#include "types.h"
/* this will parse the spot vector and prepare some numbers of list
 * equal to the number of unique motifs queried. The file will contain the
 * info of the matches in below given format. This file is to be processed 
 * further to get the desired result.
 * <query_file_name> <ref_file_name> <start_sample> <end_sample> <score>
 */

//declaring functions
int spot_comp_func(const void* a, const void* b);
int doub_cmp_fun (const void *pa, const void *pb );

// my main
int main(int argc, char* argv[]) {
  
  if (argc != 2 +1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: %s <result_file> <top_n>", argv[0]);
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
  

  // loop over spot vector and get the unique refs.
  printf("Finding unique refs and queries... \n");
  char refs[100][500];
  char queries[100][500];
  int ref_ctr = 0;  // ref counter
  int query_ctr =0; // query counter 
  {
    // for ref
    int i;
    for (i = 0; i < spot_vec_ptr -> size; i++) {
      int j = 0;
      for (j = 0; j < ref_ctr; j++)
	if (strcmp(spot_vec_ptr->vec[i].reference_path, refs[j]) == 0)
	  break;
      // if current ref is not found in queries then add it to the list.
      if (j == ref_ctr) {
	strcpy(refs[ref_ctr], spot_vec_ptr->vec[i].reference_path);
	++ref_ctr;
      }
    }
    printf("%d unique refs found namely\n", ref_ctr);
    char dummy[100];
    strcpy(dummy, argv[1]);
    strcat(dummy, "_ref.map");
    FILEPTR refmap = fopen(dummy, "w");
    for (i = 0; i <  ref_ctr; i++){
      printf("%d)\t%s\n", i, refs[i]);
      fprintf(refmap, "%d %s\n", i, refs[i]);
    }


    // for query 

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

    printf("%d unique refs found namely\n", ref_ctr);
    //char dummy[100];
    strcpy(dummy, argv[1]);
    strcat(dummy, "_query.map");
    FILEPTR querymap = fopen(dummy, "w");

 
    for (i = 0; i <  query_ctr; i++){ 
      printf("%d)\t%s\n", i, queries[i]); 
      fprintf(querymap, "%d %s\n", i, queries[i]); 
    }

    fclose(querymap);
    fclose(refmap);
    
  }


  // loop over and do something
  int i,j,k;  

  /* // first creating file pointers for each of the query output. */
  /* FILEPTR fp[query_ctr]; */
  /* int fp_ctr[query_ctr]; */
  /* for (i = 0; i < query_ctr; i++){ */
  /*   char dummy[500]; char dummy2[100]; */
  /*   strcpy(dummy, argv[4]); */
  /*   strcat(dummy, "_motif_"); */
  /*   sprintf(dummy2, "%d", i+1); */
  /*   strcat(dummy, dummy2); */
  /*   fp[i] = fopen(dummy, "wb"); */
  /*   fp_ctr[i] = 0; */
  /* } */

  FILEPTR fout, fout2;
  {
    char dummy[500];
    strcpy(dummy, argv[1]);
    strcat(dummy, ".parsed");
    fout = fopen(dummy, "w");

    strcpy(dummy, argv[1]);
    strcat(dummy, ".spotdump");
    fout2 = fopen(dummy, "w");
  }

  
  

  for (i = 0; i < ref_ctr; i++) {
    // for the five ragas
    double bhairavi[3000][2]; int br1 = 0;
    double kamboji[3000][2]; int kam1 = 0;
    double kalyani[3000][2]; int kal1 = 0;
    double shankarabharanam[3000][2]; int shank1 = 0;
    double varali[3000][2]; int var1 = 0;

    for (j = 0; j < spot_vec_ptr -> size; j++) {
      
      if (spot_vec_ptr->vec[j].group_end == 0) continue;
      char raganame[500], dum[500], dum2[500];
      find_dir_in_path(spot_vec_ptr->vec[j].query_path, dum);
      find_filename_in_path(dum, dum2);
      find_dir_in_path(dum, dum);
      find_filename_in_path(dum, raganame);
      
      // some filtering
      if (strcmp(raganame, "Bhairavi") == 0 && strcmp(dum2, "12")==0) continue;
      if (strcmp(raganame, "Bhairavi") == 0 && strcmp(dum2, "13")==0) continue;
      if (strcmp(raganame, "Shankarabharanam") == 0 && strcmp(dum2, "6") == 0) continue;
      if (strcmp(raganame, "Shankarabharanam") == 0 && strcmp(dum2, "7") == 0) continue;
      if (strcmp(raganame, "Shankarabharanam") == 0 && strcmp(dum2, "8") == 0) continue;
      //      if (strcmp(raganame, "Shankarabharanam") == 0 && strcmp(dum2, "1") == 0) continue;

      //      if it is ref
      if (strcmp(spot_vec_ptr->vec[j].reference_path, refs[i])==0) {
	//	then get the raga info of query
	//	printf("----%s---\n", dum);
	//	printf("Raga is: %s\n", raganame);
	
	if (strcmp(raganame, "Bhairavi") == 0){
	  bhairavi[br1][0] = j;
          bhairavi[br1][1] = spot_vec_ptr->vec[j].best_member_ptr->score;
	  //	  printf("%f\n",bhairavi[br1][1]);
	  ++br1;
	} else if (strcmp(raganame, "Kamboji") == 0){
	  kamboji[kam1][0] = j;
          kamboji[kam1][1] = spot_vec_ptr->vec[j].best_member_ptr->score;
	  ++kam1;
	} else if (strcmp(raganame, "Kalyani") == 0){
	  kalyani[kal1][0] = j;
          kalyani[kal1][1] = spot_vec_ptr->vec[j].best_member_ptr->score;
	  ++kal1;
	} else if (strcmp(raganame, "Shankarabharanam") == 0){
	  shankarabharanam[shank1][0] = j;
          shankarabharanam[shank1][1] = spot_vec_ptr->vec[j].best_member_ptr->score;
	  ++shank1;
	} else if (strcmp(raganame, "Varali") == 0){
	  varali[var1][0] = j;
          varali[var1][1] = spot_vec_ptr->vec[j].best_member_ptr->score;
	  ++var1;
	}
      }
    }

    qsort(bhairavi, br1, sizeof(bhairavi[0]), doub_cmp_fun);
    qsort(kamboji, kam1, sizeof(kamboji[0]), doub_cmp_fun);
    qsort(kalyani, kal1, sizeof(kalyani[0]), doub_cmp_fun);
    qsort(varali, var1, sizeof(varali[0]), doub_cmp_fun);
    qsort(shankarabharanam, shank1, sizeof(shankarabharanam[0]), doub_cmp_fun);
    
    for (k = 0; k < top_n; k++) {
      if (k >= br1) {
	bhairavi[k][1] = 0;
	bhairavi[k][0] = -1;
      }
      fprintf(fout, "%f ", bhairavi[k][1]);
      if (bhairavi[k][0] > -1)
	printf_SpotPTR(fout2, &(spot_vec_ptr->vec[(int)bhairavi[k][0]]));
    }
    for (k = 0; k < top_n; k++) {
       if (k >= kam1) {
	kamboji[k][1] = 0;
	kamboji[k][0] = -1;
      }
      fprintf(fout, "%f ", kamboji[k][1]);
      if (kamboji[k][0] > -1)
	printf_SpotPTR(fout2, &(spot_vec_ptr->vec[(int)kamboji[k][0]]));
    }
    for (k = 0; k < top_n; k++) {
      if (k >= kal1) {
	kalyani[k][1] = 0;
	kalyani[k][0] = -1;
      }
      fprintf(fout, "%f ", kalyani[k][1]);
      if (kalyani[k][0] > -1)
	printf_SpotPTR(fout2, &(spot_vec_ptr->vec[(int)kalyani[k][0]]));
    }

    for (k = 0; k < top_n; k++) {
      if (k >= shank1) {
	shankarabharanam[k][1] = 0;
	shankarabharanam[k][0] = -1;
      }

      fprintf(fout, "%f ", shankarabharanam[k][1]);
      if (shankarabharanam[k][0] > -1)
	printf_SpotPTR(fout2, &(spot_vec_ptr->vec[(int)shankarabharanam[k][0]]));
    }
    for (k = 0; k < top_n; k++) {
      if (k >= var1) {
	varali[k][1] = 0;
	varali[k][0] = -1;
      }
      
      fprintf(fout, "%f ", varali[k][1]);
      if (varali[k][0] > -1)
	printf_SpotPTR(fout2, &(spot_vec_ptr->vec[(int)varali[k][0]]));
    }
    fprintf(fout, "\n");
  }

  fclose(fout2);
  fclose(fout);
  // free and destroy
  free_members_SpotVectorPTR(spot_vec_ptr);
  destroy_SpotVectorPTR(spot_vec_ptr);

  return 0;
}



// define functions
int spot_comp_func(const void* a, const void* b) {
  if (strcmp(((Spot*) a) -> query_path, ((Spot*) b) -> query_path) == 0) {
    if (((Spot*) a) -> group_end_q > ((Spot*) b) -> group_start_q  && ((Spot*) a) -> group_start_q < ((Spot*) b) -> group_end_q) {
      if (((Spot*) a) -> best_member_ptr -> score == ((Spot*) b) -> best_member_ptr -> score) 
	return 0;
      else if (((Spot*) a) -> best_member_ptr -> score > ((Spot*) b) -> best_member_ptr -> score) 
	return -1; // this is done for desceding order
      else 
	return 1; // this is done for descending order
    } else  if (((Spot*) a) -> group_end_q < ((Spot*) b) -> group_start_q)
      return -1;
    else 
      return 1;
  } else {
    return strcmp(((Spot*) a) -> query_path, ((Spot*) b) -> query_path);
  }
}


int doub_cmp_fun (const void *pa, const void *pb )
{
    const double *a = pa;
    const double *b = pb;
    if (a[1] <= b[1]) 
        return 1;
    return -1;
}
