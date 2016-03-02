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
int doub_cmp_fun2 (const void *pa, const void *pb );

// my main
int main(int argc, char* argv[]) {
  
  if (argc != 3 +1) {
    printf("!!!!!!Error!!!!!\n");
    printf("Usage: %s <result_file> <top_n> <avg-flag(1/0)>", argv[0]);
    exit(-1);
  }
  
  int top_n = atoi(argv[2]);
  int avg_flag = atoi(argv[3]);
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
  char refs[200][500];
  char queries[200][500];
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
  int j,k;  

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
    strcat(dummy, "_");
    strcat(dummy, argv[2]);
    strcat(dummy, "_");
    strcat(dummy, argv[3]);
    strcat(dummy, ".parsedv3");
    fout = fopen(dummy, "w");

    strcpy(dummy, argv[1]);
    strcat(dummy, "_");
    strcat(dummy, argv[2]);
    strcat(dummy, "_");
    strcat(dummy, argv[3]);
    strcat(dummy, ".spotdumpv3");
    fout2 = fopen(dummy, "w");
  }

  
  

  int qi;
  for (qi = 0; qi < query_ctr; qi++) {
    // for the five ragas

    double bhairavi[3000]; int br1 = 0;
    double kamboji[3000]; int kam1 = 0;
    double kalyani[3000]; int kal1 = 0;
    double shankarabharanam[3000]; int shank1 = 0;
    double varali[3000]; int var1 = 0;

    
    int ri;
    double scorelist[5000]; int slisti = 0;
    for (ri = 0; ri < ref_ctr; ri++) {
      
      double tempragascore[3000][2]; int tempi = 0;      
      
      char refraganame[500], dum[500], dum2[500];
      find_dir_in_path(refs[ri], dum);
      //find_filename_in_path(dum, dum2);
      //find_dir_in_path(dum, dum);
      find_filename_in_path(dum, refraganame);
      
      char queryraganame[500]; //, dum[500], dum2[500];
      find_dir_in_path(queries[qi], dum);
      find_filename_in_path(dum, dum2);
      find_dir_in_path(dum, dum);
      find_filename_in_path(dum, queryraganame);

      for (j = 0; j < spot_vec_ptr -> size; j++) {
      
	if (spot_vec_ptr->vec[j].group_end == 0) continue;


	// some filtering
	// if (strcmp(raganame, "Bhairavi") == 0 && strcmp(dum2, "12")==0) continue;
	//if (strcmp(raganame, "Bhairavi") == 0 && strcmp(dum2, "13")==0) continue;
	//if (strcmp(raganame, "Shankarabharanam") == 0 && strcmp(dum2, "6") == 0) continue;
	//if (strcmp(raganame, "Shankarabharanam") == 0 && strcmp(dum2, "7") == 0) continue;
	//if (strcmp(raganame, "Shankarabharanam") == 0 && strcmp(dum2, "8") == 0) continue;
	//      if (strcmp(raganame, "Shankarabharanam") == 0 && strcmp(dum2, "1") == 0) continue;

	//      if it is ref
	if (strcmp(spot_vec_ptr->vec[j].query_path, queries[qi])==0 && strcmp(spot_vec_ptr->vec[j].reference_path, refs[ri])==0) {
	  //	then get the raga info of query
	  //	printf("----%s---\n", dum);
	  //	  printf("Raga is: %s, %s\n", queryraganame, refraganame);
	
	  tempragascore[tempi][0] = j;
	  tempragascore[tempi][1] = spot_vec_ptr->vec[j].best_member_ptr->score;
	  ++tempi;
	}
      }
    
      qsort(tempragascore, tempi, sizeof(tempragascore[0]), doub_cmp_fun);
	
      for (k = 0; k < top_n; k++) {
	if (k >= tempi) {
	  tempragascore[k][1] = 0;
	  tempragascore[k][0] = 0;
	}
	if (tempragascore[k][0] > 0)
	  printf_SpotPTR(fout2, &(spot_vec_ptr->vec[(int)tempragascore[k][0]]));
	
	if (strcmp(refraganame, "Bhairavi") == 0){
	  bhairavi[br1++] = tempragascore[k][1];
	} else if (strcmp(refraganame, "Kamboji") == 0){
	  kamboji[kam1++] = tempragascore[k][1];
	} else if (strcmp(refraganame, "Kalyani") == 0){
	  kalyani[kal1++] = tempragascore[k][1];
	} else if (strcmp(refraganame, "Shankarabharanam") == 0){
	  shankarabharanam[shank1++] = tempragascore[k][1];
	} else if (strcmp(refraganame, "Varali") == 0){
	  varali[var1++] = tempragascore[k][1];
	}
      }
    }

    
    qsort(bhairavi, br1, sizeof(double), doub_cmp_fun2);
    qsort(kamboji, kam1, sizeof(double), doub_cmp_fun2);
    qsort(kalyani, kal1, sizeof(double), doub_cmp_fun2);
    qsort(shankarabharanam, shank1, sizeof(double), doub_cmp_fun2);
    qsort(varali, var1, sizeof(double), doub_cmp_fun2);
    
    /*
    int zi;
    printf("***********************************\n");
    for (zi = 0; zi < br1; zi++)
      printf("%f\n", bhairavi[zi]);
    printf("***********************************\n");
    */

    if (avg_flag) {
      // taking mean of the score list
      double tempmean = 0;
      int ti =0;
      for (ti = 0; ti <  top_n; ti++) {
	if (ti >= br1) {
	  bhairavi[ti] = 0.0;
	}
	tempmean += bhairavi[ti];
      }
      tempmean /= top_n;
   
      fprintf(fout, "%f ", tempmean);
    

      tempmean = 0;
      for (ti = 0; ti <  top_n; ti++){
	if (ti >= kam1) {
	  kamboji[ti] = 0.0;
	}
	tempmean += kamboji[ti];
      }
      tempmean /= top_n;
      fprintf(fout, "%f ", tempmean);

      tempmean = 0;
      for (ti = 0; ti <  top_n; ti++) {
	if (ti >= kal1) {
	  kalyani[ti] = 0.0;
	}
	tempmean += kalyani[ti];
      }
      tempmean /= top_n;
      fprintf(fout, "%f ", tempmean);

      tempmean = 0;
      for (ti = 0; ti <  top_n; ti++) {
	if (ti >= shank1) {
	  shankarabharanam[ti] = 0.0;
	}
	tempmean += shankarabharanam[ti];
      }
      tempmean /= top_n;
      fprintf(fout, "%f ", tempmean);

      tempmean = 0;
      for (ti = 0; ti <  top_n; ti++) {
	if (ti >= var1) {
	  varali[ti] = 0.0;
	}
	tempmean += varali[ti];
      }
      tempmean /= top_n;
      fprintf(fout, "%f ", tempmean);
    } else {
      int ti;
      for (ti = 0; ti < top_n; ti++) {
	fprintf(fout, "%f %f %f %f %f", bhairavi[ti], kamboji[ti], kalyani[ti], shankarabharanam[ti], varali[ti]);
      }
    }


    //////////////////////////////////////////////
    /*    double tempmean = 0;
    int ti =0;
    for (ti = 0; ti <  br1; ti++)
      tempmean = tempmean > bhairavi[ti] ? tempmean : bhairavi[ti];

    fprintf(fout, "%f ", tempmean);

    tempmean = 0;
    for (ti = 0; ti <  kam1; ti++)
      tempmean = tempmean >  kamboji[ti] ? tempmean :  kamboji[ti];

    fprintf(fout, "%f ", tempmean);

    tempmean = 0;
    for (ti = 0; ti <  kal1; ti++)
      tempmean = tempmean > kalyani[ti] ? tempmean : kalyani[ti];

    fprintf(fout, "%f ", tempmean);

    tempmean = 0;
    for (ti = 0; ti <  shank1; ti++)
      tempmean = tempmean > shankarabharanam[ti] ? tempmean : shankarabharanam[ti];

    fprintf(fout, "%f ", tempmean);

    tempmean = 0;
    for (ti = 0; ti <  var1; ti++)
      tempmean = tempmean > varali[ti] ? tempmean : varali[ti];

    fprintf(fout, "%f ", tempmean);
    //////////////////////////////////
    */

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
  const double *a = (double*)pa;
  const double *b = (double*)pb;
  if (a[1] <= b[1]) 
    return 1;
  return -1;
}


int doub_cmp_fun2 (const void *pa, const void *pb )
{
  const double *a = (double*)pa;
  const double *b = (double*)pb;
  if (*a <= *b) 
    return 1;
  return -1;
}
