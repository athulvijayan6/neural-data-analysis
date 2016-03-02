/* Description: Contains wrapper modules to read yaml file and use it as 
 * a dictionary as it's done in python.
 * 
 * NOTE: Currently this version of code only work for one level yaml 
 * dictionaries of type- 
 * #################################
 * ## yaml file example ############ type1 int2int and type2 int2string
 * <key 1 int>:
 * - <value 1 int or string>
 * - <... int or string>
 * - <value 2 int string>
 * <...>
 * <key n int>:
 * - <value 1 int or string>
 * - <... int or string>
 * - <value 2 int or string>
 * ## yaml file example ends here###
 * #################################
 * THIS CODE WILL ONLY WORK FOR THIS SPECIFIC TYPE OF YAML FILE.
 *
 * Creation Date: 26 March 2014
 * Author: Shrey Dutta
 */


#include "yamlwrapper.h"


int is_key(const char line[MAX_LIM]) {
  int i =0;
  for (i=0; i < MAX_LIM; ++i){
    if (line[i]==':')
      return 1;
    if (line[i]=='\0')
      return 0;
  }
  return 0;
}

int is_value(const char line[MAX_LIM]) {

  if (line[0] == '-')
    return 1;

  return 0;
}

void get_key(const char line[MAX_LIM], char* key) {
  int i =0;
  for (i=0; i < MAX_LIM; i++){
    if (line[i]==':') {
      key[i] = '\0';
      break;
    }
    key[i] = line[i];     
  }
}

void get_value(const char line[MAX_LIM], char* value) {
  int i =2;
  while(i < MAX_LIM){
    if (line[i]=='\0') {
      value[i-2] = line[i];
     
      // just to remove the the next line character
      // removes it from the string
      // helps in type 2
      if (value[i-3] == '\n') 
	value[i-3] = '\0';
      
      break;
    }

    value[i-2] = line[i];

    ++i;     
  }
}





IntVector* yaml_read_file_type1(FILE* yaml_file) {
  
  char line[MAX_LIM];
  int i;

  IntVector* ret = (IntVector*) malloc(MAX_LIM * sizeof(IntVector));

  for (i=0; i < MAX_LIM; i++)
    alloc_members_IntVectorPTR(&ret[i], MAX_LIM);
  
  int i_key;
  
  while (fgets(line, 500, yaml_file) != NULL) {
    if (is_key(line)) {
      char * key = (char*) malloc(MAX_LIM*sizeof(char));
      get_key(line, key);
      i_key = atoi(key);
      free(key);
    } else if (is_value(line)){
      char * value = (char*) malloc(MAX_LIM* sizeof(char));
      get_value(line, value);
      int i_value = atoi(value);
      ret[i_key].vec[ret[i_key].size] = i_value;
      ret[i_key].size += 1;
      free(value);
    } else {
      printf("Error[Shrey]: Neither Key Nor Value\n");
      exit(-1);
    }
  }
  return ret;
}


//TODO: compile them in a make file and test if they are working properly

StringVector* yaml_read_file_type2(FILE* yaml_file) {
  
  char line[MAX_LIM];
  int i;

  StringVector* ret = (StringVector*) malloc(MAX_LIM * sizeof(StringVector));

  for (i=0; i < MAX_LIM; i++)
    alloc_members_StringVectorPTR(&ret[i], MAX_LIM);
  
  int i_key;
  
  while (fgets(line, MAX_LIM, yaml_file) != NULL) {
    if (is_key(line)) {
      char * key = (char*) malloc(MAX_LIM*sizeof(char));
      get_key(line, key);
      i_key = atoi(key);
      free(key);
    } else if (is_value(line)){
      char * value = (char*) malloc(MAX_LIM* sizeof(char));
      get_value(line, value);
      
      strcpy(ret[i_key].vec[ret[i_key].size].data, value);
      ret[i_key].size += 1;
      free(value);
    } else {
      printf("Error[Shrey]: Neither Key Nor Value\n");
      exit(-1);
    }
  }
  return ret;
}


void yaml_delete_type1_ds(IntVector * x) {
  int i;
  for (i = 0; i < MAX_LIM; i++) 
    free_members_IntVectorPTR(&x[i]);
  free(x);
  
}


void yaml_delete_type2_ds(StringVector * x) {
  int i;
  for (i = 0; i < MAX_LIM; i++) 
    free_members_StringVectorPTR(&x[i]);
  free(x);
  
}
