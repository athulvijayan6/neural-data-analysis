/*
 TO test yamlwrapper
 */


#include "yamlwrapper.h"
#include "types.h"


int main(int argc, char * argv[]) {
  
  if (argc!= 2+1) {
    printf("USAGE: yamlwrappertest <yamlfile.yaml> <filetype>\n");
    exit(-1);
  }  
  

  char * yaml_filename = argv[1];
  int type = atoi(argv[2]);
  

  FILE * yaml_fp = fopen(yaml_filename, "r");

  
  printf("I am in\n");

  if (type == 1){
    printf("I am type 1\n");
    IntVector * yaml_data =   yaml_read_file_type1(yaml_fp);
    int i;
    for(i = 0; i < MAX_LIM; i++) {
      if (yaml_data[i].size > 0) { //ignoring keys with 0 size
	printf("%d:\t", i);
        int j;
        for (j = 0; j < yaml_data[i].size; j++) {
	  printf("%d ", yaml_data[i].vec[j]);
        }
        printf("\n");
      }
    }

  } else {
    printf("I am type 2\n");
    StringVector * yaml_data = yaml_read_file_type2(yaml_fp);
    int i;
    for(i = 0; i < MAX_LIM; i++) {
      if (yaml_data[i].size > 0) { //ignoring keys with 0 size
	printf("%d:\n", i);
        int j;
        for (j = 0; j < yaml_data[i].size; j++) {
	  printf("%s\n", yaml_data[i].vec[j].data);
        }
        printf("\n");
      }
    }


  }
    

  return 0;
}
