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

#ifndef YAMLWRAPPER_H
#define YAMLWRAPPER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "types.h"

#define MAX_LIM 500 


IntVector* yaml_read_file_type1(FILE* yaml_file);
StringVector* yaml_read_file_type2(FILE* yaml_file);


#endif


