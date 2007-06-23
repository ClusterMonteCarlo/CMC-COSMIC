#include "cmc.h"
#include "cmc_vars.h"
#ifdef DEBUGGING

void load_id_table(GHashTable* ids, char *filename) {
  FILE * id_file;
  long id, i;
  gpointer id_pointer;
  int read_ok;
  
  id_file= fopen(filename, "r");
  if (id_file== NULL) {
    printf("%s does not exist!\n", filename);
    exit(1);
  }
  read_ok=fscanf(id_file, "%ld\n", &id);
  id_array= g_array_new(FALSE, FALSE, sizeof(long));
  while (read_ok>0) {
    g_array_append_val(id_array, id);
    read_ok= fscanf(id_file, "%ld\n", &id);
  };
  for (i=0; i< id_array->len; i++) {
    id_pointer= &g_array_index(id_array, long, i);
    g_hash_table_insert(ids, id_pointer, id_pointer);
  }
}

#endif
