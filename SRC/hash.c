/* hash tables
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>
#include "memory.h"

/*

#define NUMPRIMES      30

static int      prime[NUMPRIMES] = {29, 229, 883, 1671, 2791, 4801, 8629,
   15289, 25303, 34843, 65269, 99709, 129403, 147673, 166669, 201403,
   222163, 242729, 261431, 303491, 320237, 402761, 501131, 602309,
 701507, 800999, 900551, 1000619, 1100837, 1200359};
*/

// new list of  primes compiled 5/11/04 - SJW
// based on table from Knuth, v2, 3rd ed, p.408 (courtesy of Eric Bach)

#define NUMPRIMES 18
static int prime[NUMPRIMES] = {29, 229, 883, 1669, 2791, 4801, 8629,
			       15289, 32749, 65521, 131071, 262139,
			       524287, 1048573, 2097143, 4194301, 
			       8388593, 16777213};


typedef struct node *ListPtr;

typedef struct node {
  int             index;
  char           *entry;
  ListPtr         next;
}               List;

typedef struct {
  ListPtr        *list;
  int             size;
}               HashTable;

HashTable      *NewHashTable();

HashTable      *NewHashTable(size)
  int             size;
{
  int             i;
  HashTable      *table;

  /* allocate a hash table with size equal to a prime greater than size */
  table = (HashTable *) Malloc(sizeof(HashTable), "table");

  if (size > prime[NUMPRIMES - 1]) {
    printf("The size %d requested for the hash table exceeds the max, which is %d.\n",
	   size, prime[NUMPRIMES-1]);
	printf("Setting table size to the max allowed\n");
	// printf("The size requested for the hash table is too large: %d.\n", size);
	//    printf("Either add larger primes to the file 'hash.c' or request\n");
	//    printf("a smaller table.\n");
	//    OutOfSpace();
	size = prime[NUMPRIMES-1];
  } else {
    for (i = 0; i < NUMPRIMES; i++)
      if (size < prime[i]) {
	size = prime[i];
	break;
      }
  }
  table->size = size;
  table->list = (ListPtr *) Malloc(size * sizeof(ListPtr), "table->list");
  for (i = 0; i < size; i++)
    table->list[i] = NULL;

  return (table);
}

/*******************************************************************/

int             hash(table, string)
  HashTable      *table;
  char           *string;
{

  unsigned        number, scale;
  char           *s;

  /* Based on the size of the hash table, "hash" converts the string to a
   * number for indexing into the table.  */

  /* 0.618.... = (sqrt(5) - 1) / 2   based on Knuth, v.3, p. 510 */

  scale = (unsigned) (0.6180339887 * table->size);

  number = 0;
  for (s = string; *s != '\0'; s++)
    number = scale * number + *s;
  return (number % table->size);
}

/*******************************************************************/

int             GetIndex(table, name)
  HashTable      *table;
  char           *name;
{
  /* Given a name, go through the hash table (down the linked list if
   * necessary) and find the index for the name.  */

  List           *ptr;
  int             index, match, i;

  /* lookup entry */
  i = hash(table, name);

  match = -1;
  for (ptr = table->list[i]; ptr != NULL; ptr = ptr->next)
    if (strcmp(ptr->entry, name) == 0) {
      match = ptr->index;
      break;
    }
  return (match);
}

/* Insert makes an entry in the hash table.  The entry is indexed on name and
 * also stores an index value.  If name has already been entered into the
 * table, the routine returns a value 1.  */

Insert(table, name, index)
  HashTable      *table;
  char           *name;
  int             index;
{
  List           *ptr;
  int             i;

  /* lookup entry */

  i = hash(table, name);

  for (ptr = table->list[i]; ptr != NULL; ptr = ptr->next)
    if (strcmp(ptr->entry, name) == 0)
      break;

  if (ptr == NULL) {		/* no entry with "name" was found */
    ptr = (List *) Malloc(sizeof(List), "ptr in Hash Table");
    ptr->entry = StrDup(name, "entry");
    ptr->index = index;

    /* put this entry first in the list */
    ptr->next = table->list[i];
    table->list[i] = ptr;
    return (0);			/* normal */
  } else
    return (1);			/* name was already found in table */
}

int PrintHashTable(table)
  HashTable      *table;
{
  int             i;
  ListPtr         ptr;

  for (i = 0; i < table->size; i++) {
    printf("%d:\n", i);
    for (ptr = table->list[i]; ptr != NULL; ptr = ptr->next)
      printf(" %d '%s'\n", ptr->index, ptr->entry);
  }
  return 0;
}

int DeleteHashTable(table)
  HashTable      *table;
{
  int             i;
  ListPtr         ptr;

  for (i = 0; i < table->size; i++) {
    ptr = table->list[i];
    while (ptr != NULL) {
      ListPtr next = ptr->next;
      Free((char *) ptr->entry);
      Free((char *) ptr);
      ptr = next;
    }
  }
  Free((char *) table->list);
  Free((char *) table);
  return 0;
}

int PrintHashTableStats(table)
  HashTable      *table;
{
  int             i, count, max = 0;
  ListPtr         ptr;

  for (i = 0; i < table->size; i++) {
    count = 0;
    for (ptr = table->list[i]; ptr != NULL; ptr = ptr->next)
      count++;

    if (count > max)
      max = count;
  }

  printf("Hash Table max size = %d\n", max);
  return 0;
}
