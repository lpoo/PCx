/* hash table definitions
 *
 * PCx beta-2.0  10/31/96.
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Steve Wright.
 *
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#ifndef HashFile
#define HashFile

typedef struct node *ListPtr;

typedef struct node {
  int     index;
  char   *entry;
  ListPtr next;
} List;

typedef struct {
  ListPtr *list;
  int      size;
} HashTable;

HashTable  *NewHashTable();

#endif
