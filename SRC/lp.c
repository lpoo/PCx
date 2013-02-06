/* print LPtype data structure
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "main.h"

int PrintLP(LP)
  LPtype   *LP;
{
  
  int row, col, ent;

  printf("LP: Cols = %d Rows = %d\n", LP->Cols, LP->Rows);

  printf("Objective:\n");
  for (col = 0; col < LP->Cols; col++)
    printf(" c[%d] = %f  (%f)\n", col, LP->c[col], LP->UpBound[col]);

  printf("\nMatrix A:\n");
  for (col = 0; col < LP->Cols; col++) {
    printf(" Col %d:\n", col);
    for (ent = LP->A.pBeginRow[col]-1; ent <= LP->A.pEndRow[col]-1; ent++)
      printf("  Row %d = %f\n", LP->A.Row[ent]-1, LP->A.Value[ent]);
  }

  printf("\nUpbound:\n");
  printf(" NumberBounds = %d\n", LP->NumberBounds);
  for (ent = 0; ent < LP->NumberBounds; ent++) {
    col = LP->BoundIndex[ent];
    printf(" Bound[%d] = %f\n", col, LP->UpBound[col]);
  }

  printf("\nRHS:\n");
  for (row = 0; row < LP->Rows; row++)
    if (LP->b[row] != 0.0)
      printf(" b[%d] = %f\n", row, LP->b[row]);

  return 0;
}
