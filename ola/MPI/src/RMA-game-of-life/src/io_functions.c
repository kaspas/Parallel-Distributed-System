#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include <game-of-life.h>
#define NTHREADS 8
/* print the life board */

void print (int8_t *board, int N,int M) {
  int   i, j;
  #pragma omp parallel for private(i,j) num_threads(NTHREADS) 
  /* for each row */
  for (j=0; j<M; j++) {

    /* print each column position... */
    for (i=0; i<N; i++) {
      printf ("%c", Board(i,j) ? 'x' : ' ');
    }

    /* followed by a carriage return */
    printf ("\n");
  }
}



/* display the table with delay and clear console */

void display_table(int8_t *board, int N,int M) {
  print (board, N,M);
  usleep(100000);  
  /* clear the screen using VT100 escape codes */
  puts ("\033[H\033[J");
}
