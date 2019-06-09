#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MT.h"

double generate_random_trial(double xmin, double xmax);

int main(){
  double xmin = 1.;
  double xmax = 10.;
  
  int i;
  for (i=0;i<10;i++){
    printf("%12.9e \n",generate_random_trial(xmin,xmax));
  }

  return 0;
}

double generate_random_trial(double xmin, double xmax)
{
  return xmin+genrand_real3()*(xmax-xmin);
}
