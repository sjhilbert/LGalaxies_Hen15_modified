#include <stdlib.h>

/*------------------------------------------------------------------------------
  Random double in interval [0., 1.)
*/
double drandom()
{
    return((double)random() / ((double)RAND_MAX + 1.));
}
