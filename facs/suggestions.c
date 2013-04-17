#define MB 1048576

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "bloom.h"
#include "prob.h"
#include "suggestions.h"

double
get_mu (BIGNUM num_hit, double prob)
{
  return ((double) num_hit) * prob;
}

double
get_sigma (BIGNUM num_hit, double prob)
{
  return (double) num_hit *prob * (1 - prob);
}

double
get_evalue (BIGNUM number, double mu, double sigma)
{
  return cdf (number, mu, sigma);
}

int
get_suggestion(struct bloomstat *stats, BIGNUM capacity, double err_rate)
{
  stats->capacity = capacity;
  stats->e = err_rate;

  /*
   * Given the desired capacity and error rate, calculate the appropriate value
   * for number of hash functions and size of array
   */

  /* assuming perfect number of cells, k directly depends on e */
  stats->ideal_hashes = (int) log (stats->e) / log (0.5);
  stats->elements = find_close_prime ((BIGNUM) 13 * stats->capacity
                   *
                   (BIGNUM) stats->ideal_hashes / (BIGNUM) 9);
  /*
     recalculate k with the actual m, not the ideal 
     wouldn't need to if it wasn't prime, but that causes problems
     for hash algs
   */

  stats->ideal_hashes = 9 * stats->elements / (13 * stats->capacity);
  
  return 0;
}

int
kmer_suggestion (BIGCAST size)
{
  if (size < 10 * MB)
    {
      //bl->k_mer = 15;
      //bl->mcf = 0.4;
      return 15;
    }
  else if (size < 20 * MB)
    {
      //bl->k_mer = 17;
      //bl->mcf = 0.4;
      return 16;
    }
  else if (size < 50 * MB)
    {
      //bl->k_mer = 17;
      //bl->mcf = 0.4;
      return 17;
    }
  else if (size < 100 * MB)
    {
      //bl->k_mer = 18;
      //bl->mcf = 0.3;
      return 18;
    }
  else if (size < 500 * MB)
    {
      return 19;
    }
  else
    {
      //bl->k_mer = 20;
      //bl->mcf = 0.3;
      return 20;
    }
}

float
mco_suggestion (int k_mer)
{
  switch (k_mer)
    {
    case 15:
      return 0.4;
    case 16:
      return 0.3;
    case 17:
      return 0.3;
    case 18:
      return 0.3;
    case 19:
      return 0.4;
    case 20:
      return 0.3;
    default:
      return 0.4;
    }
}

BIGNUM
find_close_prime (BIGNUM m)
{
  if ((m % 2) == 0)
    m += 1;

  while (!is_prime (m))
    {
      m += 2;
    }
  return m;
}

int
is_prime (BIGNUM m)
{
  BIGNUM a = (BIGNUM) sqrtl ((long double) m);
  BIGNUM currval = 3;
  if (m % 2 == 0 && m != 2)
    return 0;
  while (m % currval != 0 && currval < a)
    {
      if (m % (currval + 2) == 0)
	return 0;
      if (m % (currval + 4) == 0)
	return 0;
      currval += 8;
    }
  return (int) m % currval;
}
