#ifndef _SUGGEST
#define _SUGGEST

double get_mu (BIGNUM num_hit, double prob);
double get_sigma (BIGNUM num_hit, double prob);
double get_evalue (BIGNUM number, double mu, double sigma);
int get_suggestion(struct bloomstat *stats, double err_rate);
int kmer_suggestion (BIGCAST size);
float mco_suggestion (int k_mer);
BIGNUM find_close_prime (BIGNUM m);
int is_prime (BIGNUM m);
#endif
