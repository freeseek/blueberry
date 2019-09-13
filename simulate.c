/* The MIT License

   Copyright (C) 2019 Giulio Genovese

   Author: Giulio Genovese <giulio.genovese@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

// gcc -g -O2 simulate.c -lm -lgsl -lgslcblas -o simulate && ./simulate

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <gsl/gsl_randist.h> // it requires libgsl-dev to be installed

#define VERSION "2019-09-13"

#define HMM_MM 0
#define HMM_MF 1
#define HMM_FF 2
#define HMM_M  3
#define HMM_F  4

static const char *hmm_str[] = {"MM", "MF", "FF", "M", "F"};
static inline double sq(double x) { return x*x; }
static inline double log_mean_exp(double x, double y) { return x == y ? x : ( x > y ? x + log( 1 + expf(y-x) ) : y + log( 1 + expf(x-y) ) ) - M_LN2; }

/*********************************
 * CODE FOR BETA BINOMIAL        *
 *********************************/

typedef struct
{
    double alpha;
    double beta;
    double *gamma_alpha;
    double *gamma_beta;
    double *gamma_alpha_beta;
    double *log_gamma_alpha;
    double *log_gamma_beta;
    double *log_gamma_alpha_beta;
} beta_binom_t;

static beta_binom_t *beta_binom_init(double p, double rho, int n)
{
    beta_binom_t *self = (beta_binom_t *)calloc(1, sizeof(beta_binom_t));

    if ( rho == 0 )
    {
        self->alpha = p;
        self->beta = 1.0 - p;
    } else {
        double s = ( 1.0 - rho ) / rho;
        self->alpha = p * s;
        self->beta = ( 1.0 - p ) * s;
    }

    self->gamma_alpha = (double *)malloc(n * sizeof(double));
    self->gamma_beta = (double *)malloc(n * sizeof(double));
    self->gamma_alpha_beta = (double *)malloc(n * sizeof(double));
    self->gamma_alpha[0] = 1.0;
    self->gamma_beta[0] = 1.0;
    self->gamma_alpha_beta[0] = 1.0;

    self->log_gamma_alpha = (double *)malloc(n * sizeof(double));
    self->log_gamma_beta = (double *)malloc(n * sizeof(double));
    self->log_gamma_alpha_beta = (double *)malloc(n * sizeof(double));
    self->log_gamma_alpha[0] = 0.0;
    self->log_gamma_beta[0] = 0.0;
    self->log_gamma_alpha_beta[0] = 0.0;

    if ( rho == 0 ) // binomial distribution case (no overdispersion)
    {
        double log_alpha = log( p );
        double log_beta = log( 1.0 - p );

        for(int i=1; i<n; i++)
        {
            self->gamma_alpha[i] = self->gamma_alpha[i - 1] * ( p / (double)i );
            self->gamma_beta[i] = self->gamma_beta[i - 1] * ( ( 1.0 - p ) / (double)i );
            self->gamma_alpha_beta[i] = self->gamma_alpha_beta[i - 1] / (double)i;

            double log_i = log( i );
            self->log_gamma_alpha[i] = self->log_gamma_alpha[i - 1] + log_alpha - log_i;
            self->log_gamma_beta[i] = self->log_gamma_beta[i - 1] + log_beta - log_i;
            self->log_gamma_alpha_beta[i] = self->log_gamma_alpha_beta[i - 1] - log_i;
        }
    } else {
        for(int i=1; i<n; i++)
        {
            self->gamma_alpha[i] = self->gamma_alpha[i - 1] * ( (self->alpha + (double)i - 1.0) / (double)i );
            self->gamma_beta[i] = self->gamma_beta[i - 1] * ( (self->beta + (double)i - 1.0) / (double)i );
            self->gamma_alpha_beta[i] = self->gamma_alpha_beta[i - 1] * ( (self->alpha + self->beta + (double)i - 1.0) / (double)i );

            self->log_gamma_alpha[i] = self->log_gamma_alpha[i - 1] + log( (self->alpha + (double)i - 1.0) / (double)i );
            self->log_gamma_beta[i] = self->log_gamma_beta[i - 1] + log( (self->beta + (double)i - 1.0) / (double)i );
            self->log_gamma_alpha_beta[i] = self->log_gamma_alpha_beta[i - 1] + log( (self->alpha + self->beta + (double)i - 1.0) / (double)i );
        }
    }

    return self;
}

static void beta_binom_destroy(beta_binom_t *self)
{
    free(self->gamma_alpha);
    free(self->gamma_beta);
    free(self->gamma_alpha_beta);
    free(self->log_gamma_alpha);
    free(self->log_gamma_beta);
    free(self->log_gamma_alpha_beta);
    free(self);
}

static inline double dbetabinom(const beta_binom_t *self, int y, int size)
{
    return self->gamma_alpha[y] * self->gamma_beta[size - y] / self->gamma_alpha_beta[size];
}

static inline double log_dbetabinom(const beta_binom_t *self, int y, int size)
{
    return self->log_gamma_alpha[y] + self->log_gamma_beta[size - y] - self->log_gamma_alpha_beta[size];
}

/*********************************
 * CODE FOR DEBUG                *
 *********************************/

static void debug_beta_binom(const beta_binom_t *self, int size)
{
    for (int i=0; i<size; i++)
    {
        fprintf(stderr, "%d %d - d=%g log_d=%g\n", i, size - i, exp( log_dbetabinom(self, i, size) ), log_dbetabinom(self, i, size));
    }
}

static void debug_solution(const int *hidden_path,
                           const int *a,
                           const int *cov,
                           const double *log_trisomy_lkl,
                           const double *log_euploid_lkl,
                           const int *trisomy_path,
                           const int *euploid_path,
                           int size)
{
    for (int i=0; i<size; i++)
    {
        fprintf(stderr, "%d\t%s\t%d\t%d\tMM=%.4f\tMF=%.4f\tFF%.4f\tM=%.4f\tF=%.4f\t%s\t%s\n", i,
            hmm_str[hidden_path[i]], a[i], cov[i],
            log_trisomy_lkl[3 * i + 0],
            log_trisomy_lkl[3 * i + 1],
            log_trisomy_lkl[3 * i + 2],
            log_euploid_lkl[2 * i + 0],
            log_euploid_lkl[2 * i + 1],
            hmm_str[trisomy_path[i]],
            hmm_str[3 + euploid_path[i]]);
    }
}

/*********************************
 * CODE TO COMPUTE VITERBI PATH  *
 *********************************/

static void log_viterbi_run(const double *log_lkl,
                            const double *log_trans,
                            int T,
                            int N,
                            int8_t *path)
{
    if (T<1) return;
    double log_prb[N];
    double new_log_prb[N];
    // allocate memory necessary for running the algorithm
    int8_t *ptr = (int8_t *)malloc(N * (T-1) * sizeof(int8_t));

    // initialize and rescale the first state
    for (int i=0; i<N; i++) log_prb[i] = log_lkl[i];

    // compute best probabilities at each position
    for (int t=1; t<T; t++)
    {
        for (int i=0; i<N; i++)
        {
            new_log_prb[i] = log_prb[i];
            ptr[N * (t-1) + i] = (int8_t)i;
        }

        // compute whether a hidden state switch would be used
        for (int prev=0; prev<N; prev++)
        {
            for (int next=0; next<N; next++)
            {
                double log_trans_prb = log_trans[N * prev + next];
                if ( log_prb[ prev ] + log_trans_prb > new_log_prb[ next ] )
                {
                    new_log_prb[ next ] = log_prb[ prev ] + log_trans_prb;
                    ptr[N * (t-1) + next] = (int8_t)prev;
                }
            }
        }

        // update and rescale the current state
        double max = -INFINITY;
        for (int i=0; i<N; i++)
        {
            log_prb[i] = new_log_prb[i] + log_lkl[N * t + i];
            if ( max < log_prb[i] ) max = log_prb[i];
        }

        // rescale Viterbi log probabilities to avoid underflow issues
        for (int i=0; i<N; i++) log_prb[ i ] -= max;
    }

    // retrace Viterbi path from probabilities
    path[T-1] = (int8_t)0;
    for (int i=1; i<N; i++) if (log_prb[(int)path[T-1]] < log_prb[i]) path[T-1] = (int8_t)i;

    // compute best path by tracing back the Markov chain
    for (int t=T-1; t>0; t--)
        path[t-1] = ptr[(t-1) * N + (int)path[t]];

    // free memory
    free(ptr);
}

static double get_log10_odds(const double *log_lkl,
                           const double *log_trans,
                           int N,
                           const int8_t *path,
                           int a,
                           int b)
{
    double log_odds = 0.0;
    for (int t=a; t<b; t++)
    {
        log_odds += log_lkl[N * t + path[t]];
        if ( t > a && path[t] != path[t-1] ) log_odds += log_trans[N * path[t-1] + path[t]];
    }
    return log_odds / M_LN10;
}

static int trans_count(const int8_t *path,
                       int a,
                       int b)
{
    int trans = 0;
    for (int t=a+1; t<b; t++)
    {
        if ( path[t-1] != path[t] ) trans++;
    }
    return trans;
}

/*********************************
 * MAIN FUNCTION                 *
 *********************************/

static double log_trisomy_trans[3 * 3] = {0};
static double log_euploid_trans[2 * 2] = {0};

static const double ff_dflt = 0.034;
static const double rho_dflt = .001;
static const double lambda_dflt = 2000.0;
static const int n_sites_dflt = 1500;
static const int init_state_dflt = HMM_MF;
static const int n_cross_dflt = 2;
static const int n_flips_dflt = 0;
static const int n_repeats_dflt = 1000;

static void usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Simulate HMM to distinguish fetal trisomy from fetal euploid model. ("VERSION")\n");
    fprintf(stderr, "Usage:   simulate [options]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    --ff <float>           fetal fraction [%f]\n", ff_dflt);
    fprintf(stderr, "    --rho <float>          beta-binomial intraclass correlation [%f]\n", rho_dflt);
    fprintf(stderr, "    --lambda <float>       average coverage [%f]\n", lambda_dflt);
    fprintf(stderr, "    --sites <int>          number of heterozygous sites [%d]\n", n_sites_dflt);
    fprintf(stderr, "    --init <MM|MF|FF|M|F>  initial chromosome state [%s]\n", hmm_str[init_state_dflt]);
    fprintf(stderr, "    --cross <int>          number of crossovers [%d]\n", n_cross_dflt);
    fprintf(stderr, "    --flips <int>          number of switch errors, -1 for not using phase [%d]\n", n_flips_dflt);
    fprintf(stderr, "    --repeats <int>        number of simulations to run [%d]\n", n_repeats_dflt);
    fprintf(stderr, "    --full-report          report all log likelihood ratios rather than a summary\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    double ff = ff_dflt;
    double rho = rho_dflt;
    double lambda = lambda_dflt;
    int n_sites = n_sites_dflt;
    int init_state = init_state_dflt;
    int n_cross = n_cross_dflt;
    int n_flips = n_flips_dflt;
    int n_repeats = n_repeats_dflt;
    int full_report = 0;

    // GSL's Taus generator:
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
    // Initialize the GSL generator with time:
    gsl_rng_set(rng, time(NULL)); // Seed with time

    int c;
    char *tmp = NULL;

    static struct option loptions[] =
    {
        {"ff", required_argument, NULL, 1},
        {"rho", required_argument, NULL, 2},
        {"lambda", required_argument, NULL, 3},
        {"sites", required_argument, NULL, 4},
        {"init", required_argument, NULL, 5},
        {"cross", required_argument, NULL, 6},
        {"flips", required_argument, NULL, 7},
        {"repeats", required_argument, NULL, 8},
        {"full-report", no_argument, NULL, 9},
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, "h?",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 1:
                ff = strtod(optarg, &tmp);
                if ( *tmp ) { fprintf(stderr, "Could not parse: --ff %s\n", optarg); exit(-1); }
                break;
            case 2:
                rho = strtod(optarg, &tmp);
                if ( *tmp ) { fprintf(stderr, "Could not parse: --rho %s\n", optarg); exit(-1); }
                break;
            case 3:
                lambda = strtod(optarg, &tmp);
                if ( *tmp ) { fprintf(stderr, "Could not parse: --lambda %s\n", optarg); exit(-1); }
                break;
            case 4:
                n_sites = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) { fprintf(stderr, "Could not parse: --sites %s\n", optarg); exit(-1); }
                break;
            case 5:
                if ( strcmp(optarg, "MM") == 0 ) init_state = HMM_MM;
                else if ( strcmp(optarg, "MF") == 0 ) init_state = HMM_MF;
                else if ( strcmp(optarg, "FF") == 0 ) init_state = HMM_FF;
                else if ( strcmp(optarg, "M") == 0 ) init_state = HMM_M;
                else if ( strcmp(optarg, "F") == 0 ) init_state = HMM_F;
                else { fprintf(stderr, "Could not parse: --init %s\n", optarg); exit(-1); }
                break;
            case 6:
                n_cross = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) { fprintf(stderr, "Could not parse: --cross %s\n", optarg); exit(-1); }
                break;
            case 7:
                n_flips = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) { fprintf(stderr, "Could not parse: --flips %s\n", optarg); exit(-1); }
                break;
            case 8:
                n_repeats = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) { fprintf(stderr, "Could not parse: --repeats %s\n", optarg); exit(-1); }
                break;
            case 9:
                full_report = 1;
                break;
            case 'h':
            case '?':
                usage();
            default:
                fprintf(stderr, "Unknown argument: %s\n", optarg);
                exit(-1);
        }
    }

    double phi = ( 1.0 - ff ) / ff;
    double cross_prb = (double)( 1 + n_cross ) / (double)n_sites; // even if there are no crossovers, a switch should be allowed
    double flip_prb = (double)( n_flips ) / (double)n_sites;
    double log10_odds_sum = 0.0;
    double log10_odds2_sum = 0.0;
    int MAX = 65535;

    // initialize transition matrices costs
    log_trisomy_trans[3 * HMM_MF + HMM_MM] = log(cross_prb);
    log_trisomy_trans[3 * HMM_FF + HMM_MM] = cross_prb * cross_prb < flip_prb ? log(flip_prb) : 2 * log(cross_prb);
    log_trisomy_trans[3 * HMM_MM + HMM_MF] = log(cross_prb);
    log_trisomy_trans[3 * HMM_FF + HMM_MF] = log(cross_prb);
    log_trisomy_trans[3 * HMM_MM + HMM_FF] = cross_prb * cross_prb < flip_prb ? log(flip_prb) : 2 * log(cross_prb);
    log_trisomy_trans[3 * HMM_MF + HMM_FF] = log(cross_prb);
    log_euploid_trans[2 * 0 + 1] = cross_prb < flip_prb ? log(flip_prb) : log(cross_prb);
    log_euploid_trans[2 * 1 + 0] = cross_prb < flip_prb ? log(flip_prb) : log(cross_prb);

    // initialize the beta binomials for likelihood computations and other arrays
    beta_binom_t *beta_binom[5];
    if ( ff > 0.0 )
    {
        double phi = ( 1.0 - ff ) / ff;
        beta_binom[HMM_MM] = beta_binom_init( ( 2.5 + phi ) / ( 3.0 + 2.0 * phi ), rho, MAX);
        beta_binom[HMM_M] = beta_binom_init( ( 1.5 + phi ) / ( 2.0 + 2.0 * phi ), rho, MAX);
        beta_binom[HMM_MF] = beta_binom_init( 0.5, rho, MAX);
        beta_binom[HMM_F] = beta_binom_init( ( 0.5 + phi ) / ( 2.0 + 2.0 * phi ), rho, MAX);
        beta_binom[HMM_FF] = beta_binom_init( ( 0.5 + phi ) / ( 3.0 + 2.0 * phi ), rho, MAX);
    } else {
        beta_binom[HMM_MM] = beta_binom_init( 0.5, rho, MAX);
        beta_binom[HMM_M] = beta_binom_init( 0.5, rho, MAX);
        beta_binom[HMM_MF] = beta_binom_init( 0.5, rho, MAX);
        beta_binom[HMM_F] = beta_binom_init( 0.5, rho, MAX);
        beta_binom[HMM_FF] = beta_binom_init( 0.5, rho, MAX);
    }
    int *hidden_path = (int *)malloc(n_sites * sizeof(int));
    int *cov = (int *)malloc(n_sites * sizeof(int));
    int *a = (int *)malloc(n_sites * sizeof(int));
    double *log_trisomy_lkl = (double *)malloc(3 * n_sites * sizeof(double));
    double *log_euploid_lkl = (double *)malloc(2 * n_sites * sizeof(double));
    int8_t *trisomy_path = (int8_t *)malloc(n_sites * sizeof(int8_t));
    int8_t *euploid_path = (int8_t *)malloc(n_sites * sizeof(int8_t));

    for (int i=0; i<n_repeats; i++)
    {
        // generate a random hidden path with the correct amount of crossover and phase switch errors
        hidden_path[0] = init_state;
        int tmp_cross = n_cross;
        int tmp_flips = n_flips;
        for (int i=1; i<n_sites; i++)
        {
            hidden_path[i] = hidden_path[i - 1];
            // check whether a crossover happened
            if ( tmp_cross > 0 && gsl_rng_uniform(rng) < ( (double)tmp_cross / (double)( n_sites - i ) ) )
            {
                switch ( hidden_path[i] )
                {
                    case HMM_MM:
                        hidden_path[i] = HMM_MF;
                        break;
                    case HMM_M:
                        hidden_path[i] = HMM_F;
                        break;
                    case HMM_MF:
                        if ( gsl_rng_uniform(rng) < 0.5 )
                            hidden_path[i] = HMM_MM;
                        else
                            hidden_path[i] = HMM_FF;
                        break;
                    case HMM_F:
                        hidden_path[i] = HMM_M;
                        break;
                    case HMM_FF:
                        hidden_path[i] = HMM_MF;
                        break;
                    default:
                        fprintf(stderr, "hidden state does not match allowed states: %d\n", hidden_path[i]);
                        exit(-1);
                }
                tmp_cross--;
            }
            // check whether a phase switch error happened
            if ( tmp_flips > 0 && gsl_rng_uniform(rng) < ( (double)tmp_flips / (double)( n_sites - i ) ) )
            {
                switch ( hidden_path[i] )
                {
                    case HMM_MM:
                        hidden_path[i] = HMM_FF;
                        break;
                    case HMM_M:
                        hidden_path[i] = HMM_F;
                        break;
                    case HMM_MF:
                        break;
                    case HMM_F:
                        hidden_path[i] = HMM_M;
                        break;
                    case HMM_FF:
                        hidden_path[i] = HMM_MM;
                        break;
                    default:
                        fprintf(stderr, "hidden state does not match allowed states: %d\n", hidden_path[i]);
                        exit(-1);
                }
                tmp_flips--;
            }
        }

        // generate random Poisson numbers (total coverage)
        for (int i=0; i<n_sites; i++) cov[i] = gsl_ran_poisson(rng, lambda);

        // generate random beta-binomial numbers (reference alternate coverage)
        for (int i=0; i<n_sites; i++)
        {
            if ( cov[i] < 15 )
            {
                double x = gsl_rng_uniform(rng);
                for (a[i]=0; a[i]<cov[i]; a[i]++)
                {
                    x -= dbetabinom(beta_binom[hidden_path[i]], a[i], cov[i]);
                    if ( x < 0 ) break;
                }
            } else {
                if ( rho == 0 )
                {
                    double p = beta_binom[hidden_path[i]]->alpha;
                    a[i] = gsl_ran_binomial(rng, p, cov[i]);
                } else {
                    double p = gsl_ran_beta(rng, beta_binom[hidden_path[i]]->alpha, beta_binom[hidden_path[i]]->beta);
                    a[i] = gsl_ran_binomial(rng, p, cov[i]);
                }
            }
        }

        // compute likelihoods according to the five models
        for (int i=0; i<n_sites; i++)
        {
            if ( n_flips < 0 )
            {
                log_trisomy_lkl[3 * i + HMM_MM] = log_mean_exp( log_dbetabinom(beta_binom[HMM_MM], a[i], cov[i]),
                                                                log_dbetabinom(beta_binom[HMM_MM], cov[i] - a[i], cov[i]) );
                log_trisomy_lkl[3 * i + HMM_MF] = log_dbetabinom(beta_binom[HMM_MF], a[i], cov[i]);
                log_trisomy_lkl[3 * i + HMM_FF] = log_trisomy_lkl[3 * i + HMM_MM];
                log_euploid_lkl[2 * i + HMM_M - 3] = log_mean_exp( log_dbetabinom(beta_binom[HMM_M], a[i], cov[i]),
                                                                   log_dbetabinom(beta_binom[HMM_M], cov[i] - a[i], cov[i]) );
                log_euploid_lkl[2 * i + HMM_F - 3] = log_euploid_lkl[2 * i + HMM_M - 3];
            } else {
                log_trisomy_lkl[3 * i + HMM_MM] = log_dbetabinom(beta_binom[HMM_MM], a[i], cov[i]);
                log_trisomy_lkl[3 * i + HMM_MF] = log_dbetabinom(beta_binom[HMM_MF], a[i], cov[i]);
                log_trisomy_lkl[3 * i + HMM_FF] = log_dbetabinom(beta_binom[HMM_FF], a[i], cov[i]);
                log_euploid_lkl[2 * i + HMM_M - 3] = log_dbetabinom(beta_binom[HMM_M], a[i], cov[i]);
                log_euploid_lkl[2 * i + HMM_F - 3] = log_dbetabinom(beta_binom[HMM_F], a[i], cov[i]);
            }
        }

        // determine best trisomy path
        log_viterbi_run(log_trisomy_lkl, log_trisomy_trans, n_sites, 3, trisomy_path);
        double trisomy_log10_odds = get_log10_odds(log_trisomy_lkl, log_trisomy_trans, 3, trisomy_path, 0, n_sites);
        int trisomy_trans = trans_count(trisomy_path, 0, n_sites);

        // determine best euploid path
        log_viterbi_run(log_euploid_lkl, log_euploid_trans, n_sites, 2, euploid_path);
        double euploid_log10_odds = get_log10_odds(log_euploid_lkl, log_euploid_trans, 2, euploid_path, 0, n_sites);
        int euploid_trans = trans_count(euploid_path, 0, n_sites);

        double log10_odds = trisomy_log10_odds - euploid_log10_odds;
        log10_odds_sum += log10_odds;
        log10_odds2_sum += sq(log10_odds);
        if ( full_report ) fprintf(stdout, "%.4f %.4f %.4f %d %s %d %d %.4f\n",
            ff, rho, lambda, n_sites, hmm_str[init_state], n_cross, n_flips, log10_odds);
    }

    // report summary across simulation repeats
    if ( !full_report )
    {
        double mean_log10_odds = log10_odds_sum / (double)n_repeats;
        double sd_log10_odds = sqrt( ( log10_odds2_sum - sq(log10_odds_sum) / (double)n_repeats ) / (double)( n_repeats - 1 ) );
        fprintf(stdout, "%.4f %.4f %.4f %d %s %d %d %d %.4f %.4f\n",
            ff, rho, lambda, n_sites, hmm_str[init_state], n_cross, n_flips, n_repeats, mean_log10_odds, sd_log10_odds );
    }

    // clean up
    free(trisomy_path);
    free(euploid_path);
    free(hidden_path);
    free(cov);
    free(a);
    free(log_trisomy_lkl);
    free(log_euploid_lkl);
    beta_binom_destroy(beta_binom[HMM_MM]);
    beta_binom_destroy(beta_binom[HMM_M]);
    beta_binom_destroy(beta_binom[HMM_MF]);
    beta_binom_destroy(beta_binom[HMM_F]);
    beta_binom_destroy(beta_binom[HMM_FF]);
    gsl_rng_free(rng);
}
