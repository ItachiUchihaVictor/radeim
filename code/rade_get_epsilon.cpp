/**
 * Compute the best epsilon possible to obtain an (eps,delta)-approximation to
 * (top-k) FIs from a sample of a given dataset with a given sample size and
 * fixed delta, using the stopping condition(s) based on the Rademacher
 * averages.
 *
 * Copyright 2015 Matteo Riondato <matteo@cs.brown.edu>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <climits>
#include <iostream>
#include <iomanip>
#include <map>
#include <random>
#include <unordered_map>
#include <vector>

#include <unistd.h>
extern char *optarg;
extern int optind;
extern int optopt;
extern int opterr;
extern int optreset;

#include "rade.h"
#include "rade_bounds.h"
#include "rade_sample.h"
#include "stats.h"
#include "util.h"

using namespace std;
using namespace std::chrono;

unsigned long sampling_time = 0ul;
unsigned long stopcond_time = 0ul;
unsigned int dataset_size = 0u;
unsigned int sample_size = 0u;
unsigned int max_transaction_length = 0u;
// The following maps store the support of the items, one as if the sample
// was a bag of transactions, the other as if the sample was a set (no
// duplicates).
map<unsigned int, unsigned int> sample_items_sample_as_bag;
map<unsigned int, unsigned int> sample_items_sample_as_set;

static const char* LENGTHS_CI_METH_STR = "lengths";
static const char* ITEMS_CI_METH_STR = "items";
static const char* SIZE_CI_METH_STR = "size";
static const char* EVC_CI_METH_STR = "evc";
static const char* REFINED_CI_METH_STR = "refined";

static const char* STANDARD_STOP_COND_STR = "standard";
static const char* NEW_STOP_COND_STR = "new";
static const char* REFINED_STOP_COND_STR = "refined";
static const char* OPTIMIZE_STOP_COND_STR = "optimize";

bool compute_for_topk = false; // when true, compute a sample for the top-k frequent itemsets
// the following variables control which methods to use to bound the number of closed itemsets
bool use_lengths_CI_method = false;
bool use_items_CI_method = false;
bool use_size_CI_method = false; 
bool use_eVC_CI_method = false;
bool use_refined_CI_method = false;
// the following variables control which stopping conditions to check
bool use_standard_stop_cond = false;
bool use_new_stop_cond = false;
bool use_refined_stop_cond = false;
bool use_optimize_stop_cond = false;
bool verbose = false; // verbose output
char *dataset = NULL; // dataset filename
double delta = 0.0; // confidence parameter
long double opt_param = 1000.0; // parameter to optimize the bound to the Rademacher average. Usually called 's'
long double opt_minimizer = opt_param; // the minimizer that optimizes the bound to the Rademacher average

// Initialize random number generator
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
mt19937_64 generator(seed);

/**
 * Print usage on stderr.
 */
void usage(const char *binary_name) {
	cerr << binary_name << 
		": given a sample size, a dataset, and a delta, compute the best epsilon for an (eps,delta)-approximation obtainable from a sample of that size"
		<< endl;
	cerr << "USAGE: " << binary_name << 
		" [-c closed_itemsets_bound_method1[,closed_itemsets_bound_method2,...] [-h] [-k] [-s stop_condition1[,stop_condition_2]] [-t optimization_parameter] [-v] sample_size delta dataset" 
		<< endl;
	cerr << "\t-c closed_itemsets_bound_method1[,closed_itemsets_bound_method2,...]: specify the method(s) to use to bound the number of closed itemsets. Valid methods are: " 
		<< LENGTHS_CI_METH_STR << ", " << ITEMS_CI_METH_STR << ", " <<
		SIZE_CI_METH_STR << ", " << EVC_CI_METH_STR << ", " << REFINED_CI_METH_STR << endl;
	cerr << "\t-h: print this help message and exit" << endl;
	cerr << "\t-k: compute the epsilon for an (eps,delta)-approx for the top-k FIs" << endl;
	cerr << "\t-s stop_condition1[,stop_condition2]: specify which stopping conditions to check. Valid conditions are: " <<
		STANDARD_STOP_COND_STR << ", " << NEW_STOP_COND_STR << ", " <<
		REFINED_STOP_COND_STR << ", " << OPTIMIZE_STOP_COND_STR << endl;
	cerr << "\t-t optimization_parameter : parameter to optimize the bound (normally called 's')" << endl;
	cerr << "\t-v: verbose output" << endl;
	cerr << "\tsample_size: size of the sample (0 < sample_size)" << endl;
	cerr << "\tdelta: confidence (0 < delta < 1)" << endl;
	cerr << "\tdataset: dataset file" << endl;
}

/**
 * Parse command line options.
 * Return -1 if everything went well, 0 if -h was specified, 1 if there were errors.
 */
int parse_command_line(int& argc, char *argv[]) {
	bool c_specified = false;
	bool s_specified = false;
	int opt;
	while ((opt = getopt(argc, argv, "c:hks:t:v")) != -1) {
		switch (opt) {
		case 'c':
			c_specified = true;
			{ 
				char *token = NULL;
				while ((token = strsep(&optarg, ",")) != NULL) {
					if (strcasecmp(token, LENGTHS_CI_METH_STR) == 0) {
						use_lengths_CI_method = true;
					} else if (strcasecmp(token, ITEMS_CI_METH_STR) == 0) {
						use_items_CI_method = true;
					} else if (strcasecmp(token, SIZE_CI_METH_STR) == 0) {
						use_size_CI_method = true;
					} else if (strcasecmp(token, EVC_CI_METH_STR) == 0) {
						use_eVC_CI_method = true;
					} else if (strcasecmp(token, REFINED_CI_METH_STR) == 0) {
						use_refined_CI_method = true;
					} else {
						cerr << ERROR_HEADER <<
							"closed itemsets bound method string '" 
							<< token 
							<< "' not recognized. It must be one of: " 
							<< LENGTHS_CI_METH_STR << ", " << ITEMS_CI_METH_STR
							<< ", " << SIZE_CI_METH_STR << ", " <<
							EVC_CI_METH_STR << ", " << REFINED_CI_METH_STR << endl;
						return 1;
					}
				}
			}
			break;
		case 'h':
			usage(argv[0]);
			return 0; 
			break;
		case 'k':
			compute_for_topk = true;
			break;
		case 's':
			s_specified = true;
			{
				char *token = NULL;
				while ((token = strsep(&optarg, ",")) != NULL) {
					if (strcasecmp(token, STANDARD_STOP_COND_STR) == 0) {
						use_standard_stop_cond = true;
					} else if (strcasecmp(token, NEW_STOP_COND_STR) == 0) {
						use_new_stop_cond = true;
					} else if (strcasecmp(token, REFINED_STOP_COND_STR) == 0) {
						use_refined_stop_cond = true;
					} else if (strcasecmp(token, OPTIMIZE_STOP_COND_STR) == 0) {
						use_optimize_stop_cond = true;
					} else {
						cerr << ERROR_HEADER << "stop condition string '" <<
							token << "' not recognized. It must be one of: " <<
							STANDARD_STOP_COND_STR << ", " << NEW_STOP_COND_STR
							<< ", " << REFINED_STOP_COND_STR << OPTIMIZE_STOP_COND_STR << endl;
						return 1;
					}
				}
			}
			break;
		case 't':
			opt_param = strtold(optarg, NULL);
			opt_minimizer = opt_param;
			break;
		case 'v':
			verbose = true;
			break;
		} // end switch
	} // end while getopt

	// If -c was not specified, enable the following methods to bound the
	// number of closed itemsets
	if (! c_specified) {
		use_lengths_CI_method = true;
		use_items_CI_method = true;
		use_refined_CI_method = true;
	}

	// If -s was not specify, check only the 'refined' stopping condition
	if (! s_specified) {
		use_refined_stop_cond = true;
	}

	if (optind != argc - 3) {
		cerr << ERROR_HEADER << "wrong number of arguments" << endl;
		return 1;
	} else {
		int signed_sample_size = atoi(argv[argc-3]);
		if (errno == ERANGE || signed_sample_size <= 0) {
			cerr << ERROR_HEADER << "sample_size should be greater than 0" <<
				endl;
			return 1;
		}
		sample_size = (unsigned int) signed_sample_size;
		delta = strtod(argv[argc-2], NULL);
		if (errno == ERANGE || delta >= 1.0 || delta <= 0.0) {
			cerr << ERROR_HEADER << 
				"delta should be greater than 0 and smaller than 1" 
				<< endl;
			return 1;
		}
		dataset = argv[argc-1];
	}

	return -1;
}

/**
 * Main
 */
int main(int argc, char* argv[]) {
	// Parse command line options
	int opt_ret = parse_command_line(argc, argv);
	if (opt_ret == 1) {
		return EXIT_FAILURE;
	} else if (opt_ret == 0) {
		return EXIT_SUCCESS;
	}

	// Get dataset size (doing it now because we assume to have it)
	dataset_size = get_size(dataset);

	// Open dataset file
	FILE *ds_FILE = fopen(dataset, "r");
	if (ds_FILE == NULL) {
		perror("Error opening dataset file");
		return errno;
	}

	sampling_time = sample(ds_FILE, sample_size);

	auto stopcond_start_time = steady_clock::now();
	// Compute bounds to the empirical Rademacher average
	long double rademacher_bound_standard = (use_standard_stop_cond) ?
		get_rademacher_bound_standard() : LDBL_MAX;
	long double rademacher_bound_new = (use_new_stop_cond) ?
		get_rademacher_bound_new() : LDBL_MAX;
	long double rademacher_bound_refined = (use_refined_stop_cond) ?
		get_rademacher_bound_refined() : LDBL_MAX;
	long double rademacher_bound_optimize = (use_optimize_stop_cond) ?
		get_rademacher_bound_optimize() : LDBL_MAX;

	// Compute epsilon(s)
	double multiplier = 2.0;
	if (compute_for_topk) {
		multiplier = 4.0;
	}
	double epsilon_standard = (use_standard_stop_cond) ?
		multiplier * get_max_deviation_bound(rademacher_bound_standard) :
		DBL_MAX;
	long double epsilon_new = (use_new_stop_cond) ?
		multiplier * get_max_deviation_bound(rademacher_bound_new) : LDBL_MAX;
	long double epsilon_refined = (use_refined_stop_cond) ?
		multiplier * get_max_deviation_bound(rademacher_bound_refined) :
		LDBL_MAX;
	long double epsilon_optimize = (use_optimize_stop_cond) ?
		multiplier * get_max_deviation_bound(rademacher_bound_optimize) : LDBL_MAX;

	auto stopcond_end_time = steady_clock::now();
	stopcond_time = duration_cast<milliseconds>(stopcond_end_time - stopcond_start_time).count();
	if (verbose) {
		if (use_standard_stop_cond) {
			cerr << VERBOSE_HEADER << "standard: rademacher_bound: " <<
				rademacher_bound_standard << ", epsilon: " <<
				epsilon_standard << endl;
		}
		if (use_new_stop_cond) {
			cerr << VERBOSE_HEADER << "new: rademacher_bound: " <<
				rademacher_bound_new << ", epsilon: " <<
				epsilon_new << endl;
		}
		if (use_refined_stop_cond) {
			cerr << VERBOSE_HEADER << "refined: rademacher_bound: " <<
				rademacher_bound_refined << ", epsilon: " <<
				epsilon_refined << endl;
		}
		if (use_optimize_stop_cond) {
			cerr << VERBOSE_HEADER << "optimize: rademacher_bound: " <<
				rademacher_bound_optimize << ", epsilon: " <<
				epsilon_optimize << endl;
		}
		if (use_standard_stop_cond + use_new_stop_cond + use_refined_stop_cond + use_optimize_stop_cond > 1) {
			cerr << VERBOSE_HEADER <<"minimum: ";
			long double minimum = min(rademacher_bound_standard,
					min(rademacher_bound_new, min(rademacher_bound_refined, rademacher_bound_optimize)));
			long double min_epsilon = 0.0;
			if (minimum == rademacher_bound_standard) {
				cerr << "rademacher_bound_standard"; 
				min_epsilon = epsilon_standard;
			} else if (minimum == rademacher_bound_new) {
				cerr << "rademacher_bound_new"; 
				min_epsilon = epsilon_new;
			} else if (minimum == rademacher_bound_refined) {
				cerr << "rademacher_bound_refined";
				min_epsilon = epsilon_refined;
			} else {
				cerr << "rademacher_bound_optimize";
				min_epsilon = epsilon_optimize;
			}
			cerr << " (" << minimum << ", " << min_epsilon << ")" << endl;
		}
	} // if (verbose)

	// Close dataset file
	fclose(ds_FILE);

	// Print infos and epsilon values
	cerr << LOG_HEADER << "creating the sample took " << sampling_time + stopcond_time <<
		" ms" << endl;
	if (use_standard_stop_cond)
		cerr << LOG_HEADER << "epsilon_standard : " << epsilon_standard << endl;
	if (use_new_stop_cond)
		cerr << LOG_HEADER << "epsilon_new : " << epsilon_new << endl;
	if (use_refined_stop_cond)
		cerr << LOG_HEADER << "epsilon_refined : " << epsilon_refined << endl;
	if (use_optimize_stop_cond)
		cerr << LOG_HEADER << "epsilon_optimize : " << epsilon_optimize << endl;

	// Write info as comma-separated values to stderr
	cerr << "dataset,delta,topk,sample_size,use_lengths_CI_method,use_items_CI_method,use_size_CI_method,use_eVC_CI_method,use_refined_CI_method,use_standard_stop_cond,use_new_stop_cond,use_refined_stop_cond,use_optimize_stop_cond,sampling_time,stopcond_time,runtime,opt_param,opt_minimizer,epsilon_standard,epsilon_new,epsilon_refined,epsilon_optimize" << endl;
	cerr << (strrchr(dataset, '/') + 1) << COMMA << delta << COMMA <<
		compute_for_topk << COMMA << sample_size << COMMA <<
		use_lengths_CI_method << COMMA << use_items_CI_method << COMMA <<
		use_size_CI_method << COMMA << use_eVC_CI_method << COMMA <<
		use_refined_CI_method << COMMA << use_standard_stop_cond << COMMA <<
		use_new_stop_cond << COMMA << use_refined_stop_cond << COMMA <<
		use_optimize_stop_cond << COMMA << sampling_time << COMMA <<
		stopcond_time << COMMA << sampling_time + stopcond_time << COMMA <<
		opt_param << COMMA << opt_minimizer << COMMA << epsilon_standard <<
		COMMA << epsilon_new << COMMA << epsilon_refined << COMMA <<
		epsilon_optimize << endl; 
	return EXIT_SUCCESS;
}
