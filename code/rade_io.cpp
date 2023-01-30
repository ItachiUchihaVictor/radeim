/**
 *
 * Copyright 2014 Matteo Riondato <matteo@cs.brown.edu>
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

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <unistd.h>
extern char *optarg;
extern int optind;
extern int optopt;
extern int opterr;
extern int optreset;

#include "rade_io.h"
#include "util.h"


using namespace std;

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
bool automatic_schedule = false; // when true, use automatic schedule
bool geometric_schedule = true; // when true, use geometric schedule
bool linear_schedule = false; // when true, use linear schedule
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
bool wait_for_all = false; // wait for all stopping conditions to be satisfied before returning
char *dataset = NULL; // dataset filename
double delta = 0.0; // confidence parameter
double epsilon = 0.0; // accuracy parameter
long double opt_param = 1000.0; // parameter to optimize the bound to the Rademacher average. Usually called 's'
long double opt_minimizer = opt_param; // the minimizer that optimizes the bound to the Rademacher average
double geometric_schedule_multiplier = 2.0; // control the growth of the sample size when using the geometric schedule
unsigned int linear_schedule_addend = 1000; // control the growth of the sample size when using the linear schedule
unsigned int initial_sample_size = 1000; // initial sample size;

/**
 * Print usage on stderr.
 */
void __attribute__((cold)) usage(const char *binary_name) {
	cerr << binary_name << ": progressive sampling algorithm to extract an (eps,delta)-approximation" << endl;
	cerr << "USAGE: " << binary_name << 
		" [-a | -g geometric_schedule_multiplier | -l linear_schedule_addend] [-c closed_itemsets_bound_method1[,closed_itemsets_bound_method2,...] [-h] [-i initial_sample_size] [-k] [-s stop_condition1[,stop_condition_2]] [-t optimization_parameter] [-vw] epsilon delta dataset" 
		<< endl;
	cerr << "\t-a: use the automatic sample schedule" << endl;
	cerr << "\t-c closed_itemsets_bound_method1[,closed_itemsets_bound_method2,...]: specify the method(s) to use to bound the number of closed itemsets. Valid methods are: " 
		<< LENGTHS_CI_METH_STR << ", " << ITEMS_CI_METH_STR << ", " <<
		SIZE_CI_METH_STR << ", " << EVC_CI_METH_STR << ", " << REFINED_CI_METH_STR << endl;
	cerr << "\t-g geometric_schedule_multiplier: multiplier for the geometric schedule" << endl;
	cerr << "\t-h: print this help message and exit" << endl;
	cerr << "\t-i initial_sample_size: specify initial sample size (default: "
		<< initial_sample_size << ")" << endl;
	cerr << "\t-k: compute the sample for an (eps,delta)-approx for the top-k FIs" << endl;
	cerr << "\t-l linear_schedule_addend: addend for the linear schedule" << endl;
	cerr << "\t-s stop_condition1[,stop_condition2]: specify which stopping conditions to check. Valid conditions are: " <<
		STANDARD_STOP_COND_STR << ", " << NEW_STOP_COND_STR << ", " <<
		REFINED_STOP_COND_STR << ", " << OPTIMIZE_STOP_COND_STR << endl;
	cerr << "\t-t optimization_parameter : parameter to optimize the bound (normally called 's')" << endl;
	cerr << "\t-v: verbose output" << endl;
	cerr << "\t-w: wait until all stopping conditions are satisfied before returning" << endl;
	cerr << "\tepsilon: accuracy (0 < epsilon < 1)" << endl;
	cerr << "\tdelta: confidence (0 < delta < 1)" << endl;
	cerr << "\tdataset: dataset file" << endl;
}

/**
 * Parse command line options.
 * Return -1 if everything went well, 0 if -h was specified, 1 if there were errors.
 */
int __attribute__((cold)) parse_command_line(int& argc, char *argv[]) {
	bool a_specified = false;
	bool c_specified = false;
	bool g_specified = false;
	bool l_specified = false;
	bool s_specified = false;
	int opt;
	while ((opt = getopt(argc, argv, "ac:g:hki:l:s:t:vw")) != -1) {
		switch (opt) {
		case 'a':
			automatic_schedule = true;
			geometric_schedule = false;
			linear_schedule = false;
			a_specified = true;
			break;
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
		case 'g':
			g_specified = true;
			geometric_schedule_multiplier = strtod(optarg, NULL);
			if (errno == ERANGE || geometric_schedule_multiplier <= 1.0) {
				cerr << ERROR_HEADER <<
					"geometric schedule multiplier must be greater than 1"
					<< endl;
				return 1;
			}
			geometric_schedule = true;
			automatic_schedule = false;
			linear_schedule = false;
			break;
		case 'h':
			usage(argv[0]);
			return 0;
			break;
		case 'i':
			{
				int tmp_initial_sample_size = atoi(optarg);
				if (tmp_initial_sample_size <= 0) {
					cerr << ERROR_HEADER <<
						"initial sample size must be positive"
						<< endl;
					return 1;
				}
				initial_sample_size = tmp_initial_sample_size;
			}
			break;
		case 'k':
			compute_for_topk = true;
			break;
		case 'l':
			l_specified = true;
			{
				int tmp_linear_schedule_addend = atoi(optarg);
				if (errno == ERANGE || tmp_linear_schedule_addend <= 0) {
					cerr << ERROR_HEADER <<
						"linear schedule addend should be positive"
						<< endl;
					return 1;
				}
				linear_schedule_addend = tmp_linear_schedule_addend;
			}
			linear_schedule = true;
			automatic_schedule = false;
			geometric_schedule = false;
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
		case 'w':
			wait_for_all = true;
			break;
		} // end switch
	} // end while getopt

	if (a_specified + g_specified + l_specified > 1) {
		cerr << ERROR_HEADER << 
			"only one between -a, -g, and -l can be specified" << endl;
		return 1;
	}

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
		epsilon = strtod(argv[argc-3], NULL);
		if (errno == ERANGE || epsilon >= 1.0 || epsilon <= 0.0) {
			cerr << ERROR_HEADER <<
				"epsilon should be greater than 0 and smaller than 1"
				<< endl;
			return 1;
		}
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
