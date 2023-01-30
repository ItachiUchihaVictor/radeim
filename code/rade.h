/**
 *
 * copyright 2014 matteo riondato <matteo@cs.brown.edu>
 *
 * licensed under the apache license, version 2.0 (the "license");
 * you may not use this file except in compliance with the license.
 * you may obtain a copy of the license at
 *
 *      http://www.apache.org/licenses/license-2.0
 *
 * unless required by applicable law or agreed to in writing, software
 * distributed under the license is distributed on an "as is" basis,
 * without warranties or conditions of any kind, either express or implied.
 * see the license for the specific language governing permissions and
 * limitations under the license.
 * 
 */


#include <map>
#include <set>

using namespace std;

// the following variables control which methods to use to bound the number of closed itemsets
extern bool use_lengths_CI_method;
extern bool use_items_CI_method;
extern bool use_size_CI_method; 
extern bool use_eVC_CI_method;
extern bool use_refined_CI_method;
// the following variables control which stopping conditions to check
extern bool use_standard_stop_cond;
extern bool use_new_stop_cond;
extern bool use_refined_stop_cond;
extern bool use_optimize_stop_cond;
extern bool verbose; // verbose output
extern double delta; // confidence parameter
extern double epsilon; // accuracy parameter
extern long double opt_param; // parameter to optimize the bound to the Rademacher average. Usually called 's'
long double opt_minimizer; // the minimizer that optimizes the bound to the Rademacher average
// The key is a transaction, the value the number of times it was sampled
extern map<set<unsigned int>, unsigned int> all_sampled_transactions;
extern unsigned int dataset_size; // Dataset size
unsigned int sample_size; // Current sample size
extern unsigned int max_transaction_length; // Current max transaction length in the sample
// The following maps store the support of the items, one as if the sample
// was a bag of transactions, the other as if the sample was a set (no
// duplicates).
extern map<unsigned int, unsigned int> sample_items_sample_as_bag;
extern map<unsigned int, unsigned int> sample_items_sample_as_set;

