#ifndef RADE_BOUNDS_CI_H
#define RADE_BOUNDS_CI_H

#include <map>
using namespace std;

long double get_log_CI_bound_lengths(const map<unsigned int, unsigned int>& items);
long double get_log_CI_bound_refined(const map<unsigned int, unsigned int>& items);
long double get_log_CI_bound_items(const map<unsigned int, unsigned int>& items);
long double get_log_CI_bound_eVC(const map<unsigned int, unsigned int>& items);
long double get_log_CI_bound(const map<unsigned int, unsigned int>& items);
long double get_log_CI_bound_size(const map<unsigned int, unsigned int>& items);
#endif
