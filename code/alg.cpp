#include "stdafx.h"
#include <nlopt.h>
#include "util.h"
#include "rade.h"
#include "rade_bounds.h"
#include "rade_bounds_ci.h"

double VS;
//int sample_size;
double ita = 0.1;

int graphNumV = 0;
int seedSizeK = 0;
static const string COMMA = ",";
static const string SPACE = " ";
static const string VERBOSE_HEADER = "INFO: ";
static const string LOG_HEADER = "LOG: ";
static const string ERROR_HEADER = "ERROR: ";

// A function to find the factorial.
double factorial(int n)
{
 	int i;
	int y = logl(n);
//	for(i = n-1; i > 1; i--)y += logl(i);

 	return (double) n;
}

int combination(int n, int r){
//	int y = logl(n);
//	for(int i=n-1;i>n-r;i--){
//		y += logl(i);
//	}
//	for(int i=1;i<r+1;i++) y -=logl(i);
//	return y;
	return factorial(n)- (factorial(r)+factorial(n-r));
}
inline long double logsum(const long double log_a, const long double log_b) {
#ifdef DEBUG
	long double to_return = max(log_a, log_b) + log1pl(expl(fminl(log_a, log_b) - max(log_a, log_b)));
	assert(to_return > max(log_a, log_b) || to_return < nexttoward(max(log_a,log_b), LDBL_MAX));
	return to_return;
#else
	return max(log_a, log_b) + log1pl(expl(min(log_a, log_b) - max(log_a, log_b)));
#endif
}

/**
 *
 * Get the quantity that bounds the maximum deviation
 */
long double  __attribute__((pure)) get_max_deviation_bound(const long double rademacher_bound, double inf) {
	static double fixed_part = sqrt(log(2.0 / ita));
//	double poss = log(5.0/ita)/sample_size;
//	double beta = inf;
//	double gamma = beta + 2.0/3.0*poss + sqrt(pow2(poss/sqrt(3.0))+2.0*beta*poss);
//	double pho = rademacher_bound + 2.0/3.0*poss + sqrt(4.0*beta*poss);
//	double rr = pho + poss/3.0 + sqrt(pow2(poss/2.0/sqrt(3.0))+pho*poss);
	return 2.0 * rademacher_bound + log(2.0/ita)/ (2.0 * sample_size);// + log(2.0/ita)/(alphaomega*sample_size);
	double alphaomega = log(2.0/ita)/(log(2.0/ita)+sqrt(log(2.0/ita)*(2*sample_size*rademacher_bound+log(2.0/ita))));
//	double orgerror = 2.0 * rademacher_bound + fixed_part / sqrt(2 * sample_size) + log(2.0/ita)/(alphaomega*sample_size);
  //      cout << "rademacher bound: " << rademacher_bound << ", poss: " << poss << ", beta: " << beta << ", gamma: " << gamma << ", pho: " << pho << ", rr: " << rr << " sample size: " << sample_size <<endl;
	
//	double error = 2.0 * rademacher_bound + poss/3 + sqrt(2.0*(gamma+4*rr)*poss);
//	cout << "orgerror: " << orgerror << ", error: " << error << endl;

//	return error;

        cout << "rademacher bound: " << rademacher_bound << ", fixed part: " << fixed_part << " sample size: " << sample_size << ", second term: " << log(2.0/ita)/(alphaomega*sample_size) << ", third term: " << fixed_part/sqrt(2 * sample_size) <<endl;
//	return 2.0 * rademacher_bound + fixed_part/ (2 * sample_size) + log(2.0/ita)/(alphaomega*sample_size);
	return 2.0 * rademacher_bound + log(2.0/ita)/ (2.0 * sample_size) + log(2.0/ita)/(alphaomega*sample_size);
	return 2.0 * rademacher_bound + fixed_part / sqrt(2 * sample_size) + log(2.0/ita)/(alphaomega*sample_size);
}



long double  __attribute__((pure)) get_updated_sample_size(const long double rademacher_bound, double epsilon) {
                // The commented out method to compute the next sample size is
                // for the "vanilla" bound to the supremum of the deviation
                //sampleSize *= (supDeviationBound * supDeviationBound) / (epsilon * epsilon);
                // The following method to compute the next sample size is
                // suitable for the refined bound to the supremum of the
                // deviations from Oneto et al., "An improved analysis of
                // the Rademacher data-dependent bound using its self bounding
                // property", Neural Networks 44, 107-111, 2013.
                // TODO MR Think about the penaltyFactor for the refined estimator
                assert(rademacher_bound < epsilon);
                const double mylog = log(2 / ita);
           //     const double epsilon_square = epsilon * epsilon;
                const double epsilon_square = epsilon;
                const double radebound_square = rademacher_bound * rademacher_bound;
                const double twentyseven_radebound_square = 27 * radebound_square;
                const double next_multiplier = mylog / (6 * (epsilon_square - 2 * epsilon * rademacher_bound + radebound_square));
                const double next_1st_term = 2 + 8 * epsilon;
                const double next_2nd_term_multiplier = sqrt(48 * rademacher_bound + 1 + 8 * epsilon + 16 * epsilon_square);
                const double atan2_1st_arg = 12 * sqrt(3) * fabs(-1 + 2 * (rademacher_bound + epsilon)) * sqrt(- (twentyseven_radebound_square - epsilon_square * (1 + 16 * epsilon) - rademacher_bound * (1 + 18 * epsilon)));
                const double atan2_2nd_arg = - ( -1 - 12 * epsilon + 8 * ( twentyseven_radebound_square + (21 - 8 * epsilon) * epsilon_square + 18 * rademacher_bound * (1 + epsilon)));
                const double atan2_res = atan2(atan2_1st_arg, atan2_2nd_arg);
                const double theta = atan2_res / 3;
                const double cos_theta = cos(theta);
                const double sin_theta = sin(theta);
                const std::array<int, 3> roots = {
                    (int) ceil(next_multiplier * (next_1st_term - 2 * next_2nd_term_multiplier * cos_theta)),
                    (int) ceil(next_multiplier * (next_1st_term + next_2nd_term_multiplier * (cos_theta + sqrt(3) * sin_theta))),
                    (int) ceil(next_multiplier * (next_1st_term + next_2nd_term_multiplier * (cos_theta - sqrt(3) * sin_theta)))};
                double sampleSize = *(std::max_element(roots.begin(), roots.end()));
                for (auto root : roots) {
                    if (root <= sample_size || root == sampleSize) {
                    //    TRACE("NEXT, root: ", root, " skipping as non informative");
                        continue;
                    }
                    const double alpha = mylog / (mylog + sqrt(mylog * (2 * root * rademacher_bound + mylog)));
                    const double first_term = sqrt(mylog / (2 * root));
                    const double second_term = rademacher_bound / (1 - alpha);
                    const double third_term = mylog / (2 * root * alpha * (1 - alpha));
                    const double next_deviation = first_term + second_term + third_term;
          //          TRACE("NEXT, root: ", root, ", alpha: ",
            //                alpha, ", first_term: ", first_term,
              //              ", second_term: ", second_term, ", third_term: ",
                //            third_term, ", deviation: ", next_deviation);
                    if (next_deviation <= epsilon && root < sampleSize) {
                  //      TRACE("Smaller next sample size: ", root, " (old: ", sampleSize, ")");
                        sampleSize = root;
                    }
                }
      //          TRACE("next sample size: ", sampleSize);
     //           assert(sampleSize > sample_size);
    //        sample_size = sampleSize - sample_size;
	return sampleSize;// log(2.0/ita)*(1+2*epsilon/pow2(alphaomega)+sqrt(1+4*epsilon/pow2(alphaomega)))/(4*epsilon);
}


/*
 * Compute the optimization function and the gradient
 */
double opt_fun(unsigned, const double *x, double *grad, void *supps_logquants) {
//	vector<long double> *supps_logquants = (vector<long double> *) supps_quants_pointer;

	long double log_sum_square_none_none = -1.0;
	long double log_sum_square_none_lin = -1.0;
	double *iter = (double *) supps_logquants;
	long double s_square = powl(x[0], 2.0);
	long double exponent = s_square * sqrt(*iter)/2.0/powl(sample_size, 2.0);
	log_sum_square_none_none = exponent;
//	log_sum_square_none_lin = exponent + logl(sqrt(*iter));
//	++iter;
//	while (iter != (*supps_logquants).end()) {
//        for(int i=0;i<5000;i++){
//		exponent = s_square * sqrt(*iter)/2.0/powl(sample_size,2.0);
//		log_sum_square_none_none = logsum(log_sum_square_none_none, exponent);
//		log_sum_square_none_lin = logsum(log_sum_square_none_lin, exponent + logl(sqrt(*iter)));
//		++iter;
//	}
	// XXX NLOpt only uses doubles, not long doubles :(
//	cout << "graph num V: " << graphNumV << " seedSizeK: " << seedSizeK << " combination(graphNumV, seedSizeK):" << combination(graphNumV, seedSizeK) << std::endl;
	if (grad != NULL) {
		//grad[0] = (double) (expl(log_sum_square_none_lin + logl(2.0) - log_sum_square_none_none) - (log_sum_square_none_none / powl(x[0], 2.0)));
		//grad[0] = (double) (sqrt(*iter)/2.0/powl(sample_size, 2.0) - logl((double) (graphNumV)) /  powl(x[0], 2.0));
		grad[0] = (double) (sqrt(*iter)/2.0/powl(sample_size, 2.0) - ((double) (combination(graphNumV, seedSizeK))) /  powl(x[0], 2.0));
	}
	//if (verbose) {
	//cerr << "x: " << x[0] << ", val: " << (double) (log_sum_square_none_none / x[0]) << " grad: " << ((grad) ? grad[0] : -1) << endl;
	//}
	//return (double) (log_sum_square_none_none / x[0] + logl((double) (graphNumV)) / x[0]);
	return (double) (log_sum_square_none_none / x[0] + ((double) (combination(graphNumV, seedSizeK))) / x[0]);
}

/*
 * Compute a bound to the empirical Rademacher average using a new refined
 * procedure that uses optimization.
 */
long double get_rademacher_bound_optimize() {
	double supps_logquants = VS;
	long double factor = 2.0 * powl(sample_size, 2.0);

//	for (auto element : VS) {
//		supps_logquants.push_back(element.second);
//	}
//	cout << "supps: ";
//	for (auto element : supps_logquants) {
//		if(element > 0) cout << element << " ";
//	}
//	cout << endl;
	// Setup and solve the optimization problem
	// FIXME We are using a derivative-free algorithm despite knowing how to compute
	// the gradient because there seem to be some problems in NLopt
	// implementation or in its interpretation of our gradient or in our
	// computation of the gradient (see opt_fun).
	// We would even know how to compute the Hessian (second derivative), but
	// we're not using that either.
	// Note that we are also using a "local" optimization algorithm because the
	// function is convex.
	nlopt_opt opt_prob = nlopt_create(NLOPT_LN_COBYLA, 1);
	nlopt_set_min_objective(opt_prob, opt_fun, &supps_logquants);
	double lb[1] = {1.0};
	nlopt_set_lower_bounds(opt_prob, lb);
	// XXX It is weird that we have to define a const variable and then use a
	// pointer to it. NLopt bug? It's not what's documented.
	const double xtol_abs = 1e-4;
	nlopt_set_xtol_abs(opt_prob, &xtol_abs);
	nlopt_set_ftol_abs(opt_prob, 1e-6);
	// Set initialization point
	double x[1];
#ifdef DEBUG
	x[0] = opt_param;
#else
	x[0] = (double) opt_minimizer;
#endif
	double min_value;
	x[0] = 1;

	nlopt_result result = nlopt_optimize(opt_prob, x, &min_value);
//#ifdef DEBUG
	cerr << "mathematica minimizer: " << opt_minimizer << ", nlopt minimizer: " << x[0] << endl;
//#endif
	if (result >= 0) {
		opt_minimizer = x[0];
//		if (verbose) {
			cerr << VERBOSE_HEADER << "get_rademacher_bound_optimize: found minimum at "
				<< x[0] << ": " << min_value << ", return code: " << result <<
				endl;
//		}
	} else  {
		cerr << "Optimization failed! return code: " << result << endl;
		exit(1);
	}
	nlopt_destroy(opt_prob);

	return min_value;
}


// code from OPIM
double Alg::MaxCoverVanilla(const int targetSize)
{
    // optimization with minimum upper bound among all rounds [Default].
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    // degMap: map degree to the nodes with this degree
    RRsets degMap(maxDeg + 1);

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                // Find the top-k marginal coverage
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                // Top-k influential nodes constructed
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

// In the case of identical marginal，the node with largest out degree is chosen
double Alg::MaxCoverOutDegPrority(const int targetSize)
{
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

// max cover used in the IM-sentinel Phase
double Alg::MaxCoverIMSentinel(std::vector<uint32_t> &seedSet, const int targetSize)
{
    // seedSet: the sentinel set obtained in the Sentinel Set Selection Phase
    std::unordered_set<uint32_t> subSeedSet(seedSet.begin(), seedSet.end());
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto node : seedSet)
    {
        sumInf += coverage[node];
        _vecSeed.push_back(node);
        coverage[node] = 0;

        for (auto edgeIdx : _hyperGraph._FRsets[node])
        {
            if (edgeMark[edgeIdx]) continue;

            edgeMark[edgeIdx] = true;

            for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
            {
                if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                coverage[nodeIdx]--;
            }
        }
    }

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

double Alg::MaxCoverSentinelSet(const int targetSize, const int totalTargetSize)
{
    //targetSize: the size of the sentinel set
    //totalTargetSize: the total number of the seed set
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = totalTargetSize;
                auto degBound = deg;
                FRset vecBound(totalTargetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                // Find the top-k marginal coverage
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                goto afterGreedy;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            _vecVldtInf.push_back(sumInf);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

afterGreedy:
    const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
    // std::cout << "  >>>[greedy-lazy] influence: " << finalInf << ", seed set: " << _vecSeed.size() << ", min-bound: " << _boundMin <<
    //           ", last-bound: " << _boundLast << std::endl;

    if (_vecSeed.size() == targetSize)
    {
        // if the sample size is sufficiently large, the result is reliable
        if (_numRRsets > 1000)
        {
            return finalInf;
        }
    }

    // if covering all the RR-sets, the sentinel set may include some nodes which cover only a small number of RR-sets.
    // such nodes should not be included.
    uint32_t threshold = 0.9 * sumInf;
    for (int i = _vecSeed.size() - 1; i > 0; i--)
    {
        if (_vecVldtInf[i - 1] >= threshold)
        {
            degMap[0].push_back(_vecSeed[i]);
            _vecSeed.pop_back();
        }
        else
        {
            break;
        }
    }

    // std::cout << "seedset size reaching 0.9 coverage: " << _vecSeed.size() << std::endl;
    // the following code is to select the nodes with large out-degree. 
    int seedSetSize = (degMap[0].size() > targetSize) ? targetSize : degMap[0].size();
    std::vector<std::pair<uint32_t, uint32_t>> vecHeap;

    for (int i = 0; i < seedSetSize; i++)
    {
        auto &node = degMap[0][i];
        vecHeap.push_back(std::make_pair(_vecOutDegree[node], node));
    }

    std::make_heap(vecHeap.begin(), vecHeap.end(), GreaterPair);
    const auto nodeNum = degMap[0].size();


    for (int i = seedSetSize; i < nodeNum; i++)
    {
        uint32_t node = degMap[0][i];
        uint32_t currDeg = _vecOutDegree[node];

        if (currDeg > vecHeap[0].first)
        {
            std::pop_heap(vecHeap.begin(), vecHeap.end());
            vecHeap.pop_back();
            vecHeap.push_back(std::make_pair(_vecOutDegree[node], node));
            std::push_heap(vecHeap.begin(), vecHeap.end());
        }
    }

    std::sort_heap(vecHeap.begin(), vecHeap.end(), GreaterPair);
    std::unordered_set<uint32_t> seedHashSet(_vecSeed.begin(), _vecSeed.end());

    for (auto &node : vecHeap)
    {
        if (_vecSeed.size() >= targetSize)
        {
            break;
        }

        if (seedHashSet.find(node.second) != seedHashSet.end())
        {
            std::cout << "node exist" << std::endl;
            continue;
        }

        _vecSeed.push_back(node.second);
        seedHashSet.insert(node.second);

        if (_vecVldtInf.size() < _vecSeed.size())
        {
            _vecVldtInf.push_back(sumInf);
        }
    }

    return 1.0 * _numV; // All RR sets are covered.
}

double Alg::MaxCoverTopK(const int targetSize)
{
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        //if (coverage[i] == 0) continue;
        degMap[coverage[i]].push_back(i);
    }

    Nodelist sortedNode(_numV); // sortedNode: record the sorted nodes in ascending order of degree
    Nodelist nodePosition(_numV); // nodePosition: record the position of each node in the sortedNode
    Nodelist degreePosition(maxDeg + 2); // degreePosition: the start position of each degree in sortedNode
    uint32_t idxSort = 0;
    size_t idxDegree = 0;

    for (auto& nodes : degMap)
    {
        degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
        idxDegree++;

        for (auto& node : nodes)
        {
            nodePosition[node] = idxSort;
            sortedNode[idxSort++] = node;
        }
    }

    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    // record the total of top-k marginal gains
    size_t sumTopk = 0;

    for (auto deg = maxDeg + 1; deg--;)
    {
        if (degreePosition[deg] <= _numV - targetSize)
        {
            sumTopk += deg * (degreePosition[deg + 1] - (_numV - targetSize));
            break;
        }

        sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
    }

    _boundMin = 1.0 * sumTopk;
    _vecSeed.clear();
    size_t sumInf = 0;

    /*
    * sortedNode: position -> node
    * nodePosition: node -> position
    * degreePosition: degree -> position (start position of this degree)
    * coverage: node -> degree
    * e.g., swap the position of a node with the start position of its degree
    * swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
    */
    for (auto k = targetSize; k--;)
    {
        const auto seed = sortedNode.back();
        sortedNode.pop_back();
        const auto newNumV = sortedNode.size();
        sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];
        sumInf += coverage[seed];
        _vecSeed.push_back(seed);
        coverage[seed] = 0;

        for (auto edgeIdx : _hyperGraph._FRsets[seed])
        {
            if (edgeMark[edgeIdx]) continue;

            edgeMark[edgeIdx] = true;

            for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
            {
                if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                const auto currPos = nodePosition[nodeIdx]; // The current position
                const auto currDeg = coverage[nodeIdx]; // The current degree
                const auto startPos = degreePosition[currDeg]; // The start position of this degree
                const auto startNode = sortedNode[startPos]; // The node with the start position
                // Swap this node to the start position with the same degree, and update their positions in nodePosition
                std::swap(sortedNode[currPos], sortedNode[startPos]);
                nodePosition[nodeIdx] = startPos;
                nodePosition[startNode] = currPos;
                // Increase the start position of this degree by 1, and decrease the degree of this node by 1
                degreePosition[currDeg]++;
                coverage[nodeIdx]--;

                // If the start position of this degree is in top-k, reduce topk by 1
                if (startPos >= newNumV - targetSize) sumTopk--;
            }
        }

        _boundLast = 1.0 * (sumInf + sumTopk);

        if (_boundMin > _boundLast) _boundMin = _boundLast;
    }

    _boundMin *= 1.0 * _numV / _numRRsets;
    _boundLast *= 1.0 * _numV / _numRRsets;
    const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
    std::cout << "  >>>[greedy-topk] influence: " << finalInf << ", min-bound: " << _boundMin <<
              ", last-bound: " << _boundLast << '\n';
    return finalInf;
}

double Alg::MaxCover(const int targetSize)
{
    if (targetSize >= 1000) return MaxCoverTopK(targetSize);

    return MaxCoverVanilla(targetSize);
}

void Alg::set_prob_dist(const ProbDist dist)
{
    _probDist = dist;
    _hyperGraph.set_prob_dist(dist);
    _hyperGraphVldt.set_prob_dist(dist);
}

void Alg::set_vanilla_sample(const bool isVanilla)
{
    if (isVanilla)
    {
        std::cout << "Vanilla sampling method is used" << std::endl;
    }

    _hyperGraph.set_vanilla_sample(isVanilla);
    _hyperGraphVldt.set_vanilla_sample(isVanilla);
}

double Alg::EfficInfVldtAlg()
{
    return EfficInfVldtAlg(_vecSeed);
}

double Alg::EfficInfVldtAlg(const Nodelist vecSeed)
{
    Timer EvalTimer("Inf. Eval.");
    std::cout << "  >>>Evaluating influence in [0.99,1.01]*EPT with prob. 99.9% ...\n";
    const auto inf = _hyperGraph.EfficInfVldtAlg(vecSeed);
    //const auto inf = _hyperGraphVldt.EfficInfVldtAlg(vecSeed);
    std::cout << "  >>>Down! influence: " << inf << ", time used (sec): " << EvalTimer.get_total_time() << '\n';
    return inf;
}

double Alg::estimateRRSize()
{
    const int sampleNum = 100;
    _hyperGraph.BuildRRsets(sampleNum);
    double avg= _hyperGraph.HyperedgeAvg();
    _hyperGraph.RefreshHypergraph();
    return avg;
}

double Alg::subsimOnly(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    ita = 1.0/_numV;
    graphNumV = _numV;
    seedSizeK = targetSize;
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
    //double numRbase = log(2/ita)/2*((pow2(1.0-1.0/e)+1.0-1.0/e)/(epsilon)-1.0+1.0/e); //size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
 //   double numRbase = log(2/ita)/epsilon;//2*((pow2(1.0-1.0/e)+1.0-1.0/e)/(epsilon)-1.0+1.0/e); //size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    double numRbase = log(2/ita)/((epsilon)); //size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = 100000 * (size_t)log2(maxNumR / numRbase) + 1;
    double numAdd = 0;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    auto numR = numRbase; // << (idx-1);

//    std::cout << std::endl;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        numR += numAdd; // << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsets(numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        const auto infSelf = MaxCover(targetSize);
        time2 += timerSubsim.get_operation_time();

        auto degVldt = infSelf; // * _numRRsets / _numV;
	auto infVldt = degVldt;
	cout << "degVldt: " << degVldt << endl;
        auto upperBound = _boundMin;
	sample_size = _numRRsets;
	VS = degVldt;

//	cout << "infVldt: " << infVldt << endl;
//	cout << "infVldt/_numRRset: " << infVldt/_numRRsets << endl;
//	cout << "infVldt/numR: " << infVldt/numR << endl;
//	cout << "infVldt/_numV: " << infVldt/_numV << endl;
//		        for(int i=0; i<sims.size(); i++) std::cout << i << " " << sims[i] << endl;
        double deviation_bound_optimize = (get_rademacher_bound_optimize());
        double empirical_error = get_max_deviation_bound(deviation_bound_optimize, infVldt/_numV);
    //    double empirical_error = (logl((seedSizeK+1))/2.0 + 1.0 ) * (1+30*(epsilon-0.01)) * get_max_deviation_bound(deviation_bound_optimize, infVldt/_numV);
	
//       double empirical_error = (logl(seedSizeK+1) ) * (1+30*(epsilon-0.001)) * get_max_deviation_bound(deviation_bound_optimize, infVldt/_numV);
	double epsilonplus = infVldt / _numV /( (pow2(1.0 - 1.0/e) + 1.0 - 1.0/e)/epsilon - 1 + 1.0/e );
        double updated_sample_size = log(2/ita) / ( 2 * infVldt / ( pow2(1.0-1.0/e) +1.0-1.0/e )/epsilon -4*deviation_bound_optimize ); // get_updated_sample_size(deviation_bound_optimize, epsilonplus);
	numAdd = updated_sample_size - numR;
//	cout << "total sample size: " << sample_size << ", updated_sample_size: " << updated_sample_size << ", empirical error: " << empirical_error << ", error parameter: " << epsilon << ", epsilonplus: " << epsilonplus << endl;
        auto lowerSelect = infVldt / _numV - empirical_error; //(pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
        auto upperOPT = infVldt / _numV/(1.0 - 1.0/e) + empirical_error; //pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        auto currApprox = lowerSelect / upperOPT;

 //       std::cout << "lower bound: " << (lowerSelect * _numV) << ", upperBound: " << (upperOPT * _numV) << std::endl;
        std::cout << "-->RADESIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox << " target appro. ratio: " << approx - epsilon << 
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';

//        double avgSize = _hyperGraph.HyperedgeAvg();
	if(numAdd <=0 ) numAdd = sample_size;

        if (currApprox >= approx - epsilon)
        {
//            numR += numAdd; // << (idx-1);
//            std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
//            timerSubsim.get_operation_time();
//            _hyperGraph.BuildRRsets(numR); // R1
//            _numRRsets = _hyperGraph.get_RR_sets_size();
//            time1 += timerSubsim.get_operation_time();
//            auto infSelf = MaxCover(targetSize);
//            time2 += timerSubsim.get_operation_time();

        //    _hyperGraph.BuildRRsets(10*numR); // R1
        //    _numRRsets = _hyperGraph.get_RR_sets_size();
        //    double infSelf = MaxCover(targetSize);
          //  infVldt = _hyperGraphVldt.CalculateInf(_vecSeed);
           // auto upperBound = _boundMin;
          //  currApprox = estinf / infVldt * (1.0 -1.0/e); // upperBound;
            //currApprox =  infVldt / upperBound;
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _hyperGraphVldt.BuildRRsets(numR); // R2
            double estinf = _hyperGraphVldt.CalculateInf(_vecSeed);
            _res.set_influence(infSelf);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(numR);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << _res.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
            return 0;
        }
    }

    return 0.0;
}

int decideMultiple(int ratio, int numRRsets)
{
    int multiple = 1;

    if (numRRsets < 100)
    {
        return multiple;
    }

    if (ratio >= 32)
    {
        multiple = 8;
    }
    else if (ratio >= 16)
    {
        multiple = 4;
    }
    else if (ratio >= 4)
    {
        multiple = 2;
    }
    else
    {
        multiple = 1;
    }

    return multiple;
}

double Alg::subsimWithTrunc(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
    const auto numRbase = size_t(3 * log(1 / delta));
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    double infVldt = 0.0;
    int multiple = 1;

    std::cout << std::endl;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsets(numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        const auto infSelf = MaxCoverOutDegPrority(targetSize);
        time2 += timerSubsim.get_operation_time();
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        infVldt = _hyperGraphVldt.EvalSeedSetInf(connSet, _numRRsets * multiple);
        time4 += timerSubsim.get_operation_time();
        const auto degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        const auto upperDegOPT = upperBound * _numRRsets / _numV;
        const auto lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;

        std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.EvalHyperedgeAvg();
        // if truncRRset is more efficient, increase the size of R2 in next iteration
        int ratio = fullRRSize / truncRRSize;
        multiple = decideMultiple(ratio, _numRRsets);

        if (currApprox >= approx - epsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Time for full RR sets: " << time1  << std::endl;
            std::cout << "==>Time for truncated RR set: " << time4 << std::endl;
            std::cout << "==>Time for greedy: " << time2 << std::endl;
            return 0;
        }
    }

    return 0.0;
}

double Alg::IncreaseR2(std::unordered_set<uint32_t> &connSet, double a, double upperOPT, double targetAppr)
{
    size_t vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
    size_t R1RRsets = _hyperGraph.get_RR_sets_size();
    int multiple = 4;
    double estimateAppr = 0.0;
    double lowerSelect = 0;
    int maxMultiple = 3;
    double infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
    double degVldt = infVldt * vldtRRsets / _numV;
    lowerSelect = (pow2(sqrt(degVldt * multiple + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0) ;
    estimateAppr = (lowerSelect / (multiple * vldtRRsets)) / (upperOPT / R1RRsets);

    if (estimateAppr < targetAppr)
    {
        return 0.0;
    }

    _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, vldtRRsets * multiple);
    infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
    vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
    degVldt = infVldt *  vldtRRsets / _numV;
    lowerSelect = (pow2(sqrt(degVldt + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0);
    double newAppr = (lowerSelect / vldtRRsets) / (upperOPT / R1RRsets);
    return (newAppr > targetAppr) ? newAppr : 0.0;
}

double Alg::FindRemSet(const int targetSize, const double epsilon, const double targetEpsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    size_t subSeedSetSize = _vecSeed.size();
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    // delta for upper bound on the number of RR sets
    const double delta_upper = delta / 3.0;
    const double alpha = sqrt(log(3.0 / delta_upper));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize - subSeedSetSize) + log(3.0 / delta_upper)));
    //const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    //const auto numRbase = size_t(_baseNumRRsets);
    //double numRbase = log(2/ita)*(2/pow2(epsilon)); //size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    double numRbase = log(2/ita)/(4*(epsilon)); //size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto maxNumR = size_t(2.0 * _numV * pow2(alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    double infVldt = 0.0;
    double currApprox = 0.0;
    double infSelf = 0.0;
    int multiple = 1;
    auto numR = numRbase;
    double numAdd = 0;
    std::unordered_set<uint32_t> subSeedSet(_vecSeed.begin(), _vecSeed.end());
    std::vector<uint32_t> vecSubSeed(_vecSeed.begin(), _vecSeed.end());

    for (auto idx = 1; idx <= numIter; idx++)
    {
	numR += numAdd;
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsetsEarlyStop(subSeedSet, numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        infSelf = MaxCoverIMSentinel(vecSubSeed, targetSize);
        time2 += timerSubsim.get_operation_time();
    //    std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
    //    infVldt = _hyperGraph.EvalSeedSetInf(connSet);
        time4 += timerSubsim.get_operation_time();
        const auto degVldt = infSelf; // * multiple * _numRRsets / _numV;
	infVldt = degVldt;
	cout << "degVldt: " << degVldt << endl;
        auto upperBound = _boundMin;
	sample_size = _numRRsets;
	VS = degVldt;

	cout << "infVldt: " << infVldt << endl;
	cout << "infVldt/_numRRset: " << infVldt/_numRRsets << endl;
	cout << "infVldt/_numV: " << infVldt/_numV << endl;
//		        for(int i=0; i<sims.size(); i++) std::cout << i << " " << sims[i] << endl;
        double deviation_bound_optimize = (get_rademacher_bound_optimize());
        double empirical_error = get_max_deviation_bound(deviation_bound_optimize, infVldt/_numV);
	double epsilonplus = infVldt / _numV /( (pow2(1.0 - 1.0/e) + 1.0 - 1.0/e)/epsilon - 1 + 1.0/e );
        double updated_sample_size = get_updated_sample_size(deviation_bound_optimize, epsilonplus);
	numAdd = updated_sample_size - numR;
	cout << "total sample size: " << sample_size << "updated_sample_size: " << updated_sample_size << ", empirical error: " << empirical_error << ", error parameter: " << epsilon << ", epsilonplus: " << epsilonplus << endl;
        const auto lowerSelect = infVldt / _numV - empirical_error; //(pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
        const auto upperOPT = infVldt / _numV/(1.0 - 1.0/e) + empirical_error; //pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;

        std::cout << "lower bound: " << (lowerSelect * _numV) << ", upperBound: " << (upperOPT * _numV) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
		   ", targetEpsilon: " << targetEpsilon << 
		   ", target appro. ratio: " << approx - targetEpsilon <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
 //       double fullRRSize = _hyperGraph.HyperedgeAvg();
	if(numAdd <=0 ) numAdd = 0.5 * sample_size;
//        double truncRRSize = _hyperGraphVldt.EvalHyperedgeAvg();

        // if truncRRset is more efficient, increase the size of R2 in next iteration
//        int ratio = fullRRSize / truncRRSize;
//        multiple = decideMultiple(ratio, _numRRsets);

        if (currApprox >= approx - targetEpsilon)
        {
//            numR += numAdd; // << (idx-1);
  //          std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
    //        timerSubsim.get_operation_time();
  //          _hyperGraph.BuildRRsets(numR); // R1
 //           _numRRsets = _hyperGraph.get_RR_sets_size();
//            time1 += timerSubsim.get_operation_time();
//            auto infSelf = MaxCover(targetSize);
//            time2 += timerSubsim.get_operation_time();

            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets);
            std::cout << "==>Time for full RR in IM-Sentinel phase: " << time1  << std::endl;
            std::cout << "==>Time for truncated RR in IM-Sentinel phase: " << time4 << std::endl;
            std::cout << "==>Time for greedy in IM-Sentinel phase: " << time2 << std::endl;
            std::cout << "==>Influence via R2 in IM-Sentinel phase: " << infVldt << ", time: " << _res.get_running_time() << '\n';
    //        _hyperGraphVldt.BuildRRsets(numR); // R2
    //        double estinf = _hyperGraphVldt.CalculateInf(_vecSeed);
            _res.set_influence(infSelf);
            return 0;
        }
    }

    return 0.0;
}

/*
double Alg::FindDynamSub(const int totalTargetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double x = (1.0 - 1.0 / totalTargetSize);
    const int minSubSize = ceil(log(1 - epsilon) / log(x));
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, totalTargetSize) + log(6.0 / delta)));
    const auto numRbasePrevious = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize);
    const auto numRbase = size_t(_baseNumRRsets);
    
    // the successful probability of at least 1-delta/3
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 6.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    int multiple = 1;
    double infVldt = 0.0;
    bool firstRound = true;


    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        //build R1
        _hyperGraph.BuildRRsets(numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        _vecVldtInf.clear();
        int targetSize = firstRound ? (totalTargetSize / 4) : (totalTargetSize / 8);
        firstRound = false;
        const auto infSelf = MaxCoverSentinelSet(targetSize, totalTargetSize);
        time2 += timerSubsim.get_operation_time();
        std::vector<double> vecAppro(_vecSeed.size());
        int lastPos = 0;
        bool found = false;

        if (_vecSeed.size() < targetSize)
        {
            //low influence
            continue;
        }

        double calcAppr = 0.0;

        for (int i = _vecSeed.size() - 1; i >= 0; i--)
        {
            infVldt = _vecVldtInf[i];
            double lowerDeg = infVldt;
            double upperDeg = _boundMin;
            double a = log(numIter * 6.0  / delta);
            double lower = pow2(sqrt(lowerDeg + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0;
            double upper = pow2(sqrt(upperDeg + a1 / 2.0) + sqrt(a1 / 2.0));
            upper = (upper > numR) ? numR : upper;
            vecAppro[i] = lower / upper;
            calcAppr = (1 - pow(x, i + 1) - epsilon) * 1.2;

            if (vecAppro[i] > calcAppr)
            {
                found = true;
                lastPos = i;
                break;
            }
        }


        if (!found)
        {
            lastPos = 0;
        }

        size_t setSize = (lastPos + 1);
        setSize = setSize > 10 ? setSize : 10;
        setSize = (setSize > targetSize) ? targetSize : setSize;
        setSize = (setSize < minSubSize) ? minSubSize : setSize;
        std::vector<uint32_t> dynSeedSet(_vecSeed.begin(), _vecSeed.begin() + setSize);
        _vecSeed.clear();
        _vecSeed.assign(dynSeedSet.begin(), dynSeedSet.end());
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        _hyperGraphVldt.RefreshHypergraph();
        _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, _numRRsets * multiple);
        infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
        time4 += timerSubsim.get_operation_time();
        double degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        double upperDegOPT = upperBound * _numRRsets / _numV;
        double lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;

        if (lowerSelect < 0)
        {
            lowerSelect = 1.0 * _vecSeed.size() / _numV * _numRRsets * multiple;
        }

        double upperOPT = pow2(sqrt(upperDegOPT + a1 / 2.0) + sqrt(a1 / 2.0));
        upperOPT = (upperOPT > _numRRsets) ? _numRRsets : upperOPT;
        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect * _numV / (_numRRsets)) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        const double approx = 1 - pow(x, _vecSeed.size());
        double targetAppr = approx - epsilon;

        if (currApprox >= targetAppr)
        {
            goto succ;
        }

        if (_numRRsets < 100)
        {
            continue;
        }

        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.HyperedgeAvg();

        if (fullRRSize / truncRRSize < 2)
        {
            continue;
        }

        double lowerThreshold = (upperOPT * _numV / _numRRsets) * targetAppr;

        if ((1.0 * infVldt / multiple) > lowerThreshold && lowerThreshold > 0)
        {
            double newAppr = IncreaseR2(connSet, a2, upperOPT, targetAppr);
            time4 += timerSubsim.get_operation_time();

            if (newAppr > targetAppr)
            {
                std::cout << "increase R2 successfully" << std::endl;
                infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
                goto succ;
            }
        }
    }

succ:
    std::cout << "==>Time for full RR in SentinelSet phase: " << time1  << std::endl;
    std::cout << "==>Time for truncated RR in SentinelSet phase: " << time4 << std::endl;
    std::cout << "==>Time for greedy in SentinelSet phase: " << time2 << std::endl;
    std::cout << "==>size of sentinel set: " << _vecSeed.size() << ", inf: " << infVldt << std::endl;
    std::cout << "==>total time for SentinelSet phase: " << timerSubsim.get_total_time() << std::endl;
    return 0.0;
}
*/
double Alg::FindDynamSub(const int totalTargetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double x = (1.0 - 1.0 / totalTargetSize);
    const int minSubSize = ceil(log(1 - epsilon) / log(x));
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, totalTargetSize) + log(6.0 / delta)));
    const auto numRbasePrevious = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize);
    const auto numRbase = size_t(_baseNumRRsets);
    
    // the successful probability of at least 1-delta/3
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 6.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    int multiple = 1;
    double infVldt = 0.0;
    bool firstRound = true;


    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        //build R1
        _hyperGraph.BuildRRsets(numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        _vecVldtInf.clear();
        int targetSize = firstRound ? (totalTargetSize / 4) : (totalTargetSize / 8);
        firstRound = false;
        const auto infSelf = MaxCoverSentinelSet(targetSize, totalTargetSize);
        time2 += timerSubsim.get_operation_time();
        std::vector<double> vecAppro(_vecSeed.size());
        int lastPos = 0;
        bool found = false;

        if (_vecSeed.size() < targetSize)
        {
            //low influence
            continue;
        }

        double calcAppr = 0.0;

        for (int i = _vecSeed.size() - 1; i >= 0; i--)
        {
            infVldt = _vecVldtInf[i];
            double lowerDeg = infVldt;
            double upperDeg = _boundMin;
            double a = log(numIter * 6.0  / delta);
            double lower = pow2(sqrt(lowerDeg + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0;
            double upper = pow2(sqrt(upperDeg + a1 / 2.0) + sqrt(a1 / 2.0));
            upper = (upper > numR) ? numR : upper;
            vecAppro[i] = lower / upper;
            calcAppr = (1 - pow(x, i + 1) - epsilon) * 1.2;

            if (vecAppro[i] > calcAppr)
            {
                found = true;
                lastPos = i;
                break;
            }
        }


        if (!found)
        {
            lastPos = 0;
        }

        size_t setSize = (lastPos + 1);
        setSize = setSize > 10 ? setSize : 10;
        setSize = (setSize > targetSize) ? targetSize : setSize;
        setSize = (setSize < minSubSize) ? minSubSize : setSize;
        std::vector<uint32_t> dynSeedSet(_vecSeed.begin(), _vecSeed.begin() + setSize);
        _vecSeed.clear();
        _vecSeed.assign(dynSeedSet.begin(), dynSeedSet.end());
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
//        _hyperGraphVldt.RefreshHypergraph();
//        _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, _numRRsets * multiple);
//        time4 += timerSubsim.get_operation_time();
        double degVldt = infSelf * _numRRsets / _numV;
        auto upperBound = _boundMin;
        infVldt = degVldt; //_hyperGraphVldt.CalculateInfEarlyStop();
        //infVldt = _hyperGraphVldt.CalculateInfEarlyStop();

//        double upperDegOPT = upperBound * _numRRsets / _numV;
 //       double lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
	sample_size = _numRRsets;
	VS = degVldt;

//		        for(int i=0; i<sims.size(); i++) std::cout << i << " " << sims[i] << endl;
        double deviation_bound_optimize = (get_rademacher_bound_optimize());
        double empirical_error = get_max_deviation_bound(deviation_bound_optimize, infVldt);
	cout << "infVldt: " << infVldt << endl;
	cout << "total sample size: " << sample_size << ", empirical error: " << empirical_error << ", error parameter: " << epsilon << endl;
        double lowerSelect = infVldt / _numRRsets - empirical_error; //(pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
        double upperOPT = infVldt / _numRRsets/(1 - 1/e) + empirical_error; //pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));


        if (lowerSelect < 0)
        {
            lowerSelect = 1.0 * _vecSeed.size() / _numV * _numRRsets * multiple;
        }

     //   double upperOPT = pow2(sqrt(upperDegOPT + a1 / 2.0) + sqrt(a1 / 2.0));
        upperOPT = (upperOPT > _numRRsets) ? _numRRsets : upperOPT;
        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect * _numV ) << ", upperBound: " << (upperOPT * _numV ) << std::endl;
        std::cout << "-->SUBSIM (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        const double approx = 1 - pow(x, _vecSeed.size());
        double targetAppr = approx - epsilon;

        if (currApprox >= targetAppr)
        {
            goto succ;
        }

        if (_numRRsets < 100)
        {
            continue;
        }

        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.HyperedgeAvg();

        if (fullRRSize / truncRRSize < 2)
        {
            continue;
        }

        double lowerThreshold = (upperOPT * _numV / _numRRsets) * targetAppr;

//        if ((1.0 * infVldt / multiple) > lowerThreshold && lowerThreshold > 0)
//        {
        //    double newAppr = IncreaseR2(connSet, a2, upperOPT, targetAppr);
        //    time4 += timerSubsim.get_operation_time();

//            if (newAppr > targetAppr)
  //          {
      //          std::cout << "increase R2 successfully" << std::endl;
      //          infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
//                goto succ;
//            }
//        }
    }

succ:
    std::cout << "==>Time for full RR in SentinelSet phase: " << time1  << std::endl;
    std::cout << "==>Time for truncated RR in SentinelSet phase: " << time4 << std::endl;
    std::cout << "==>Time for greedy in SentinelSet phase: " << time2 << std::endl;
    std::cout << "==>size of sentinel set: " << _vecSeed.size() << ", inf: " << infVldt << std::endl;
    std::cout << "==>total time for SentinelSet phase: " << timerSubsim.get_total_time() << std::endl;
    return 0.0;
}


double Alg::subsimWithHIST(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    _baseNumRRsets = 3 * log(1 / delta);
    ita = 1.0/_numV;

    graphNumV = _numV;
    seedSizeK = targetSize;
    std::cout << std::endl;
    std::cout << "Sentinel Set Selection Phase" << std::endl;
    FindDynamSub(targetSize, epsilon / 2, delta / 2);
    _hyperGraph.RefreshHypergraph();
    _hyperGraphVldt.RefreshHypergraph();

    std::cout << std::endl;
    std::cout << "IM-Sentinel Phase" << std::endl;
    FindRemSet(targetSize, epsilon / 2, epsilon, delta / 2);
    _res.set_running_time(timerSubsim.get_total_time());
    return 0.0;
}
