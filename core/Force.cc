#include "core/Force.h"

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>

using namespace Glucose;

Force::Force(const ClauseAllocator& ca, const vec<CRef>& cs)
	: _ca(ca), _cs(cs)
{
}

void Force::execute(vec<Var>& order)
{
	int num_vars = 0;
	// Map solver variables to FORCE variables
	std::unordered_map<Var, Var> renumber;
    for (int i = 0; i < _cs.size(); i++) {
        const Clause& c = _ca[_cs[i]];
        for (int j = 0; j < c.size(); j++) {
        	if (renumber.find(var(c[j])) == renumber.end()) {
        		renumber.emplace(var(c[j]), num_vars);
        		num_vars++;
        	}
        }
    }

    std::vector<int> idx_of_var;
    idx_of_var.reserve(num_vars);
    std::vector<Var> var_of_idx;
    var_of_idx.reserve(num_vars);
    for (int i = 0; i < num_vars; i++) {
    	idx_of_var.push_back(i);
    	var_of_idx.push_back(i);
    }

    // Random initial order
    unsigned int seed = 1;
    srand(seed);
    random_shuffle(idx_of_var.begin(), idx_of_var.end());

    long best_sum_span = LONG_MAX;
    std::vector<Var> best_idx_or_var(idx_of_var);

    std::vector<VarInfo> info;
    int num = 5 * std::log(num_vars);
    for (int i = 0; i < num; i++) {
    	info.assign(num_vars, {0.0, 0});

    	for (int j = 0; j < _cs.size(); j++) {
    		const Clause& c = _ca[_cs[j]];
    		int sum = 0;
    		for (int k = 0; k < c.size(); k++) {
    			sum += idx_of_var[renumber[var(c[k])]];
    		}
    		double cog_of_clause = (double) sum / c.size();
    		for (int k = 0; k < c.size(); k++) {
    			VarInfo& vi = info[renumber[var(c[k])]];
    			vi.sum_cog += cog_of_clause;
    			vi.num_occurs++;
    		}
    	}

    	// Center of gravity
    	std::vector<double> cog_of_var;
    	cog_of_var.reserve(num_vars);
    	for (int j = 0; j < num_vars; j++) {
    		VarInfo& vi = info[j];
    		cog_of_var.push_back(vi.sum_cog / vi.num_occurs);
    	}

    	std::sort(var_of_idx.begin(), var_of_idx.end(), VarLessThan<double>(cog_of_var));

    	for (int j = 0; j < num_vars; j++) {
    		idx_of_var[var_of_idx[j]] = j;
    	}

    	// Compute the sum of span
    	long sum_span = 0;
    	for (int j = 0; j < _cs.size(); j++) {
    		const Clause& c = _ca[_cs[j]];
    		int bottom = idx_of_var[renumber[var(c[0])]];
    		int top = bottom;
    		for (int k = 1; k < c.size(); k++) {
    			int idx = idx_of_var[renumber[var(c[k])]];
    			if (bottom > idx) {
    				bottom = idx;
    			}
    			if (top < idx) {
    				top = idx;
    			}
    		}
    		sum_span += (top - bottom);
    	}
    	if (sum_span < best_sum_span) {
    		best_sum_span = sum_span;
    		best_idx_or_var.assign(idx_of_var.begin(), idx_of_var.end());
    	}
    }

    // Map FORCE variables to solver variables
    std::vector<Var> rev_renumber(num_vars, 0);
    for (const auto& it : renumber) {
    	assert(it.second < rev_renumber.size());
    	rev_renumber[it.second] = it.first;
    }

    order.capacity(num_vars);
    for (const auto& var : var_of_idx) {
    	order.push(rev_renumber[var]);
    }
}
