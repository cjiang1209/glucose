#ifndef Force_h
#define Force_h

#include <vector>
#include <cstdlib>

#include "mtl/Vec.h"
#include "core/SolverTypes.h"

namespace Glucose {

class Force
{
private:
	struct VarInfo
	{
		double sum_cog;
		int num_occurs;
	};

	template<typename T>
	struct VarLessThan
	{
		const std::vector<T> _weights;

		VarLessThan(const std::vector<T> weights)
			: _weights(weights)
		{
		}

		bool operator()(int x, int y)
		{
			return _weights[x] < _weights[y];
		}
	};

	const ClauseAllocator& _ca;
	const vec<CRef>& _cs;

	template<typename RandomIt>
	void random_shuffle(RandomIt first, RandomIt last)
	{
	    typename std::iterator_traits<RandomIt>::difference_type i, n;
	    n = last - first;
	    for (i = n-1; i > 0; --i) {
	        std::swap(first[i], first[std::rand() % (i+1)]);
	    }
	}

public:
	Force(const ClauseAllocator& ca, const vec<CRef>& cs);
	void execute(vec<Var>& order);
};

}

#endif
