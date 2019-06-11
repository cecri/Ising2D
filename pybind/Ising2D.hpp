#ifndef ISING2D_HPP
#define ISING2D_HPP

#include <cstdint>
#include <vector>
#include <random>
#include <unordered_set>
#include <deque>
/**
 */
class Ising2D
{
private:
	uint32_t nrows_;
	uint32_t ncols_;

public:
	Ising2D(uint32_t nrows, uint32_t ncols)
		: nrows_(nrows), ncols_(ncols)
	{
	}
	uint32_t get_nrows() const
	{
		return nrows_;
	}
	uint32_t get_ncols() const
	{
		return ncols_;
	}
	uint32_t size() const
	{
		return nrows_*ncols_;
	}

	double energy(const std::vector<int8_t>& conf) const
	{
		double res = 0.0;
		for(auto m: all_neighbors())
		{
			res += -conf[m.first]*conf[m.second];
		}
		return res;
	}

	inline uint32_t to_idx(uint32_t row, uint32_t col) const
	{
		return col*nrows_ + row;
	}

	inline std::pair<uint32_t, uint32_t> to_coord(uint32_t idx) const
	{
		return std::make_pair(idx%nrows_, idx/nrows_);
	}

	std::vector<uint32_t> neighbors(uint32_t idx) const
	{
		std::vector<uint32_t> res;
		auto m = to_coord(idx);
		if(m.first >= 1)
			res.emplace_back(to_idx(m.first-1, m.second));
		if(m.first + 1 < nrows_)
			res.emplace_back(to_idx(m.first+1, m.second));
		if(m.second >= 1)
			res.emplace_back(to_idx(m.first, m.second-1));
		if(m.second + 1 < ncols_)
			res.emplace_back(to_idx(m.first, m.second+1));
		return res;
	}

	double magnetization(const std::vector<int8_t>& conf) const
	{
		return std::accumulate(conf.begin(), conf.end(), 0.0);
	}
	
	std::vector<std::pair<uint32_t, uint32_t> > all_neighbors() const
	{
		std::vector<std::pair<uint32_t, uint32_t> > res;
		for(uint32_t i = 0; i < nrows_; i++)
		{
			for(uint32_t j = 0; j < ncols_-1; j++)
			{
				auto n1 = to_idx(i,j);
				auto n2 = to_idx(i,j+1);
				res.emplace_back(n1,n2);
			}
		}
		for(uint32_t i = 0; i < nrows_-1; i++)
		{
			for(uint32_t j = 0; j < ncols_; j++)
			{
				auto n1 = to_idx(i,j);
				auto n2 = to_idx(i+1,j);
				res.emplace_back(n1,n2);
			}
		}
		return res;
	}
};

class WolffSampler
{
private:
	const Ising2D& model_;
	const double beta_;
	std::vector<int8_t> conf_;
	std::default_random_engine re_;

public:
	WolffSampler(const Ising2D& model, const double beta)
		: model_(model), beta_(beta)
	{
	}
	void set_seed(uint32_t seed)
	{
		re_.seed(seed);
	}

	void randomize_conf()
	{
		std::uniform_int_distribution<> uid(0, 1);
		conf_.resize(model_.size(), 0);
		for(auto& x : conf_)
		{
			x = 1-2*uid(re_);
		}
	}

	std::vector<int8_t> get_conf() const
	{
		return conf_;
	}

	void sweep()
	{
		std::unordered_set<uint32_t> cluster;
		std::deque<uint32_t> queue;

		auto n = model_.size();
		std::uniform_int_distribution<> uid(0, n-1);
		std::uniform_real_distribution<> urd(0.0, 1.0);

		double prob = 1.0 - std::exp(-2.0*beta_);
		uint32_t ini = uid(re_);
		queue.push_back(ini);

		while(!queue.empty())
		{
			uint32_t k = queue.front();
			queue.pop_front();
			cluster.emplace(k);
			for(auto x: model_.neighbors(k))
			{
				if(conf_[x] == conf_[k])
				{
					if((cluster.count(x) == 0) &&
							(urd(re_) < prob))
					{
						queue.push_back(x);
					}
				}
			}
		}

		for(auto x: cluster)
		{
			conf_[x] *= -1;
		}
	}
	
};
#endif//ISING2D_HPP
