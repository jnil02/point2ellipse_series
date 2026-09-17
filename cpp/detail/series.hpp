#pragma once

#include <functional>               // std::function for generating functions.
#include <vector>                   // std::vector for cached terms.
#include <memory>                   // std::shared_ptr for polymorphism.
#include <optional>                 // std::optional for clean sentinel in TSeries.
#include <stdexcept>                // std::runtime_error.
#include <gmpxx.h>                  // mpq_class for scalar multiply in TSeriesBase.

// ---------------------------------------------------------------------------
// Generic template series classes keyed on value type T.
// T must support: T operator+(T, const T&), T operator*(const T&, const T&),
//                 T operator*(T, const mpq_class&), and default-construct to zero.
// ---------------------------------------------------------------------------

template<typename T>
class TSeriesBase {
public:
	virtual T getItem(int n) = 0;
	virtual std::shared_ptr<TSeriesBase<T>> operator*(const mpq_class &s) = 0;
	virtual std::shared_ptr<TSeriesBase<T>> operator*(const std::shared_ptr<TSeriesBase<T>> &other) = 0;
	virtual std::shared_ptr<TSeriesBase<T>> operator+(const std::shared_ptr<TSeriesBase<T>> &other) = 0;
	virtual ~TSeriesBase() = default;
};

template<typename T>
class TSeries : public TSeriesBase<T>, public std::enable_shared_from_this<TSeries<T>> {
	std::function<T(int)> gen;
	std::vector<std::optional<T>> items;
public:
	explicit TSeries(std::function<T(int)> gen)
			: gen(std::move(gen)), items(16, std::nullopt) {}

	T getItem(int n) override {
		while (items.size() <= static_cast<size_t>(n))
			items.resize(items.size() * 2, std::nullopt);
		if (!items[n])
			items[n] = gen(n);
		return *items[n];
	}

	std::shared_ptr<TSeriesBase<T>> operator*(const mpq_class &s) override {
		auto self = this->shared_from_this();
		return std::make_shared<TSeries<T>>([self, s](int n) { return self->getItem(n) * s; });
	}

	std::shared_ptr<TSeriesBase<T>> operator*(const std::shared_ptr<TSeriesBase<T>> &other) override {
		auto self = this->shared_from_this();
		return std::make_shared<TSeries<T>>([self, other](int n) {
			T sum{};
			for (int l = 0; l <= n; ++l)
				sum = sum + self->getItem(l) * other->getItem(n - l);
			return sum;
		});
	}

	std::shared_ptr<TSeriesBase<T>> operator+(const std::shared_ptr<TSeriesBase<T>> &other) override {
		auto self = this->shared_from_this();
		return std::make_shared<TSeries<T>>([self, other](int n) {
			return self->getItem(n) + other->getItem(n);
		});
	}
};

template<typename T>
class TSeriesFactor : public TSeriesBase<T>, public std::enable_shared_from_this<TSeriesFactor<T>> {
	T val;
public:
	explicit TSeriesFactor(T val) : val(std::move(val)) {}

	T getItem(int n) override { return (n == 0) ? val : T{}; }

	std::shared_ptr<TSeriesBase<T>> operator*(const mpq_class &s) override {
		return std::make_shared<TSeriesFactor<T>>(val * s);
	}

	std::shared_ptr<TSeriesBase<T>> operator*(const std::shared_ptr<TSeriesBase<T>> &other) override {
		std::shared_ptr<TSeriesBase<T>> self = this->shared_from_this();
		return std::make_shared<TSeries<T>>([self, other](int n) {
			T sum{};
			for (int l = 0; l <= n; ++l)
				sum = sum + self->getItem(l) * other->getItem(n - l);
			return sum;
		});
	}

	std::shared_ptr<TSeriesBase<T>> operator+(const std::shared_ptr<TSeriesBase<T>> &) override {
		throw std::runtime_error("Cannot add to TSeriesFactor.");
	}
};

template<typename T>
class TSeriesEmpty : public TSeriesBase<T> {
public:
	T getItem(int) override { throw std::runtime_error("Empty series: no item."); }

	std::shared_ptr<TSeriesBase<T>> operator*(const mpq_class &) override {
		throw std::runtime_error("Empty series: no multiplication.");
	}

	std::shared_ptr<TSeriesBase<T>> operator*(const std::shared_ptr<TSeriesBase<T>> &) override {
		throw std::runtime_error("Empty series: no multiplication.");
	}

	std::shared_ptr<TSeriesBase<T>> operator+(const std::shared_ptr<TSeriesBase<T>> &other) override {
		return other;
	}
};
