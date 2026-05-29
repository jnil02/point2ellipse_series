#pragma once

#include <symengine/expression.h>   // SymEngine symbolic expressions.
#include <functional>               // std::function for generating functions.
#include <vector>                   // std::vector for cached terms.
#include <memory>                   // std::shared_ptr for polymorphism.
#include <optional>                 // std::optional for clean sentinel in TSeries.
#include <stdexcept>                // std::runtime_error.
#include <gmpxx.h>                  // mpq_class for scalar multiply in TSeriesBase.

/**
 * Abstract base class for lazy series arithmetics.
 */
class SeriesBase {
public:
	virtual SymEngine::Expression getItem(int n) = 0;
	/*
	 * The reason we use shared_ptr return values is that, due to lazy
	 * evaluation, the new series need to keep a reference to the original
	 * series so series arithmetic essentially builds up a graph of underlying
	 * series. Without the shared_ptr, we would have to make deep copies in
	 * every operation.
	 *
	 * Use const reference for arguments such that we can accept temporaries
	 * and still not get a reference increment until we actually store it.
	 */
	virtual std::shared_ptr<SeriesBase> operator*(const SymEngine::Expression& other) = 0;
	virtual std::shared_ptr<SeriesBase> operator*(const std::shared_ptr<SeriesBase>& other) = 0;
	virtual std::shared_ptr<SeriesBase> operator+(const std::shared_ptr<SeriesBase>& other) = 0;
	virtual ~SeriesBase() = default;
};

/**
 * Regular series with a generating function for the terms.
 * Lazy evaluation and caching of results.
 * enable_shared_from_this since a reference to self is passed to resulting
 * series from series arithmetics.
 */
class Series : public SeriesBase, public std::enable_shared_from_this<Series> {
private:
	SymEngine::RCP<const SymEngine::Basic> m_basic;  // Placeholder for uninitialized terms.
	std::function<SymEngine::Expression(int)> gen;   // Generating function.
	std::vector<SymEngine::Expression> items;        // Cache for computed terms.

public:
	explicit Series(std::function<SymEngine::Expression(int)> gen);
	SymEngine::Expression getItem(int n) override;
	std::shared_ptr<SeriesBase> operator*(const SymEngine::Expression& other) override;
	std::shared_ptr<SeriesBase> operator*(const std::shared_ptr<SeriesBase>& other) override;
	std::shared_ptr<SeriesBase> operator+(const std::shared_ptr<SeriesBase>& other) override;
};

/**
 * Series representation of a factor, i.e. only zeroth item non-zero.
 */
class SeriesFactor : public SeriesBase {
private:
	SymEngine::Expression a; // Constant factor value.

public:
	explicit SeriesFactor(SymEngine::Expression a);
	SymEngine::Expression getItem(int n) override;
	std::shared_ptr<SeriesBase> operator*(const SymEngine::Expression& other) override;
	std::shared_ptr<SeriesBase> operator*(const std::shared_ptr<SeriesBase>& other) override;
	std::shared_ptr<SeriesBase> operator+(const std::shared_ptr<SeriesBase>& other) override;
};

/**
 * Empty series which can be used as an initial value when series are added.
 * Only addition is supported. Multiplication throws runtime exception.
 */
class SeriesEmpty : public SeriesBase {
public:
	SymEngine::Expression getItem(int) override;
	std::shared_ptr<SeriesBase> operator*(const SymEngine::Expression& other) override;
	std::shared_ptr<SeriesBase> operator*(const std::shared_ptr<SeriesBase>& other) override;
	std::shared_ptr<SeriesBase> operator+(const std::shared_ptr<SeriesBase>& other) override;
};

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
