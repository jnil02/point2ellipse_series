#pragma once

#include <symengine/expression.h>   // SymEngine symbolic expressions.
#include <functional>               // std::function for generating functions.
#include <vector>                   // std::vector for cached terms.
#include <memory>                   // std::shared_ptr for polymorphism.

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
