
#include <stdexcept>                // std::runtime_error for exceptions.

#include "series.hpp"

// ---------------------------------- Series -----------------------------------

Series::Series(std::function<SymEngine::Expression(int)> gen)
		: m_basic(SymEngine::integer(123456789)), gen(std::move(gen)),
		items(16, m_basic) {}

SymEngine::Expression Series::getItem(int n) {
	while (items.size() <= static_cast<size_t>(n))
		items.resize(items.size() * 2, m_basic);
	if (items[n] == m_basic)
		items[n] = gen(n);
	return items[n];
}

std::shared_ptr<SeriesBase> Series::operator*(const SymEngine::Expression& other) {
	auto self = shared_from_this();
	return std::make_shared<Series>([self, other](int n) {
		return self->getItem(n) * other;
	});
}

std::shared_ptr<SeriesBase> Series::operator*(const std::shared_ptr<SeriesBase>& other) {
	auto self = shared_from_this();
	if (auto otherSeries = std::dynamic_pointer_cast<Series>(other)) {
		return std::make_shared<Series>([self, otherSeries](int n) {
			SymEngine::Expression sum = SymEngine::Expression(0);
			for (int l = 0; l <= n; ++l)
				sum = sum + self->getItem(l) * otherSeries->getItem(n - l);
			return sum;
		});
	}
	else if (auto factor = std::dynamic_pointer_cast<SeriesFactor>(other)) {
		return std::make_shared<Series>([self, factor](int n) {
			return self->getItem(n) * factor->getItem(0);
		});
	}
	else {
		throw std::runtime_error("Unsupported multiplication with non-Series object.");
	}
}

std::shared_ptr<SeriesBase> Series::operator+(const std::shared_ptr<SeriesBase>& other) {
	auto self = shared_from_this();
	return std::make_shared<Series>([self, other](int n) {
		return self->getItem(n) + other->getItem(n);
	});
}

// ------------------------------- SeriesFactor --------------------------------

SeriesFactor::SeriesFactor(SymEngine::Expression a) : a(std::move(a)) {}

SymEngine::Expression SeriesFactor::getItem(int n) {
	return (n == 0) ? a : SymEngine::Expression(0);
}

std::shared_ptr<SeriesBase> SeriesFactor::operator*(const SymEngine::Expression& other) {
	return std::make_shared<SeriesFactor>(a * other);
}

std::shared_ptr<SeriesBase> SeriesFactor::operator*(const std::shared_ptr<SeriesBase>& other) {
	if (auto otherSeries = std::dynamic_pointer_cast<Series>(other)) {
		return (*otherSeries) * a;
	}
	else if (auto otherFactor = std::dynamic_pointer_cast<SeriesFactor>(other)) {
		return std::make_shared<SeriesFactor>(a * otherFactor->a);
	}
	else {
		throw std::runtime_error("Unsupported multiplication with non-Series object.");
	}
}

std::shared_ptr<SeriesBase> SeriesFactor::operator+(const std::shared_ptr<SeriesBase>&) {
	throw std::runtime_error("Cannot add to SeriesFactor.");
}

// -------------------------------- SeriesEmpty --------------------------------

SymEngine::Expression SeriesEmpty::getItem(int) {
	throw std::runtime_error("No series item available for empty series.");
}

std::shared_ptr<SeriesBase> SeriesEmpty::operator*(const SymEngine::Expression&) {
	throw std::runtime_error("Multiplication not defined for empty series.");
}

std::shared_ptr<SeriesBase> SeriesEmpty::operator*(const std::shared_ptr<SeriesBase>&) {
	throw std::runtime_error("Multiplication not defined for empty series.");
}

std::shared_ptr<SeriesBase> SeriesEmpty::operator+(const std::shared_ptr<SeriesBase>& other) {
	return other;
}
