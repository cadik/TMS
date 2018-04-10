// Pallas Solver
// Copyright 2015. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: ryan.latture@gmail.com (Ryan Latture)
//
// Adapted from Value Semantics and Concept-based Polymorphism by Sean Parent
// Copyright 2013 Adobe Systems Incorporated
// Distributed under the MIT License (see license at http://stlab.adobe.com/licenses.html)
// original code can be found at:
// https://raw.githubusercontent.com/wiki/sean-parent/sean-parent.github.io/presentations/2013-09-06-inheritance/value-semantics-unique.cpp

#ifndef PALLAS_HISTORY_CONCEPT_H
#define PALLAS_HISTORY_CONCEPT_H

#include <memory>
#include "pallas/types.h"

namespace pallas {

    class HistoryConcept {
    public:
        template <typename T>
        HistoryConcept(T x) : self_(new model<T>(std::move(x)))
        { }

        HistoryConcept(const HistoryConcept& x) : self_(x.self_->copy_())
        { }
        HistoryConcept(HistoryConcept&&) noexcept = default;

        HistoryConcept& operator=(const HistoryConcept& x)
        { HistoryConcept tmp(x); *this = std::move(tmp); return *this; }
        HistoryConcept& operator=(HistoryConcept&&) noexcept = default;

        friend void dump(const HistoryConcept &x, HistoryWriter& writer)
        { x.self_->dump_(writer); }

    private:
        struct concept_t {
            virtual ~concept_t() = default;
            virtual concept_t* copy_() const = 0;
            virtual void dump_(HistoryWriter&) const = 0;
        };
        template <typename T>
        struct model : concept_t {
            model(T x) : data_(std::move(x)) { }
            concept_t* copy_() const { return new model(*this); }
            void dump_(HistoryWriter &writer) const
            { dump(data_, writer); }

            T data_;
        };

        std::unique_ptr<const concept_t> self_;
    };

    /**
     * @brief Stores series history outputs that characterize the progression of optimization.
     */
    using HistorySeries = std::vector<HistoryConcept>;

    /**
     * @brief Writes the optimization history series to the stream contained in the `HistoryWriter`.
     *
     * @param history HistorySeries. Series of history data to write.
     * @param writer HistoryWriter. Dumps the history series to the contained stream.
     */
    void dump(const HistorySeries& history, HistoryWriter& writer);
} // namespace pallas

#endif //PALLAS_HISTORY_CONCEPT_H
