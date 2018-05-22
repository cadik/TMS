//
// Created by ryan on 5/9/15.
//

#ifndef PALLAS_TYPES_H
#define PALLAS_TYPES_H

#include "ceres/ceres.h"
#include "ceres/types.h"
#include "Eigen/Core"
#include "rapidjson/writer.h"

namespace pallas {

    using GradientCostFunction = ceres::FirstOrderFunction;
    using LeastSquaresCostFunction = ceres::CostFunction;

    using LeastSquaresProblem = ceres::Problem;
    using GradientProblem = ceres::GradientProblem;

    using GradientLocalMinimizer = ceres::GradientProblemSolver;
    using LeastSquaresLocalMinimizer = ceres::Solver;

    using Vector = ceres::Vector;
    using VectorRef = ceres::VectorRef;
    using ConstVectorRef = ceres::ConstVectorRef;

    using LineSearchDirectionType = ceres::LineSearchDirectionType;

    using TerminationType = ceres::TerminationType;

    using Vector2d = Eigen::Vector2d;

    using HistoryWriter = rapidjson::Writer<rapidjson::StringBuffer>;

    enum MutationStrategyType {
        BEST_1,
        RAND_1,
        RAND_TO_BEST_1,
        BEST_2,
        RAND_2,
    };

    enum CrossoverStrategyType {
        BINOMIAL,
        EXPONENTIAL,
    };

    enum CostFunctionType {
        GRADIENT,
        LEAST_SQUARES,
    };

    enum PopulationInitializationType {
        LATIN_HYPERCUBE,
        RANDOM,
    };

    enum CoolingScheduleType {
        BOLTZMANN,
        CAUCHY,
        FAST,
    };

    const char* CoolingScheduleTypeToString(CoolingScheduleType type);
    bool StringToCoolingScheduleType(std::string value, CoolingScheduleType* type);
    
    const char* MutationStrategyTypeToString(MutationStrategyType type);
    bool StringToMutationStrategyType(std::string value, MutationStrategyType* type);

    const char* CrossoverStrategyTypeToString(CrossoverStrategyType type);
    bool StringToCrossoverStrategyType(std::string value, CrossoverStrategyType* type);

    const char* PopulationInitializationTypeToString(PopulationInitializationType type);
    bool StringToPopulationInitializationType(std::string value, PopulationInitializationType* type);

} // namespace pallas

#endif //PALLAS_TYPES_H
