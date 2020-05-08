
/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_PRECONDITIONERFACTORY_HEADER
#define OPM_PRECONDITIONERFACTORY_HEADER

#include <opm/simulators/linalg/OwningBlockPreconditioner.hpp>
#include <opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/amgcpr.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/kamg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/preconditioners.hh>

#include <boost/property_tree/ptree.hpp>

#include <map>
#include <memory>

namespace Opm
{

/// This is an object factory for creating preconditioners.  The
/// user need only interact with the factory through the static
/// methods addStandardPreconditioners() and create(). In addition
/// a user can call the addCreator() static method to add further
/// preconditioners.
template <class Operator, class Comm>
class PreconditionerFactory
{
public:
    /// Linear algebra types.
    using Matrix = typename Operator::matrix_type;
    using Vector = typename Operator::domain_type; // Assuming symmetry: that domain and range types are the same.

    /// The type of pointer returned by create().
    using PrecPtr = std::shared_ptr<Dune::PreconditionerWithUpdate<Vector, Vector>>;

    /// The type of creator functions passed to addCreator().
    using Creator = std::function<PrecPtr(const Operator&, const boost::property_tree::ptree&, const std::function<Vector()>&)>;
    using ParCreator = std::function<PrecPtr(const Operator&, const boost::property_tree::ptree&, const std::function<Vector()>&, const Comm&)>;

    /// Create a new serial preconditioner and return a pointer to it.
    /// \param op    operator to be preconditioned.
    /// \param prm   parameters for the preconditioner, in particular its type.
    /// \param weightsCalculator Calculator for weights used in CPR.
    /// \return      (smart) pointer to the created preconditioner.
    static PrecPtr create(const Operator& op, const boost::property_tree::ptree& prm,
                          const std::function<Vector()>& weightsCalculator = std::function<Vector()>())
    {
        return instance().doCreate(op, prm, weightsCalculator);
    }

    /// Create a new parallel preconditioner and return a pointer to it.
    /// \param op    operator to be preconditioned.
    /// \param prm   parameters for the preconditioner, in particular its type.
    /// \param comm  communication object (typically OwnerOverlapCopyCommunication).
    /// \param weightsCalculator Calculator for weights used in CPR.
    /// \return      (smart) pointer to the created preconditioner.
    static PrecPtr create(const Operator& op, const boost::property_tree::ptree& prm,
                          const std::function<Vector()>& weightsCalculator, const Comm& comm)
    {
        return instance().doCreate(op, prm, weightsCalculator, comm);
    }

    /// Create a new parallel preconditioner and return a pointer to it.
    /// \param op    operator to be preconditioned.
    /// \param prm   parameters for the preconditioner, in particular its type.
    /// \param comm  communication object (typically OwnerOverlapCopyCommunication).
    /// \return      (smart) pointer to the created preconditioner.
    static PrecPtr create(const Operator& op, const boost::property_tree::ptree& prm, const Comm& comm)
    {
        return instance().doCreate(op, prm, std::function<Vector()>(), comm);
    }
    /// Add a creator for a serial preconditioner to the PreconditionerFactory.
    /// After the call, the user may obtain a preconditioner by
    /// calling create() with the given type string as a parameter
    /// contained in the property_tree.
    /// \param type     the type string we want the PreconditionerFactory to
    ///                 associate with the preconditioner.
    /// \param creator  a function or lambda creating a preconditioner.
    static void addCreator(const std::string& type, Creator creator)
    {
        instance().doAddCreator(type, creator);
    }

    /// Add a creator for a parallel preconditioner to the PreconditionerFactory.
    /// After the call, the user may obtain a preconditioner by
    /// calling create() with the given type string as a parameter
    /// contained in the property_tree.
    /// \param type     the type string we want the PreconditionerFactory to
    ///                 associate with the preconditioner.
    /// \param creator  a function or lambda creating a preconditioner.
    static void addCreator(const std::string& type, ParCreator creator)
    {
        instance().doAddCreator(type, creator);
    }

private:
    using CriterionBase
        = Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricMatrixDependency<Matrix, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;

    // Helpers for creation of AMG preconditioner.
    static Criterion amgCriterion(const boost::property_tree::ptree& prm)
    {
        Criterion criterion(15, prm.get<int>("coarsenTarget", 1200));
        criterion.setDefaultValuesIsotropic(2);
        criterion.setAlpha(prm.get<double>("alpha", 0.33));
        criterion.setBeta(prm.get<double>("beta", 1e-5));
        criterion.setMaxLevel(prm.get<int>("maxlevel", 15));
        criterion.setSkipIsolated(prm.get<bool>("skip_isolated", false));
        criterion.setNoPreSmoothSteps(prm.get<int>("pre_smooth", 1));
        criterion.setNoPostSmoothSteps(prm.get<int>("post_smooth", 1));
        criterion.setDebugLevel(prm.get<int>("verbosity", 0));
        return criterion;
    }

    template <typename Smoother>
    static auto amgSmootherArgs(const boost::property_tree::ptree& prm)
    {
        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = prm.get<int>("iterations", 1);
	//int iluwitdh = prm.get<int>("iluwidth", 0);
	//MILU_VARIANT milu = convertString2Milu(prm.get<std::string>("milutype",std::string("ilu")));
	//smootherArgs.setMilu(milu);
	//smootherArgs.setN(iluwitdh);
        // smootherArgs.overlap=SmootherArgs::vertex;
        // smootherArgs.overlap=SmootherArgs::none;
        // smootherArgs.overlap=SmootherArgs::aggregate;
        smootherArgs.relaxationFactor = prm.get<double>("relaxation", 1.0);
        return smootherArgs;
    }
    // template <>
    // static auto amgSmootherArgs< Opm::ParallelOverlappingILU0< Matrix, Vector, Vector, Comm>  >(const boost::property_tree::ptree& prm)
    // {
    // 	using Smoother = Opm::ParallelOverlappingILU0<Matrix, Vector, Vector, Comm>;
    //     using SmootherArgs = typename Dune::Amg::SmootherTraits<  Smoother >::Arguments;
    //     SmootherArgs smootherArgs;
    //     smootherArgs.iterations = prm.get<int>("iterations", 1);
    // 	int iluwitdh = prm.get<int>("iluwidth", 0);
    // 	MILU_VARIANT milu = convertString2Milu( prm.get<std::string>("milutype",std::string("ilu") ) );
    // 	smootherArgs.setMilu(milu);
    // 	smootherArgs.setN(iluwitdh);
    //     // smootherArgs.overlap=SmootherArgs::vertex;
    //     // smootherArgs.overlap=SmootherArgs::none;
    //     // smootherArgs.overlap=SmootherArgs::aggregate;
    //     smootherArgs.relaxationFactor = prm.get<double>("relaxation", 1.0);
    //     return smootherArgs;
    // }

    template <class Smoother>
    static PrecPtr makeAmgPreconditioner(const Operator& op, const boost::property_tree::ptree& prm, bool useKamg = false)
    {
        auto crit = amgCriterion(prm);
        auto sargs = amgSmootherArgs<Smoother>(prm);
	if(useKamg){
	    return std::make_shared<
		Dune::DummyUpdatePreconditioner<
		    Dune::Amg::KAMG< Operator, Vector, Smoother>
		    >
		>(op, crit, sargs,
		  prm.get<size_t>("max_krylov", 1),
		  prm.get<double>("min_reduction", 1e-1)  );
	}else{
            return std::make_shared<Dune::Amg::AMGCPR<Operator, Vector, Smoother>>(op, crit, sargs);
        }
    }

    // Add a useful default set of preconditioners to the factory.
    // This is the default template, used for parallel preconditioners.
    // (Serial specialization below).
    template <class CommArg>
    void addStandardPreconditioners(const CommArg*)
    {
        using namespace Dune;
        using O = Operator;
        using M = Matrix;
        using V = Vector;
        using P = boost::property_tree::ptree;
        using C = Comm;
        doAddCreator("ILU0", [](const O& op, const P& prm, const std::function<Vector()>&, const C& comm) {
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), comm, 0, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("ParOverILU0", [](const O& op, const P& prm, const std::function<Vector()>&, const C& comm) {
            const double w = prm.get<double>("relaxation", 1.0);
            const int n = prm.get<int>("ilulevel", 0);
            // Already a parallel preconditioner. Need to pass comm, but no need to wrap it in a BlockPreconditioner.
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), comm, n, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("ILUn", [](const O& op, const P& prm, const std::function<Vector()>&, const C& comm) {
            const int n = prm.get<int>("ilulevel", 0);
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), comm, n, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("Jac", [](const O& op, const P& prm, const std::function<Vector()>&,
                               const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqJac<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("GS", [](const O& op, const P& prm, const std::function<Vector()>&, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqGS<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("SOR", [](const O& op, const P& prm, const std::function<Vector()>&, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("SSOR", [](const O& op, const P& prm, const std::function<Vector()>&, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("amg", [](const O& op, const P& prm, const std::function<Vector()>&, const C& comm) {
            const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
	    auto crit = amgCriterion(prm);
            if (smoother == "ParOverILU0") {
                using Smoother = Opm::ParallelOverlappingILU0<M, V, V, C>;
                auto sargs = amgSmootherArgs<Smoother>(prm);
                return std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
	    // } else if (smoother == "ILU0") {
	    // 	using Smoother = Dune::BlockPreconditioner<V, V, C, Dune::SeqILU0<M, V, V> >;
            //     auto sargs = amgSmootherArgs<Smoother>(prm);
            //     return std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
	    // } else if (smoother == "Jac") {
	    // 	using Smoother = Dune::BlockPreconditioner<V, V, C, Dune::SeqJac<M, V, V>>;
            //     auto sargs = amgSmootherArgs<Smoother>(prm);
            //     return std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
             } else {
                OPM_THROW(std::invalid_argument, "Properties: No smoother with name " << smoother <<".");
            }
	});
        doAddCreator("cpr", [](const O& op, const P& prm, const std::function<Vector()> weightsCalculator, const C& comm) {
            assert(weightsCalculator);
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, false, Comm>>(op, prm, weightsCalculator, comm);
        });
        doAddCreator("cprt", [](const O& op, const P& prm, const std::function<Vector()> weightsCalculator, const C& comm) {
            assert(weightsCalculator);
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, true, Comm>>(op, prm, weightsCalculator, comm);
        });
    }

    // Add a useful default set of preconditioners to the factory.
    // This is the specialization for the serial case.
    void addStandardPreconditioners(const Dune::Amg::SequentialInformation*)
    {
        using namespace Dune;
        using O = Operator;
        using M = Matrix;
        using V = Vector;
        using P = boost::property_tree::ptree;
        doAddCreator("ILU0", [](const O& op, const P& prm, const std::function<Vector()>&) {
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V>>(
                op.getmat(), 0, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("ParOverILU0", [](const O& op, const P& prm, const std::function<Vector()>&) {
            const double w = prm.get<double>("relaxation", 1.0);
            const int n = prm.get<int>("ilulevel", 0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V>>(
                op.getmat(), n, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("ILUn", [](const O& op, const P& prm, const std::function<Vector()>&) {
            const int n = prm.get<int>("ilulevel", 0);
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V>>(
                op.getmat(), n, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("Jac", [](const O& op, const P& prm, const std::function<Vector()>&) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqJac<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("GS", [](const O& op, const P& prm, const std::function<Vector()>&) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqGS<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("SOR", [](const O& op, const P& prm, const std::function<Vector()>&) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqSOR<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("SSOR", [](const O& op, const P& prm, const std::function<Vector()>&) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqSSOR<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("amg", [](const O& op, const P& prm, const std::function<Vector()>&) {
            const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
            if (smoother == "ILU0" || smoother == "ParOverILU0") {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                using Smoother = SeqILU<M, V, V>;
#else
                using Smoother = SeqILU0<M, V, V>;
#endif
                return makeAmgPreconditioner<Smoother>(op, prm);
            } else if (smoother == "Jac") {
                using Smoother = SeqJac<M, V, V>;
                return makeAmgPreconditioner<Smoother>(op, prm);
            } else if (smoother == "SOR") {
                using Smoother = SeqSOR<M, V, V>;
                return makeAmgPreconditioner<Smoother>(op, prm);
            } else if (smoother == "SSOR") {
                using Smoother = SeqSSOR<M, V, V>;
                return makeAmgPreconditioner<Smoother>(op, prm);
            } else if (smoother == "ILUn") {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                using Smoother = SeqILU<M, V, V>;
#else
                using Smoother = SeqILUn<M, V, V>;
#endif
                return makeAmgPreconditioner<Smoother>(op, prm);
            } else {
                OPM_THROW(std::invalid_argument, "Properties: No smoother with name " << smoother <<".");
            }
        });
	doAddCreator("kamg", [](const O& op, const P& prm, const std::function<Vector()>&) {
            const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
            if (smoother == "ILU0" || smoother == "ParOverILU0") {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                using Smoother = SeqILU<M, V, V>;
#else
                using Smoother = SeqILU0<M, V, V>;
#endif
                return makeAmgPreconditioner<Smoother>(op, prm, true);
            } else if (smoother == "Jac") {
                using Smoother = SeqJac<M, V, V>;
                return makeAmgPreconditioner<Smoother>(op, prm, true);
            } else if (smoother == "SOR") {
                using Smoother = SeqSOR<M, V, V>;
                return makeAmgPreconditioner<Smoother>(op, prm, true);
	    // } else if (smoother == "GS") {
            //     using Smoother = SeqGS<M, V, V>;
            //     return makeAmgPreconditioner<Smoother>(op, prm, true);
            } else if (smoother == "SSOR") {
                using Smoother = SeqSSOR<M, V, V>;
                return makeAmgPreconditioner<Smoother>(op, prm, true);
            } else if (smoother == "ILUn") {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                using Smoother = SeqILU<M, V, V>;
#else
                using Smoother = SeqILUn<M, V, V>;
#endif
                return makeAmgPreconditioner<Smoother>(op, prm, true);
            } else {
                OPM_THROW(std::invalid_argument, "Properties: No smoother with name " << smoother <<".");
            }
        });
        doAddCreator("famg", [](const O& op, const P& prm, const std::function<Vector()>&) {
            auto crit = amgCriterion(prm);
            Dune::Amg::Parameters parms;
            parms.setNoPreSmoothSteps(1);
            parms.setNoPostSmoothSteps(1);
            return wrapPreconditioner<Dune::Amg::FastAMG<O, V>>(op, crit, parms);
        });
        doAddCreator("cpr", [](const O& op, const P& prm, const std::function<Vector()>& weightsCalculator) {
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, false>>(op, prm, weightsCalculator);
        });
        doAddCreator("cprt", [](const O& op, const P& prm, const std::function<Vector()>& weightsCalculator) {
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, true>>(op, prm, weightsCalculator);
        });
    }


    // The method that implements the singleton pattern,
    // using the Meyers singleton technique.
    static PreconditionerFactory& instance()
    {
        static PreconditionerFactory singleton;
        return singleton;
    }

    // Private constructor, to keep users from creating a PreconditionerFactory.
    PreconditionerFactory()
    {
        Comm* dummy = nullptr;
        addStandardPreconditioners(dummy);
    }

    // Actually creates the product object.
    PrecPtr doCreate(const Operator& op, const boost::property_tree::ptree& prm,
                     const std::function<Vector()> weightsCalculator)
    {
        const std::string& type = prm.get<std::string>("type", "ParOverILU0");
        auto it = creators_.find(type);
        if (it == creators_.end()) {
            std::ostringstream msg;
            msg << "Preconditioner type " << type << " is not registered in the factory. Available types are: ";
            for (const auto& prec : creators_) {
                msg << prec.first << ' ';
            }
            msg << std::endl;
            OPM_THROW(std::invalid_argument, msg.str());
        }
        return it->second(op, prm, weightsCalculator);
    }

    PrecPtr doCreate(const Operator& op, const boost::property_tree::ptree& prm,
                     const std::function<Vector()> weightsCalculator, const Comm& comm)
    {
        const std::string& type = prm.get<std::string>("type", "ParOverILU0");
        auto it = parallel_creators_.find(type);
        if (it == parallel_creators_.end()) {
            std::ostringstream msg;
            msg << "Parallel preconditioner type " << type
                << " is not registered in the factory. Available types are: ";
            for (const auto& prec : parallel_creators_) {
                msg << prec.first << ' ';
            }
            msg << std::endl;
            OPM_THROW(std::invalid_argument, msg.str());
        }
        return it->second(op, prm, weightsCalculator, comm);
    }

    // Actually adds the creator.
    void doAddCreator(const std::string& type, Creator c)
    {
        creators_[type] = c;
    }

    // Actually adds the creator.
    void doAddCreator(const std::string& type, ParCreator c)
    {
        parallel_creators_[type] = c;
    }

    // This map contains the whole factory, i.e. all the Creators.
    std::map<std::string, Creator> creators_;
    std::map<std::string, ParCreator> parallel_creators_;
};

} // namespace Dune

#endif // OPM_PRECONDITIONERFACTORY_HEADER
