/*--------------------------------------------------------------------------*/
/*-------------------- File MasterProblemBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class MasterProblemBlock, which derives from Block to
 * implement the Master Problem of a generic (Generalized) Bundle algorithm
 * within the SMS++ framework.
 *
 * MasterProblemBlock represents and solves the Master Problem (MP) of a
 * Bundle method on a sum-function
 *
 *     f( x ) = b * x + \sum_{ k \in K } f^k( x )
 *
 * where each component f^k is accessed either via a black-box oracle
 * ("hard" component) or has a known compact convex description ("easy"
 * component, cf. [Frangioni, Gorgone, MP 2014]).
 *
 * Two complementary formulations of the MP are supported, both directly
 * encoded as a Block structure (and therefore solvable by any Solver
 * registered to MasterProblemBlock, typically a [MILP]Solver):
 *
 * - the **primal** form (cf. eq. (24)/(28) of the reference paper) reads
 *
 *      min   b * d + \sum_k v^k + (1/2t) || d ||^2_2
 *      s.t.  v^k >= g^k_i * d + alpha^k_i        for each i in B^k       (P)
 *            v^k >= LB^k
 *
 *   where d is the step from the current stability center x_bar, v^k is
 *   the epigraph variable for component k, and (g^k_i, alpha^k_i) are the
 *   linearizations stored in the bundle B^k;
 *
 * - the **dual** form (cf. eq. (25)/(26)/(31) of the reference paper) reads
 *
 *      max   \sum_k \sum_i theta^k_i (c^k - x_bar A^k) u^k_i
 *            + x_bar z - (t/2) || z ||^2_2 - \sum_k \sum_i theta^k_i alpha^k_i
 *      s.t.  \sum_i theta^k_i + gamma^k = lambda          for each k     (D)
 *            z = b - \sum_k A^k ( \sum_i theta^k_i u^k_i )
 *            theta^k_i >= 0 ,  gamma^k >= 0 ,  lambda >= 0
 *
 *   plus the additional native u^k variables and constraints of any "easy"
 *   component k (its compact polyhedral description being directly inserted
 *   in the MP, rather than inner-approximated by extreme points).
 *
 * In both formulations:
 *
 * - one PolyhedralFunctionBlock is registered as a sub-Block of
 *   MasterProblemBlock for each "hard" component, holding its bundle B^k
 *   (either in primal or in dual representation, consistently with the
 *   chosen MP form);
 *
 * - each "easy" component is registered as a sub-Block containing the
 *   compact convex description of f^k (typically a linear/convex program
 *   inherited from a LagBFunction);
 *
 * - the MasterProblemBlock itself owns the *coupling* part of the MP, i.e.,
 *   the static Variable and Constraint that glue the per-component sub-MPs
 *   together (lambda, r, omega, gamma^k, the coupling constraint on z,
 *   the level/global-LB rows, ...), and an FRealObjective with the
 *   stabilizing quadratic term plus the linear part.
 *
 * Three stabilization schemes are supported, selected via the
 * #stabilization_type enum:
 *
 * - **Proximal**: quadratic term (1/2t) || d ||^2_2 in (P), equivalently
 *   - (t/2) || z ||^2_2 in (D);
 *
 * - **Level**: an additional level constraint v <= f_lev together with a
 *   non-negative multiplier omega in the dual;
 *
 * - **Doubly-Stabilized**: combines both, see [Astorino, Frangioni,
 *   Fuduli, Gorgone, MP 2017].
 *
 * MasterProblemBlock is meant to be driven by a (Generalized) BundleSolver
 * which is responsible for keeping the bundles B^k updated as the algorithm
 * proceeds; the BundleSolver does *not* directly call any MILP backend, it
 * only manipulates the MP at the Block/Modification level and triggers
 * compute() on the [MILP]Solver attached to the MasterProblemBlock.
 *
 * For the underlying theory and notation, the reader is referred to:
 *
 *  A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
 *  with ``Easy'' Components: Applications to Multicommodity Network Design"
 *  Mathematical Programming 145(1), 133 - 161, 2014
 *
 * available at
 *
 *  \link http://www.di.unipi.it/~frangio/abstracts.html#MP11c \endlink
 *
 *  A. Frangioni "Generalized Bundle Methods"
 *  SIAM Journal on Optimization 13(1), 117 - 156, 2002
 *
 * available at
 *
 *  \link http://www.di.unipi.it/~frangio/abstracts.html#SIOPT02 \endlink
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Calandrini, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MasterProblemBlock
 #define __MasterProblemBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "BlockSolverConfig.h"
#include "ColVariable.h"
#include "FRowConstraint.h"

#include <string>
#include <vector>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MasterProblemBlock_CLASSES Classes in MasterProblemBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS MasterProblemBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A Block representing the Master Problem of a (Generalized) Bundle Method
/** MasterProblemBlock implements, as an SMS++ Block, the Master Problem (MP)
 * of a Generalized Bundle Method (cf. \link MasterProblemBlock.h \endlink for
 * the mathematical description of the supported primal and dual forms, and
 * for the underlying references).
 *
 * The Block exposes the coupling part of the MP (the static Variable and
 * Constraint linking together the per-component bundles, plus the stabilizing
 * quadratic Objective), while each component f^k is represented by a
 * dedicated sub-Block:
 *
 * - "hard" components -> a PolyhedralFunctionBlock holding the bundle B^k
 *   (in the same primal/dual representation as the master);
 *
 * - "easy" components -> a sub-Block containing the compact convex
 *   description of f^k (e.g. inherited from a LagBFunction).
 *
 * A regular Solver (typically a [MILP]Solver from the SMS++ MILPSolver
 * module) is attached to MasterProblemBlock through register_Solver(), and
 * is then asked to solve the MP at every Bundle iteration. */

class MasterProblemBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

 /// stabilization scheme used by the Master Problem
 /** The enum #stabilization_type lists the stabilization mechanisms that
  * MasterProblemBlock can insert into the MP. They differ in the term that
  * is added to the cutting-plane model to dampen the oscillations of the
  * unstabilized Kelley method:
  *
  * - #kProximal: proximal quadratic term (1/2t) || d ||^2_2;
  *
  * - #kLevel: level constraint v <= f_lev with dual multiplier omega >= 0;
  *
  * - #kDoublyStabilized: both terms combined. */

 enum stabilization_type {
  kProximal         = 0 ,  ///< proximal stabilization
  kLevel            = 1 ,  ///< level stabilization
  kDoublyStabilized = 2    ///< doubly-stabilized bundle method
  };

/*----------------------------- CONSTANTS ----------------------------------*/

/** @} ---------------------------------------------------------------------*/
/*----------- CONSTRUCTING AND DESTRUCTING MasterProblemBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MasterProblemBlock
 *  @{ */

 /// constructor: initializes every algorithmic field to a safe default
 /** All "size" fields are set to 0; the default stabilization is
  * #kDoublyStabilized and the default form is the dual one (IsPrimal ==
  * false). The actual size of the MP is then established by SetDim() and
  * the formulation is selected by CreateEmptyMP(). */

 explicit MasterProblemBlock( Block * father = nullptr )
  : Block( father ) , IsPrimal( false ) , StblType( kDoublyStabilized ) ,
    MaxBSize( 0 ) , MaxSGLen( 0 ) , NumVars( 0 ) ,
    NoTotCmps( 0 ) , NoEasyCmps( 0 ) , NoHardCmps( 0 ) , DoEasy( 0 ) { }

/*--------------------------------------------------------------------------*/
 /// destructor: releases all the resources owned by MasterProblemBlock
 /** The destructor releases the dynamic Variable/Constraint that
  * MasterProblemBlock owns via the abstract representation, and detaches
  * any registered Solver. The actual cleanup is delegated to clear(). */

 ~MasterProblemBlock() override { clear(); }

/*--------------------------------------------------------------------------*/
 /// release the abstract representation and any per-component state
 /** clear() resets MasterProblemBlock to an "empty" state: all per-component
  * sub-Blocks pointers are forgotten (the sub-Blocks themselves are
  * dismissed via the Block destructor of the base class), every size field
  * goes back to 0 and the MP type information is reset. A subsequent
  * SetDim() + CreateEmptyMP() is then needed to rebuild the MP. */

 void clear();

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// provide MasterProblemBlock with the basic dimensions of the MP
 /** Provides MasterProblemBlock with the four numbers fully describing the
  * "size" of the MP it has to solve:
  *
  * - \p MxBSz is the maximum number of items (subgradients + feasibility
  *   cuts) that can be kept *per component*; it is stored in #MaxBSize;
  *
  * - \p NVars is the number of optimization variables in the Primal MP
  *   (i.e., the dimension of the step d); it is stored in #NumVars;
  *
  * - \p NrFi is the total number of components of the sum-function; it is
  *   stored in #NoTotCmps;
  *
  * - \p NrFiEasy is the number of "easy" components, which receive a
  *   compact in-place representation rather than being inner-approximated
  *   by a bundle of linearizations; it is stored in #NoEasyCmps and the
  *   number of "hard" components is derived as NoTotCmps - NoEasyCmps.
  *
  * Calling SetDim() resets any previous state: existing sub-Blocks are
  * dropped and the abstract representation is destroyed; CreateEmptyMP()
  * must be called afterwards to actually populate the new MP. */

 void SetDim( int MxBSz , int NVars , int NrFi , int NrFiEasy );

/*--------------------------------------------------------------------------*/
 /// register the inner Solver of MasterProblemBlock
 /** Attaches a Solver to MasterProblemBlock, used to (re-)solve the MP at
  * every iteration of the surrounding Bundle algorithm. \p solv_cfg_filename
  * is the path of a file containing a BlockSolverConfig (in either text or
  * netCDF format): the file is deserialized via Configuration::deserialize()
  * and apply()-ed to MasterProblemBlock.
  *
  * Two kinds of error are reported via std::invalid_argument:
  *
  * - \p solv_cfg_filename is empty: there is *no* default backend hard-wired
  *   into MasterProblemBlock; the choice of the actual [MILP]Solver must
  *   always be expressed by the caller through a BlockSolverConfig (and
  *   resolved by the SMS++ Solver factory at apply() time), exactly like
  *   for any other Block in the framework;
  *
  * - the file does not contain a BlockSolverConfig.
  *
  * After a successful call the registered Solver is ready to be
  * compute()-ed. */

 void register_Solver( std::string && solv_cfg_filename );

/*--------------------------------------------------------------------------*/
 /// initialize an *empty* Master Problem with the given stabilization
 /** Populates the abstract representation of MasterProblemBlock with the
  * coupling part of the Master Problem (the static Variable / Constraint
  * and the Objective), in either the primal or the dual formulation
  * depending on the chosen stabilization scheme and on the presence of
  * "easy" components.
  *
  * At call time the per-component bundles are *empty*: the dynamic cuts
  * of each "hard" component are added on the fly by the surrounding Bundle
  * algorithm via the Modification interface of the corresponding
  * PolyhedralFunctionBlock sub-Block.
  *
  * As a general rule the primal form is preferred when no "easy" component
  * is present, while the dual form is the only viable choice when
  * \p NoEasy > 0 *and* \p DoEasyCmp != 0 (because easy components are
  * naturally expressed in the dual MP).
  *
  * @param Stbl       the stabilization scheme, see #stabilization_type;
  * @param NoCmps     total number of components (NoTotCmps);
  * @param DoEasyCmp  bit-wise flag controlling the easy-component handling
  *                   (it is forwarded to each LagBFunction sub-Block);
  * @param NoEasy     number of "easy" components;
  * @param IsEasy     boolean vector of length \p NoCmps such that
  *                   <tt>IsEasy[k] == true</tt> iff component \p k is
  *                   "easy"; a copy is kept inside the class. */

 void CreateEmptyMP( stabilization_type Stbl , int NoCmps , int DoEasyCmp ,
                     int NoEasy , std::vector< bool > IsEasy );

/*--------------------------------------------------------------------------*/
 /// initialize the primal version of the Master Problem
 /** Builds the static Variable, the static Constraint and the Objective of
  * the primal MP (cf. eq. (P) in the file documentation). This method is
  * called internally by CreateEmptyMP() when the dual form is not strictly
  * required, and is exposed publicly so that derived classes can override
  * the construction (e.g. for problem-specific stabilizing terms).
  *
  * The current implementation is a work-in-progress placeholder; calling
  * it on a non-empty Block has unspecified effects. */

 void CreatePrimalMP( stabilization_type Stbl );

/*--------------------------------------------------------------------------*/
 /// initialize the dual version of the Master Problem
 /** Builds the static Variable, the static Constraint and the Objective of
  * the dual MP (cf. eq. (D) in the file documentation). This method is
  * called internally by CreateEmptyMP() whenever "easy" components are
  * present, since their compact description is naturally expressed in the
  * dual form; it is exposed publicly so that derived classes can override
  * the construction.
  *
  * The current implementation is a work-in-progress placeholder; calling
  * it on a non-empty Block has unspecified effects. */

 void CreateDualMP( stabilization_type Stbl );

/*--------------------------------------------------------------------------*/
 /// hand the abstract representation of the MP to the registered Solver
 /** load_problem() iterates over the Solver registered to
  * MasterProblemBlock and instructs each of them to (re-)load the abstract
  * representation of the MP, including any extra information coming from
  * the per-component sub-Blocks (e.g. the compact constraints of every
  * "easy" component).
  *
  * It must be called after SetDim() + CreateEmptyMP() + register_Solver()
  * and before the first compute(); subsequent compute() calls do *not*
  * require a fresh load_problem() unless the MP type/size changes. */

 void load_problem( void );

/*--------------------------------------------------------------------------*/
 /// load a MasterProblemBlock out of an istream
 /** Required by the abstract interface of Block, but currently not
  * implemented (a MasterProblemBlock is always built programmatically by
  * its driving BundleSolver): always throws an std::logic_error. */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error(
       "MasterProblemBlock::load not implemented yet" ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // - - - - - - - - -  algorithmic / structural parameters - - - - - - - - - -

 bool IsPrimal;     ///< whether the MP is in its primal or dual form

 stabilization_type StblType;  ///< type of stabilization, see #stabilization_type

 int MaxBSize;      ///< maximum number of items kept per bundle B^k

 int MaxSGLen;      ///< maximum length of a subgradient (currently == NumVars)

 int NumVars;       ///< number of optimization variables in the Primal MP

 int NoTotCmps;     ///< total number of components in the sum-function

 int NoEasyCmps;    ///< number of "easy" components

 int NoHardCmps;    ///< number of "hard" components ( == NoTotCmps - NoEasyCmps)

 int DoEasy;        ///< bit-wise flag controlling the easy-component handling

 std::vector< bool > IsEasyCmp;  ///< IsEasyCmp[k] == true iff component k is easy

 // - - - - - - - - -  pointers to the per-component sub-Blocks - - - - - - - -

 std::vector< Block * > EasyCmps;  ///< sub-Blocks of the "easy" components

 std::vector< Block * > HardCmps;  ///< sub-Blocks of the "hard" components

 // - - - - - - - - - - - - - -  static MP entities  - - - - - - - - - - - - -

 ColVariable Var_lambda;           ///< dual multiplier of the global v >= ...
                                   ///< row, "lambda" in (D)

 ColVariable Var_r;                ///< dual multiplier of the global LB row

 ColVariable Var_omega;            ///< dual multiplier of the level / X row

 std::vector< ColVariable > Var_gamma;
                                   ///< per-hard-component LB^k multipliers

 FRowConstraint NormalizationCns;  ///< global normalization row lambda+r-omega=1

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class MasterProblemBlock )

/** @} end( group( MasterProblemBlock_CLASSES ) ) ----------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* __MasterProblemBlock */

/*--------------------------------------------------------------------------*/
/*--------------------- End File MasterProblemBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
