/*--------------------------------------------------------------------------*/
/*-------------------- File MasterProblemBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class MasterProblemBlock, which derives from Block to
 * implement the Master Problem of a generic Bundle algorithm
 * within the SMS++ framework.
 *
 * MasterProblemBlock represents and solves the Master Problem (MP) of a Bundle
 * method on a sum-function
 *
 *     f( x ) = b * x + \sum_{ k \in K } f^k( x )
 *
 * where each component f^k is accessed either via a black-box oracle ("hard"
 * component) or has a known compact convex description ("easy" component, cf.
 * [Frangioni, Gorgone, MP 2014]).
 *
 * Two complementary formulations of the MP are supported, both directly encoded
 * as a Block structure (and therefore solvable by any Solver registered to
 * MasterProblemBlock, typically a [MILP]Solver):
 *
 * - the **primal** form reads
 *
 *      min   b * d + \sum_k v^k + (1/2t) || d ||^2_2
 *      s.t.  v^k >= g^k_i * d + alpha^k_i        for each i in B^k       (P)
 *            v^k >= LB^k
 *
 *   where d is the step from the current stability center x_bar, v^k is
 *   the epigraph variable for component k, and (g^k_i, alpha^k_i) are the
 *   linearizations stored in the bundle B^k;
 *
 * - the **dual** form reads
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
 * MasterProblemBlock is meant to be driven by a bundle driver
 * which is responsible for keeping the bundles B^k updated as the algorithm
 * proceeds; the driver does *not* directly call any MILP backend, it only
 * manipulates the MP at the Block/Modification level and triggers compute() on
 * the [MILP]Solver attached to the MasterProblemBlock.
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
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Calandrini \n
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
#include "C05Function.h"
#include "ColVariable.h"
#include "FRowConstraint.h"
#include "LinearFunction.h"
#include "OneVarConstraint.h"

#include <iosfwd>
#include <list>
#include <string>
#include <unordered_set>
#include <vector>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------- FORWARD-DECLARED FRIENDS ---------------------------*/
/*--------------------------------------------------------------------------*/

class BendersBFunction;
class LagBFunction;

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
/// A Block representing the Master Problem of a Bundle Method
/** MasterProblemBlock implements, as an SMS++ Block, the Master Problem (MP)
 * of a Generalized Bundle Method (cf. \link MasterProblemBlock.h \endlink for
 * the mathematical description of the supported primal and dual forms, and for
 * the underlying references).
 *
 * The Block exposes the coupling part of the MP (the static Variable and
 * Constraint linking together the per-component bundles, plus the stabilizing
 * quadratic Objective), while each component f^k is represented by a dedicated
 * sub-Block:
 *
 * - "hard" components -> a PolyhedralFunctionBlock holding the bundle B^k
 *   (in the same primal/dual representation as the master);
 *
 * - "easy" components -> a sub-Block containing the compact convex
 *   description of f^k (e.g. inherited from a LagBFunction).
 *
 * A regular Solver (typically a [MILP]Solver from the SMS++ MILPSolver module)
 * is attached to MasterProblemBlock through register_Solver(), and is then
 * asked to solve the MP at every Bundle iteration. */

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
  * MasterProblemBlock can insert into the MP. They differ in the term that is
  * added to the cutting-plane model to dampen the oscillations of the
  * unstabilized Kelley method:
  *
  * - #kNone: no stabilization (pure Kelley cutting-plane);
  *
  * - #kProximal: proximal quadratic term (1/2t) || d ||^2_2;
  *
  * - #kLevel: level constraint v <= f_lev with dual multiplier omega >= 0;
  *
  * - #kTrustRegion: hard trust region || d ||_inf <= t;
  *
  * - #kDoublyStabilized: kProximal + kLevel combined;
  *
  * - #kUpperLower: two-sided upper/lower bundle [u_bar, l_bar] with a
  *   proximal term anchored at the upper bundle u_bar; #kProximal,
  *   #kLevel and #kDoublyStabilized are recovered as special cases by
  *   degenerating one of the two bundles.
  *
  * Only #kProximal, #kLevel and #kDoublyStabilized are currently fully
  * wired in CreatePrimalMP / CreateDualMP; #kNone, #kTrustRegion and
  * #kUpperLower are reserved enumerators and CreatePrimalMP /
  * CreateDualMP throw `std::logic_error` if invoked with one of them. */

 enum stabilization_type {
  kProximal         = 0 ,  ///< proximal stabilization
  kLevel            = 1 ,  ///< level stabilization
  kDoublyStabilized = 2 ,  ///< doubly-stabilized bundle method
  kNone             = 3 ,  ///< no stabilization (pure cutting plane)
  kTrustRegion      = 4 ,  ///< trust region || d ||_inf <= t
  kUpperLower       = 5    ///< upper / lower bundle pair
  };

/*----------------------------- CONSTANTS ----------------------------------*/

/** @} ---------------------------------------------------------------------*/
/*----------- CONSTRUCTING AND DESTRUCTING MasterProblemBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MasterProblemBlock
 *  @{ */

 /// constructor: initializes every algorithmic field to a safe default
 /** All "size" fields are set to 0; the default stabilization is
  * #kDoublyStabilized and the default form is the dual one (IsPrimal == false).
  * The actual size of the MP is then established by SetDim() and the
  * formulation is selected by CreateEmptyMP(). */

 explicit MasterProblemBlock( Block * father = nullptr )
  : Block( father ) , IsPrimal( false ) , IsConvex( true ) ,
    StblType( kDoublyStabilized ) ,
    MaxBSize( 0 ) , MaxSGLen( 0 ) , NumVars( 0 ) ,
    NoTotCmps( 0 ) , NoEasyCmps( 0 ) , NoHardCmps( 0 ) , DoEasy( 0 ) ,
    t_stab( 1.0 ) , f_lev( 0.0 ) ,
    z_obj_idx( -1 ) , r_obj_idx( -1 ) , omega_obj_idx( -1 ) { }

/*--------------------------------------------------------------------------*/
 /// destructor: releases all the resources owned by MasterProblemBlock
 /** The destructor releases the dynamic Variable/Constraint that
  * MasterProblemBlock owns via the abstract representation, and detaches any
  * registered Solver. The actual cleanup is delegated to clear(). */

 ~MasterProblemBlock() override { clear(); }

/*--------------------------------------------------------------------------*/
 /// release the abstract representation and any per-component state
 /** clear() resets MasterProblemBlock to an "empty" state: all per-component
  * sub-Blocks pointers are forgotten (the sub-Blocks themselves are dismissed
  * via the Block destructor of the base class), every size field goes back to 0
  * and the MP type information is reset. A subsequent SetDim() +
  * CreateEmptyMP() is then needed to rebuild the MP. */

 void clear();

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// one-shot configuration of MasterProblemBlock
 /** Single-call configuration that bundles dimensioning, primal/dual
  * variant selection and easy-component registration into one entry
  * point. This is the API the surrounding driver /
  * BendersDecompositionSolver uses to populate MasterProblemBlock with
  * everything it needs.
  *
  * MasterProblemBlock never sees the "hard" components directly: only
  * their number (\p num_hard_cmps) is fed in, since their bundles B^k
  * are built up internally as empty PolyhedralFunctionBlock(s) under
  * *this* and populated by the surrounding bundle driver
  * through the add_cut() / remove_cut() interface. The "easy" components
  * are passed in via \p easy_components, type-dispatched (LagBFunction,
  * BendersBFunction, ...) and absorbed into the master in the form
  * dictated by the requested primal / dual MP. The dualization details
  * stay inside MasterProblemBlock, so callers only need to know which
  * components are easy.
  *
  * The Block containing the model variables and constraints is
  * *stolen* (i.e. its ownership is transferred via
  * Block::transfer_ownership_to(this)) into the master tree, so that
  * the inner :MILPSolver naturally sees those variables/constraints
  * when it loads the MP. The optional \p ignored_blocks list is then
  * forwarded to the inner Solver via Solver::set_excluded_blocks(),
  * telling it which sub-Block subtrees of the stolen Block (and of the
  * easy components) have to be skipped — for instance the original
  * sub-Block that carried the C05Function whose linearizations are
  * now represented by a PolyhedralFunctionBlock attached here.
  *
  * @param primal            if true the MP is built in its primal form,
  *                          else in its dual form. The primal form
  *                          currently supports only BendersBFunction as
  *                          easy components; the dual form supports
  *                          only LagBFunction. Any other C05Function
  *                          type, or a mismatch with \p primal, makes
  *                          configure() throw.
  * @param max_bundle_size   maximum number of items (subgradients +
  *                          feasibility cuts) kept per hard component;
  *                          stored in #MaxBSize.
  * @param num_vars          number of optimization variables in the
  *                          underlying sum-function space (the dimension
  *                          of the step d, of the dual z, ...); stored
  *                          in #NumVars.
  * @param num_hard_cmps     number of "hard" components, i.e. components
  *                          whose bundle is built up incrementally by the
  *                          driving bundle driver through
  *                          add_cut() / remove_cut(); stored in
  *                          #NoHardCmps.
  * @param easy_components   pointers to the easy C05Function components.
  *                          Each pointer must be non-null and must be of
  *                          a supported type (LagBFunction in the dual
  *                          MP, BendersBFunction in the primal MP);
  *                          its size is stored in #NoEasyCmps.
  * @param original_block    the Block holding the model variables and
  *                          constraints to be plugged into the master.
  *                          If not nullptr, it is stolen via
  *                          Block::transfer_ownership_to(this).
  * @param ignored_blocks    set of sub-Blocks the inner Solver must
  *                          ignore (forwarded to BlockSolverConfig::apply
  *                          when register_Solver() is called, which in
  *                          turn calls Solver::set_excluded_blocks() on
  *                          the freshly-created Solver). Default = empty.
  * @param reg               stabilization scheme, see
  *                          #stabilization_type.
  * @param convex            true if the underlying C05Function is convex
  *                          (the master is minimised), false if it is
  *                          concave (the master is maximised). */

 void configure( bool primal ,
                 int max_bundle_size ,
                 int num_vars ,
                 int num_hard_cmps ,
                 const std::vector< C05Function * > & easy_components ,
                 Block * original_block ,
                 std::unordered_set< Block * > ignored_blocks = {} ,
                 stabilization_type reg = kDoublyStabilized ,
                 bool convex = true );

/*--------------------------------------------------------------------------*/
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
  * Calling SetDim() resets any previous state: existing sub-Blocks are dropped
  * and the abstract representation is destroyed; CreateEmptyMP() must be called
  * afterwards to actually populate the new MP. */

 void SetDim( int MxBSz , int NVars , int NrFi , int NrFiEasy );

/*--------------------------------------------------------------------------*/
 /// register the inner Solver of MasterProblemBlock
 /** Attaches a Solver to MasterProblemBlock, used to (re-)solve the MP at
  * every iteration of the surrounding Bundle algorithm. \p solv_cfg_filename is
  * the path of a file containing a BlockSolverConfig (in either text or netCDF
  * format): the file is deserialized via Configuration::deserialize() and
  * apply()-ed to MasterProblemBlock.
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
  * After a successful call the registered Solver is ready to be compute()-ed.
  */

 void register_Solver( std::string && solv_cfg_filename );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// register the inner Solver + install an exclusion list on it
 /** Same as register_Solver(filename) but, after the inner Solver is
  * attached, also forwards \p ignored_blocks to it via
  * Solver::set_excluded_blocks() (through BlockSolverConfig::apply()'s
  * second argument). Use this overload when configure() has stolen
  * sub-Block subtrees whose Variable / Constraint should be invisible to
  * the inner Solver (e.g. the original sub-Block whose C05Function is now
  * represented by an internally-managed PolyhedralFunctionBlock). The
  * eager descendant expansion is performed Solver-side; the caller only
  * has to name the roots of the sub-trees to skip. */

 void register_Solver( std::string && solv_cfg_filename ,
                       std::unordered_set< Block * > && ignored_blocks );

/*--------------------------------------------------------------------------*/
 /// initialize an *empty* Master Problem with the given stabilization
 /** Populates the abstract representation of MasterProblemBlock with the
  * coupling part of the Master Problem (the static Variable / Constraint and
  * the Objective), in either the primal or the dual formulation depending on
  * the chosen stabilization scheme and on the presence of "easy" components.
  *
  * At call time the per-component bundles are *empty*: the dynamic cuts of each
  * "hard" component are added on the fly by the surrounding Bundle algorithm
  * via the Modification interface of the corresponding PolyhedralFunctionBlock
  * sub-Block.
  *
  * As a general rule the primal form is preferred when no "easy" component is
  * present, while the dual form is the only viable choice when \p NoEasy > 0
  * *and* \p DoEasyCmp != 0 (because easy components are naturally expressed in
  * the dual MP).
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
  * the primal MP (cf. eq. (P) in the file documentation). This method is called
  * internally by CreateEmptyMP() when the dual form is not strictly required,
  * and is exposed publicly so that derived classes can override the
  * construction (e.g. for problem-specific stabilizing terms).
  *
  * The current implementation is a work-in-progress placeholder; calling it on
  * a non-empty Block has unspecified effects. */

 void CreatePrimalMP( stabilization_type Stbl );

/*--------------------------------------------------------------------------*/
 /// initialize the dual version of the Master Problem
 /** Builds the static Variable, the static Constraint and the Objective of
  * the dual MP (cf. eq. (D) in the file documentation). This method is called
  * internally by CreateEmptyMP() whenever "easy" components are present, since
  * their compact description is naturally expressed in the dual form; it is
  * exposed publicly so that derived classes can override the construction.
  *
  * The current implementation is a work-in-progress placeholder; calling it on
  * a non-empty Block has unspecified effects. */

 void CreateDualMP( stabilization_type Stbl );

/*--------------------------------------------------------------------------*/
 /// absorb the row-mapping of a BendersBFunction into the primal MP
 /** Embeds an easy BendersBFunction into the primal Master Problem. The
  * BendersBFunction is associated with the optimization problem
  *
  *      phi( x ) = min { c( y ) : ( g - F x ) <= E( y ) <= ( h - F x ) ,
  *                                y in Y }
  *
  * and its inner Block carries the constraints in fixed-rhs form
  *
  *      bar{w}_i <= E_i( y ) <= bar{z}_i
  *
  * (the current values returned by BendersBFunction::map_f_value() at
  * the previous stability centre). To absorb the BendersBFunction into
  * the primal master, every active mapping row i (read from
  * BendersBFunction::get_A() / get_b() / get_constraints() /
  * get_sides()) is relaxed on the inner Block (LHS = -INF, RHS = +INF
  * on the corresponding RowConstraint, depending on its
  * #ConstraintSide) and re-installed on *this* as an FRowConstraint
  *
  *      E_i( y ) - ( A x )_i  [<=, =, >=]  b_i
  *
  * built by appending the linear coupling -A_i . x to the function of
  * the original RowConstraint. The cast on the function of the original
  * RowConstraint is type-specific (LinearFunction, DQuadFunction,
  * QuadFunction); unsupported types currently trigger an exception.
  *
  * @param bbf  pointer to the BendersBFunction to absorb; must be non-null. */

 void absorb_BBF_into_primal_MP( BendersBFunction * bbf );

/*--------------------------------------------------------------------------*/
 /// absorb a LagBFunction's inner constraints into the dual MP
 /** Embeds an easy LagBFunction into the dual Master Problem. The
  * polyhedral feasible set U^k of the LagBFunction is described by the
  * RowConstraints living in its inner Block in the fixed-rhs form
  *
  *      lhs_i <= ( E^k_i u^k ) <= rhs_i
  *
  * (with constants lhs_i / rhs_i). The stationarity condition of the
  * dual MP at the pi^k variables of the easy linear program (cf.
  * ) reads instead
  *
  *      E^k_i u^k + lambda * e^k_i = 0    for every row i,
  *
  * where lambda is the global multiplier owned by *this* (Var_lambda)
  * and e^k_i is taken from the side of the original RowConstraint:
  *
  *  - if S_i = eRHS  -> e^k_i := rhs_i ;
  *  - if S_i = eLHS  -> e^k_i := lhs_i ;
  *  - if S_i = eBoth -> e^k_i := rhs_i (= lhs_i, equality constraint).
  *
  * For every RowConstraint of the inner Block this method snapshots its
  * (lhs_i, rhs_i, function), relaxes it (LHS = -INF, RHS = +INF) and
  * installs on *this* a fresh FRowConstraint with the cloned
  * LinearFunction augmented by the +e^k_i coefficient on Var_lambda and
  * with LHS = RHS = 0. Quadratic / non-linear function types are not
  * supported and trigger an exception, since the underlying primal LP of
  * the easy component is linear by assumption.
  *
  * @param lbf  pointer to the LagBFunction to absorb; must be non-null
  *             and its inner Block must have been already transferred to
  *             *this* by configure(). */

 void absorb_LBF_into_dual_MP( LagBFunction * lbf );

/*--------------------------------------------------------------------------*/
 /// hand the abstract representation of the MP to the registered Solver
 /** load_problem() iterates over the Solver registered to
  * MasterProblemBlock and instructs each of them to (re-)load the abstract
  * representation of the MP, including any extra information coming from the
  * per-component sub-Blocks (e.g. the compact constraints of every "easy"
  * component).
  *
  * It must be called after SetDim() + CreateEmptyMP() + register_Solver() and
  * before the first compute(); subsequent compute() calls do *not* require a
  * fresh load_problem() unless the MP type/size changes. */

 void load_problem( void );

/*--------------------------------------------------------------------------*/
 /// returns the k-th hard-component sub-Block
 /** Returns the PolyhedralFunctionBlock used to model the k-th "hard"
  * component of the sum-function, for k in [0, NoHardCmps). The pointer stays
  * valid as long as the MP is not torn down by clear() or SetDim(); it is meant
  * to be used by the driver to feed cuts into the underlying
  * PolyhedralFunction via the Modification interface. */

 [[nodiscard]] Block * get_hard_component( int k ) const;

/*--------------------------------------------------------------------------*/
 /// returns the k-th easy-component sub-Block
 /** Returns the Block registered as the k-th "easy" component, for k in
  * [0, EasyCmps.size()), or nullptr if it has not been registered yet. */

 [[nodiscard]] Block * get_easy_component( int k ) const;

/*--------------------------------------------------------------------------*/
 /// returns ||z*||^2, the squared 2-norm of the aggregated dual subgradient
 /** In the dual MP the auxiliary variables z encode the aggregated
  * subgradient z* = b - sum_k A^k ( sum_i theta^k_i u^k_i ); this method
  * returns the current squared 2-norm sum_j z_j^2 read off the ColVariable
  * values, which is meaningful only after the master Solver has solved the MP
  * at least once.
  *
  * In the primal MP the same quantity has to be reconstructed from the per-PFB
  * optimal multipliers and is not available here; the method returns 0 in that
  * case. */

 [[nodiscard]] double get_dual_norm_squared( void ) const;

/*--------------------------------------------------------------------------*/
 /// add a new linearization (g, alpha_raw) at slot \p slot of bundle B^k
 /** Appends one new linearization to the bundle B^k of the k-th "hard"
  * component, occupying the persistent NDOFi-style "slot" \p slot of that
  * component. \p slot must lie in [0, MaxBSize) and must currently be empty (an
  * std::logic_error is thrown otherwise). \p g is the new subgradient (size
  * NumVars; moved into the PolyhedralFunction); \p alpha is the *raw* cut
  * constant (the b such that  F_k( x ) >= b + g . x  in convex-min sign,
  * un-translated). The caller is *not* expected to pre-promote alpha to
  * the linearization-error-at-x_bar form: the master keeps the raw value
  * in its bundle and rebuilds the linearization-error form on demand from
  * the cached F_k( x_bar ) installed via set_reference().
  *
  * The (k, slot) pair lets the surrounding driver keep its
  * ItemVcblr / InvItemVcblr name-management *unchanged* across subsequent
  * remove_cut() calls: MasterPB internally tracks the slot -> local-row mapping
  * into the PolyhedralFunction of HardCmps[k] (rows are renumbered on
  * delete_row(), the slot is not).
  *
  * The call issues the appropriate PolyhedralFunctionModAddd through the usual
  * Modification interface, so the [MILP]Solver attached to MasterProblemBlock
  * picks it up on the next compute().
  *
  * Before committing the new row, the incoming (g, alpha) is compared
  * coefficient-by-coefficient against every linearization already in the
  * bundle of HardCmps[k]: if a match is found (within an absolute
  * tolerance scaled by the magnitudes involved) no insertion takes place
  * and the return value is the (k, slot) of the duplicate. Otherwise the
  * cut is installed as described above and the return value is the
  * sentinel #kCutInserted (== -1), signalling that the master problem has
  * effectively changed. The surrounding driver consumes this return code
  * to decide whether the new cut warrants treating the master as
  * "changed" for the purposes of the SS/NS bookkeeping. */

 static constexpr int kCutInserted = -1;

 int add_cut( int k , int slot , std::vector< double > && g , double alpha ,
              bool is_vert = false );

/*--------------------------------------------------------------------------*/
 /// add a new linearization (g, alpha) to bundle B^k, auto-allocating a slot
 /** Like add_cut(k, slot, g, alpha, is_vert) but the slot is chosen by
  * MasterPB itself; returns the chosen slot (or -1 if no slot is free, i.e. all
  * MaxBSize positions of B^k are occupied). \p is_vert tags the cut as a
  * vertical (feasibility) one so that the PolyhedralFunction back-end records
  * it as such (no v_k epigraph variable on the LHS of the row). */

 int add_cut( int k , std::vector< double > && g , double alpha ,
              bool is_vert = false );

/*--------------------------------------------------------------------------*/
 /// returns the number of cuts currently in the bundle of the k-th hard cmp
 /** Returns the size of the bundle B^k, i.e. the number of subgradient /
  * feasibility cuts currently stored in the underlying PolyhedralFunction of
  * HardCmps[k]. \p k must lie in [0, NoHardCmps); returns 0 if the k-th hard
  * component is not registered yet. */

 [[nodiscard]] int get_bundle_size( int k ) const;

/*--------------------------------------------------------------------------*/
 /// returns the highest "name" currently in use in the k-th hard bundle
 /** In NDOFi-style master problems every item carries a persistent
  * integer "name" used to map slots to items; MasterPB does not maintain such a
  * mapping (PolyhedralFunction rows are simply indexed [0, get_nrows())
  * sequentially), so this method returns get_bundle_size(k). Provided as a
  * drop-in for the surrounding driver code that iterates with Index j =
  * Master->MaxName() ; j-- ; ... */

 [[nodiscard]] int get_max_cut_name( int k ) const
  { return( get_bundle_size( k ) ); }

/*--------------------------------------------------------------------------*/
 /// returns whether every hard-component bundle is empty
 /** Returns true iff get_bundle_size(k) == 0 for every k in
  * [0, NoHardCmps); a convenience for the typical "empty MP" check performed by
  * the surrounding Bundle algorithm at the very first iteration (or after a
  * complete bundle reset). */

 [[nodiscard]] bool is_bundle_empty( void ) const;

/*--------------------------------------------------------------------------*/
 /// remove the cut occupying slot \p slot from bundle B^k
 /** Deletes the cut stored at the NDOFi-style slot \p slot of the
  * k-th hard component; \p slot must lie in [0, MaxBSize) and must currently be
  * occupied (no-op if empty). The corresponding row of the PolyhedralFunction
  * is delete_row()-ed, and the slot->local-row mapping of MasterPB is updated
  * coherently (every other slot whose local row was past the removed one is
  * shifted down by one). */

 void remove_cut( int k , int slot );

/*--------------------------------------------------------------------------*/
 /// replace the cut at slot \p slot of B^k with (g, alpha)
 /** Overwrites both the subgradient and the linearization error of the
  * cut at slot \p slot of B^k. \p g and \p alpha follow the same sign
  * convention as add_cut(); the slot must currently be occupied. */

 void modify_cut( int k , int slot ,
                  std::vector< double > && g , double alpha );

/*--------------------------------------------------------------------------*/
 /// replace only the raw constant alpha_raw at slot \p slot of B^k
 /** Like modify_cut() but only the constant term is touched (a cheaper
  * Modification, C05FunctionModLin instead of full row replacement).
  *
  * Stored cut constants are otherwise *immutable* from the driver's point
  * of view: this entry point exists exclusively to keep the master in
  * sync with C05Function modifications that change the native cut
  * constant (e.g. when the inner Block evolves and the linearization is
  * still in the global pool); the legacy "shift after a stability-centre
  * move" use is gone, replaced by set_reference(). */

 void modify_alpha( int k , int slot , double alpha );

/*--------------------------------------------------------------------------*/
 /// returns the optimal multiplier theta at slot \p slot of B^k
 /** Returns the current value of the dual multiplier theta^k stored at
  * NDOFi-style slot \p slot of HardCmps[k]; the lookup goes through the slot ->
  * local-row mapping. Meaningful only after solve_master() and only in the dual
  * MP form; returns 0 if the slot is empty. */

 [[nodiscard]] double get_theta( int k , int slot ) const;

/*--------------------------------------------------------------------------*/
 /// returns the linearization error alpha at slot \p slot of B^k
 /** Returns the constant term b_i of the row of the PolyhedralFunction
  * of HardCmps[k] that currently lives in slot \p slot. Returns 0 if the slot
  * is empty. */

 [[nodiscard]] double get_alpha( int k , int slot ) const;

/*--------------------------------------------------------------------------*/
 /// returns the subgradient at slot \p slot of B^k
 /** Returns a const reference to the row vector of the
  * PolyhedralFunction of HardCmps[k] that currently lives in slot \p slot.
  * Returns an empty vector if the slot is empty. */

 [[nodiscard]] const std::vector< double > &
                            get_subgradient( int k , int slot ) const;

/*--------------------------------------------------------------------------*/
 /// is the cut at slot \p slot of B^k a "true" subgradient (vs feasibility)?
 /** Returns true iff the cut occupying slot \p slot of HardCmps[k] is a
  * diagonal linearization (a "true" subgradient), false if it is a vertical /
  * feasibility cut (in the PolyhedralFunction's is_row_vertical sense). Returns
  * false for empty slots; the NDOFi counterpart is Master->IsSubG( name ). */

 [[nodiscard]] bool is_subgradient( int k , int slot ) const;

/*--------------------------------------------------------------------------*/
 /// returns the total number of vertical (feasibility) cuts across all B^k
 /** NDOFi counterpart of Master->BCSize(): scans every HardCmps[k] and
  * counts the rows i for which PolyhedralFunction::is_row_vertical(i) is true,
  * returning the global tally. O(sum_k get_bundle_size(k)). */

 [[nodiscard]] int get_vertical_count( void ) const;

/*--------------------------------------------------------------------------*/
 /// returns Sigma* = aggregated linearization error of the dual bundle
 /** Computes the bundle-method "Sigma" = aggregated linearization error at
  * the current stability centre x_bar.
  *
  * Each cut is stored in its native (un-translated) form, i.e. as a pair
  * (g, alpha_raw) that bounds the function as
  *     F_k( x ) >= alpha_raw + g . x      (convex case)
  *     F_k( x ) <= alpha_raw + g . x      (concave case)
  * The per-cut linearization error at x_bar is therefore
  *     e_i = | F_k( x_bar ) - alpha_raw_i - g_i . x_bar |
  * (positive by construction). This method assembles
  *     Sigma_k = sum_i theta^k_i * e_i
  * on the fly from the bundle, the cached F_k( x_bar ) installed via
  * set_reference / set_F_at_x_bar, and the cached stability centre f_x_bar.
  *
  * \p k == -1 (default) sums over every hard component (the global aggregated
  * linearization error); \p k in [0, NoHardCmps) restricts the sum to component
  * k. Meaningful only after solve_master() and only in the dual MP form;
  * returns 0 otherwise. */

 [[nodiscard]] double get_aggregated_alpha( int k = -1 ) const;

/*--------------------------------------------------------------------------*/
 /// returns sum_i theta^k_i * alpha_raw_i (the un-translated convex combination)
 /** Companion of #get_aggregated_alpha that exposes the *raw* convex
  * combination of the native cut constants stored in the bundle, without
  * applying the F_k( x_bar ) / g . x_bar translation. Useful for callers
  * that want to perform the translation themselves or just need the raw
  * aggregate (e.g. diagnostics). \p k == -1 sums over every hard component;
  * returns 0 outside the dual MP form. */

 [[nodiscard]] double get_raw_aggregated_alpha( int k = -1 ) const;

/*--------------------------------------------------------------------------*/
 /// returns the model value v*[k] of the cutting-plane model at the step
 /** NDOFi counterpart of Master->ReadFiBLambda() / ReadFiBLambda(k+1).
  * In the primal MP this is just Var_v_hard[k].get_value() (per-cmp, \p k >= 0)
  * or the sum over every k (\p k == -1). In the dual MP the same quantity is
  * read off the per-PFB Objective contribution, which by LP duality equals
  * sum_i theta^k_i alpha^k_i (modulo the gamma^k * LB^k term when a global LB
  * is active), hence the call collapses onto get_aggregated_alpha(k).
  * Meaningful only after solve_master(). */

 [[nodiscard]] double get_FiBLambda( int k = -1 ) const;

/*--------------------------------------------------------------------------*/
 /// returns the aggregated subgradient z* in a fresh std::vector
 /** NDOFi counterpart of Master->ReadZ() (the dense, global, no-bse
  * variant). In the dual MP each z_j is a ColVariable and the call is just a
  * `for j: out[j] = Var_z[j].get_value()` materialisation; in the primal MP the
  * aggregated subgradient has to be reconstructed from the per-PFB theta
  * multipliers (currently not implemented: returns an empty vector). Meaningful
  * only after solve_master(). */

 [[nodiscard]] std::vector< double > get_z_vector( void ) const;

/*--------------------------------------------------------------------------*/
 /// returns Z[k] = sum_i theta^k_i * g^k_i for the k-th hard component
 /** Materialises in a fresh std::vector<double> of length NumVars the
  * aggregated subgradient of the k-th hard component, i.e. the convex
  * combination of the bundle rows g^k_i weighted by their optimal theta^k_i
  * multipliers. NDOFi counterpart of Master->ReadZ(..., k+1). Meaningful only
  * after solve_master() in the dual MP form; returns an all-zero vector in the
  * primal MP or when the slot is unavailable. */

 [[nodiscard]] std::vector< double > get_aggregated_subgradient( int k ) const;

/*--------------------------------------------------------------------------*/
 /// returns the primal direction d in a fresh std::vector
 /** NDOFi counterpart of Master->Readd(). In the primal MP the step
  * direction d is a static-variable group (Var_d); the call materialises its
  * current value. In the dual MP under proximal stabilization the primal
  * direction is d* = -t * z*, hence the result is -t_stab times get_z_vector().
  * Meaningful only after solve_master(). */

 [[nodiscard]] std::vector< double > get_d_vector( void ) const;

/*--------------------------------------------------------------------------*/
 /// returns z* . d, the scalar product of the aggregated subgradient and d
 /** NDOFi counterpart of Master->ReadGid() (the global, no-name variant).
  * Under proximal stabilization in the dual MP this is -t * || z* ||^2; the
  * implementation computes it directly via get_dual_norm_squared() and t_stab.
  * Meaningful only after solve_master(). */

 [[nodiscard]] double get_Gid_aggregate( void ) const;

/*--------------------------------------------------------------------------*/
 /// returns g^k_slot . d, the scalar product of the cut subgradient and d
 /** NDOFi counterpart of Master->ReadGid( name ) (the per-name variant).
  * \p k and \p slot identify the cut in the bundle of the k-th hard component;
  * the call returns sum_j g^k_slot[j] * d[j] using get_subgradient(k, slot) and
  * get_d_vector(). Returns 0 if the slot is empty or if d is not available. */

 [[nodiscard]] double get_Gid( int k , int slot ) const;

/*--------------------------------------------------------------------------*/
 /// invalidate every cut of B^k (NBModification on the underlying f_polyf)
 /** NDOFi counterpart of Master->ChgSubG(...) when the whole bundle B^k
  * has to be marked "stale" because the underlying subgradient values changed
  * (e.g. the active Variable of the C05Function were re-bound). MasterPB issues
  * an NBModification on the PolyhedralFunctionBlock of HardCmps[k], which is
  * the SMS++ "everything in this Block has just changed" signal: the
  * [MILP]Solver attached to MasterProblemBlock will then re-load the full
  * PolyhedralFunction on the next compute(). Returns silently if the slot is
  * empty / not a PFB. */

 void invalidate_subgradients( int k );

/*--------------------------------------------------------------------------*/
 /// bulk replacement of the linearization errors of B^k
 /** NDOFi counterpart of Master->ChgAlfa(Alfa.data(), k+1). Replaces
  * the constant terms of the bundle B^k with the values in \p alphas; \p alphas
  * must have size at least equal to MaxBSize and is read slot-by-slot through
  * slot_to_local[k], skipping empty slots. The underlying PolyhedralFunction
  * issues one C05FunctionModLin per touched row (a future optimisation could
  * collapse them with a single modify_constants call). */

 void set_alphas_bulk( int k , const std::vector< double > & alphas );

/*--------------------------------------------------------------------------*/
 /// linear-in-t sensitivity of v*(t) around the current t_stab
 /** NDOFi counterpart of Master->SensitAnals(vl, vc). The
  * The driver uses these two scalars to predict
  *     v( t_new ) >= vc + t_new * vl
  * in the long-term "hard" t-strategy. With the proximal dual MP
  *   v*(t) = sum_k sum_i theta^k_i alpha^k_i + gamma^k LB^k + x_bar . z*
  *           - (t/2) || z* ||^2
  * and assuming theta^k_i / z* approximately constant in a neighbourhood of the
  * current t_stab, the implementation returns
  *   vl = - || z* ||^2 / 2
  *   vc = v*(t_stab) - vl * t_stab           (linear extrapolation)
  * In the primal MP both are set to 0 (the sensitivity is not yet implemented
  * there). */

 void sensitivity_analysis( double & vl , double & vc ) const;

/*--------------------------------------------------------------------------*/
 /// returns the whole vector of linearization errors of B^k
 /** Returns a const reference to the full RealVector b of the
  * PolyhedralFunction of HardCmps[k] (i.e. every alpha^k_i in one shot, indexed
  * by the *local* position i in B^k). The reference stays valid until the next
  * bundle modification (add_cut / remove_cut / modify_cut / modify_alpha) on
  * the same component. */

 [[nodiscard]] const std::vector< double > & get_alphas( int k ) const;

/*--------------------------------------------------------------------------*/
 /// returns a fresh vector of optimal multipliers theta^k_i for B^k
 /** Materialises the std::list<ColVariable> f_theta of HardCmps[k] into
  * a std::vector<double> of length get_bundle_size(k), with entry i
  * containing the optimal theta^k_i (the call is O(n)). Meaningful only
  * after solve_master() and only in the dual MP form; returns an empty
  * vector otherwise. */

 [[nodiscard]] std::vector< double > get_thetas( int k ) const;

/*--------------------------------------------------------------------------*/
 /// set the stability-centre shift C used by the LB / Lvl rows
 /** Updates the additive shift
  *     C = -f^0( x_bar ) - sum_{k in H} f^{k,H}( x_bar )
  *  that turns LB into LB_xbar = LB + C and
  * Lvl into Lvl_xbar = Lvl + C. The shift is depositied as a member and
  * applied on top of the value passed to set_global_LB() / set_f_lev() /
  * set_LB() when their coefficients are committed to the dual Objective.
  * Should be called by the driver whenever the stability
  * centre x_bar moves; default is 0 (= no shift, suitable for kProximal
  * without explicit global LB / level). */

 void set_C( double C );

/*--------------------------------------------------------------------------*/
 /// set the global lower bound LB on the value of the sum-function
 /** In the dual MP the global lower bound enters as the linear coefficient
  * of the r multiplier in the master Objective (+ r * LB_xbar term, with
  * LB_xbar = LB + C, cf. set_C). This API forwards \p LB to that coefficient
  * via modify_term on the FRealObjective LinearFunction. No-op under the
  * primal MP (where the global LB has to be expressed differently). */

 void set_global_LB( double LB );

/*--------------------------------------------------------------------------*/
 /// set the stability centre x_bar of the master problem
 /** In the dual MP the stability centre x_bar enters as the linear
  * coefficient on every z_j (the term + x_bar * z in the Objective); this API
  * forwards every \p x_bar[j] to the LinearFunction coefficient of Var_z[j]. \p
  * x_bar must have size NumVars. No-op under the primal MP (where the step
  * direction d is the natural primal variable and x_bar contributes to the
  * constant gradient b through set_b). */

 void set_x_bar( const std::vector< double > & x_bar );

/*--------------------------------------------------------------------------*/
 /// atomically refresh the master reference (x_bar + per-component F_k)
 /** Single-call companion of set_x_bar that, in one shot, installs the new
  * stability centre and the (hard) per-component reference values F_k( x_bar ).
  * The atomic form exists so callers can not forget to update either side:
  * the master's notion of "current reference" is one indivisible piece of
  * state.
  *
  * \p x_bar must have size NumVars (same contract as set_x_bar).
  *
  * \p F_at_x_bar must have one entry per *hard* component, in the same order
  * used by HardCmps (i.e. F_at_x_bar[ k ] is F_k( x_bar ) for the k-th hard
  * component). Internally each value is cached as a member and used to
  * translate, on demand, the raw α stored in the bundle (each cut keeps the
  * native b such that  F_k( x ) >= b + g . x  ) into the linearization-error
  *   e_i = F_k( x_bar ) - b_i - g_i . x_bar
  * which the bundle-method stop test consumes; the cuts themselves stay
  * untouched. Default for every entry is 0.
  *
  * The two sides may still be poked individually via the lower-level
  * set_x_bar() / set_F_at_x_bar() if a partial refresh is genuinely
  * meaningful (e.g. the very first call when F_k( x_bar ) is not yet known),
  * but the supported pattern is the atomic one. */

 void set_reference( const std::vector< double > & x_bar ,
                     const std::vector< double > & F_at_x_bar );

/*--------------------------------------------------------------------------*/
 /// lower-level: set the per-component reference value F_k( x_bar )
 /** Per-component poke used internally by set_reference and exposed as a
  * fallback. The driver is normally expected to call set_reference()
  * instead. Default for every entry is 0.
  *
  * \p k is the index in HardCmps; must satisfy 0 <= k < NoHardCmps. */

 void set_F_at_x_bar( int k , double value );

/*--------------------------------------------------------------------------*/
 /// read back the cached reference value F_k( x_bar )
 /** Returns the last value passed to set_F_at_x_bar / set_reference for \p k ;
  * 0 if never set or if \p k is out of range. */

 [[nodiscard]] double get_F_at_x_bar( int k ) const;

/*--------------------------------------------------------------------------*/
 /// set the per-coordinate box  L <= x <= U  on the optimization variables
 /** Installs the lower / upper bounds on the optimization variable x that
  * the surrounding driver is minimising on. In terms of the step
  * variable  d = x - x_bar  the box reads
  *     L - x_bar <= d <= U - x_bar
  * which on the dual MP enters as
  *     + s^+_j ( L_j - x_bar_j ) - s^-_j ( U_j - x_bar_j )
  * in the objective with non-negative slacks s^+, s^- .
  *
  * For each coordinate j the corresponding slack(s) are unfixed only if
  * the matching side is finite; otherwise the slack stays fixed to 0,
  * meaning the bound does not really apply.
  *
  * \p L and \p U must both have size NumVars; passing an empty vector is
  * equivalent to "no bound on that side" (= all slacks of that side
  * stay / become fixed to 0). No-op under the primal MP. */

 void set_box( const std::vector< double > & L ,
               const std::vector< double > & U );

/*--------------------------------------------------------------------------*/
 /// set the upper bundle u_bar of an #kUpperLower stabilization
 /** Companion of #set_x_bar to be used under #kUpperLower stabilization
  * (which is currently a reserved-but-not-implemented type, see
  * #stabilization_type): \p u_bar is the upper bundle centre, around
  * which the proximal term is anchored on the "upper" side. Throws when
  * called under any other stabilization. \p u_bar must have size
  * NumVars.
  *
  * Not yet implemented: the method exists as an API skeleton so callers
  * can be written against the final interface today, but it currently
  * throws std::logic_error unconditionally. */

 void set_u_bar( const std::vector< double > & u_bar );

/*--------------------------------------------------------------------------*/
 /// set the lower bundle l_bar of an #kUpperLower stabilization
 /** Companion of #set_u_bar: \p l_bar is the lower bundle centre, around
  * which the proximal term is anchored on the "lower" side. Throws when
  * called under any stabilization other than #kUpperLower. \p l_bar
  * must have size NumVars.
  *
  * Not yet implemented: the method exists as an API skeleton so callers
  * can be written against the final interface today, but it currently
  * throws std::logic_error unconditionally. */

 void set_l_bar( const std::vector< double > & l_bar );

/*--------------------------------------------------------------------------*/
 /// install the constant linear part \p b of the original sum-function
 /** The original objective minimized by the bundle method has the structure
  * \f[
  *    \min \; b^\top \lambda \,+\, \sum_{k=1}^K F_k( \lambda )
  * \f]
  * where every \f$F_k\f$ is a generic convex component (handled through a
  * PolyhedralFunctionBlock sub-Block) and \f$b\f$ is the gradient of the
  * single linear "0-th" component, an affine term of the sum-function whose
  * value at \f$\lambda\f$ is exactly \f$b^\top \lambda\f$.
  *
  * In the dual master problem the bundle variables \f$\theta\f$ and
  * \f$\gamma\f$ are tied to the step direction \f$z\f$ by the coupling
  * equations
  * \f[
  *    z_j \;=\; b_j \;-\; \sum_{k=1}^K \sum_{i \in B_k}
  *                          \theta^k_i \, A^k_{i,j}
  *    \quad \forall j = 0, \dots, NumVars - 1
  * \f]
  * (CreateDualMP builds them with rhs zero, and this method installs the
  * actual rhs \f$b_j\f$ for every \f$j\f$). The \f$b_j\f$ term acts as the
  * "free" drift of \f$z\f$ that exists irrespective of the bundle content:
  * if all \f$\theta^k_i\f$ vanish (empty bundle or a degenerate dual
  * solution) the equation collapses to \f$z = b\f$, so \f$z\f$ inherits the
  * direction of steepest ascent of the linear component.
  *
  * \p b must have size NumVars. Under the primal MP the linear term is
  * absorbed directly in the master Objective via set_b and this method is
  * a no-op. */

 void set_linear_part( const std::vector< double > & b );

/*--------------------------------------------------------------------------*/
 /// append \p n new optimization variables to the Master Problem
 /** Drop-in for Master->AddVars(n). Extends the master problem from
  * NumVars to NumVars + n coordinates by appending n new entries to Var_d /
  * Var_v_hard / Var_z and growing CouplingCns accordingly. This is a structural
  * change and forces a fresh load_problem() of the registered [MILP]Solver on
  * the next compute().
  *
  * NOT YET IMPLEMENTED -- throws std::logic_error. Adding NumVars on the fly
  * requires rebuilding the diagonal-quadratic part of the Objective (the per-d
  * / per-z triples) and re-wiring every PolyhedralFunction- Block sub-Block via
  * set_variables() / set_conjugate_constraint(). */

 void add_vars( int n );

/*--------------------------------------------------------------------------*/
 /// remove a subset of optimization variables from the Master Problem
 /** Drop-in for Master->RmvVars(subset, sz). Removes the \p sz coordinates
  * listed in \p subset (or *all* coordinates if \p subset == nullptr) from
  * Var_d / Var_v_hard / Var_z, and patches CouplingCns / every
  * PolyhedralFunctionBlock sub-Block accordingly.
  *
  * NOT YET IMPLEMENTED -- throws std::logic_error. Same caveats as add_vars: it
  * is a structural change that forces a fresh load_problem(). */

 void remove_vars( const int * subset , int sz );

/*--------------------------------------------------------------------------*/
 /// forward an ostream to the [MILP]Solver attached to MasterPB
 /** Sets the log ostream on the first Solver registered to
  * MasterProblemBlock. The verbosity itself is left to the Solver's own
  * ComputeConfig (typically the `intLogVerb` knob in the MPBSolverCfg
  * file), so the surrounding driver does not need to
  * carry a duplicate "MP log verbosity" parameter. No-op if no Solver
  * is registered. */

 void forward_log( std::ostream * log_stream );

/*--------------------------------------------------------------------------*/
 /// cap the master Solver running time to \p t (seconds)
 /** Forwards `dblMaxTime = t` to the first Solver registered to
  * MasterProblemBlock via set_par. No-op if no Solver is registered or if the
  * Solver does not advertise the dblMaxTime parameter. */

 void set_max_time( double t );

/*--------------------------------------------------------------------------*/
 /// solve the Master Problem by triggering compute() on the inner Solver
 /** Forwards to the compute() method of the first Solver registered to
  * MasterProblemBlock through register_Solver(); throws std::logic_error if no
  * Solver is attached. Returns the status returned by compute(). */

 int solve_master( void );

/*--------------------------------------------------------------------------*/
 /// returns the current optimal value of the global multiplier lambda
 /** lambda is the master-side non-negative dual multiplier paired with the
  * model-value equation of the lower model (
  * stationarity (i): lambda + r - omega = 1), shared across every hard
  * component (cf. PolyhedralFunctionBlock::set_lambda). Meaningful only
  * after solve_master(). */

 [[nodiscard]] double get_lambda( void ) const
  { return( Var_lambda.get_value() ); }

/*--------------------------------------------------------------------------*/
 /// returns the current optimal value of r
 /** r is the dual multiplier of the global LB row. Meaningful only after
  * solve_master(). */

 [[nodiscard]] double get_r( void ) const { return( Var_r.get_value() ); }

/*--------------------------------------------------------------------------*/
 /// returns the current optimal value of omega
 /** omega is the dual multiplier of the level/X row. Fixed to 0 under
  * #kProximal. Meaningful only after solve_master(). */

 [[nodiscard]] double get_omega( void ) const
  { return( Var_omega.get_value() ); }

/*--------------------------------------------------------------------------*/
 /// returns the current optimal value of z_j
 /** z_j is the j-th coordinate of the aggregated subgradient z*. \p j
  * must lie in [0, NumVars). Meaningful only after solve_master(). */

 [[nodiscard]] double get_z( int j ) const { return( Var_z[ j ].get_value() ); }

/*--------------------------------------------------------------------------*/
 /// update the proximal stabilization parameter t
 /** Sets the proximal stabilization parameter t used in the quadratic /
  * linear term of the MP Objective. The new value is stored in #t_stab; the
  * coefficients of the Objective are updated only if the MP has already been
  * built (i.e. after a CreateEmptyMP() call) so that the same value can also be
  * configured before the MP is generated.
  *
  * \note this method only updates the proximal coefficient. The Bundle
  *       algorithm is responsible for issuing the call at every t-change
  *       and for triggering a re-solve of the MP. */

 void set_t( double t );

/*--------------------------------------------------------------------------*/
 /// returns the current value of the proximal stabilization parameter t
 /** Returns the proximal stabilization parameter t currently used in the
  * quadratic / linear term of the MP Objective. Needed by the caller (e.g.
  * the driver) to convert the master-side step
  * d* = -t * z* back to a t-independent displacement when applying a
  * different Tau (as in t-strategies). */

 [[nodiscard]] double get_t() const { return( t_stab ); }

/*--------------------------------------------------------------------------*/
 /// update the linear coefficient vector b of the primal Objective
 /** Sets the linear coefficient on every d_i of the primal MP Objective
  * to the value \p b[i] (i.e. the "b" of the paper's primal MP b * d + sum_k
  * v^k + (1/(2t)) || d ||^2_2). \p b must have size #NumVars. The call is a
  * no-op under the dual MP. */

 void set_b( const std::vector< double > & b );

/*--------------------------------------------------------------------------*/
 /// set the lower bound LB^k of the k-th hard component
 /** Sets the constraint v^k >= LB^k of the primal MP (or, equivalently,
  * the linear coefficient on gamma^k of the dual one). \p k must lie in [0,
  * NoHardCmps). LB^k = -INF turns the bound off. */

 void set_LB( int k , double LB );

/*--------------------------------------------------------------------------*/
 /// update the level stabilization parameter f_lev
 /** Sets the level target f_lev used by the #kLevel / #kDoublyStabilized
  * stabilization schemes. The value enters the dual MP objective as the linear
  * coefficient of the omega multiplier (+ omega * f_lev). It is a no-op when
  * #StblType == #kProximal (omega is fixed to 0 in that case).
  *
  * \note the surrounding Bundle algorithm is responsible for choosing
  *       f_lev consistently with its level-management strategy. */

 void set_f_lev( double f );

/*--------------------------------------------------------------------------*/
 /// load a MasterProblemBlock out of an istream
 /** Required by the abstract interface of Block, but currently not
  * implemented (a MasterProblemBlock is always built programmatically by its
  * driving solver): always throws an std::logic_error. */

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
 bool IsConvex;     ///< whether the C05Function is convex (true) or concave

 stabilization_type StblType;
                          ///< type of stabilization, see #stabilization_type

 int MaxBSize;      ///< maximum number of items kept per bundle B^k

 int MaxSGLen;      ///< maximum length of a subgradient (currently == NumVars)

 int NumVars;       ///< number of optimization variables in the Primal MP

 int NoTotCmps;     ///< total number of components in the sum-function

 int NoEasyCmps;    ///< number of "easy" components

 int NoHardCmps;    ///< number of "hard" components
                    ///< ( == NoTotCmps - NoEasyCmps)

 int DoEasy;        ///< bit-wise flag controlling the easy-component handling

 std::vector< bool > IsEasyCmp;
                              ///< IsEasyCmp[k] == true iff component k is easy

 // - - - - - - - - -  pointers to the per-component sub-Blocks - - - - - - - -

 std::vector< Block * > EasyCmps;  ///< sub-Blocks of the "easy" components

 std::vector< Block * > HardCmps;  ///< sub-Blocks of the "hard" components

 /// metadata of one absorbed BendersBFunction row in the primal MP
 /** When absorb_BBF_into_primal_MP() processes the i-th mapping row of a
  * BendersBFunction it (i) relaxes the original RowConstraint C_i in the
  * inner Block and (ii) re-installs an equivalent FRowConstraint
  *
  *     g_i( y ) - A_i . d  [<=, =, >=]  ( A_i . x_bar + b_i )
  *
  * on *this*, where A_i = bbf->get_A()[ i ], b_i = bbf->get_b()[ i ] and
  * the d on the master side stands for the step variables Var_d (the
  * BBF active variables x are split as x = x_bar + d). EasyBBFRow stores
  * the bookkeeping needed by set_x_bar() to refresh the absorbed
  * RowConstraint's right-hand side(s) whenever the stability centre
  * changes. */

 struct EasyBBFRow {
  FRowConstraint * cns;                ///< absorbed RowConstraint on *this*
  std::vector< double > A_row;         ///< i-th row of bbf->get_A()
  double b_i;                          ///< i-th entry of bbf->get_b()
  char side;                           ///< BendersBFunction::ConstraintSide
  double orig_lhs;                     ///< C_i->get_lhs() at absorption time
  double orig_rhs;                     ///< C_i->get_rhs() at absorption time
  };

 std::vector< EasyBBFRow > EasyBBFRows;

 // - - - - - - - - - - -  static MP entities (primal form)  - - - - - - - - -

 std::vector< ColVariable > Var_d;
                                ///< the step variables d (free, size NumVars)

 std::vector< ColVariable > Var_v_hard;
                                ///< the epigraph variables v^k
                                ///< per hard component

 std::vector< BoxConstraint > Bounds_v_hard;
                                ///< per-hard-component LB^k: v^k >= LB^k
                                ///< (rhs = +INF by default)

 FRowConstraint LevelCns;                ///< primal level constraint
                                         ///< sum_k v^k <= f_lev
                                         ///< (kLevel / kDoublyStabilized only)

 // - - - - - - - - - - - -  static MP entities (dual form)  - - - - - - - - -

 ColVariable Var_lambda;
                                   ///< global non-negative dual multiplier
                                   ///< of the model-value equation of the
                                   ///< lower model. The *same* Var_lambda
                                   ///< enters with coefficient +1 the
                                   ///< simplex (= normalization) row of
                                   ///< every hard-component PFB sub-Block
                                   ///< (cf. PolyhedralFunctionBlock::set_-
                                   ///< lambda) and with coefficient +1 the
                                   ///< master-side NormalizationCns
                                   ///< (lambda + r - omega = 1).

 std::vector< ColVariable > Var_s_plus;
                                   ///< non-negative slack multipliers s^+
                                   ///< paired with the lower side of the
                                   ///< box  L - x_bar <= d  (cf.
                                   ///< 
                                   ///< / (48)). One entry per coordinate
                                   ///< (size NumVars); coordinates without
                                   ///< a finite L are kept fixed to 0,
                                   ///< i.e. the corresponding slack does
                                   ///< not really exist.

 std::vector< ColVariable > Var_s_minus;
                                   ///< non-negative slack multipliers s^-
                                   ///< paired with the upper side of the
                                   ///< box  d <= U - x_bar  (cf.
                                   ///< 
                                   ///< / (48)). One entry per coordinate
                                   ///< (size NumVars); coordinates without
                                   ///< a finite U are kept fixed to 0.

 ColVariable Var_r;                ///< dual multiplier of the global LB row

 ColVariable Var_omega;            ///< dual multiplier of the level / X row

 std::vector< ColVariable > Var_z;
                                   ///< auxiliary dual variables z (one per
                                   ///< coordinate of the original sum-function
                                   ///< variable space; size NumVars)

 FRowConstraint NormalizationCns;
                                   ///< global normalization row
                                   ///< sum_k lambda_k + r - omega = NoHardCmps
                                   ///< (one +1 lambda_k coefficient per hard
                                   ///< component, see Var_lambdas)

 std::list< FRowConstraint > CouplingCns;
                                   ///< coupling rows z_j = b_j
                                   ///< (j = 0 .. NumVars-1); populated by each
                                   ///< hard-cmp sub-Block via PolyhedralFunc-
                                   ///< tionBlock::set_conjugate_constraint

 // - - - - - - - - - - - - - - -  stabilization parameters  - - - - - - - -

 double t_stab;     ///< current value of the proximal parameter t

 double f_C = 0.0;  ///< additive shift used to express LB_xbar = LB + C
                    ///< and Lvl_xbar = Lvl + C, with C = -f^0( x_bar )
                    ///< - sum_H f^{k,H}( x_bar ); refreshed by set_C() at
                    ///< every stability-centre move

 std::vector< double > f_L;
                    ///< per-coordinate lower bound L on x for the box
                    ///< constraint L - x_bar <= d <= U - x_bar; a
                    ///< non-finite L_t means "no lower bound on
                    ///< coordinate t", in which case s^+_t stays fixed
                    ///< to 0 (= the slack does not really exist)

 std::vector< double > f_U;
                    ///< per-coordinate upper bound U on x; mirrors f_L
                    ///< with s^-_t playing the role of s^+_t

 std::vector< double > f_x_bar;
                    ///< snapshot of the latest stability centre passed
                    ///< to set_x_bar(), used by set_box() and the box
                    ///< refresh logic to compute (L - x_bar) and
                    ///< (U - x_bar) without re-asking the caller

 std::vector< double > f_F_at_x_bar;
                    ///< per-hard-component cache of  F_k( x_bar ) ; sized
                    ///< NoHardCmps in configure(), refreshed by the driver
                    ///< via set_F_at_x_bar() at every stability-centre
                    ///< change. Used (or will be used) to derive Sigma in
                    ///< linearization-error form without forcing the
                    ///< driver to pre-promote each cut's alpha

 double f_lev;      ///< current value of the level f_lev (level / doubly only)

 int z_obj_idx;     ///< index of the first z_j entry in the DQuadFunction
                    ///< triples (the NumVars entries z_0..z_{NumVars-1} are
                    ///< laid out contiguously), or -1 if absent

 int r_obj_idx;     ///< index of the r multiplier in the DQuadFunction
                    ///< triples; carries the (+ r * LB) global lower
                    ///< bound contribution, refreshed by set_global_LB

 int omega_obj_idx; ///< index of omega in the DQuadFunction triples, or -1
                    ///< if omega does not contribute to the master Objective
                    ///< (i.e. under #kProximal)

 int s_plus_obj_idx  = -1;
                    ///< index of the first s^+_j entry in the DQuadFunction
                    ///< triples (the NumVars s^+ entries are laid out
                    ///< contiguously); carries the +sgn*(L_j - x_bar_j)
                    ///< coefficient updated by set_x_bar / set_box

 int s_minus_obj_idx = -1;
                    ///< index of the first s^-_j entry in the DQuadFunction
                    ///< triples; carries the -sgn*(U_j - x_bar_j) coefficient

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// slot -> local-row mapping into PolyhedralFunction::v_A of HardCmps[k]
 /** slot_to_local[k] has size MaxBSize: slot_to_local[k][s] is the local
  * index of the cut occupying slot s in the PolyhedralFunction of the k-th hard
  * component, or -1 if the slot is empty. PolyhedralFunction rows shift on
  * delete_row(); this map absorbs the shift so that the surrounding
  * the driver can keep its NDOFi-style persistent name (ItemVcblr) intact
  * across remove_cut() calls. */
 std::vector< std::vector< int > > slot_to_local;

 /// model Block stolen by configure() (or nullptr)
 Block * f_original_block = nullptr;

 /// set of sub-Blocks the inner Solver must ignore; populated by
 /// configure() (and/or by the two-argument register_Solver() overload)
 /// and forwarded to the inner Solver via BlockSolverConfig::apply()'s
 /// second argument when register_Solver() is called (or already if a
 /// Solver is registered). Stored as eagerly-expandable set of typed
 /// Block * (the descendants are added by Solver::set_excluded_blocks
 /// itself); see Solver.h.
 std::unordered_set< Block * > f_ignored_blocks;

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
