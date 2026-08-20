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
 * Following the SMS++ distinction between physical and abstract representations,
 * the physical representation of a MasterProblemBlock is the compact state that
 * characterizes the current master problem: the registered sub-Blocks, the
 * current stability center, and the stabilization parameters (the proximal
 * parameter t, the level value Lvl and the selected stabilization type). The
 * Variable, Constraint and Objective objects forming the primal or dual MP are
 * the corresponding abstract representation generated from that physical state.
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

#include <algorithm>
#include <numeric>
#include <iosfwd>
#include <list>
#include <stdexcept>
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
class DQuadFunction;
class PolyhedralFunctionBlock;

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

 /** \param is_easy
  *        per-global-component flag distinguishing easy from hard. When
  *        provided (and of size num_hard_cmps + easy_components.size())
  *        it is copied into MasterProblemBlock::IsEasyCmp, which the
  *        per-component getters (\see get_FiBLambda(k), and similar)
  *        consult to dispatch the EASY branch. Empty by default for
  *        backward compatibility with callers that do not yet maintain
  *        such a flag vector; in that case IsEasyCmp stays at the
  *        all-false default initialised by SetDim() and only the legacy
  *        all-hard path is taken.
  *  \param hard_cmp_scaling
  *        bit-wise numerical scaling requested for every hard-component
  *        PolyhedralFunctionBlock: bit 0 enables local row scaling and bit 1
  *        enables global epigraph scaling. Thus 0 = none, 1 = local only,
  *        2 = global only, 3 = both. */
 void configure( bool primal ,
                 int max_bundle_size ,
                 int num_vars ,
                 int num_hard_cmps ,
                 const std::vector< C05Function * > & easy_components ,
                 Block * original_block ,
                 std::unordered_set< Block * > ignored_blocks = {} ,
                 stabilization_type reg = kDoublyStabilized ,
                 bool convex = true ,
                 const std::vector< bool > & is_easy = {} ,
                 int hard_cmp_scaling = 0 ,
                 const std::vector< std::vector< Index > > & easy_local2global = {} );

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
 /** Populates the physical representation of MasterProblemBlock with the
  * compact ingredients of the Master Problem, in either the primal or the dual
  * formulation depending on the chosen stabilization scheme and on the presence
  * of "easy" components. The Variable, Constraint and Objective objects are
  * left to the generate_abstract_*() methods.
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
 /// initialize the primal physical representation of the Master Problem
 /** Builds only the physical representation of the primal MP: the registered
  * PolyhedralFunctionBlock sub-Blocks, the current stability-center data and
  * the stabilization parameters. The corresponding Variable, Constraint and
  * Objective objects are generated later by generate_abstract_variables(),
  * generate_abstract_constraints() and generate_objective().
  *
  * This method is called internally by CreateEmptyMP() when the dual form is not
  * strictly required, and is exposed publicly so that derived classes can
  * override the physical construction. */

 void CreatePrimalMP( stabilization_type Stbl );

/*--------------------------------------------------------------------------*/
 /// initialize the dual physical representation of the Master Problem
 /** Builds only the physical representation of the dual MP: the registered
  * PolyhedralFunctionBlock sub-Blocks, the current stability-center data and
  * the stabilization parameters. The corresponding Variable, Constraint and
  * Objective objects are generated later by generate_abstract_variables(),
  * generate_abstract_constraints() and generate_objective().
  *
  * This method is called internally by CreateEmptyMP() whenever "easy"
  * components are present, since their compact description is naturally
  * expressed in the dual form; it is exposed publicly so that derived classes
  * can override the physical construction. */

 void CreateDualMP( stabilization_type Stbl );

/*--------------------------------------------------------------------------*/
 /// generate the Variable groups in the abstract representation
 /** Dispatches to the primal or dual generator according to the currently
  * selected MP form. The optional Configuration is interpreted as in
  * PolyhedralFunctionBlock: a SimpleConfiguration<int> may carry already-built
  * abstract-representation bits, so repeated calls can be harmless no-ops. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the Constraint groups in the abstract representation
 /** Dispatches to the primal or dual constraint generator. It requires that
  * generate_abstract_variables() has already materialized the Variable side and
  * records the generated-constraint stage to make repeated calls idempotent. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the Objective in the abstract representation
 /** Dispatches to the primal or dual objective generator. It requires the
  * abstract Variable side to exist and records the generated-objective stage
  * to make repeated calls idempotent. */

 void generate_objective( Configuration * objc = nullptr ) override;

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
 /// Append the contribution of all easy LagBFunction components to the dual
 /// coupling rows.
 /** In the dual MP, each easy LagBFunction contributes to the stationarity
  * rows associated with the master-space variables. For every active variable
  * y_i of an easy LagBFunction, this method retrieves the corresponding
  * Lagrangian term
  *
  *     g_i( u ) = sum_h a_ih u_h + c_i
  *
  * through LagBFunction::get_Lagrangian_term( i ), maps the local active
  * index i to the global master coordinate j, and appends the
  * objective-sense-adjusted easy subgradient to CouplingCns[ j ].
  * Thus the dual coupling row becomes
  *
  *     z_j - lambda b_j + s^+_j - s^-_j
  *          - sign( F_internal ) sum_E g^k_j( u^k )
  *          - hard-component terms = 0.
  *
  * This method is meaningful only in the dual MP and only after CouplingCns
  * has been initialized.
  */

 void add_LBF_to_coupling_rows(
  std::vector< LinearFunction::v_coeff_pair > & vp_Cns );

/*--------------------------------------------------------------------------*/
 /// Map a local active-variable index of an easy component to the global
 /// master-space coordinate.
 /** For dense Lambda representations, no explicit map is stored and the
  * local index is already the global index, so this method returns local_i.
  *
  * For sparse Lambda representations, EasyLocal2Global[ easy_id ][ local_i ]
  * stores the global master coordinate corresponding to the local active
  * variable local_i of easy component easy_id.
  *
  * @param easy_id  local index of the easy component among EasyCmps.
  * @param local_i  local active-variable index inside that easy component.
  *
  * @return the global coordinate j such that the corresponding contribution
  *         must be inserted into CouplingCns[ j ].
  */

 [[nodiscard]] Index easy_local_to_global( Index easy_id ,
                                           Index local_i ) const;

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
 /// returns ||z*||^2, the squared 2-norm of the aggregate subgradient
 /** Reads the normalized aggregate returned by get_z_vector() and squares it.
  * This matters in pure-level dual form because Var_z stores eta z*, not z*
  * itself; callers of this method always get the BundleSolver z* norm. */

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
  * picks it up on the next compute(). Each bundle item owns a physical row,
  * even when its coefficients coincide with another item: duplicate
  * management is an algorithmic responsibility of BundleSolver, not of the
  * master representation. The return value is #kCutInserted. */

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
 /// returns whether any hard component is in the "stuck simplex" state
 /** Returns true iff there exists a hard component k whose dual simplex
  * row
  *     sum_i theta^k_i + gamma^k = lambda
  * cannot be made satisfiable by the inner :MILPSolver, because:
  *  - gamma^k is fixed to 0 (no global LB has been installed on
  *    HardCmps[k], so the gamma^k * LB^k contribution is structurally
  *    absent), AND
  *  - the bundle of theta^k_i is empty (no cuts pushed via add_cut()
  *    have survived for this component).
  * Under these conditions the row collapses to lambda = 0, which
  * contradicts K * lambda = 1 (always true with Var_r and Var_omega
  * pinned). The surrounding driver can use this hook to detect a
  * transient warm-start state - typically after an oracle mutation
  * has invalidated every cut of some component - and short-circuit
  * the master solve with a trivial direction, letting the next
  * Fi(.) evaluation refill the affected bundles before the master
  * is asked to solve a non-degenerate dual MP. */

 [[nodiscard]] bool has_pinned_empty_cmp( void ) const;

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
 /// find an item with exactly the same coefficient vector
 /** Returns the slot of a cut in component \p k whose coefficient vector is
  * element-wise identical to \p g and whose diagonal/vertical type matches
  * \p is_vert. The comparison is exact, matching the historical
  * MPSolver::CheckBCopy semantics. Returns -1 if no copy exists. */

 [[nodiscard]] int find_identical_cut(
                         int k , const std::vector< double > & g ,
                         bool is_vert = false ) const;

/*--------------------------------------------------------------------------*/
 /// convert an incoming raw constant to the representation stored by MPB
 /** Applies the same reference-frame and objective-sense transformation used
  * by add_cut() without modifying the master. */

 [[nodiscard]] double get_stored_constant(
                         int k , const std::vector< double > & g ,
                         double alpha , bool is_vert = false ) const;

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
  * k. In pure-level dual form the row masses carry eta, so the returned Sigma
  * is divided by eta to stay in the normalized bundle-method units. Meaningful
  * only after solve_master() and only in the dual MP form; returns 0 otherwise. */

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
 /// returns sum_i theta_i alpha_i + gamma LB for the k-th hard component
 /** This is the un-normalized raw constant used when building the aggregate
  * V2 cut. The theta part is the same raw aggregate returned by
  * get_raw_aggregated_alpha(k); the gamma contribution is the horizontal
  * lower-bound cut, if a genuine finite LB^k is installed. */

 [[nodiscard]] double get_raw_aggregated_alpha_with_LB( int k ) const;

/*--------------------------------------------------------------------------*/
 /// returns the model value v*[k] of the cutting-plane model at the step
 /** NDOFi counterpart of Master->ReadFiBLambda() / ReadFiBLambda(k+1).
  * In the primal MP this is just Var_v_hard[k].get_value() (per-cmp, \p k >= 0)
  * or the sum over every k (\p k == -1). In the dual MP the same quantity is
  * read from the dual optimality identities: proximal uses d* = -t z*, while
  * pure level uses d* = -eta z* because the level row multiplier carries the
  * aggregate mass. Meaningful only after solve_master(). */

 [[nodiscard]] double get_FiBLambda( int k = -1 ) const;

/*--------------------------------------------------------------------------*/
 /// returns the aggregated subgradient z* in a fresh std::vector
 /** NDOFi counterpart of Master->ReadZ() (the dense, global, no-bse
  * variant). In the dual MP Var_z is returned directly for proximal/doubly
  * stabilization, but in pure level it is divided by eta because the dual
  * stationarity vector stores eta z*. In the primal MP z* is reconstructed from
  * the solved displacement and the active stabilization. Meaningful only after
  * solve_master(). */

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
 /** NDOFi counterpart of Master->Readd(). In the translated primal MP the step
  * d is stored directly in Var_d. In the raw primal MP Var_d stores absolute x
  * and the call returns x - x_bar. In the dual MP the proximal identity is
  * d* = -t z*; in pure level Var_z stores eta z*, so d* = -Var_z. Meaningful
  * only after solve_master(). */

 [[nodiscard]] std::vector< double > get_d_vector( void ) const;

/*--------------------------------------------------------------------------*/
 /// returns the multiplier of the level constraint, if available
 /** In primal form this is the dual value of the explicit level row; in
  *  dual form this is the value of the explicit omega variable. Meaningful
  *  after solve_master() in level / doubly-stabilized modes. */

 [[nodiscard]] double get_level_multiplier( void ) const;

/*--------------------------------------------------------------------------*/
 /// returns z* . d, the scalar product of the aggregated subgradient and d
 /** NDOFi counterpart of Master->ReadGid() (the global, no-name variant).
  * get_z_vector() and get_d_vector() first map the solved primal/dual master
  * representation to the common physical quantities; this method then computes
  * the same scalar product for every stabilization. Meaningful only after
  * solve_master(). */

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
 /// select the master storage frame ( translated vs raw/iterate )
 /** Select which algebraically-equivalent storage frame the master uses.
  * The choice is INTERNAL to MasterProblemBlock and invisible to the driver,
  * which always passes / reads raw alpha.
  *
  * - 0 = displacement form ( default, production ): cuts are stored in the
  *   linearization-error frame  b = F_k( x_bar ) - alpha + g . x_bar; the
  *   x_bar dependence is "baked" into the cut constants ( which stay ~ 0 at
  *   tight cuts, well scaled ), the linear coefficient on Var_z is 0, and
  *   set_reference() re-shifts the constants as the centre moves. This path
  *   is left exactly as-is.
  *
  * - 1 = raw/iterate form: in the primal MP the optimization variable is the
  *   absolute point x and the proximal term is ||x-x_bar||^2/(2t). In the
  *   dual MP the x_bar dependence is carried by the explicit linear
  *   coefficient sgn * x_bar on Var_z and the cut constant is stored
  *   function-value-relative as b = -alpha + F_k( x_bar ).
  *
  * Set before configure(); switching mid-solve is not supported.
  * Out-of-range -> 0. */

 void set_v2_form( int form = 1 ) noexcept
  { f_v2_form = ( form < 0 || form > 1 ) ? 0 : form; }

/*--------------------------------------------------------------------------*/
 /// query the master storage frame ( 0 = translated, 1 = raw/iterate )

 [[ nodiscard ]] int get_v2_form() const noexcept { return( f_v2_form ); }

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

 /// read back the current stability centre

 [[nodiscard]] const std::vector< double > & get_x_bar( void ) const
  { return( f_x_bar ); }

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
 /// set the "lazy reference" tolerance for the displacement form
 /** The diagonal cuts store the linearization error at a reference x_ref. By
  * default x_ref tracks the stability centre x_bar at every serious step
  * ( tol == 0 ), i.e. the plain displacement form that re-bakes the per-cut
  * g . delta dot product on every move. Setting tol > 0 lets x_ref lag behind
  * x_bar: the residual g . ( x_bar - x_ref ) is carried by the explicit lin-z
  * coefficient on Var_z ( O( NumVars ) per move, no per-cut dot product ), and
  * x_ref is re-aligned to x_bar ( paying the per-cut shift once ) only when
  * || x_bar - x_ref ||_inf > tol. The master Problem is UNCHANGED for any tol
  * ( exact identity ); tol only trades the per-serious-step cost against the
  * per-cut constant magnitude. tol == 0 is a strict no-op. */

 void set_xref_tol( double tol ) noexcept
  { f_xref_tol = ( tol > 0.0 ) ? tol : 0.0; }

/*--------------------------------------------------------------------------*/
 /// read back the cached reference value F_k( x_bar )
 /** Returns the last value passed to set_F_at_x_bar / set_reference for \p k ;
  * 0 if never set or if \p k is out of range. */

 [[nodiscard]] double get_F_at_x_bar( int k ) const;

 /// read back all cached component values of the current reference

 [[nodiscard]] const std::vector< double > & get_F_at_x_bar( void ) const
  { return( f_F_at_x_bar ); }

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
  * In the primal MP the same box is installed directly on the absolute x
  * variables in raw form, or translated to [ L - x_bar , U - x_bar ] on the
  * displacement variables.
  *
  * \p L and \p U must both have size NumVars; passing an empty vector is
  * equivalent to "no bound on that side" (= all slacks of that side
  * stay / become fixed to 0). */

 void set_box( const std::vector< double > & L ,
               const std::vector< double > & U );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// modify the box on a contiguous range of coordinates
 /** Installs the values in \p L and \p U on the left-closed, right-open
  * interval [ range.first , range.second ). Each nonempty vector must have
  * range.second - range.first entries. An empty \p L (respectively \p U)
  * removes that side of the box throughout the range.
  *
  * A changed lower side issues a MasterProblemRngdMod of type
  * #MasterProblemMod::BoxLowerChanged carrying the new lower values; a
  * changed upper side analogously issues #MasterProblemMod::BoxUpperChanged.
  * Each Modification therefore carries one direct vector of new values. */

 void set_box( std::vector< double > L , std::vector< double > U ,
               Range range );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// modify the box on an arbitrary subset of coordinates
 /** Installs L[ i ] and U[ i ] on coordinate subset[ i ]. Each nonempty
  * value vector must have subset.size() entries. An empty \p L
  * (respectively \p U) removes that side of the box on the entire subset.
  * Repeated or out-of-range indices are rejected.
  *
  * If \p ordered is false, the subset is sorted and both value vectors are
  * reordered with it before the change is applied. Each changed side issues
  * its own MasterProblemSbstMod of type BoxLowerChanged or BoxUpperChanged;
  * each Modification exposes a strictly increasing subset and one aligned
  * vector of new values. */

 void set_box( std::vector< double > L , std::vector< double > U ,
               Subset subset , bool ordered = false );

 /// read back the lower-bound vector of the current box

 [[nodiscard]] const std::vector< double > & get_box_lower( void ) const
  { return( f_L ); }

 /// read back the upper-bound vector of the current box

 [[nodiscard]] const std::vector< double > & get_box_upper( void ) const
  { return( f_U ); }

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
  * absorbed directly in the master Objective via the lower-level set_b()
  * helper; this method still owns the physical-state update and issues the
  * corresponding Modification. */

 void set_linear_part( const std::vector< double > & b );

 /// read back the current linear part b

 [[nodiscard]] const std::vector< double > & get_linear_part( void ) const
  { return( f_linear_part ); }

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
 /// manage the one-shot proximal objective used to seed pure level

 [[nodiscard]] bool has_initial_level_objective( void ) const;

 void remove_initial_level_objective( void );

 [[nodiscard]] bool uses_pure_level_aggregation( void ) const {
  return( StblType == kLevel && ! has_initial_level_objective() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current optimal value of the global multiplier lambda
 /** lambda is the master-side non-negative dual multiplier paired with the
  * model-value equation of the lower model (
  * stationarity (i): lambda + r - omega = 1), shared across every hard
  * component (cf. PolyhedralFunctionBlock::set_lambda). Meaningful only
  * after solve_master(). */

 [[nodiscard]] double get_lambda( void ) const
  { return( IsPrimal ? 1.0 : Var_lambda.get_value() ); }

/*--------------------------------------------------------------------------*/
 /// returns the physical gamma multiplier of the k-th hard component
 /** The PFB stores gamma in the same locally-scaled units as its internal
  * normalization row. Multiplying it by get_v_scale() maps it back to the
  * physical mass where sum_i theta_i + gamma = lambda. */

 [[nodiscard]] double get_gamma( int k ) const;

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
  * linear term of the MP Objective. The new value is stored in #t_stab. In the
  * primal MP the abstract Objective update is deferred and batched with any
  * intervening centre / linear-part changes immediately before the next actual
  * solve; in the dual MP the quadratic coefficient is updated immediately.
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
  * no-op under the dual MP. The abstract primal Objective is synchronized
  * immediately before the next actual master solve, together with any pending
  * t / centre change. */

 void set_b( const std::vector< double > & b );

/*--------------------------------------------------------------------------*/
 /// set the lower bound LB^k of the k-th hard component
 /** Sets the constraint v^k >= LB^k of the primal MP (or, equivalently,
  * the linear coefficient on gamma^k of the dual one). \p k must lie in [0,
  * NoHardCmps). LB^k = -INF turns the bound off. */

 void set_LB( int k , double LB );

/*--------------------------------------------------------------------------*/
 /// install (\p on == true) or remove a fictitious model LB v^k >= 0
 /** Installs (\p on == true) or removes (\p on == false) a *fictitious*
  * lower bound  v^k >= 0  on the model value of the k-th hard component.
  *
  * It is meant for the transient state in which the bundle of component
  * k is empty: an empty PolyhedralFunction models F_k as the max over an
  * empty set of affine pieces, i.e. the improper constant -infinity, so
  * the master would be primal-unbounded / dual-infeasible. The fictitious
  * v^k >= 0 (the historical QP-bundle device) makes the empty component
  * "disappear": v^k is held at 0 and gamma^k becomes a free, zero-cost
  * multiplier that absorbs the per-PFB simplex mass lambda.
  *
  * The bound is expressed directly in the d-space / linearization-error
  * frame (no F_k(x_bar) translation, unlike a genuine native bound passed
  * to set_LB), because the model value v^k already lives in that frame.
  * The driver must remove it (\p on == false) as soon as the component
  * gets a real cut, restoring the gamma^k-fixed "no bound" state. \p k
  * must lie in [0, NoHardCmps). */

 void set_fictitious_LB( int k , bool on );

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

 std::vector< LagBFunction * > EasyCmps;  ///< the "easy" components

 std::vector< std::vector< Index > > EasyLocal2Global;
                    ///< map local index of easy components into MP global ones

 // - - - - - - - - -  pointers to the per-component sub-Blocks - - - - - - - -

 std::vector< Block * > EasyCmps_SB;  ///< sub-Blocks of the "easy" components

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

 std::vector< std::vector< FRowConstraint > > EasyLBFCns;
                                   ///< generated rows indexed by original
                                   ///< component k and row

 // - - - - - - - - - - -  static MP entities (primal form)  - - - - - - - - -

 std::vector< ColVariable > Var_d;
                                ///< d in translated primal form, absolute x
                                ///< in raw primal form (free, size NumVars)

 std::vector< BoxConstraint > Bounds_d;
                                ///< per-coordinate primal box on x (raw form)
                                ///< or on d (translated form)

 std::vector< ColVariable > Var_v_hard;
                                ///< the epigraph variables v^k
                                ///< per hard component

 std::vector< BoxConstraint > Bounds_v_hard;
                                ///< per-hard-component LB^k: v^k >= LB^k
                                ///< (rhs = +INF by default)

 FRowConstraint LevelCns;                ///< primal level constraint
                                         ///< b*d + sum_k v^k <= f_lev
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

 bool f_primal_objective_dirty = false;
                    ///< whether the primal objective must be synchronized
                    ///< before the next actual master solve

 double f_C = 0.0;  ///< additive shift used to express LB_xbar = LB + C
                    ///< and Lvl_xbar = Lvl + C, with C = -f^0( x_bar )
                    ///< - sum_H f^{k,H}( x_bar ); refreshed by set_C() at
                    ///< every stability-centre move

 double f_global_LB_xbar = - Inf< Function::FunctionValue >();
                    ///< translated global lower bound LB + f_C last
                    ///< installed by set_global_LB(); -INF means "no
                    ///< global LB". Read by set_fictitious_LB() to place
                    ///< the fictitious per-component bound strictly below
                    ///< the aggregate one (legacy OSIMPSolver device)

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

 std::vector< double > f_linear_part;
                    ///< gradient b of the linear 0-th component installed
                    ///< in the dual coupling rows by set_linear_part()

 std::vector< double > f_F_at_x_bar;
                    ///< per-hard-component cache of  F_k( x_bar ) ; sized
                    ///< NoHardCmps in configure(), refreshed by the driver
                    ///< via set_F_at_x_bar() at every stability-centre
                    ///< change. Used (or will be used) to derive Sigma in
                    ///< linearization-error form without forcing the
                    ///< driver to pre-promote each cut's alpha

 double f_xref_tol = 0.0;
                    ///< "lazy reference" tolerance (see set_xref_tol). 0 =
                    ///< x_ref tracks x_bar at every serious step (plain
                    ///< displacement, strict no-op); > 0 defers the per-cut
                    ///< g . ( x_bar - x_ref ) shift to the lin-z and re-aligns
                    ///< x_ref only when || x_bar - x_ref ||_inf exceeds it

 std::vector< double > f_x_ref;
                    ///< the lazy storage reference: the diagonal cut constants
                    ///< bake g . x_ref ( not g . x_bar ); the residual
                    ///< g . ( x_bar - x_ref ) rides in the Var_z lin-z. Equal
                    ///< to f_x_bar whenever f_xref_tol == 0

 // the reference vector against which the diagonal cut constants are baked:
 // x_ref under an active lazy reference ( f_xref_tol > 0, sized ), x_bar
 // otherwise. At f_xref_tol == 0 this is x_bar, so every cut-baking site
 // ( add_cut / modify_cut / modify_alpha / get_raw_aggregated_alpha )
 // reduces to the plain displacement frame
 const std::vector< double > & cut_ref() const {
  return ( ( f_xref_tol > 0.0 ) &&
           ( f_x_ref.size() == f_x_bar.size() ) ) ? f_x_ref : f_x_bar;
  }

 [[nodiscard]] double get_master_objective_value() const;

 std::vector< double > f_LB_raw;
                    ///< per-hard-component cache of the raw native lower
                    ///< bound LB^k installed via set_LB(); -INFshift means
                    ///< "no bound" (gamma^k stays fixed to 0). MPB feeds
                    ///< the dual master gamma * LB^k contribution in
                    ///< linearization-error form ( LB^k - F_k( x_bar ) )
                    ///< and refreshes it on every set_reference() so the
                    ///< driver does not have to re-call set_LB

 double f_lev;      ///< current value of the level f_lev (level / doubly only)

 int f_v2_form = 0; ///< master storage frame: 0 = translated/displacement
                    ///< (default, production; cuts in the lin-error frame,
                    ///< x_bar baked into the constants, 0 linear z term),
                    ///< 1 = raw/iterate form. The primal variable is absolute
                    ///< x; in the dual MP x_bar is in the explicit +x_bar^T R
                    ///< z term and cut constants are function-value-relative.
                    ///< See set_v2_form

 int HardCmpScaling = 0;
                    ///< bit-wise PFB scaling for hard components:
                    ///< bit 0 = local rows, bit 1 = global epigraph

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

 int easy_obj_idx = -1;
                    ///< index of the first displacement-form easy objective
                    ///< correction in the root DQuadFunction, or -1

 int level_model_obj_idx = -1;
                    ///< first v^k term in the primal one-shot level probe
                    ///< objective

 bool f_dual_level_probe_active = false;
                    ///< true while pure-level dual form is temporarily solved
                    ///< as a proximal seed before switching to eta/omega level

 char f_abs_rep = 0;
                    ///< built-stage bits for MPB's abstract representation

 std::vector< ColVariable * > EasyObjVars;
                    ///< unique easy inner variables receiving x_bar * g(u)

 std::vector< std::vector< std::pair< Index , double > > > EasyObjCoeffs;
                    ///< per easy variable, the (global j, g coefficient)
                    ///< contributions determining its objective coefficient

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 void generate_primal_abstract_variables();
                    ///< materialize the primal master variables and wire the
                    ///< hard-component PolyhedralFunctionBlock active variables

 void generate_primal_abstract_constraints();
                    ///< materialize the primal box, level row and PFB rows

 void generate_primal_objective();
                    ///< materialize the primal master objective and PFB
                    ///< objective pieces

 void generate_dual_abstract_variables();
                    ///< materialize the dual master variables and wire the
                    ///< hard-component PolyhedralFunctionBlock active variables

 void generate_dual_abstract_constraints();
                    ///< materialize the dual normalization/coupling rows and
                    ///< the hard-component PFB rows

 void generate_dual_objective();
                    ///< materialize the dual master objective and PFB
                    ///< objective pieces

 void refresh_primal_objective();
                    ///< emit one batched objective Modification from the
                    ///< current t, x_bar and linear-part state

 void refresh_primal_level_linear_part();
                    ///< refresh the b coefficients in the primal level row
                    ///< without changing its v^k terms

 void refresh_box_coordinate( Index j , DQuadFunction * dqf );
                    ///< synchronize one cached box coordinate with the
                    ///< generated primal/dual abstract representation

 static PolyhedralFunctionBlock *
 pfb_at( const std::vector< Block * > & HardCmps , int k , const char * fn );

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

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS MasterProblemMod ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class to describe physical modifications specific to MasterProblemBlock
/** MasterProblemMod describes changes in the physical representation of a
 * MasterProblemBlock. In the SMS++ terminology this is a physical
 * Modification: it is issued after the MasterProblemBlock state has already
 * changed, and therefore it does not concern the Block again
 * (concerns_Block() remains false, as in the base Modification class).
 *
 * The affected physical data are the MPB-level ingredients of the master
 * problem, such as the stability centre, the stabilization parameters and the
 * MPB-owned coefficients / bounds. Modifications of the registered
 * PolyhedralFunctionBlock sub-Blocks are not duplicated here: those sub-Blocks
 * issue their own PolyhedralFunctionMod objects. */

class MasterProblemMod : public Modification
{

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/

 /// type of physical change in the MasterProblemBlock
 /** This enum follows the convention of C05FunctionMod,
  * PolyhedralFunctionMod and LagBFunctionMod: the actual field is stored as
  * an int so that derived classes may extend the set of supported
  * modifications if more detailed payloads are needed later. */

 enum master_problem_mod_type {
  ReferenceChanged ,       ///< the complete cached reference has changed
  TChanged ,               ///< the proximal stabilization parameter t changed
  LevelChanged ,           ///< the level stabilization value changed
  BoxLowerChanged ,        ///< the lower side of the variable box changed
  BoxUpperChanged ,        ///< the upper side of the variable box changed
  LinearPartChanged ,      ///< the MPB-owned linear part b changed
  LowerBoundChanged ,      ///< a component or global lower bound changed
  MasterProblemModLastParam
  ///< first value available to derived classes
  };

/*---------------------------- CONSTRUCTOR ---------------------------------*/
 /// constructor: takes the affected MasterProblemBlock and the change type

 explicit MasterProblemMod( MasterProblemBlock * block , int type )
  : Modification() , f_Block( block ) , f_type( type ) {}

/*------------------------------ DESTRUCTOR --------------------------------*/

 ~MasterProblemMod() override = default;  ///< destructor: does nothing

/*---------------------------- ACCESSORS -----------------------------------*/
 /// returns the Block this Modification was originated from

 [[nodiscard]] Block * get_Block( void ) const override { return( f_Block ); }

 /// accessor to the type of physical change

 [[nodiscard]] int type( void ) const { return( f_type ); }

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream & output ) const override {
  output << "MasterProblemMod on MasterProblemBlock [" << f_Block
         << "]: type = " << f_type << std::endl;
  }

/*-------------------------- PROTECTED FIELDS ------------------------------*/

 MasterProblemBlock * f_Block;  ///< affected MasterProblemBlock

 int f_type;                    ///< type of physical change

/*--------------------------------------------------------------------------*/

 };  // end( class MasterProblemMod )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS MasterProblemParamMod -------------------------*/
/*--------------------------------------------------------------------------*/
/// physical Modification carrying a new scalar master parameter value
/** MasterProblemParamMod specializes MasterProblemMod for the inexpensive
 * scalar payloads associated with #MasterProblemMod::TChanged and
 * #MasterProblemMod::LevelChanged. Storing the value that was current when
 * the Modification was issued preserves that information even if the same
 * parameter changes again before an asynchronous Solver processes it. */

class MasterProblemParamMod : public MasterProblemMod
{

 public:

 /// constructor: takes the affected Block, parameter type and new value

 MasterProblemParamMod( MasterProblemBlock * block , int type ,
                        double new_value )
  : MasterProblemMod( block , type ) , f_new_value( new_value ) {
  if( type != TChanged && type != LevelChanged )
   throw( std::invalid_argument(
    "MasterProblemParamMod: only TChanged and LevelChanged are supported" ) );
  }

 ~MasterProblemParamMod() override = default;  ///< destructor: does nothing

 /// returns the parameter value installed when the Modification was issued

 [[nodiscard]] double new_value( void ) const { return( f_new_value ); }

 protected:

 /// print the MasterProblemParamMod

 void print( std::ostream & output ) const override {
  output << "MasterProblemParamMod on MasterProblemBlock [" << f_Block
         << "]: type = " << f_type << ", new value = " << f_new_value
         << std::endl;
  }

 double f_new_value;  ///< value installed by the corresponding setter

 };  // end( class MasterProblemParamMod )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS MasterProblemRngdMod --------------------------*/
/*--------------------------------------------------------------------------*/
/// physical master Modification localized to a contiguous range
/** Besides the affected range, MasterProblemRngdMod owns the new values that
 * were installed when the Modification was issued. Its type identifies the
 * affected master datum: lower bounds, upper bounds, or (when partial linear
 * updates are added) linear-part coefficients. Values are snapshots, not
 * deltas. */

class MasterProblemRngdMod : public MasterProblemMod
{

 public:

 using Range = Block::Range;
 using Values = std::vector< double >;

 /// constructor: takes ownership of the new values

 MasterProblemRngdMod( MasterProblemBlock * block , int type , Range range ,
                       Values && new_values )
  : MasterProblemMod( block , type ) , f_range( range ) ,
    f_new_values( std::move( new_values ) ) {
  if( f_range.second < f_range.first )
   throw( std::invalid_argument(
    "MasterProblemRngdMod: invalid range" ) );
  const auto size = f_range.second - f_range.first;
  if( f_new_values.size() != size )
   throw( std::invalid_argument(
    "MasterProblemRngdMod: values and range sizes do not match" ) );
  }

 ~MasterProblemRngdMod() override = default;

 /// returns the affected left-closed, right-open range

 [[nodiscard]] const Range & range( void ) const { return( f_range ); }

 /// returns the new values aligned with range()

 [[nodiscard]] const Values & new_values( void ) const
  { return( f_new_values ); }

 protected:

 void print( std::ostream & output ) const override {
  output << "MasterProblemRngdMod on MasterProblemBlock [" << f_Block
         << "]: type = " << f_type << ", range = [ " << f_range.first
         << " , " << f_range.second << " ), values = "
         << f_new_values.size() << std::endl;
  }

 Range f_range;       ///< affected coordinates at issue time

 Values f_new_values; ///< new values aligned with f_range

 };  // end( class MasterProblemRngdMod )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS MasterProblemSbstMod --------------------------*/
/*--------------------------------------------------------------------------*/
/// physical master Modification localized to an arbitrary subset
/** MasterProblemSbstMod is the subset counterpart of
 * MasterProblemRngdMod. It owns a strictly increasing subset and one vector
 * of new values positionally aligned with that subset. If the incoming subset
 * is not ordered, the constructor sorts it and reorders the values in exactly
 * the same way. */

class MasterProblemSbstMod : public MasterProblemMod
{

 public:

 using Index = Block::Index;
 using Subset = Block::Subset;
 using Values = std::vector< double >;

 /// constructor: takes ownership of the subset and new values

 MasterProblemSbstMod( MasterProblemBlock * block , int type ,
                       Subset && subset , Values && new_values ,
                       bool ordered = false )
  : MasterProblemMod( block , type ) , f_subset( std::move( subset ) ) ,
    f_new_values( std::move( new_values ) ) {
  if( f_new_values.size() != f_subset.size() )
   throw( std::invalid_argument(
    "MasterProblemSbstMod: values and subset sizes do not match" ) );

  if( ( ! ordered ) && f_subset.size() > 1 ) {
   std::vector< std::size_t > order( f_subset.size() );
   std::iota( order.begin() , order.end() , std::size_t( 0 ) );
   std::sort( order.begin() , order.end() ,
              [ this ]( auto i , auto j )
              { return( f_subset[ i ] < f_subset[ j ] ); } );

   Subset sorted_subset( f_subset.size() );
   Values sorted_values( f_new_values.size() );
   for( std::size_t i = 0 ; i < order.size() ; ++i ) {
    sorted_subset[ i ] = f_subset[ order[ i ] ];
    sorted_values[ i ] = f_new_values[ order[ i ] ];
    }
   f_subset = std::move( sorted_subset );
   f_new_values = std::move( sorted_values );
   }

  for( std::size_t i = 1 ; i < f_subset.size() ; ++i )
   if( f_subset[ i - 1 ] >= f_subset[ i ] )
    throw( std::invalid_argument(
     "MasterProblemSbstMod: unordered or repeated subset" ) );
  }

 ~MasterProblemSbstMod() override = default;

 /// returns the strictly increasing affected subset

 [[nodiscard]] const Subset & subset( void ) const { return( f_subset ); }

 /// returns the new values aligned with subset()

 [[nodiscard]] const Values & new_values( void ) const
  { return( f_new_values ); }


 protected:
 void print( std::ostream & output ) const override {
  output << "MasterProblemSbstMod on MasterProblemBlock [" << f_Block
         << "]: type = " << f_type << ", subset size = "
         << f_subset.size() << ", values = "
         << f_new_values.size() << std::endl;
  }

 Subset f_subset;     ///< affected coordinates at issue time

 Values f_new_values; ///< new values aligned with f_subset

 };  // end( class MasterProblemSbstMod )

/*--------------------------------------------------------------------------*/
/*--------------- CLASS MasterProblemLowerBoundMod ------------------------*/
/*--------------------------------------------------------------------------*/
/// physical Modification carrying a new master lower-bound value
/** The target is either one hard component or the aggregate global lower
 * bound. Both identifiers and values are small enough to preserve each
 * individual change for an asynchronously processing Solver. */

class MasterProblemLowerBoundMod : public MasterProblemMod
{

 public:

 /// sentinel identifying the aggregate global lower bound

 static constexpr int GlobalLowerBound = -1;

 /// constructor: component index (or GlobalLowerBound) and new value

 MasterProblemLowerBoundMod( MasterProblemBlock * block , int component ,
                             double new_value )
  : MasterProblemMod( block , LowerBoundChanged ) ,
    f_component( component ) , f_new_value( new_value ) {
  if( component < GlobalLowerBound )
   throw( std::invalid_argument(
    "MasterProblemLowerBoundMod: invalid lower-bound target" ) );
  }

 ~MasterProblemLowerBoundMod() override = default;

 /// returns the hard-component index, or GlobalLowerBound

 [[nodiscard]] int component( void ) const { return( f_component ); }

 /// returns true when this Modification concerns the global lower bound

 [[nodiscard]] bool is_global( void ) const
  { return( f_component == GlobalLowerBound ); }

 /// returns the lower-bound value installed by the setter

 [[nodiscard]] double new_value( void ) const { return( f_new_value ); }

 protected:

 void print( std::ostream & output ) const override {
  output << "MasterProblemLowerBoundMod on MasterProblemBlock [" << f_Block
         << "]: target = " << f_component
         << " (-1 = global), new value = " << f_new_value << std::endl;
  }

 int f_component;       ///< hard-component index, or GlobalLowerBound

 double f_new_value;    ///< value installed by the corresponding setter

 };  // end( class MasterProblemLowerBoundMod )

/** @} end( group( MasterProblemBlock_CLASSES ) ) ----------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* __MasterProblemBlock */

/*--------------------------------------------------------------------------*/
/*--------------------- End File MasterProblemBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
