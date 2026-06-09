/*--------------------------------------------------------------------------*/
/*------------------- File MasterProblemBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MasterProblemBlock class.
 *
 * MasterProblemBlock generates, on demand, either the primal or the dual
 * form of the Bundle master problem as an SMS++ Block tree, so that any
 * [MILP]Solver can be attached on top of it through a standard
 * BlockSolverConfig. See CreatePrimalMP() and CreateDualMP() for the
 * actual construction of the static Variable, Constraint and Objective.
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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MasterProblemBlock.h"

#include <cstdio>
#include <iostream>

#include "BendersBFunction.h"
#include "BlockSolverConfig.h"
#include "Configuration.h"
#include "DQuadFunction.h"
#include "FRealObjective.h"
#include "LagBFunction.h"
#include "LinearFunction.h"
#include "MILPSolver.h"
#include "Modification.h"
#include "PolyhedralFunctionBlock.h"
#include "QuadFunction.h"
#include "Solver.h"

#include <algorithm>
#include <functional>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY REGISTRATION ---------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( MasterProblemBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- CLEAR / REINITIALIZE -----------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::clear()
{
 // forget any per-component sub-Block (the sub-Block objects themselves
 // are owned by the base Block, which will dispose of them in due time)
 EasyCmps.clear();
 HardCmps.clear();

 // absorbed BendersBFunction RowConstraints live as static_constraint
 // on *this* and are owned by the base Block; just drop the bookkeeping
 EasyBBFRows.clear();
 EasyLBFCns.clear();

 // reset every size / structural field
 MaxBSize   = 0;
 MaxSGLen   = 0;
 NumVars    = 0;
 NoTotCmps  = 0;
 NoEasyCmps = 0;
 NoHardCmps = 0;
 DoEasy     = 0;
 IsEasyCmp.clear();

 // back to the "default" MP type
 IsPrimal      = false;
 StblType      = kDoublyStabilized;
 HardCmpScaling = 0;
 t_stab        = 1.0;
 f_lev         = 0.0;
 z_obj_idx     = -1;
 r_obj_idx     = -1;
 omega_obj_idx = -1;

 // drop the dynamically-sized variable / constraint groups (primal d /
 // v^k / Bounds_v_hard, dual z / CouplingCns); the scalar members
 // (Var_lambda, Var_r, Var_omega, NormalizationCns, LevelCns) keep
 // their default state and are released by their own destructors when
 // *this is destroyed
 Var_d.clear();
 Var_v_hard.clear();
 Bounds_v_hard.clear();
 Var_z.clear();
 Var_lambda.set_value( 0.0 );
 Var_lambda.is_fixed( false , eNoMod );
 Var_s_plus.clear();
 Var_s_minus.clear();
 CouplingCns.clear();
 slot_to_local.clear();

 }  // end( MasterProblemBlock::clear )

/*--------------------------------------------------------------------------*/
/*-------------------------- DIMENSION SETUP -------------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::SetDim( int MxBSz , int NVars ,
                                 int NrFi , int NrFiEasy )
{
 if( ( MxBSz < 0 ) || ( NVars < 0 ) || ( NrFi < 0 ) || ( NrFiEasy < 0 ) ||
     ( NrFiEasy > NrFi ) )
  throw( std::invalid_argument(
       "MasterProblemBlock::SetDim: invalid dimensions" ) );

 // a SetDim() call always starts from a clean slate, since the structure of
 // the master problem may change from one call to the next (different
 // number of components, different number of variables, ...). The actual
 // building of the abstract representation is delegated to a subsequent
 // call to CreateEmptyMP().
 clear();

 MaxBSize   = MxBSz;
 MaxSGLen   = NVars;
 NumVars    = NVars;
 NoTotCmps  = NrFi;
 NoEasyCmps = NrFiEasy;
 NoHardCmps = NrFi - NrFiEasy;

 // by default no component is "easy"; the actual map is established by
 // CreateEmptyMP()
 IsEasyCmp.assign( NoTotCmps , false );

 // reserve the slots for the per-component sub-Blocks; the actual
 // PolyhedralFunctionBlock / easy-cmp Block objects are allocated by
 // CreateEmptyMP() once the formulation is known. The per-hard-cmp LB
 // multipliers gamma^k live inside each PolyhedralFunctionBlock sub-Block
 // (its own f_gamma) and are therefore *not* materialized here.
 EasyCmps.reserve( NoEasyCmps );
 HardCmps.reserve( NoHardCmps );

 // pre-allocate the slot -> local-row mapping: MaxBSize slots per hard
 // component, all initially empty (encoded as -1)
 slot_to_local.assign( NoHardCmps , std::vector< int >( MaxBSize , -1 ) );

 // per-hard-component cache of F_k( x_bar ); default 0 until the
 // driver feeds the actual reference values via set_F_at_x_bar
 f_F_at_x_bar.assign( NoHardCmps , 0.0 );

 // gradient of the linear 0-th component; default 0 until the driver
 // installs the actual coefficients via set_linear_part()
 f_linear_part.assign( NumVars , 0.0 );

 // The initial master is built around the zero stability centre. Keeping
 // the correctly-sized vector from the start also lets the first
 // set_reference() translate cuts that were inserted before F(x_bar) became
 // known, instead of mistaking that update for an uninitialised reference.
 f_x_bar.assign( NumVars , 0.0 );

 // per-hard-component cache of the raw native lower bound LB^k;
 // -INFshift = "no bound" by default until set_LB() is called
 f_LB_raw.assign( NoHardCmps , - Inf< double >() );

 }  // end( MasterProblemBlock::SetDim )

/*--------------------------------------------------------------------------*/
/*--------------------------- SOLVER REGISTRATION --------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::register_Solver( std::string && solv_cfg_filename )
{
 // no default backend is hard-wired into MasterProblemBlock: the choice of
 // the actual [MILP]Solver must always be expressed by the caller through
 // a BlockSolverConfig and resolved by the SMS++ Solver factory
 if( solv_cfg_filename.empty() )
  throw( std::invalid_argument(
       "MasterProblemBlock::register_Solver: empty configuration filename; a "
       "BlockSolverConfig is required to attach a Solver to the MasterProblemBlock" ) );

 auto cfg = Configuration::deserialize( solv_cfg_filename );
 auto MPBSC = dynamic_cast< BlockSolverConfig * >( cfg );
 if( ! MPBSC ) {
  delete cfg;
  throw( std::invalid_argument(
       "MasterProblemBlock::register_Solver: the provided Configuration "
       "file is not a BlockSolverConfig" ) );
  }

 // forward the exclusion list (if any) as the second argument of apply():
 // BlockSolverConfig::apply will install it on the freshly-created Solver
 // via Solver::set_excluded_blocks() BEFORE registering it to *this*, so
 // load_problem() already sees the right exclusion set
 const std::unordered_set< Block * > * ignored = f_ignored_blocks.empty()
                                                 ? nullptr : & f_ignored_blocks;
 MPBSC->apply( this , ignored );
 MPBSC->clear();
 delete MPBSC;

 }  // end( MasterProblemBlock::register_Solver )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::register_Solver(
                          std::string && solv_cfg_filename ,
                          std::unordered_set< Block * > && ignored_blocks )
{
 // remember the exclusion list so that register_Solver() and configure()
 // can be called in any order: whichever lands last picks up the list
 // and forwards it to the inner Solver via BlockSolverConfig::apply()
 f_ignored_blocks = std::move( ignored_blocks );
 register_Solver( std::move( solv_cfg_filename ) );

 }  // end( MasterProblemBlock::register_Solver( + ignored_blocks ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CONFIGURE ----------------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::configure(
                          bool primal ,
                          int max_bundle_size ,
                          int num_vars ,
                          int num_hard_cmps ,
                          const std::vector< C05Function * > & easy_components ,
                          Block * original_block ,
                          std::unordered_set< Block * > ignored_blocks ,
                          stabilization_type reg ,
                          bool convex ,
                          const std::vector< bool > & is_easy ,
                          int hard_cmp_scaling )
{
 // - - - sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( max_bundle_size < 0 || num_vars < 0 || num_hard_cmps < 0 )
  throw( std::invalid_argument(
       "MasterProblemBlock::configure: negative size" ) );
 if( hard_cmp_scaling < 0 || hard_cmp_scaling > 3 )
  throw( std::invalid_argument(
       "MasterProblemBlock::configure: hard_cmp_scaling must be between 0 and 3" ) );

 const int n_easy  = int( easy_components.size() );
 const int n_total = num_hard_cmps + n_easy;

 // - - - sizes / form / stabilisation - - - - - - - - - - - - - - - - - - -
 // SetDim() starts from a clean slate (it calls clear()) and re-initialises
 // MaxBSize / NumVars / NoTotCmps / NoEasyCmps / NoHardCmps
 SetDim( max_bundle_size , num_vars , n_total , n_easy );
 EasyLBFCns.resize( n_total );

 // record the per-global-k easy/hard flag so the per-cmp getters
 // (get_FiBLambda, get_aggregated_alpha, ...) can dispatch the EASY
 // branch. SetDim() above already initialised IsEasyCmp to all-false;
 // we overwrite it here when the caller supplies a flag vector of the
 // right size. The default empty vector preserves the legacy all-hard
 // behaviour for callers that have not yet adopted the new parameter
 if( ! is_easy.empty() && int( is_easy.size() ) == n_total )
  IsEasyCmp = is_easy;

 IsPrimal = primal;
 IsConvex = convex;
 StblType = reg;
 HardCmpScaling = hard_cmp_scaling;

 // - - - steal original_block into the master - - - - - - - - - - - - - -
 // detach the model Block from its previous parent and reattach it under
 // *this*. The Variables / Constraints living there will be seen by the
 // inner Solver as part of the master, except for sub-Block subtrees
 // listed in ignored_blocks (which the inner Solver will skip via
 // Solver::set_excluded_blocks(), installed by BlockSolverConfig::apply
 // inside register_Solver)
 f_original_block = original_block;
 if( f_original_block )
  f_original_block->transfer_ownership_to( this );

 // - - - dispatch each easy component into the master - - - - - - - - - - -
 // The MP-side embedding of an easy component is asymmetric on two axes:
 //
 //   axis 1: kind of C05Function (LagBFunction vs BendersBFunction);
 //   axis 2: form of the master problem (primal vs dual).
 //
 // Only two of the four combinations have a "natural" closed-form
 // embedding that does not require dualizing a Block:
 //
 //  - LagBFunction in the *dual* MP: take the inner Block (the primal of
 //    the dualized problem) and the linear map A coupling its variables
 //    to the LagBFunction's active y. The coupling is absorbed into the
 //    z-row equations of the dual MP via set_conjugate_constraint().
 //
 //  - BendersBFunction in the *primal* MP: the inner Block carries the
 //    fixed-rhs constraints  g_i( y ) [<=, =, >=] bar{d}_i. The embedding
 //    requires relaxing those constraints in the inner Block and
 //    re-installing them in the master as  g_i( y ) [<=, =, >=] ( A x +
 //    b )_i, where x are the active variables of the BendersBFunction
 //    and A, b come from get_A() / get_b(). This step requires casting
 //    g_i to a concrete Function type (LinearFunction, DQuadFunction,
 //    QuadFunction) and so is type-specific.
 //
 // The other two combinations (LagBFunction in primal, BendersBFunction
 // in dual) would require an explicit Block dualization, which is not
 // currently available in SMS++; they are reported as "unsupported".
 // Any C05Function type other than LagBFunction / BendersBFunction (for
 // instance a PolyhedralFunction, which is a hard component by design)
 // is also reported as "unsupported": such a component should never be
 // marked easy.
 int next_global_easy = 0;
 for( int k = 0 ; k < n_easy ; ++k ) {
  int component = k;
  if( int( is_easy.size() ) == n_total ) {
   while( next_global_easy < n_total && ! is_easy[ next_global_easy ] )
    ++next_global_easy;
   if( next_global_easy >= n_total )
    throw( std::invalid_argument(
         "MasterProblemBlock::configure: inconsistent easy-component flags" ) );
   component = next_global_easy++;
   }

  auto * c05 = easy_components[ k ];
  if( ! c05 )
   throw( std::invalid_argument(
        "MasterProblemBlock::configure: easy component " +
        std::to_string( k ) + " has a null C05Function" ) );

  if( auto * lbf = dynamic_cast< LagBFunction * >( c05 ) ) {
   if( primal )
    throw( std::invalid_argument(
         "MasterProblemBlock::configure: easy component " +
         std::to_string( k ) + " is a LagBFunction in the primal MP; "
         "this requires explicit dualization of the inner Block, which "
         "is not currently available in SMS++" ) );
   auto * inner = lbf->get_inner_block();
   if( ! inner )
    throw( std::invalid_argument(
         "MasterProblemBlock::configure: easy component " +
         std::to_string( k ) + " is a LagBFunction with no inner Block" ) );
   // generate the per-row stationarity constraints
   //     E^k_i u^k + lambda * e^k_i = 0
   // BEFORE moving the inner Block under *this*:
   // absorb_LBF_into_dual_MP() reaches into the inner Block via
   // lbf->get_inner_block(), which still resolves to `inner` here but
   // would return nullptr after transfer_ownership_to
   absorb_LBF_into_dual_MP( lbf , component );
   inner->transfer_ownership_to( this );
   EasyCmps.push_back( inner );
   continue;
   }

  if( auto * bbf = dynamic_cast< BendersBFunction * >( c05 ) ) {
   if( ! primal )
    throw( std::invalid_argument(
         "MasterProblemBlock::configure: easy component " +
         std::to_string( k ) + " is a BendersBFunction in the dual MP; "
         "this requires explicit dualization of the inner Block, which "
         "is not currently available in SMS++" ) );
   auto * inner = bbf->get_inner_block();
   if( ! inner )
    throw( std::invalid_argument(
         "MasterProblemBlock::configure: easy component " +
         std::to_string( k ) +
         " is a BendersBFunction with no inner Block" ) );
   absorb_BBF_into_primal_MP( bbf );
   inner->transfer_ownership_to( this );
   EasyCmps.push_back( inner );
   continue;
   }

  throw( std::invalid_argument(
       "MasterProblemBlock::configure: easy component " +
       std::to_string( k ) + " has an unsupported C05Function type; "
       "supported are LagBFunction (in the dual MP) and "
       "BendersBFunction (in the primal MP)" ) );
  }

 for( int k = 0 ; k < n_total ; ++k )
  if( ! EasyLBFCns[ k ].empty() )
   add_static_constraint( EasyLBFCns[ k ] ,
                          "MPB_LBF_easy_" + std::to_string( k ) );

 // - - - stash ignored sub-Blocks until a Solver is registered - - - - - -
 // configure() and register_Solver() can be called in any order; whichever
 // lands last forwards the set to the inner Solver via
 // Solver::set_excluded_blocks(). If a Solver is already attached, we
 // refresh its exclusion list directly here (any subsequent reload of the
 // model picks it up)
 f_ignored_blocks = std::move( ignored_blocks );
 if( ! f_ignored_blocks.empty() && ! get_registered_solvers().empty() )
  get_registered_solvers().front()->set_excluded_blocks( & f_ignored_blocks );

 // - - - build the static MP - - - - - - - - - - - - - - - - - - - - - - -
 // CreatePrimalMP / CreateDualMP populate the static abstract
 // representation of the master (the variables d / z / r / omega / v^k,
 // the coupling rows and the master Objective) according to the chosen
 // form; the per-hard-component PolyhedralFunctionBlock sub-Blocks are
 // allocated here and the surrounding driver then
 // feeds the linearizations into them via the Modification interface
 if( IsPrimal )
  CreatePrimalMP( StblType );
 else
  CreateDualMP( StblType );

 // - - - alert attached Solvers that the Block has been (re)built - - - - -
 // configure() is typically called *after* register_Solver(), which means
 // the inner :MILPSolver was attached to an empty Block and load_problem()
 // saw 0 variables / 0 rows. add_static_variable / add_static_constraint
 // do NOT emit any Modification, so the only way to force a full reload
 // from the inner Solver is the "nuclear option": an NBModification on
 // this Block, which process_modifications() will translate into a fresh
 // load_problem() call
 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( MasterProblemBlock::configure )

/*--------------------------------------------------------------------------*/
/*-------------------------- CREATE EMPTY MP -------------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::CreateEmptyMP( stabilization_type Stbl , int NoCmps ,
                                        int DoEasyCmp , int NoEasy ,
                                        std::vector< bool > IsEasy )
{
 // SetDim() must have been called first with consistent values
 if( ( NoCmps != NoTotCmps ) || ( NoEasy != NoEasyCmps ) ||
     ( int( IsEasy.size() ) != NoCmps ) )
  throw( std::logic_error(
       "MasterProblemBlock::CreateEmptyMP: dimensions inconsistent "
       "with the last SetDim() call" ) );

 StblType  = Stbl;
 DoEasy    = DoEasyCmp;
 IsEasyCmp = std::move( IsEasy );

 // the dual form is the only viable choice as soon as there is at least
 // one "easy" component to be inserted as-is; otherwise the primal form
 // is preferred since it directly minimizes on the step d, which is the
 // natural variable space of a Bundle method
 IsPrimal = ( NoEasyCmps == 0 );

 if( IsPrimal )
  CreatePrimalMP( StblType );
 else
  CreateDualMP( StblType );

 }  // end( MasterProblemBlock::CreateEmptyMP )

/*--------------------------------------------------------------------------*/
/*------------------------- PRIMAL / DUAL MP -------------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::CreatePrimalMP( stabilization_type Stbl )
{
 // kNone, kTrustRegion and kUpperLower are reserved enumerators of
 // stabilization_type; their full wiring is tracked separately and
 // until then callers must fall back to kProximal / kLevel /
 // kDoublyStabilized
 if( Stbl == kNone || Stbl == kTrustRegion || Stbl == kUpperLower )
  throw( std::logic_error(
       "MasterProblemBlock::CreatePrimalMP: stabilization type " +
       std::to_string( int( Stbl ) ) +
       " is reserved but not yet implemented" ) );

 StblType = Stbl;
 IsPrimal = true;

 // ---- static Variable: the step d ---------------------------------------
 //
 // The epigraph variable v^k of each hard component is NOT a static
 // master-side ColVariable here: it is the f_v of the per-component
 // PolyhedralFunctionBlock allocated below in the linearized-primal
 // representation (is_linearized() == true). The PFB also owns:
 //  - its f_bcv box constraint on f_v (loaded with the PolyhedralFunction
 //    global LB / UB)
 //  - its f_const dynamic FRowConstraint group, holding the
 //    v_k >= a_i^k . d + b_i^k cuts as they are pushed by
 //    PolyhedralFunction::add_row (i.e. by MasterProblemBlock::add_cut
 //    on this side); the cuts share the master-side Var_d as the
 //    "x" variables of the PolyhedralFunction.
 // The master-side Objective below then references the PFB f_v pointers
 // directly via get_v(), so no v_k coupling row is needed

 Var_d.clear();
 Var_d.resize( NumVars );          // x (raw form) or d (translated form)
 if( NumVars > 0 )
  add_static_variable( Var_d , f_v2_form ? "MPB_x" : "MPB_d" );

 Bounds_d.clear();
 Bounds_d.resize( NumVars );
 for( int j = 0 ; j < NumVars ; ++j ) {
  Bounds_d[ j ].set_variable( & Var_d[ j ] , eNoMod );
  Bounds_d[ j ].set_lhs( - Inf< double >() , eNoMod );
  Bounds_d[ j ].set_rhs( Inf< double >() , eNoMod );
  }
 if( NumVars > 0 )
  add_static_constraint( Bounds_d , "MPB_box" );

 // ---- one PolyhedralFunctionBlock sub-Block per "hard" component ---------
 //
 // The PFB is wired in its *linearized primal* representation (rep == 1).
 // Each PFB's variables are bound to Var_d. In raw form these are the
 // absolute x variables; in translated form they are the step d.
 // We do NOT call set_lambda() here: that is dual-form specific
 // and would introduce an extra slack in the bundle's normalization; the
 // linearized-primal rep does not have a lambda at all (no dual
 // normalization row).
 //
 // The PolyhedralFunction interior bundle is empty here: the
 // the driver feeds rows (g, alpha) into each f_polyf via
 // the Modification interface as new linearizations are produced

 HardCmps.clear();
 HardCmps.reserve( NoHardCmps );
 // PFB stvv bits 2 and 3 enable local and global scaling, respectively.
 // Keep bit 0 set so this remains the linearized-primal representation:
 // no scaling / local / global / both map to 1 / 5 / 9 / 13.
 const int scaling_cfg = ( ( HardCmpScaling & 1 ) ? 4 : 0 ) |
                         ( ( HardCmpScaling & 2 ) ? 8 : 0 );
 const SimpleConfiguration< int > rep_lin_primal( 1 | scaling_cfg );

 for( int k = 0 ; k < NoHardCmps ; ++k ) {
  auto * pfb = new PolyhedralFunctionBlock( this );

  // wire the "x" variables of f_polyf to the master-side Var_d, so every
  // cut v_k >= a . d + b pushed via add_row lands on the right columns
  PolyhedralFunction::VarVector vv;
  vv.reserve( NumVars );
  for( auto & di : Var_d )
   vv.push_back( & di );
  pfb->get_PolyhedralFunction().set_variables( std::move( vv ) );

  // ensure f_polyf has NO global lower bound (for the default convex
  // case, f_bound defaults to +INF which would translate the linearized
  // primal box-constraint f_bcv to "+INF <= f_v <= +INF" -> infeasible
  // master; we want f_v to be free). modify_bound is called with eNoMod
  // because we are still in setup phase and no Solver is attached yet
  pfb->get_PolyhedralFunction().modify_bound(
                  - Inf< Function::FunctionValue >() , eNoMod );

  pfb->generate_abstract_variables(
                       const_cast< SimpleConfiguration< int > * >( & rep_lin_primal ) );
  pfb->generate_abstract_constraints();
  pfb->generate_objective();

  HardCmps.push_back( pfb );
  add_nested_Block( pfb );
  }

 // ---- level constraint sum_k v^k <= f_lev (kLevel / kDoublyStabilized) ---
 //
 // The terms reference each PFB's f_v directly via get_v(); LevelCns is a
 // master-side static FRowConstraint so the inner Solver picks it up

 if( ( Stbl == kLevel || Stbl == kDoublyStabilized ) && NoHardCmps > 0 ) {
  LinearFunction::v_coeff_pair lvl_terms;
  lvl_terms.reserve( NoHardCmps );
  for( int k = 0 ; k < NoHardCmps ; ++k ) {
   auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
   if( pfb )
    lvl_terms.emplace_back( pfb->get_v() , 1.0 );
   }

  LevelCns.set_lhs( - Inf< double >() , eNoMod );
  LevelCns.set_rhs( f_lev , eNoMod );
  LevelCns.set_function(
     new LinearFunction( std::move( lvl_terms ) ) , eNoMod );
  add_static_constraint( LevelCns , "MPB_level" );
  }

 // ---- Objective ---------------------------------------------------------
 //
 // translated: min b*d + sum_k f_v[k] + (1/(2t)) ||d||^2
 // raw:        min b*x + sum_k f_v[k] + (1/(2t)) ||x-x_bar||^2
 //
 // The proximal quadratic term (1/(2t))||d||^2_2 is present only under
 // kProximal / kDoublyStabilized; under pure kLevel the master is an LP
 // and the d_i contribute only through the linear b*d part.
 //
 // The linear coefficient on d (the "b" of the paper, i.e. the constant
 // gradient of the linear part of the original sum-function) is left at 0
 // here; the driver is expected to install it through the dedicated
 // set_b() API.
 //
 // Each v^k contributes its f_v with linear coefficient +1 (the
 // PolyhedralFunctionBlock's own Objective also references f_v with
 // coefficient +1 in linearized-primal rep, but that lives in the sub-Block
 // and is summed in by the SMS++ scan of v_BFS, so we must NOT include
 // f_v[k] in the master's own Objective to avoid double-counting). We
 // therefore only build the d-side quadratic / linear terms here

 const bool has_quad =
  ( Stbl == kProximal ) || ( Stbl == kDoublyStabilized );

 DQuadFunction::v_coeff_triple triples;
 triples.reserve( NumVars );

 const double quad_coeff = has_quad ? 1.0 / ( 2.0 * t_stab ) : 0.0;
 for( int i = 0 ; i < NumVars ; ++i ) {
  const double lin_coeff = ( f_v2_form && has_quad )
                           ? - f_x_bar[ i ] / t_stab : 0.0;
  triples.emplace_back( & Var_d[ i ] , lin_coeff , quad_coeff );
  }

 FRealObjective * obj;
 if( triples.empty() ) {
  obj = new FRealObjective( this , new LinearFunction() );
  }
 else {
  obj = new FRealObjective( this , new DQuadFunction( std::move( triples ) ) );
  }
 obj->set_sense( Objective::eMin , eNoMod );
 set_objective( obj , eNoMod );
 f_primal_objective_dirty = false;

 // bookkeeping fields: Bounds_v_hard / Var_v_hard are *unused* in this
 // refactor; we keep them clean so the dual MP code path (which still
 // owns them) stays unaffected
 Bounds_v_hard.clear();
 Var_v_hard.clear();

 }  // end( MasterProblemBlock::CreatePrimalMP )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::CreateDualMP( stabilization_type Stbl )
{
 // kTrustRegion and kUpperLower are reserved enumerators with no
 // implementation yet; reject them up front. kNone is handled below
 // by the natural flow (no proximal quadratic term, no level row).
 if( Stbl == kTrustRegion || Stbl == kUpperLower )
  throw( std::logic_error(
       "MasterProblemBlock::CreateDualMP: stabilization type " +
       std::to_string( int( Stbl ) ) +
       " is reserved but not yet implemented" ) );

 StblType = Stbl;
 IsPrimal = false;

 // ---- static dual multipliers: lambda (global), r, omega -----------------
 //
 // lambda is the *single* global non-negative multiplier paired with the
 // model-value equation
 //
 //     v = d b + sum_E pi^k e^k + sum_H v_x^k
 //
 // of the lower model  and the
 // stationarity condition (ii) at v_x^k: lambda = sum_i theta^k_i +
 // gamma^k for every k in H). Hence the same lambda enters the simplex
 // (= normalization) row of *every* hard component PFB sub-Block with
 // coefficient +1 (cf. PolyhedralFunctionBlock::set_lambda, called below
 // with a single shared pointer), so each per-PFB row reads
 //
 //     sum_i theta^k_i + gamma^k + lambda = 1.
 //
 // r is the non-negative multiplier of the global lower bound v >= LB;
 // omega is the non-negative multiplier of the level row Lvl >= v,
 // live only under #kLevel / #kDoublyStabilized. The master-side
 // stationarity at v then reads
 //
 //     lambda + r - omega = 1
 //
 // which becomes the single static "normalization" constraint installed
 // below; the per-PFB rows are owned by the PFB sub-Blocks themselves.

 Var_lambda.is_positive( true , eNoMod );

 // Var_r is the dual multiplier of the global LB row. It is structurally
 // present in every stabilization type, but only carries an objective
 // contribution when set_global_LB() is called with a finite LB; until
 // then we pin it to 0 so the master normalization
 //     K * lambda + r - omega = 1
 // collapses to  K * lambda = 1 ,  lambda = 1 / K , and the per-PFB
 //     sum_i theta^k_i + gamma^k = lambda
 // sum across the K hard components to  K * lambda = 1  - i.e. the
 // classic simplex  total mass = 1  of the bundle dual. Leaving r free
 // with a 0 objective coefficient would allow r to pick any value in
 // [ 0 , 1 ], pushing lambda below 1 / K and shrinking the simplex to
 // mass < 1 ; the master would still be bounded but Sigma would be a
 // scaled version of the classical aggregated linearization error, and
 // the bundle stop test would fire prematurely on the scaled value.
 // set_global_LB() flips is_fixed dynamically once a finite LB is
 // provided (and back to 0 otherwise)
 Var_r.is_positive( true , eNoMod );
 Var_r.set_value( 0 );
 Var_r.is_fixed( true , eNoMod );

 // Var_omega is the dual multiplier of the level row. Under #kProximal
 // it is structurally absent and stays fixed to 0; under #kLevel /
 // #kDoublyStabilized it is meaningful only once set_f_lev() installs
 // a finite level target, so it starts pinned to 0 as well and gets
 // unfixed by set_f_lev() the moment a finite f_lev is provided
 // (symmetrically to how set_global_LB() handles Var_r). Without this
 // pin, the K * lambda + r - omega = 1 normalization would let omega
 // drift to +infinity, pulling lambda along the direction
 // omega -> +infinity, lambda -> +infinity / K, with the same gamma^k
 // unboundedness consequence described above for Var_r
 Var_omega.is_positive( true , eNoMod );
 Var_omega.set_value( 0 );
 Var_omega.is_fixed( true , eNoMod );
 add_static_variable( Var_lambda , "MPB_lambda" );
 add_static_variable( Var_r      , "MPB_r"      );
 add_static_variable( Var_omega  , "MPB_omega"  );

 // ---- z auxiliary variables (free, one per coordinate) -------------------

 Var_z.clear();
 Var_z.resize( NumVars );
 if( NumVars > 0 )
  add_static_variable( Var_z , "MPB_z" );

 // ---- s^+ / s^- non-negative slack multipliers (box) ---------------------
 // One per coordinate; the slack is meaningful only when the matching
 // box side is finite (L_t - (x_bar)_t for s^+, U_t - (x_bar)_t for s^-).
 // Until the box is wired in we keep both fixed to 0 so the slack does
 // not contribute to the dual problem.
 Var_s_plus.clear();
 Var_s_minus.clear();
 Var_s_plus.resize( NumVars );
 Var_s_minus.resize( NumVars );
 for( int j = 0 ; j < NumVars ; ++j ) {
  Var_s_plus[ j ].is_positive( true , eNoMod );
  Var_s_plus[ j ].set_value( 0.0 );
  Var_s_plus[ j ].is_fixed( true , eNoMod );
  Var_s_minus[ j ].is_positive( true , eNoMod );
  Var_s_minus[ j ].set_value( 0.0 );
  Var_s_minus[ j ].is_fixed( true , eNoMod );
  }
 if( NumVars > 0 ) {
  add_static_variable( Var_s_plus  , "MPB_s_plus"  );
  add_static_variable( Var_s_minus , "MPB_s_minus" );
  }

 // ---- global normalization row: lambda + r - omega = 1 -----------------
 // The disaggregated proximal dual stationarity in v_k (the per-component
 // model decrease) reads
 //     dL/dv_k = 1 - sum_i theta^k_i - gamma^k - r + omega = 0
 // hence per component sum_i theta^k_i + gamma^k = 1 - r + omega =: lambda.
 // The per-PFB simplex row enforces sum_i theta^k_i + gamma^k = lambda for
 // every hard k, so this global row only has to pin lambda + r - omega = 1
 // (coefficient 1 on lambda, NOT NoHardCmps): the dual aggregates the
 // *sum* of the K components f = sum_k f_k, and each one carries unit
 // convexity mass independently. A K factor here would instead force
 // lambda = 1 / K when r = omega = 0, i.e. it would average the components
 // rather than sum them, shrinking the aggregate subgradient z = sum theta g
 // by 1 / K and the proximal step d = - t z with it.
 {
  LinearFunction::v_coeff_pair norm_terms;
  norm_terms.reserve( 3 );
  norm_terms.emplace_back( & Var_lambda , 1.0 );
  norm_terms.emplace_back( & Var_r      ,  1.0 );
  norm_terms.emplace_back( & Var_omega  , -1.0 );

  NormalizationCns.set_lhs( 1.0 , eNoMod );
  NormalizationCns.set_rhs( 1.0 , eNoMod );
  NormalizationCns.set_function(
     new LinearFunction( std::move( norm_terms ) ) , eNoMod );
  add_static_constraint( NormalizationCns , "MPB_norm" );
  }

 // ---- coupling rows ) -----------------
 //
 //   z_j - lambda * b_j + s^+_j - s^-_j - sum_{k in E} ( u^k A^k )_j
 //         - sum_{k in H, i in beta_SG^k} theta^k_i g^k_{i,j}
 //         - sum_{k in H, h in beta_F^k}  theta^k_h g^k_{h,j}        = 0
 //
 // The lambda*b_j coefficient is filled in by set_linear_part() once
 // the driver knows the linear part b of the original sum-function;
 // the theta-side terms are appended by PolyhedralFunctionBlock::set_-
 // conjugate_constraint() called below for every hard component; the
 // u^k A^k terms (easy components in the dual MP) live inside the
 // stolen LagBFunction inner Block. CouplingCns is exposed as a *dynamic*
 // group because PolyhedralFunctionBlock::set_conjugate_constraint takes
 // a std::list< FRowConstraint > & by reference; the list size itself
 // stays constant (NumVars) throughout the algorithm.
 //
 // The fixed entries on Var_lambda / Var_s_plus[ j ] / Var_s_minus[ j ]
 // are inserted *now* (with a 0 coefficient on Var_lambda, which
 // set_linear_part will later refresh to -b_j) so that the positions
 // 1..3 of every CouplingCns[ j ] LinearFunction are stable; subsequent
 // calls to PolyhedralFunctionBlock::set_conjugate_constraint append
 // their theta-terms in tail without disturbing them.

 CouplingCns.clear();
 CouplingCns.resize( NumVars );
 {
  auto it = CouplingCns.begin();
  for( int j = 0 ; j < NumVars ; ++j , ++it ) {
   LinearFunction::v_coeff_pair vp;
   vp.reserve( 4 );
   vp.emplace_back( & Var_z[ j ]        ,  1.0 );  // pos 0
   vp.emplace_back( & Var_lambda        ,  0.0 );  // pos 1 (set by set_linear_part)
   vp.emplace_back( & Var_s_plus[ j ]   ,  1.0 );  // pos 2
   vp.emplace_back( & Var_s_minus[ j ]  , -1.0 );  // pos 3
   it->set_function( new LinearFunction( std::move( vp ) , 0.0 ) , eNoMod );
   it->set_lhs( 0.0 , eNoMod );
   it->set_rhs( 0.0 , eNoMod );
   }
  }
 if( NumVars > 0 )
  add_dynamic_constraint( CouplingCns , "MPB_coupling" );

 // ---- one PolyhedralFunctionBlock sub-Block per "hard" component ---------
 //
 // The PFB is wired in its *linearised dual* representation (rep == 3):
 //
 //  - f_gamma >= 0 is the per-component LB^k multiplier;
 //  - f_theta (dynamic) is the list of theta^k_i bundle multipliers;
 //  - the per-PFB normalization row sum_i theta^k_i + gamma^k = 1 is
 //    later augmented with the master-side lambda via set_lambda(), so
 //    that it becomes sum_i theta^k_i + gamma^k + lambda = 1, which is
 //    just the paper's normalization sum_i theta^k_i + gamma^k = lambda
 //    re-grouped (the constant 1 on the right-hand side stays, lambda
 //    appears on the LHS with sign +1);
 //  - the per-PFB Objective contributes sum_i theta^k_i * b^k_i +
 //    gamma^k * LB^k (the latter is 0 unless f_polyf carries an explicit
 //    lower bound), which is precisely the per-component piece of the
 //    dual objective in (D);
 //  - the per-PFB contribution to the master-side z coupling rows is
 //    installed via set_conjugate_constraint(CouplingCns), which augments
 //    every CouplingCns[j] with the terms +theta^k_i * a^k_{i,j}.
 //
 // The PolyhedralFunction interior bundle is empty here: the
 // the driver feeds rows (g, alpha) into each f_polyf via
 // the Modification interface as new linearizations are produced.

 HardCmps.clear();
 HardCmps.reserve( NoHardCmps );
 // Keep bits 0 and 1 set so this remains the linearized-dual
 // representation: no scaling / local / global / both map to 3 / 7 / 11 / 15.
 const int scaling_cfg = ( ( HardCmpScaling & 1 ) ? 4 : 0 ) |
                         ( ( HardCmpScaling & 2 ) ? 8 : 0 );
 const SimpleConfiguration< int > rep_dual( 3 | scaling_cfg );

 for( int k = 0 ; k < NoHardCmps ; ++k ) {
  auto * pfb = new PolyhedralFunctionBlock( this );

  // f_polyf needs an "active variables" vector of size NumVars; in the
  // dual representation these only serve as keys for the bookkeeping of
  // set_conjugate_constraint (they are not the lambda of the original
  // PolyhedralFunction; the bundle is fed in dual space). Var_z fits the
  // role: one ColVariable* per coordinate, with NumVars in total.
  PolyhedralFunction::VarVector vv;
  vv.reserve( NumVars );
  for( auto & zj : Var_z )
   vv.push_back( & zj );
  pfb->get_PolyhedralFunction().set_variables( std::move( vv ) );

  // Align the PFB's PolyhedralFunction sense to IsConvex so that
  // generate_objective() emits the same eMin/eMax of the surrounding
  // master MP and of any easy sub-Block grafted under *this* by the
  // configure() dispatch. Specifically (cf. PolyhedralFunctionBlock.cpp:215
  // — `dual_min = !convex`):
  //  - C05Function convex  ⇒ master MP is min ⇒ PFB needs eMin
  //                          ⇒ PolyhedralFunction must be "concave"
  //  - C05Function concave ⇒ master MP is max ⇒ PFB needs eMax
  //                          ⇒ PolyhedralFunction stays "convex" (default)
  pfb->get_PolyhedralFunction().set_is_convex( ! IsConvex , eNoMod );

  pfb->generate_abstract_variables(
                              const_cast< SimpleConfiguration< int > * >( & rep_dual ) );
  pfb->generate_abstract_constraints();
  pfb->generate_objective();

  pfb->set_lambda( & Var_lambda );
  pfb->set_conjugate_constraint( CouplingCns );

  HardCmps.push_back( pfb );
  add_nested_Block( pfb );
  }

 // ---- master-side FRealObjective: -(t/2) || z ||^2_2 + omega * f_lev -----
 //
 // The bundle-summing terms theta^k_i b^k_i + gamma^k LB^k of every hard
 // component already live in the sub-PFB Objectives and are accumulated by
 // the SMS++ engine when the master is solved. Here we only need to add
 // the master-side stabilization terms:
 //
 //  - the quadratic stabilization  -(t/2) || z ||^2_2  is present under
 //    #kProximal and #kDoublyStabilized; it is *absent* under pure
 //    #kLevel, where the master is an LP (no proximal term);
 //
 //  - the level/X linear coefficient on omega (i.e. + omega * f_lev) is
 //    present under #kLevel / #kDoublyStabilized; under #kProximal omega
 //    is fixed to 0 and the term vanishes.
 //
 // The other linear part x_bar * z is left at 0 here: the driver
 // injects it through LinearFunction::modify_coefficient as the stability
 // centre changes.

 // Triple layout in the DQuadFunction (always the same regardless of
 // stabilization, so that set_x_bar / set_global_LB / set_f_lev can locate
 // their target coefficient by a fixed offset):
 //
 //   triples[ 0 .. NumVars - 1 ]   :  z_j with (linear = 0, quad = -t/2
 //                                   if has_quad else 0); set_x_bar moves
 //                                   the linear coefficient to x_bar[j]
 //   triples[ NumVars ]            :  r with (linear = 0, quad = 0);
 //                                   set_global_LB moves the linear
 //                                   coefficient to LB
 //   triples[ NumVars + 1 ]        :  omega with (linear = f_lev, quad = 0)
 //                                   IF has_omega_lin; set_f_lev refreshes
 //                                   the linear coefficient
 const bool has_quad =
  ( Stbl == kProximal ) || ( Stbl == kDoublyStabilized );
 const bool has_omega_lin =
  ( Stbl == kLevel ) || ( Stbl == kDoublyStabilized );

 // The dual MP is written in the same sense as the C05Function (= eMax
 // for concave/max, eMin for convex/min) so it composes with the easy
 // sub-Blocks grafted under *this* by configure() and with the hard PFB
 // sub-Blocks (set_is_convex above) without triggering the mixed max/min
 // check in [MILP]Solver::load_problem. The textbook form is
 //   max  Σ θ α − t/2 ‖z‖² + x̄·z − f_lev·ω + LB·r
 // when the C05Function is concave; the convex case is the negation of
 // every coefficient with set_sense flipped to eMin. The sign factor
 // sgn collapses both cases into a single expression below.
 const double sgn = IsConvex ? -1.0 : 1.0;

 DQuadFunction::v_coeff_triple triples;
 triples.reserve( NumVars + 1 + ( has_omega_lin ? 1 : 0 ) );

 z_obj_idx = NumVars > 0 ? 0 : -1;
 const double quad_coeff = has_quad ? - sgn * t_stab / 2.0 : 0.0;
 for( int j = 0 ; j < NumVars ; ++j )
  triples.emplace_back( & Var_z[ j ] , 0.0 , quad_coeff );

 r_obj_idx = int( triples.size() );
 triples.emplace_back( & Var_r , 0.0 , 0.0 );

 omega_obj_idx = -1;
 if( has_omega_lin ) {
  omega_obj_idx = int( triples.size() );
  // : the omega-side contribution to the dual objective is
  //     - omega * Lvl_xbar  ,   Lvl_xbar = f_lev + f_C
  // in the textbook eMax form; the convex case flips the whole row sign
  const double omega_lin = std::isfinite( f_lev )
                            ? - sgn * ( f_lev + f_C ) : 0.0;
  triples.emplace_back( & Var_omega , omega_lin , 0.0 );
  }

 // : the box slacks s^+ / s^- contribute
 //     + s^+ ( L - x_bar ) - s^- ( U - x_bar )
 // in the textbook eMax form. Coefficients start at 0 (= no box wired
 // yet, slacks fixed); set_box() / set_x_bar() update them later. The
 // s^+ / s^- triples are laid out contiguously after omega so that the
 // refresh logic can address them by a fixed base offset
 s_plus_obj_idx  = -1;
 s_minus_obj_idx = -1;
 if( NumVars > 0 ) {
  s_plus_obj_idx = int( triples.size() );
  for( int j = 0 ; j < NumVars ; ++j )
   triples.emplace_back( & Var_s_plus[ j ] , 0.0 , 0.0 );
  s_minus_obj_idx = int( triples.size() );
  for( int j = 0 ; j < NumVars ; ++j )
   triples.emplace_back( & Var_s_minus[ j ] , 0.0 , 0.0 );
  }

 auto obj = new FRealObjective( this ,
                                new DQuadFunction( std::move( triples ) ) );
 obj->set_sense( IsConvex ? Objective::eMin : Objective::eMax , eNoMod );
 set_objective( obj , eNoMod );

 // Easy components are not allocated here: they are registered one by one
 // by the surrounding driver via register_easy_component(),
 // which adds the easy-cmp sub-Block and augments every CouplingCns[j]
 // with the +A^k_{i,j} u^k_i terms produced by that component.

 // Two coefficients are intentionally left unset by CreateDualMP and must
 // be filled in by the surrounding driver as soon as
 // the corresponding pieces of the sum-function become known:
 //  - the omega * h coefficient on the level row, whenever X is a
 //    polyhedron with explicit right-hand side h;
 //  - the constant term b_j on every CouplingCns[j] rhs, coming from the
 //    linear part of the original sum-function.

 }  // end( MasterProblemBlock::CreateDualMP )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::absorb_BBF_into_primal_MP( BendersBFunction * bbf )
{
 if( ! bbf )
  throw( std::invalid_argument(
       "MasterProblemBlock::absorb_BBF_into_primal_MP: null BendersBFunction" ) );

 // sanity: the BendersBFunction must have as many active variables as
 // the primal MP step (one entry of Var_d per active x of the BBF), so
 // that the coupling -A_i . d makes sense
 if( int( bbf->get_num_active_var() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::absorb_BBF_into_primal_MP: BendersBFunction "
       "has " + std::to_string( bbf->get_num_active_var() ) +
       " active variables but the primal MP has " +
       std::to_string( NumVars ) + " step variables" ) );

 const auto & A = bbf->get_A();
 const auto & b = bbf->get_b();
 const auto & C = bbf->get_constraints();
 const auto & S = bbf->get_sides();

 const std::size_t m = A.size();
 if( b.size() != m || C.size() != m || S.size() != m )
  throw( std::logic_error(
       "MasterProblemBlock::absorb_BBF_into_primal_MP: inconsistent "
       "BendersBFunction mapping sizes" ) );

 EasyBBFRows.reserve( EasyBBFRows.size() + m );

 for( std::size_t i = 0 ; i < m ; ++i ) {
  RowConstraint * ci = C[ i ];
  if( ! ci )
   throw( std::invalid_argument(
        "MasterProblemBlock::absorb_BBF_into_primal_MP: null RowConstraint "
        "in BendersBFunction mapping at row " + std::to_string( i ) ) );

  auto * fci = dynamic_cast< FRowConstraint * >( ci );
  if( ! fci )
   throw( std::invalid_argument(
        "MasterProblemBlock::absorb_BBF_into_primal_MP: RowConstraint at "
        "row " + std::to_string( i ) + " is not a FRowConstraint" ) );

  auto * fun = fci->get_function();
  auto * qf  = dynamic_cast< QuadFunction   * >( fun );
  // QuadFunction derives from DQuadFunction, so test the most-derived
  // type first; DQuadFunction is the diagonal-only specialisation and
  // LinearFunction is the further linear specialisation
  auto * dqf = qf ? nullptr : dynamic_cast< DQuadFunction * >( fun );
  auto * lf  = ( qf || dqf ) ? nullptr
                             : dynamic_cast< LinearFunction * >( fun );
  if( ! lf && ! dqf && ! qf )
   throw( std::logic_error(
        "MasterProblemBlock::absorb_BBF_into_primal_MP: row " +
        std::to_string( i ) + " has a Function of unsupported type; "
        "currently only LinearFunction, DQuadFunction and QuadFunction "
        "are supported" ) );

  // 1. snapshot the original (lhs, rhs) of C_i and relax it on the inner
  //    Block: from this point on the constraint lives only on *this*
  const double orig_lhs = fci->get_lhs();
  const double orig_rhs = fci->get_rhs();
  fci->set_lhs( - Inf< double >() , eNoMod );
  fci->set_rhs(   Inf< double >() , eNoMod );

  // 2. build the master-side function g_i( y ) - A_i . d as a fresh
  //    Function of the same concrete type as the original, with a copy
  //    of the y-side coefficients plus the -A_i[ j ] coupling on every
  //    Var_d[ j ] (the coupling is *linear*; its quadratic contribution
  //    is zero by construction)
  Function * new_fun = nullptr;
  if( lf ) {
   LinearFunction::v_coeff_pair pairs;
   pairs.reserve( lf->get_num_active_var() + NumVars );
   for( Function::Index j = 0 ; j < lf->get_num_active_var() ; ++j )
    pairs.emplace_back(
         static_cast< ColVariable * >( lf->get_active_var( j ) ) ,
         lf->get_coefficient( j ) );
   for( int j = 0 ; j < NumVars ; ++j )
    pairs.emplace_back( & Var_d[ j ] , - A[ i ][ j ] );
   new_fun = new LinearFunction( std::move( pairs ) );
   }
  else if( dqf ) {
   DQuadFunction::v_coeff_triple triples;
   triples.reserve( dqf->get_num_active_var() + NumVars );
   for( Function::Index j = 0 ; j < dqf->get_num_active_var() ; ++j )
    triples.emplace_back(
         static_cast< ColVariable * >( dqf->get_active_var( j ) ) ,
         dqf->get_linear_coefficient( j ) ,
         dqf->get_quadratic_coefficient( j ) );
   for( int j = 0 ; j < NumVars ; ++j )
    triples.emplace_back( & Var_d[ j ] , - A[ i ][ j ] , 0.0 );
   new_fun = new DQuadFunction( std::move( triples ) );
   }
  else {  // QuadFunction
   // copy the diagonal triples ( y_j, linear, diag-quad ) and the
   // off-diagonal terms ( i_y, j_y, off-diag-quad ) verbatim; the
   // off-diagonal matrix indices remain valid because the appended
   // Var_d[ j ] are placed *after* the original y_j entries.
   DQuadFunction::v_coeff_triple triples;
   triples.reserve( qf->get_num_active_var() + NumVars );
   for( Function::Index j = 0 ; j < qf->get_num_active_var() ; ++j )
    triples.emplace_back(
         static_cast< ColVariable * >( qf->get_active_var( j ) ) ,
         qf->get_linear_coefficient( j ) ,
         qf->DQuadFunction::get_quadratic_coefficient( j ) );
   for( int j = 0 ; j < NumVars ; ++j )
    triples.emplace_back( & Var_d[ j ] , - A[ i ][ j ] , 0.0 );
   QuadFunction::v_off_diag_term off_diag;
   qf->get_v_nd_var( off_diag );
   new_fun = new QuadFunction( std::move( triples ) , std::move( off_diag ) );
   }

  // 3. allocate the master-side FRowConstraint with the cloned function;
  //    LHS/RHS are temporarily set assuming x_bar == 0, set_x_bar() will
  //    refresh them on every change of the stability centre
  auto * new_cns = new FRowConstraint();
  new_cns->set_function( new_fun , eNoMod );
  const double base_side = b[ i ];  // == A_i . 0 + b_i
  switch( S[ i ] ) {
   case( BendersBFunction::eLHS ):
    new_cns->set_lhs( base_side , eNoMod );
    new_cns->set_rhs( orig_rhs  , eNoMod );
    break;
   case( BendersBFunction::eRHS ):
    new_cns->set_lhs( orig_lhs  , eNoMod );
    new_cns->set_rhs( base_side , eNoMod );
    break;
   case( BendersBFunction::eBoth ):
    new_cns->set_lhs( base_side , eNoMod );
    new_cns->set_rhs( base_side , eNoMod );
    break;
   default:
    throw( std::invalid_argument(
         "MasterProblemBlock::absorb_BBF_into_primal_MP: unknown "
         "BendersBFunction::ConstraintSide at row " +
         std::to_string( i ) ) );
   }

  // 4. attach the new RowConstraint to *this* (as a static one, named
  //    after the absorbed component for diagnostic purposes)
  add_static_constraint( *new_cns , "MPB_BBF_easy" );

  // 5. record metadata for set_x_bar() to keep the absorbed sides in
  //    sync with the stability centre
  EasyBBFRow row;
  row.cns      = new_cns;
  row.A_row    = A[ i ];
  row.b_i      = b[ i ];
  row.side     = static_cast< char >( S[ i ] );
  row.orig_lhs = orig_lhs;
  row.orig_rhs = orig_rhs;
  EasyBBFRows.emplace_back( std::move( row ) );
  }
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::absorb_LBF_into_dual_MP( LagBFunction * lbf ,
                                                  int component )
{
 if( ! lbf )
  throw( std::invalid_argument(
       "MasterProblemBlock::absorb_LBF_into_dual_MP: null LagBFunction" ) );
 if( component < 0 || component >= int( EasyLBFCns.size() ) )
  throw( std::invalid_argument(
       "MasterProblemBlock::absorb_LBF_into_dual_MP: invalid component" ) );

 // walk every RowConstraint of the inner Block of the LagBFunction
 // (already transferred under *this* by configure()) and re-install it
 // on the master as the stationarity row of the dual MP at pi^k
 // ():
 //
 //     E^k_i u^k - lambda * e^k_i = 0    for every row i.
 //
 // The original RowConstraint of the inner Block is relaxed in place
 // (LHS = -INF, RHS = +INF) so that the inner :MILPSolver loaded by the
 // dual MP back-end does not enforce it twice.

 auto * inner = lbf->get_inner_block();
 if( ! inner )
  throw( std::invalid_argument(
       "MasterProblemBlock::absorb_LBF_into_dual_MP: null inner Block" ) );

 auto & easy_cns = EasyLBFCns[ component ];

 auto count_any = []( const boost::any & any ) -> std::size_t {
  if( boost::any_cast< FRowConstraint >( & any ) )
   return( 1 );
  if( auto one_p = boost::any_cast< FRowConstraint * >( & any ) )
   return( *one_p ? 1 : 0 );
  if( auto vec =
            boost::any_cast< std::vector< FRowConstraint > >( & any ) )
   return( vec->size() );
  if( auto vec_p =
            boost::any_cast< std::vector< FRowConstraint > * >( & any ) )
   return( *vec_p ? ( *vec_p )->size() : 0 );
  if( auto lst = boost::any_cast< std::list< FRowConstraint > >( & any ) )
   return( lst->size() );
  if( auto lst_p =
            boost::any_cast< std::list< FRowConstraint > * >( & any ) )
   return( *lst_p ? ( *lst_p )->size() : 0 );
  return( 0 );
  };

 std::size_t num_rows = 0;
 for( const auto & any : inner->get_static_constraints() )
  num_rows += count_any( any );
 for( const auto & any : inner->get_dynamic_constraints() )
  num_rows += count_any( any );
 easy_cns.reserve( num_rows );

 auto process_one = [ this , & easy_cns ]( FRowConstraint & ci ) {
  auto * lf = dynamic_cast< LinearFunction * >( ci.get_function() );
  if( ! lf )
   throw( std::logic_error(
        "MasterProblemBlock::absorb_LBF_into_dual_MP: only LinearFunction "
        "constraints are supported on the inner Block (the easy component's "
        "feasible set must be polyhedral)" ) );

  // pick e^k_i from the side the original constraint was using: an
  // equality constraint has lhs = rhs (= eBoth in BBF jargon); an
  // upper-bounded one stores e^k_i as the rhs; a lower-bounded one as
  // the lhs. Pure two-sided inequalities are ambiguous in ,
  // since stationarity wants a single e^k_i; we currently accept them
  // by taking the finite side (rhs if finite, else lhs).
  const double orig_lhs = ci.get_lhs();
  const double orig_rhs = ci.get_rhs();
  double e_i;
  if( orig_lhs == orig_rhs )                e_i = orig_rhs;
  else if( std::isfinite( orig_rhs ) )      e_i = orig_rhs;
  else if( std::isfinite( orig_lhs ) )      e_i = orig_lhs;
  else
   throw( std::invalid_argument(
        "MasterProblemBlock::absorb_LBF_into_dual_MP: inner RowConstraint "
        "is unbounded on both sides; cannot derive e^k_i" ) );

  // build the master-side LinearFunction: original ( u^k, E^k_{i,j} )
  // pairs plus the -e^k_i coupling on Var_lambda
  LinearFunction::v_coeff_pair pairs;
  pairs.reserve( lf->get_num_active_var() + 1 );
  for( Function::Index j = 0 ; j < lf->get_num_active_var() ; ++j )
   pairs.emplace_back(
        static_cast< ColVariable * >( lf->get_active_var( j ) ) ,
        lf->get_coefficient( j ) );
  pairs.emplace_back( & Var_lambda , - e_i );

  easy_cns.emplace_back();
  auto & new_cns = easy_cns.back();
  new_cns.set_function( new LinearFunction( std::move( pairs ) ) , eNoMod );
  new_cns.set_lhs( 0.0 , eNoMod );
  new_cns.set_rhs( 0.0 , eNoMod );

  // relax the original RowConstraint on the inner Block
  ci.set_lhs( - Inf< double >() , eNoMod );
  ci.set_rhs(   Inf< double >() , eNoMod );
  };

 // static / dynamic constraints are exposed as const vector<boost::any>
 // by Block::get_static_constraints / get_dynamic_constraints; we need
 // mutable access to relax LHS/RHS of the absorbed RowConstraints, so
 // we const_cast the outer container (the inner Block now lives under
 // *this*, so there is no aliasing concern with other observers).
 // Each entry holds either a single FRowConstraint, a vector of them
 // or a std::list of them (dynamic only). Blocks may register such
 // groups either *by value* (the Block owns the storage and boost::any
 // holds a copy/container directly) or *by pointer* (the Block keeps
 // the storage as a member and registers its address into boost::any);
 // both flavours must be handled here. BoxConstraints (variable bounds
 // of the inner Block) are NOT absorbed: they live in u-space and the
 // inner :MILPSolver of the dual MP will enforce them as ordinary
 // variable bounds, so we skip them silently
 auto walk_any = [ & process_one ]( boost::any & any ) {
  if( auto * one = boost::any_cast< FRowConstraint >( & any ) )
   process_one( *one );
  else if( auto * one_p = boost::any_cast< FRowConstraint * >( & any ) ) {
   if( *one_p ) process_one( **one_p );
   }
  else if( auto * vec = boost::any_cast< std::vector< FRowConstraint > >( & any ) ) {
   for( auto & ci : *vec ) process_one( ci );
   }
  else if( auto * vec_p =
                boost::any_cast< std::vector< FRowConstraint > * >( & any ) ) {
   if( *vec_p ) for( auto & ci : **vec_p ) process_one( ci );
   }
  else if( auto * lst = boost::any_cast< std::list< FRowConstraint > >( & any ) ) {
   for( auto & ci : *lst ) process_one( ci );
   }
  else if( auto * lst_p =
                boost::any_cast< std::list< FRowConstraint > * >( & any ) ) {
   if( *lst_p ) for( auto & ci : **lst_p ) process_one( ci );
   }
  else if( boost::any_cast< std::vector< BoxConstraint > * >( & any ) )
   ;  // variable bounds, kept on the inner Block - intentionally skipped
  else if( boost::any_cast< std::vector< BoxConstraint > >( & any ) )
   ;  // ditto, by-value form
  };

 auto & s_cnst =
  const_cast< Vec_any & >( inner->get_static_constraints() );
 for( auto & any : s_cnst ) walk_any( any );

 auto & d_cnst =
  const_cast< Vec_any & >( inner->get_dynamic_constraints() );
 for( auto & any : d_cnst ) walk_any( any );
 }

/*--------------------------------------------------------------------------*/

Block * MasterProblemBlock::get_hard_component( int k ) const
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( nullptr );
 return( HardCmps[ k ] );
 }

/*--------------------------------------------------------------------------*/

Block * MasterProblemBlock::get_easy_component( int k ) const
{
 if( k < 0 || k >= int( EasyCmps.size() ) )
  return( nullptr );
 return( EasyCmps[ k ] );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_dual_norm_squared( void ) const
{
 if( IsPrimal ) {
  if( t_stab <= 0.0 ||
      ( StblType != kProximal && StblType != kDoublyStabilized ) )
   return( 0.0 );

  const auto d = get_d_vector();
  const double inv_t = 1.0 / t_stab;
  return( std::transform_reduce(
             d.cbegin() , d.cend() , 0.0 , std::plus<>() ,
             [ inv_t ]( double dj ) {
              const double zj = - inv_t * dj;
              return( zj * zj );
              } ) );
  }

 return( std::transform_reduce(
              Var_z.cbegin() , Var_z.cend() , 0.0 , std::plus<>() ,
              []( const ColVariable & v ) {
               const double x = v.get_value();
               return( x * x );
               } ) );
 }

/*--------------------------------------------------------------------------*/

int MasterProblemBlock::add_cut( int k , int slot ,
                                 std::vector< double > && g , double alpha ,
                                 bool is_vert )
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  throw( std::invalid_argument(
       "MasterProblemBlock::add_cut: hard-component index out of range" ) );
 if( slot < 0 || slot >= int( slot_to_local[ k ].size() ) )
  throw( std::invalid_argument(
       "MasterProblemBlock::add_cut: slot out of range" ) );

 auto pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  throw( std::logic_error(
       "MasterProblemBlock::add_cut: HardCmps[k] is not a "
       "PolyhedralFunctionBlock" ) );

 const double incoming_b = get_stored_constant( k , g , alpha , is_vert );
 if( ! IsPrimal )
  for( auto & gj : g )
   gj = - gj;

 // defensive eviction: the bundle treats `slot` as a *global* pool index
 // (one occupant across all k); MasterPB stores slot_to_local as a 2D
 // matrix [k][slot] for ergonomic access, so a stale entry can survive
 // in row k' != k after the bundle reassigns `slot` to a new component
 // without an explicit remove_cut(k', slot) call. Silently re-claim the
 // slot regardless of who owned it (single-occupancy write-through)
 if( slot_to_local[ k ][ slot ] >= 0 )
  remove_cut( k , slot );
 else
  for( int kk = 0 ; kk < int( slot_to_local.size() ) ; ++kk )
   if( kk != k && slot_to_local[ kk ][ slot ] >= 0 ) {
    remove_cut( kk , slot );
    break;  // at most one occupant by invariant
    }

 // the new cut is appended to the end of v_A/v_b of PolyhedralFunction;
 // record the resulting local index in slot_to_local. No NBModification
 // reload is needed in linearized-primal rep: the PolyhedralFunctionBlock
 // dispatcher reacts to PolyhedralFunctionModAddd with a single
 // add_dynamic_constraints(f_const, newc, eNoBlck) call, which surfaces
 // as a BlockModAdd<FRowConstraint> picked up by the attached :MILPSolver
 // via the standard dynamic_modification -> add_dynamic_constraint path
 //
 // The driver always passes the physical cut (g, alpha). MPB stores the row in
 // the active representation: primal raw keeps it as-is, while the dual stores
 // the row with the sign required by the coupling equations.
 const int new_local = int( pfb->get_PolyhedralFunction().get_nrows() );
 pfb->get_PolyhedralFunction().add_row( std::move( g ) , incoming_b ,
                                        eModBlck , is_vert );
 slot_to_local[ k ][ slot ] = new_local;
 return( kCutInserted );
 }

/*--------------------------------------------------------------------------*/

int MasterProblemBlock::add_cut( int k , std::vector< double > && g ,
                                 double alpha , bool is_vert )
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( -1 );

 // pick the first free slot of HardCmps[k]
 const auto & st = slot_to_local[ k ];
 for( int slot = 0 ; slot < int( st.size() ) ; ++slot )
  if( st[ slot ] < 0 ) {
   add_cut( k , slot , std::move( g ) , alpha , is_vert );
   return( slot );
   }
 return( -1 );
 }

/*--------------------------------------------------------------------------*/

int MasterProblemBlock::get_bundle_size( int k ) const
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( 0 );
 auto pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( 0 );
 return( int( pfb->get_PolyhedralFunction().get_nrows() ) );
 }

/*--------------------------------------------------------------------------*/

bool MasterProblemBlock::is_bundle_empty( void ) const
{
 return( std::all_of( HardCmps.cbegin() , HardCmps.cend() ,
                      []( const Block * b ) {
                       const auto pfb = dynamic_cast< const PolyhedralFunctionBlock * >( b );
                       return( ! pfb ||
                               const_cast< PolyhedralFunctionBlock * >( pfb )
                                   ->get_PolyhedralFunction().get_nrows() == 0 );
                       } ) );
 }

/*--------------------------------------------------------------------------*/

bool MasterProblemBlock::has_pinned_empty_cmp( void ) const
{
 // detect a hard component whose per-cmp simplex row
 //     sum_i theta^k_i + gamma^k = lambda
 // cannot be satisfied: gamma^k is fixed to 0 (no global LB installed
 // on this PFB) *and* the bundle of theta^k_i is empty. The row then
 // collapses to lambda = 0, contradicting K * lambda = 1 (which holds
 // whenever Var_r and Var_omega are pinned). A solve in this state
 // would be reported as kError / kInfeasible by the inner :MILPSolver;
 // the caller (solve_master) uses this hook to short-circuit instead
 for( auto * b : HardCmps ) {
  auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( b );
  if( ! pfb )
   continue;
  auto & poly = pfb->get_PolyhedralFunction();
  if( poly.get_nrows() > 0 )
   continue;  // bundle non-empty for this cmp, row is satisfiable

  // bundle empty for this cmp: check gamma^k
  auto * gv = pfb->get_static_variable< ColVariable >( "PolyF_gamma" );
  if( gv && gv->is_fixed() )
   return( true );
  }
 return( false );
 }

/*--------------------------------------------------------------------------*/

PolyhedralFunctionBlock *
MasterProblemBlock::pfb_at( const std::vector< Block * > & HardCmps ,
                            int k , const char * fn )
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  throw( std::invalid_argument(
       std::string( "MasterProblemBlock::" ) + fn +
       ": hard-component index out of range" ) );
 auto pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  throw( std::logic_error(
       std::string( "MasterProblemBlock::" ) + fn +
       ": HardCmps[k] is not a PolyhedralFunctionBlock" ) );
 return( pfb );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::remove_cut( int k , int slot )
{
 auto pfb = pfb_at( HardCmps , k , "remove_cut" );
 if( slot < 0 || slot >= int( slot_to_local[ k ].size() ) )
  throw( std::invalid_argument(
       "MasterProblemBlock::remove_cut: slot out of range" ) );

 const int loc = slot_to_local[ k ][ slot ];
 if( loc < 0 )
  return;  // already empty

 // remove the row from PolyhedralFunction and patch the slot->local map:
 // every other slot whose local index was past loc shifts down by one,
 // the slot itself becomes empty
 pfb->get_PolyhedralFunction().delete_row(
                                       PolyhedralFunction::Index( loc ) );
 slot_to_local[ k ][ slot ] = -1;
 for( auto & s : slot_to_local[ k ] )
  if( s > loc )
   --s;
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::modify_cut( int k , int slot ,
                                     std::vector< double > && g , double alpha )
{
 auto pfb = pfb_at( HardCmps , k , "modify_cut" );
 if( slot < 0 || slot >= int( slot_to_local[ k ].size() ) ||
     slot_to_local[ k ][ slot ] < 0 )
  throw( std::invalid_argument(
       "MasterProblemBlock::modify_cut: slot empty or out of range" ) );

 // mirror the on-insert conversion done in add_cut(): callers pass the
 // physical cut (g, alpha), MPB stores the active primal/dual representation.
 const int loc = slot_to_local[ k ][ slot ];
 auto & poly = pfb->get_PolyhedralFunction();
 const bool is_vert = poly.is_row_vertical( PolyhedralFunction::Index( loc ) );

 const double b_store = get_stored_constant( k , g , alpha , is_vert );
 if( ! IsPrimal )
  for( auto & gj : g )
   gj = - gj;

 poly.modify_row( PolyhedralFunction::Index( loc ) ,
                  std::move( g ) , b_store );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::modify_alpha( int k , int slot , double alpha )
{
 auto pfb = pfb_at( HardCmps , k , "modify_alpha" );
 if( slot < 0 || slot >= int( slot_to_local[ k ].size() ) ||
     slot_to_local[ k ][ slot ] < 0 )
  throw( std::invalid_argument(
       "MasterProblemBlock::modify_alpha: slot empty or out of range" ) );

 // mirror the on-insert conversion. The stored row may be dual-signed, so
 // recover the physical g before converting the new physical alpha.
 const int loc = slot_to_local[ k ][ slot ];
 auto & poly = pfb->get_PolyhedralFunction();
 const bool is_vert = poly.is_row_vertical( PolyhedralFunction::Index( loc ) );

 std::vector< double > physical_g;
 const auto & A = poly.get_A();
 if( loc >= 0 && loc < int( A.size() ) ) {
  physical_g = A[ loc ];
  if( ! IsPrimal )
   for( auto & gj : physical_g )
    gj = - gj;
  }
 const double b_store =
  get_stored_constant( k , physical_g , alpha , is_vert );

 poly.modify_constant( PolyhedralFunction::Index( loc ) , b_store );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_theta( int k , int slot ) const
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( 0.0 );
 if( slot < 0 || slot >= int( slot_to_local[ k ].size() ) )
  return( 0.0 );
 const int loc = slot_to_local[ k ][ slot ];
 if( loc < 0 )
  return( 0.0 );

 const auto * pfb = 
     dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( 0.0 );

 return( pfb->get_row_multiplier(
          PolyhedralFunctionBlock::Index( loc ) ) );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_gamma( int k ) const
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( 0.0 );

 const auto * pfb =
  dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( 0.0 );

 if( IsPrimal ) {
  auto & poly =
   const_cast< PolyhedralFunctionBlock * >( pfb )->get_PolyhedralFunction();
  double diagonal_mass = 0.0;
  for( PolyhedralFunction::Index i = 0 ; i < poly.get_nrows() ; ++i )
   if( ! poly.is_row_vertical( i ) )
    diagonal_mass += pfb->get_row_multiplier( i );

  const double gamma = 1.0 - diagonal_mass;
  return( gamma > 0.0 ? gamma : 0.0 );
  }

 auto * gamma =
  const_cast< PolyhedralFunctionBlock * >( pfb )
    ->get_static_variable< ColVariable >( "PolyF_gamma" );
 if( ! gamma )
  return( 0.0 );

 return( pfb->get_v_scale() * gamma->get_value() );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_alpha( int k , int slot ) const
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( 0.0 );
 if( slot < 0 || slot >= int( slot_to_local[ k ].size() ) )
  return( 0.0 );
 const int loc = slot_to_local[ k ][ slot ];
 if( loc < 0 )
  return( 0.0 );

 const auto pfb = dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( 0.0 );
 const auto & b = const_cast< PolyhedralFunctionBlock * >( pfb )
                       ->get_PolyhedralFunction().get_b();
 if( loc >= int( b.size() ) )
  return( 0.0 );
 return( b[ loc ] );
 }

/*--------------------------------------------------------------------------*/

int MasterProblemBlock::find_identical_cut(
                         int k , const std::vector< double > & g ,
                         bool is_vert ) const
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( -1 );

 const auto pfb =
  dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( -1 );

 auto & poly = const_cast< PolyhedralFunctionBlock * >( pfb )
                    ->get_PolyhedralFunction();
 const auto & A = poly.get_A();
 for( int slot = 0 ; slot < int( slot_to_local[ k ].size() ) ; ++slot ) {
  const int loc = slot_to_local[ k ][ slot ];
  if( loc < 0 || loc >= int( A.size() ) )
   continue;
  if( poly.is_row_vertical( PolyhedralFunction::Index( loc ) ) != is_vert )
   continue;
  if( A[ loc ].size() != g.size() )
   continue;
  bool same = true;
  for( std::size_t j = 0 ; j < g.size() ; ++j )
   if( A[ loc ][ j ] != ( IsPrimal ? g[ j ] : - g[ j ] ) ) {
    same = false;
    break;
    }
  if( same )
   return( slot );
  }

 return( -1 );
}

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_stored_constant(
                         int k , const std::vector< double > & g ,
                         double alpha , bool is_vert ) const
{
 double stored = alpha;
 if( is_vert )
  return( IsPrimal ? stored : - stored );

 if( IsPrimal ) {
  if( f_v2_form || k < 0 || k >= int( f_F_at_x_bar.size() ) )
   return( stored );

  const auto & xbar = f_x_bar;
  const std::size_t n = std::min( g.size() , xbar.size() );
  double dot = 0.0;
  for( std::size_t j = 0 ; j < n ; ++j )
   dot += g[ j ] * xbar[ j ];
  return( stored + dot - f_F_at_x_bar[ k ] );
  }

 if( k < 0 || k >= int( f_F_at_x_bar.size() ) )
  return( stored );

 if( f_v2_form )   // iterate form: g . x_bar lives in the linear z term
  stored = - alpha + f_F_at_x_bar[ k ];
 else {
  const auto & xref = cut_ref();
  const std::size_t n = std::min( g.size() , xref.size() );
  double dot = 0.0;
  for( std::size_t j = 0 ; j < n ; ++j )
   dot += g[ j ] * xref[ j ];
  stored = f_F_at_x_bar[ k ] - alpha - dot;
  }

 if( ! IsConvex )
  stored = - stored;

 return( stored );
}

/*--------------------------------------------------------------------------*/

const std::vector< double > &
                          MasterProblemBlock::get_subgradient( int k , int slot ) const
{
 static const std::vector< double > empty_vec;

 if( k < 0 || k >= int( HardCmps.size() ) )
  return( empty_vec );
 if( slot < 0 || slot >= int( slot_to_local[ k ].size() ) )
  return( empty_vec );
 const int loc = slot_to_local[ k ][ slot ];
 if( loc < 0 )
  return( empty_vec );

 const auto pfb = dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( empty_vec );
 const auto & A = const_cast< PolyhedralFunctionBlock * >( pfb )
                       ->get_PolyhedralFunction().get_A();
 if( loc >= int( A.size() ) )
  return( empty_vec );
 return( A[ loc ] );
 }

/*--------------------------------------------------------------------------*/

bool MasterProblemBlock::is_subgradient( int k , int slot ) const
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( false );
 if( slot < 0 || slot >= int( slot_to_local[ k ].size() ) )
  return( false );
 const int loc = slot_to_local[ k ][ slot ];
 if( loc < 0 )
  return( false );

 const auto pfb = dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( false );
 return( ! const_cast< PolyhedralFunctionBlock * >( pfb )
              ->get_PolyhedralFunction().is_row_vertical(
                                            PolyhedralFunction::Index( loc ) ) );
 }

/*--------------------------------------------------------------------------*/

int MasterProblemBlock::get_vertical_count( void ) const
{
 int total = 0;
 for( const auto * b : HardCmps ) {
  const auto pfb = dynamic_cast< const PolyhedralFunctionBlock * >( b );
  if( ! pfb )
   continue;
  auto & pf = const_cast< PolyhedralFunctionBlock * >( pfb )
                  ->get_PolyhedralFunction();
  const auto n = pf.get_nrows();
  for( PolyhedralFunction::Index i = 0 ; i < n ; ++i )
   if( pf.is_row_vertical( i ) )
    ++total;
  }
 return( total );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_FiBLambda( int k ) const
{
 if( IsPrimal ) {
  // primal linearized rep: every hard component k owns its own epigraph
  // variable f_v inside the PolyhedralFunctionBlock; the master-side
  // Var_v_hard is unused. k in [0, NoHardCmps) returns v^k* of
  // that component; k == -1 (default) sums v^k* over every hard cmp.
  // The proximal/level stabilization terms on d / v are not subtracted
  // out: the surrounding bundle solver expects the "model value" v*
  auto delta_v_of = [ this ]( int kk ) -> double {
   if( kk < 0 || kk >= int( HardCmps.size() ) )
    return( 0.0 );
   auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ kk ] );
   if( ! pfb )
    return( 0.0 );
   const auto * v = pfb->get_v();
   if( ! v )
    return( 0.0 );

   const double reference =
    ( f_v2_form && kk < int( f_F_at_x_bar.size() ) )
    ? f_F_at_x_bar[ kk ] : 0.0;
   return( v->get_value() - reference );
   };
  if( k >= 0 )
   return( delta_v_of( k ) );

  double sum = 0.0;
  for( int kk = 0 ; kk < int( HardCmps.size() ) ; ++kk )
   sum += delta_v_of( kk );

  const auto d = get_d_vector();
  const std::size_t n = std::min( d.size() , f_linear_part.size() );
  double bd = 0.0;
  for( std::size_t j = 0 ; j < n ; ++j )
   bd += f_linear_part[ j ] * d[ j ];
  sum += bd;

  return( sum );
  }

 // dual MP: v*[k] is the cutting-plane model's *predicted decrease* of
 // component k at the proximal step d* = -t z*, NOT the aggregated
 // linearization error Sigma_k (a frequent confusion: Sigma_k >= 0 is
 // the model gap at the stability centre, while v*[k] <= 0 is the gain
 // the model promises by moving to Lambda1). With g^k = sum_i theta^k_i
 // g^k_i the per-cmp aggregated (textbook) subgradient and z* the total
 // aggregated subgradient (linear part included),
 //   v*[k] = F^k_model( Lambda1 ) - F^k( Lambda )
 //         = < g^k , d* > - Sigma_k  =  -t < g^k , z* > - Sigma_k
 // and, summing the linear-part contribution < b , d* > too,
 //   v* (total) = -t || z* ||^2 - Sigma   ( = -( Sigma + 2 D*_t(z*) ) ).
 //
 // get_aggregated_subgradient(k) returns the physical aggregate g^k even
 // when the dual PFB stores the sign-flipped row -g^k_i internally, so the
 // per-component expression keeps the textbook sign:
 //   v*[k] = -Sigma_k - t < get_aggregated_subgradient(k) , z* >
 //
 // The previous implementation returned + Sigma_k, i.e. the
 // linearization error with the wrong sign and the wrong meaning: it
 // turned vStar positive, inflated the SS target
 //   UpTrgt = UpRifFi + (1 - m2) vStar ,
 // made the bundle accept spurious serious steps, and corrupted the
 // trajectory (the master appeared to "ascend" on a minimisation)
 const double t = t_stab;
 auto has_model_row = [ this ]( int kk ) -> bool {
  if( kk < 0 || kk >= int( HardCmps.size() ) )
   return( false );

  if( kk < int( f_LB_raw.size() ) && std::isfinite( f_LB_raw[ kk ] ) )
   return( true );

  const auto * pfb =
   dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ kk ] );
  if( ! pfb )
   return( false );

  auto & poly =
   const_cast< PolyhedralFunctionBlock * >( pfb )->get_PolyhedralFunction();
  const auto n = poly.get_nrows();
  for( PolyhedralFunction::Index i = 0 ; i < n ; ++i )
   if( ! poly.is_row_vertical( i ) )
    return( true );

  return( false );
  };

 if( k >= 0 ) {
  if( ! has_model_row( k ) )
   return( Inf< double >() );

  const double sigma_k = get_aggregated_alpha( k );
  const auto zk = get_aggregated_subgradient( k );
  const auto zt = get_z_vector();
  double dot = 0.0;
  const std::size_t n = std::min( zk.size() , zt.size() );
  for( std::size_t j = 0 ; j < n ; ++j )
   dot += zk[ j ] * zt[ j ];
  return( - sigma_k - t * dot );
  }

 for( int kk = 0 ; kk < int( HardCmps.size() ) ; ++kk )
  if( ! has_model_row( kk ) )
   return( Inf< double >() );

 // total v* = -( Sigma + t || z* ||^2 ), with || z* ||^2 the squared
 // norm of the *total* aggregated subgradient (linear part included)
 return( - ( get_aggregated_alpha( -1 ) + t * get_dual_norm_squared() ) );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > MasterProblemBlock::get_z_vector( void ) const
{
 std::vector< double > out;
 if( IsPrimal ) {
  // Under proximal stabilization, stationarity of the primal master gives
  // z* + d*/t = 0. This recovers the complete essential subgradient,
  // including the linear part and active-domain/box multipliers.
  out = get_d_vector();
  if( t_stab > 0.0 &&
      ( StblType == kProximal || StblType == kDoublyStabilized ) )
   for( auto & zj : out )
    zj = - zj / t_stab;
  else
   std::fill( out.begin() , out.end() , 0.0 );
  return( out );
  }
 out.reserve( Var_z.size() );
 for( const auto & zj : Var_z )
  out.push_back( zj.get_value() );
 return( out );
 }

/*--------------------------------------------------------------------------*/

std::vector< double >
MasterProblemBlock::get_aggregated_subgradient( int k ) const
{
 // The abstract multiplier is stored differently in the two linearized PFB
 // representations: as the dual value of a cut constraint in the primal one,
 // and as an explicit theta variable in the dual one. PFB::get_row_multiplier()
 // hides this distinction and also maps the value back to the units of the
 // original unscaled row. The dual PFB stores rows as -g to satisfy its
 // coupling equations; flip them back here so callers always see the physical
 // aggregate subgradient.
 std::vector< double > out( NumVars , 0.0 );

 if( k < 0 || k >= int( HardCmps.size() ) )
  return( out );

 const auto * pfb =
  dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( out );

 const auto & A =
  const_cast< PolyhedralFunctionBlock * >( pfb )
    ->get_PolyhedralFunction().get_A();
 const double row_sign = IsPrimal ? 1.0 : -1.0;

 for( std::size_t i = 0 ; i < A.size() ; ++i ) {
  const double theta =
   pfb->get_row_multiplier( PolyhedralFunctionBlock::Index( i ) );

  if( theta == 0.0 )
   continue;

  const auto & gi = A[ i ];
  const std::size_t n = std::min( gi.size() , out.size() );

  for( std::size_t j = 0 ; j < n ; ++j )
   out[ j ] += row_sign * theta * gi[ j ];
  }

 return( out );
}

/*--------------------------------------------------------------------------*/

std::vector< double > MasterProblemBlock::get_d_vector( void ) const
{
 std::vector< double > out;

 if( IsPrimal ) {
  out.reserve( Var_d.size() );
  for( std::size_t i = 0 ; i < Var_d.size() ; ++i )
   out.push_back( Var_d[ i ].get_value() -
                  ( f_v2_form ? f_x_bar[ i ] : 0.0 ) );
  return( out );
  }

 // dual MP, proximal stabilization: d* = -t * z*
 out.reserve( Var_z.size() );
 for( const auto & zj : Var_z )
  out.push_back( - t_stab * zj.get_value() );
 return( out );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_Gid_aggregate( void ) const
{
 if( IsPrimal ) {
  // primal MP: <z*, d*> = sum_k <z^k, d>, with z^k aggregated through
  // the per-PFB get_aggregated_subgradient and d directly read from
  // the master-side Var_d
  const auto z = get_z_vector();
  const auto d = get_d_vector();
  if( z.empty() || d.empty() )
   return( 0.0 );
  const std::size_t n = std::min( z.size() , d.size() );
  double s = 0.0;
  for( std::size_t j = 0 ; j < n ; ++j )
   s += z[ j ] * d[ j ];
  return( s );
  }

 // dual MP under proximal stabilization:
 //   d* = -t * z*, so z* . d* = -t * || z* ||^2
 return( - t_stab * get_dual_norm_squared() );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_Gid( int k , int slot ) const
{
 const auto & g = get_subgradient( k , slot );
 if( g.empty() )
  return( 0.0 );
 const auto d = get_d_vector();
 if( d.size() != g.size() )
  return( 0.0 );
 const double stored_dot =
  std::inner_product( g.begin() , g.end() , d.begin() , 0.0 );
 return( IsPrimal ? stored_dot : - stored_dot );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::sensitivity_analysis( double & vl ,
                                               double & vc ) const
{
 if( IsPrimal ) {
  vl = 0.0;
  vc = 0.0;
  return;
  }
 const auto nz2 = get_dual_norm_squared();
 vl = - nz2 / 2.0;
 vc = get_FiBLambda() - vl * t_stab;
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::invalidate_subgradients( int k )
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return;
 auto pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return;
 pfb->add_Modification( std::make_shared< NBModification >( pfb ) );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_alphas_bulk( int k ,
                                          const std::vector< double > & alphas )
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return;
 const auto & st = slot_to_local[ k ];
 const int n = std::min( int( alphas.size() ) , int( st.size() ) );
 for( int slot = 0 ; slot < n ; ++slot )
  if( st[ slot ] >= 0 )
   modify_alpha( k , slot , alphas[ slot ] );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_aggregated_alpha( int k ) const
{
 if( IsPrimal ) {
  // The primal master directly provides d* and the model decrease v*. For
  // every component the aggregate-cut identity is
  //     v*[k] = < z*[k] , d* > - Sigma*[k].
  // Using it here automatically includes horizontal lower-bound mass and
  // vertical/domain multipliers without exposing their abstract dual objects.
  const auto d = get_d_vector();

  if( k >= 0 ) {
   const auto zk = get_aggregated_subgradient( k );
   const std::size_t n = std::min( zk.size() , d.size() );
   double zd = 0.0;
   for( std::size_t j = 0 ; j < n ; ++j )
    zd += zk[ j ] * d[ j ];
   return( zd - get_FiBLambda( k ) );
   }

  return( get_Gid_aggregate() - get_FiBLambda() );
  }

 // b[i] is stored in the physical PolyhedralFunction units. The multiplier
 // must therefore also be expressed in physical units. get_row_multiplier()
 // hides both the active PFB representation and any internal row scaling.
 auto contrib = [ this ]( int kk ) -> double {
  if( kk < 0 || kk >= int( HardCmps.size() ) )
   return( 0.0 );

  const auto * pfb =
   dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ kk ] );
  if( ! pfb )
   return( 0.0 );

  auto & poly =
   const_cast< PolyhedralFunctionBlock * >( pfb )->get_PolyhedralFunction();
  const auto & b = poly.get_b();

  // recover the true physical linearization error Sigma_k = sum_i theta_i
  // ( the lin-error at x_bar ) from the objective-sense-signed stored b[ i ].
  //  iterate form: b omits g . x_bar entirely ( it rides in the +x_bar^T R
  //    lin-z ), so add back sum_i theta_i ( A_i . x_bar );
  //  displacement form with a lazy reference: b baked g . x_ref, so add back
  //    only the residual sum_i theta_i ( A_i . ( x_bar - x_ref ) );
  //  plain displacement ( f_xref_tol == 0 ): b is already the lin-error at
  //    x_bar, add nothing.
  const bool iterate = ( f_v2_form != 0 );
  const auto & xref = cut_ref();
  const bool lazy = ( ! iterate ) && ( f_xref_tol > 0.0 ) &&
                    ( xref.size() == f_x_bar.size() ) &&
                    ( &xref != &f_x_bar );
  const auto & A = poly.get_A();

  double s = 0.0;
  for( std::size_t i = 0 ; i < b.size() ; ++i ) {
   const double theta =
    pfb->get_row_multiplier( PolyhedralFunction::Index( i ) );
   s += theta * ( IsConvex ? b[ i ] : - b[ i ] );
   if( ( iterate || lazy ) && i < A.size() ) {
    const auto & Ai = A[ i ];
    const std::size_t n = std::min( Ai.size() , f_x_bar.size() );
    double add = 0.0;
    for( std::size_t j = 0 ; j < n ; ++j )
     add += Ai[ j ] * ( iterate ? f_x_bar[ j ]
                                : ( f_x_bar[ j ] - xref[ j ] ) );
    s += theta * add;
    }
   }

  // The component lower bound is a valid horizontal row. If gamma carries
  // mass, it contributes no subgradient, but it does contribute its
  // linearization error F_k(x_bar) - LB_k to Sigma_k.
  if( kk < int( f_LB_raw.size() ) && kk < int( f_F_at_x_bar.size() ) &&
      std::isfinite( f_LB_raw[ kk ] ) )
   s += get_gamma( kk ) * ( f_F_at_x_bar[ kk ] - f_LB_raw[ kk ] );

  return( s );
  };

 if( k >= 0 )
  return( contrib( k ) );

 const double mp_obj = get_master_objective_value();
 if( std::isfinite( mp_obj ) ) {
  // The full dual master objective contains the model gap and the quadratic
  // stabilization term. In the convex/min representation it is Sigma + D;
  // in the concave/max representation the objective sense stores the negated
  // value. Recover the common minimization value first, then subtract D.
  const auto * obj = get_objective();
  double signed_obj =
   ( obj && obj->get_sense() == Objective::eMax ) ? - mp_obj : mp_obj;

  // The iterate form keeps the hard-component cuts relative to F_k(x_bar),
  // but writes the linear 0-th component as b . x rather than b . (x-x_bar).
  // Since the coupling rows contain -lambda * b, converting the signed
  // objective back to the displacement-form model gap requires adding the
  // solution-independent constant lambda * b . x_bar.
  if( f_v2_form && f_linear_part.size() == f_x_bar.size() )
   signed_obj += get_lambda() *
    std::inner_product( f_linear_part.begin() , f_linear_part.end() ,
                        f_x_bar.begin() , 0.0 );

  return( signed_obj - 0.5 * t_stab * get_dual_norm_squared() );
  }

 double total = 0.0;
 for( int kk = 0 ; kk < int( HardCmps.size() ) ; ++kk )
  total += contrib( kk );

 return( total );
}

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_master_objective_value( void ) const
{
 // At the very first iteration, no known objective value is yet known.
 // We simply return 0.0 as the main solver should be aware of discarding
 // such value.
 if( is_bundle_empty() ){
  return( 0.0 );
 }

 const auto & solvers = get_registered_solvers();
 if( solvers.empty() )
  return( Inf< double >() );

 const auto * obj = get_objective();
 if( ! obj )
  return( Inf< double >() );

 auto * slv = solvers.front();
 double val = Inf< double >();
 if( obj->get_sense() == Objective::eMin )
  val = slv->get_ub();
 else if( obj->get_sense() == Objective::eMax )
  val = slv->get_lb();

 if( std::isfinite( val ) )
  return( val );

 // Some Solver implementations may expose only the bound side even at
 // optimality. Use it as a harmless fallback when it is finite.
 val = ( obj->get_sense() == Objective::eMin ) ? slv->get_lb()
                                               : slv->get_ub();
 return( std::isfinite( val ) ? val : Inf< double >() );
}

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_raw_aggregated_alpha( int k ) const
{
 if( IsPrimal ) {
  auto contrib = [ this ]( int kk ) -> double {
   if( kk < 0 || kk >= int( HardCmps.size() ) )
    return( 0.0 );

   const auto * pfb =
    dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ kk ] );
   if( ! pfb )
    return( 0.0 );

   auto & poly =
    const_cast< PolyhedralFunctionBlock * >( pfb )->get_PolyhedralFunction();
   const auto & A = poly.get_A();
   const auto & b = poly.get_b();
   const double F_xb = kk < int( f_F_at_x_bar.size() )
                       ? f_F_at_x_bar[ kk ] : 0.0;

   double alpha = 0.0;
   for( PolyhedralFunction::Index i = 0 ; i < poly.get_nrows() ; ++i ) {
    if( poly.is_row_vertical( i ) )
     continue;

    double raw_b = b[ i ];
    if( ! f_v2_form && i < A.size() ) {
     const std::size_t n = std::min( A[ i ].size() , f_x_bar.size() );
     double gxbar = 0.0;
     for( std::size_t j = 0 ; j < n ; ++j )
      gxbar += A[ i ][ j ] * f_x_bar[ j ];
     raw_b += F_xb - gxbar;
     }

    alpha += pfb->get_row_multiplier( i ) * raw_b;
    }

   return( alpha );
   };

  if( k >= 0 )
   return( contrib( k ) );

  double total = 0.0;
  for( int kk = 0 ; kk < int( HardCmps.size() ) ; ++kk )
   total += contrib( kk );
  return( total );
  }

 // A[i] and b[i] are stored in the physical PolyhedralFunction units.
 // get_row_multiplier() returns the corresponding physical theta_i,
 // independently of the active PFB representation and internal scaling.
 auto contrib = [ this ]( int kk ) -> double {
  if( kk < 0 || kk >= int( HardCmps.size() ) )
   return( 0.0 );

  const auto * pfb =
   dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ kk ] );
  if( ! pfb )
   return( 0.0 );

  auto & poly =
   const_cast< PolyhedralFunctionBlock * >( pfb )
     ->get_PolyhedralFunction();

  const auto & b = poly.get_b();
  const auto & A = poly.get_A();

  // recover the raw alpha aggregate sum_i theta_i alpha_i from the stored b,
  // first mapping the objective-sense-signed b back to physical error units.
  //  displacement form: b_i = F_k - alpha_i + g_i . x_bar, hence
  //    sum theta alpha = F_k sum_th - b_sum + sum theta ( A_i . x_bar ).
  //  iterate form: b_i = F_k - alpha_i ( no g . x_bar ), hence
  //    sum theta alpha = F_k sum_th - b_sum, with NO z_dot term.
  const bool iterate = ( f_v2_form != 0 );

  double b_sum  = 0.0;
  double sum_th = 0.0;
  double z_dot  = 0.0;

  for( std::size_t i = 0 ; i < b.size() ; ++i ) {
   const double theta =
    pfb->get_row_multiplier( PolyhedralFunction::Index( i ) );

   b_sum  += theta * ( IsConvex ? b[ i ] : - b[ i ] );
   sum_th += theta;

   if( ( ! iterate ) && i < A.size() ) {
    const auto & Ai = A[ i ];
    // reconstruct the raw alpha from b in its OWN stored frame: b baked
    // g . x_ref (= x_bar when not lazy), so the z . x term must use x_ref
    const auto & xref = cut_ref();
    const std::size_t n = std::min( Ai.size() , xref.size() );

    double dot = 0.0;
    for( std::size_t j = 0 ; j < n ; ++j )
     dot += Ai[ j ] * xref[ j ];

    z_dot += theta * dot;
    }
   }

  const double F_xb = ( kk < int( f_F_at_x_bar.size() ) )
                      ? f_F_at_x_bar[ kk ] : 0.0;

  return( F_xb * sum_th - b_sum + z_dot );
  };

 if( k >= 0 )
  return( contrib( k ) );

 double total = 0.0;
 for( int kk = 0 ; kk < int( HardCmps.size() ) ; ++kk )
  total += contrib( kk );

 return( total );
}

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_raw_aggregated_alpha_with_LB( int k ) const
{
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( 0.0 );

 double alpha = get_raw_aggregated_alpha( k );

 if( k < int( f_LB_raw.size() ) && std::isfinite( f_LB_raw[ k ] ) )
  alpha += get_gamma( k ) * f_LB_raw[ k ];

 return( alpha );
}

/*--------------------------------------------------------------------------*/

const std::vector< double > & MasterProblemBlock::get_alphas( int k ) const
{
 static const std::vector< double > empty_vec;

 if( k < 0 || k >= int( HardCmps.size() ) )
  return( empty_vec );
 const auto pfb = dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( empty_vec );
 return( const_cast< PolyhedralFunctionBlock * >( pfb )
              ->get_PolyhedralFunction().get_b() );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > MasterProblemBlock::get_thetas( int k ) const
{
 std::vector< double > out;
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( out );

 const auto * pfb =
  dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( out );

 const auto n = pfb->get_PolyhedralFunction().get_nrows();
 out.reserve( n );

 for( PolyhedralFunctionBlock::Index i = 0 ; i < n ; ++i )
  out.push_back( pfb->get_row_multiplier( i ) );

 return( out );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_C( double C )
{
 // store the shift; the next call to set_global_LB / set_f_lev will
 // pick it up when committing the LB / Lvl coefficient to the dual
 // Objective. Note that we do *not* re-commit the existing coefficient
 // here: the caller is expected to drive the refresh by calling
 // set_global_LB / set_f_lev again with the up-to-date LB / Lvl, which
 // is the natural pattern when x_bar moves.
 f_C = C;
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_box( const std::vector< double > & L ,
                                  const std::vector< double > & U )
{
 if( ! L.empty() && int( L.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: L must be empty or of size NumVars" ) );
 if( ! U.empty() && int( U.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: U must be empty or of size NumVars" ) );

 f_L = L;
 f_U = U;

 if( IsPrimal ) {
  if( int( Bounds_d.size() ) != NumVars )
   return;

  for( int j = 0 ; j < NumVars ; ++j ) {
   double lhs = ( ! f_L.empty() && std::isfinite( f_L[ j ] ) )
                ? f_L[ j ] : - Inf< double >();
   double rhs = ( ! f_U.empty() && std::isfinite( f_U[ j ] ) )
                ? f_U[ j ] : Inf< double >();
   if( ! f_v2_form ) {
    const double xj = j < int( f_x_bar.size() ) ? f_x_bar[ j ] : 0.0;
    if( std::isfinite( lhs ) )
     lhs -= xj;
    if( std::isfinite( rhs ) )
     rhs -= xj;
    }
   Bounds_d[ j ].set_lhs( lhs );
   Bounds_d[ j ].set_rhs( rhs );
   }
  return;
  }

 auto obj = dynamic_cast< FRealObjective * >( get_objective() );
 auto dqf = obj ? dynamic_cast< DQuadFunction * >( obj->get_function() )
                : nullptr;
 if( ! dqf || s_plus_obj_idx < 0 || s_minus_obj_idx < 0 )
  return;

 // : s^+_j gets coefficient +sgn*(L_j - x_bar_j), s^-_j
 // gets -sgn*(U_j - x_bar_j); the slack stays fixed to 0 (and the
 // coefficient stays 0) whenever the corresponding bound is non-finite.
 const double sgn = IsConvex ? -1.0 : 1.0;
 for( int j = 0 ; j < NumVars ; ++j ) {
  const double xj = ( j < int( f_x_bar.size() ) ) ? f_x_bar[ j ] : 0.0;

  const bool has_L = ! f_L.empty() && std::isfinite( f_L[ j ] );
  Var_s_plus[ j ].is_fixed( ! has_L , eNoMod );
  if( ! has_L ) Var_s_plus[ j ].set_value( 0.0 );
  dqf->modify_term( DQuadFunction::Index( s_plus_obj_idx + j ) ,
                    has_L ? sgn * ( f_L[ j ] - xj ) : 0.0 , 0.0 );

  const bool has_U = ! f_U.empty() && std::isfinite( f_U[ j ] );
  Var_s_minus[ j ].is_fixed( ! has_U , eNoMod );
  if( ! has_U ) Var_s_minus[ j ].set_value( 0.0 );
  dqf->modify_term( DQuadFunction::Index( s_minus_obj_idx + j ) ,
                    has_U ? - sgn * ( f_U[ j ] - xj ) : 0.0 , 0.0 );
  }
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_global_LB( double LB )
{
 if( IsPrimal || r_obj_idx < 0 )
  return;

 auto obj = dynamic_cast< FRealObjective * >( get_objective() );
 if( ! obj )
  return;
 auto dqf = dynamic_cast< DQuadFunction * >( obj->get_function() );
 if( ! dqf )
  return;

 // dual master sense matches IsConvex; the +LB_xbar*r term of the
 // textbook max form is negated when IsConvex (the whole objective is
 // in min form). The shifted lower bound is  LB_xbar = LB + f_C
 // . A non-finite LB (= no global bound
 // known) collapses the LB*r term to zero rather than leaking Inf
 // into the Objective.
 const bool finite = std::isfinite( LB );
 const double sgn = IsConvex ? -1.0 : 1.0;
 const double coeff = finite ? sgn * ( LB + f_C ) : 0.0;
 dqf->modify_term( DQuadFunction::Index( r_obj_idx ) , coeff , 0.0 );

 // cache the translated global lower bound LB + f_C (the value the
 // r-multiplier carries) so that set_fictitious_LB() can place the
 // fictitious per-component bound strictly *below* it, mirroring the
 // legacy OSIMPSolver device (gamma_i cost = global_LB - 1)
 f_global_LB_xbar = finite ? ( LB + f_C ) : - Inf< Function::FunctionValue >();

 // Var_r is meaningful only when LB is finite: unfix it now (it was
 // pinned to 0 in CreateDualMP), or re-pin it to 0 if LB has gone
 // back to -infinity. See the matching comment in CreateDualMP for
 // the rationale on the master normalization
 if( finite ) {
  if( Var_r.is_fixed() )
   Var_r.is_fixed( false , eNoMod );
  }
 else {
  Var_r.set_value( 0 );
  if( ! Var_r.is_fixed() )
   Var_r.is_fixed( true , eNoMod );
  }
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_F_at_x_bar( int k , double value )
{
 if( k < 0 || k >= int( f_F_at_x_bar.size() ) )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_F_at_x_bar: k out of range" ) );
 f_F_at_x_bar[ k ] = value;
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_F_at_x_bar( int k ) const
{
 if( k < 0 || k >= int( f_F_at_x_bar.size() ) )
  return( 0.0 );
 return( f_F_at_x_bar[ k ] );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_reference(
                              const std::vector< double > & x_bar ,
                              const std::vector< double > & F_at_x_bar )
{
 if( int( F_at_x_bar.size() ) != int( f_F_at_x_bar.size() ) )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_reference: F_at_x_bar size != NoHardCmps" ) );
 if( int( x_bar.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_reference: x_bar size != NumVars" ) );

 // capture the old reference *before* writing the new one, so the per-cut
 // shift below can be expressed as a delta against the old state. Diagonal
 // cuts in the master store the linearization error at the current x_bar,
 // signed to match the master objective sense, not the raw cut constant: this
 // is an internal numerical-stability choice and the caller never sees it. When
 // the reference moves to ( x_bar_new , F_new ) the stored b[ i ] must be
 // refreshed to track the new x_bar
 const std::vector< double > old_x_bar = f_x_bar;
 const std::vector< double > old_F     = f_F_at_x_bar;

 // the iterate form ( f_v2_form == 1 ) stores no g . x_bar in the cuts, so
 // their refresh is the uniform per-component dF alone ( no g-shift ); only
 // the displacement form carries a per-cut g-shift dxbar
 const bool with_gx = ( f_v2_form == 0 );

 // Decide the per-cut shift dxbar AND re-align the lazy storage reference
 // BEFORE refreshing x_bar, so set_x_bar() computes the final lin-z gap in a
 // single pass ( = 0 when re-aligning, since f_x_ref is moved to x_bar here ).
 //  - plain displacement ( f_xref_tol == 0 ): shift by g . ( x_bar -
 //    old_x_bar ) every serious step ( exactly the legacy behaviour );
 //  - lazy reference ( f_xref_tol > 0 ): shift by g . ( x_bar - x_ref ) and
 //    re-align x_ref to x_bar ONLY when || x_bar - x_ref ||_inf exceeds the
 //    tolerance; otherwise defer the whole g . ( x_bar - x_ref ) to the lin-z
 //    and shift the cut constants by the uniform dF alone.
 // ( iterate form: with_gx is false, dxbar stays empty -> dF-only refresh. )
 std::vector< double > dxbar;
 if( ( ! IsPrimal ) && with_gx && old_x_bar.size() == x_bar.size() ) {
  if( f_xref_tol > 0.0 ) {
   if( f_x_ref.size() != x_bar.size() )
    f_x_ref = old_x_bar;
   double drift = 0.0;
   for( std::size_t j = 0 ; j < x_bar.size() ; ++j )
    drift = std::max( drift , std::abs( x_bar[ j ] - f_x_ref[ j ] ) );
   if( drift > f_xref_tol ) {
    dxbar.resize( x_bar.size() );
    for( std::size_t j = 0 ; j < x_bar.size() ; ++j )
     dxbar[ j ] = x_bar[ j ] - f_x_ref[ j ];
    f_x_ref = x_bar;  // re-align now so set_x_bar's lin-z gap is 0
    }
   }
  else {
   dxbar.resize( x_bar.size() );
   for( std::size_t j = 0 ; j < x_bar.size() ; ++j )
    dxbar[ j ] = x_bar[ j ] - old_x_bar[ j ];
   }
  }

 // refresh x_bar via the regular setter (which also updates any state that
 // depends on it, e.g. the lin-z gap on Var_z and the box-slack
 // coefficients); the re-align decision is already baked into f_x_ref
 set_x_bar( x_bar );

 // first-contact init of the lazy storage reference ( harmless / unused when
 // f_xref_tol == 0, where dxbar shifts against old_x_bar )
 if( ( ! IsPrimal ) && f_x_ref.size() != f_x_bar.size() )
  f_x_ref = f_x_bar;

 // commit the new per-component reference values
 f_F_at_x_bar = F_at_x_bar;

 // shift each *diagonal* cut's stored constant to track the moving reference.
 //  displacement form ( f_v2_form == 0 ): the stored full lin-error
 //    b[ i ] = sense_sign * ( F_k - alpha + g . x_bar ) is refreshed by
 //    sense_sign * ( ( F_new - F_old ) +
 //                    A[ i ] . ( x_bar_new - x_bar_old ) ).
 //  iterate form ( f_v2_form == 1 ): the stored constant is only the
 //    F_k-relative b[ i ] = F_k - alpha ( the g . x_bar cross-term lives in
 //    the explicit +x_bar^T R lin-z ), so its refresh is the UNIFORM
 //    per-component delta dF = F_new - F_old with NO per-cut dot product.
 // (vertical cuts encode a feasibility constraint whose right-hand side is
 // reference-independent and are left alone). The math is the dual of the
 // legacy "shift after a current-point change" loop, but here the bundle is
 // walked once and the update is invisible to the surrounding driver
 if( IsPrimal )
  return;  // primal MP carries the raw lin-error elsewhere; nothing to shift

 for( int k = 0 ; k < int( HardCmps.size() ) ; ++k ) {
  if( old_x_bar.size() != f_x_bar.size() )
   break;  // first ever set_reference; old_x_bar empty, nothing to shift

  auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
  if( ! pfb )
   continue;
  auto & poly = pfb->get_PolyhedralFunction();
  const auto & A = poly.get_A();
  const auto & b = poly.get_b();

  const double dF = ( k < int( old_F.size() ) )
                    ? ( F_at_x_bar[ k ] - old_F[ k ] )
                    : F_at_x_bar[ k ];

  // walk every row of HardCmps[ k ]; vertical rows are skipped because
  // their right-hand side does not depend on the reference
  for( std::size_t i = 0 ; i < b.size() ; ++i ) {
   if( poly.is_row_vertical( PolyhedralFunction::Index( i ) ) )
    continue;
   double dx_dot = 0.0;
   if( ! dxbar.empty() ) {  // displacement g-shift ( x_bar - old_x_bar, or
    const auto & Ai = A[ i ];          // x_bar - x_ref under a lazy reference;
    const std::size_t n = std::min( Ai.size() , dxbar.size() );  // empty when
    for( std::size_t j = 0 ; j < n ; ++j )       // deferred or in iterate form )
     dx_dot += Ai[ j ] * dxbar[ j ];
    }
   const double delta = ( IsConvex ? 1.0 : -1.0 ) * ( dF + dx_dot );
   if( delta == 0.0 )
    continue;
   poly.modify_constant( PolyhedralFunction::Index( i ) , b[ i ] + delta );
   }

  // also refresh the per-cmp gamma * ( F_k(x_bar) - LB ) coefficient: the
  // raw LB^k is invariant under the reference move, only the ( F - LB )
  // linearization-error form changes (same frame as set_LB and the diagonal
  // cuts). Skip when the bound is genuinely absent (raw == -INF)
  if( k < int( f_LB_raw.size() ) ) {
   const double LB_raw = f_LB_raw[ k ];
   if( std::isfinite( LB_raw ) )
    poly.modify_bound( ( IsConvex ? 1.0 : -1.0 ) *
                       ( F_at_x_bar[ k ] - LB_raw ) );
   }
  }
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_x_bar( const std::vector< double > & x_bar )
{
 if( int( x_bar.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_x_bar: x_bar must have NumVars entries" ) );

 // record the new stability centre: every linearization-error translation
 // (add_cut / modify_cut / modify_alpha incoming_b, the set_reference cut
 // shift, get_aggregated_alpha's z . x_bar reconstruction) reads f_x_bar,
 // so it MUST track the actual centre. Without this assignment f_x_bar
 // stayed pinned at its zero-initialised value and every translation
 // silently assumed x_bar == 0: harmless while the centre never leaves the
 // origin (the unbounded instances), but as soon as a finite optimum makes
 // the centre move the cut shift dropped its geometric A . (x_bar_new -
 // x_bar_old) term, producing negative (invalid) stored linearization
 // errors and a corrupted master
 f_x_bar = x_bar;

 // refresh the absorbed BendersBFunction RowConstraints (primal MP):
 // their right-hand side(s) are  A_i . x_bar + b_i  for the side(s)
 // marked dynamic by the original ConstraintSide; the static side(s)
 // keep their snapshot from absorb_BBF_into_primal_MP()
 if( IsPrimal ) {
  // In raw form x_bar changes the linear part of
  // ||x-x_bar||^2/(2t). Defer the abstract objective update until the next
  // actual solve, where it is batched with any intervening t / b changes.
  if( f_v2_form )
   f_primal_objective_dirty = true;

  if( ! f_v2_form && int( Bounds_d.size() ) == NumVars )
   for( int j = 0 ; j < NumVars ; ++j ) {
    const double lhs = ( ! f_L.empty() && std::isfinite( f_L[ j ] ) )
                       ? f_L[ j ] - f_x_bar[ j ] : - Inf< double >();
    const double rhs = ( ! f_U.empty() && std::isfinite( f_U[ j ] ) )
                       ? f_U[ j ] - f_x_bar[ j ] : Inf< double >();
    Bounds_d[ j ].set_lhs( lhs );
    Bounds_d[ j ].set_rhs( rhs );
    }

  for( auto & row : EasyBBFRows ) {
   double coupling = row.b_i;
   for( int j = 0 ; j < NumVars ; ++j )
    coupling += row.A_row[ j ] * x_bar[ j ];
   switch( row.side ) {
    case( BendersBFunction::eLHS ):
     row.cns->set_lhs( coupling , eNoMod );
     break;
    case( BendersBFunction::eRHS ):
     row.cns->set_rhs( coupling , eNoMod );
     break;
    case( BendersBFunction::eBoth ):
     row.cns->set_lhs( coupling , eNoMod );
     row.cns->set_rhs( coupling , eNoMod );
     break;
    }
   }
  return;
  }

 if( z_obj_idx < 0 )
  return;

 auto obj = dynamic_cast< FRealObjective * >( get_objective() );
 if( ! obj )
  return;
 auto dqf = dynamic_cast< DQuadFunction * >( obj->get_function() );
 if( ! dqf )
  return;

 // refresh the quadratic coefficient ±t/2 (or 0 under kLevel) of every z_j.
 // The cuts are stored in the linearization-error frame b = F(x_bar) - alpha
 // - g . x_bar, which has already absorbed the stability centre into sigma:
 // the disaggregated proximal dual then reads max_theta [ - Sigma -
 // (t/2)||z||^2 ] with NO explicit x_bar . z term (the centre lives entirely
 // in d-space), so the linear z coefficient is normally left at 0.
 // EXCEPTIONS carrying an explicit linear z term:
 //  iterate form ( f_v2_form == 1 ): the cuts carry NO g . x_bar at all, so
 //   the master objective reproduces it as the explicit +x_bar^T z cross-term
 //   in the textbook max form, lin_coeff = sgn * x_bar_j;
 //  lazy reference ( displacement, f_xref_tol > 0 ): the cuts bake g . x_ref,
 //   so only the deferred residual rides here, lin_coeff = -sgn (x_bar_j -
 //   x_ref_j) ( = 0 when x_ref tracks x_bar, i.e. tol == 0 -> plain frame ).
 // ( for the convex / eMin master sgn = -1, giving + the gap, matching the
 // plain displacement frame bit-for-bit at the optimum ). The sign convention
 // ±t/2 follows the Objective sense (negated under IsConvex).
 const double sgn = IsConvex ? -1.0 : 1.0;
 const double quad_coeff = ( StblType == kProximal ||
                             StblType == kDoublyStabilized )
                           ? - sgn * t_stab / 2.0 : 0.0;
 const bool iterate = ( f_v2_form != 0 );
 const bool lazy = ( ! iterate ) && ( f_xref_tol > 0.0 ) &&
                   ( int( f_x_ref.size() ) == NumVars );
 for( int j = 0 ; j < NumVars ; ++j ) {
  double lin_coeff = 0.0;
  if( iterate )
   lin_coeff = sgn * f_x_bar[ j ];
  else if( lazy )
   lin_coeff = - sgn * ( f_x_bar[ j ] - f_x_ref[ j ] );
  dqf->modify_term( DQuadFunction::Index( z_obj_idx + j ) ,
                    lin_coeff , quad_coeff );
  }

 // refresh the Lambda-box slack coefficients: in the displacement form
 // s^+_j carries sgn*(L_j - x_bar_j) and s^-_j carries -sgn*(U_j - x_bar_j),
 // both moving with the stability centre; in the iterate form the box stays
 // invariant ( s^+_j = sgn*L_j, s^-_j = -sgn*U_j ) because x_bar lives in the
 // Var_z linear coefficient instead. Entries with a non-finite bound keep
 // their fixed-to-0 slack and 0 coefficient
 if( ( s_plus_obj_idx >= 0 ) && ( s_minus_obj_idx >= 0 ) &&
     ( ! f_L.empty() || ! f_U.empty() ) )
  for( int j = 0 ; j < NumVars ; ++j ) {
   const double xj = ( ( ! iterate ) && j < int( f_x_bar.size() ) )
                     ? f_x_bar[ j ] : 0.0;
   const bool has_L = ! f_L.empty() && std::isfinite( f_L[ j ] );
   dqf->modify_term( DQuadFunction::Index( s_plus_obj_idx + j ) ,
                     has_L ? sgn * ( f_L[ j ] - xj ) : 0.0 , 0.0 );
   const bool has_U = ! f_U.empty() && std::isfinite( f_U[ j ] );
   dqf->modify_term( DQuadFunction::Index( s_minus_obj_idx + j ) ,
                     has_U ? - sgn * ( f_U[ j ] - xj ) : 0.0 , 0.0 );
   }
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_u_bar( const std::vector< double > & )
{
 // skeleton: makes the upper-bundle centre API available so callers can
 // be written against the final interface today, but the kUpperLower
 // stabilization is not yet implemented
 throw( std::logic_error(
      "MasterProblemBlock::set_u_bar: kUpperLower stabilization is not "
      "yet implemented" ) );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_l_bar( const std::vector< double > & )
{
 // skeleton: see set_u_bar; kUpperLower stabilization is not yet
 // implemented
 throw( std::logic_error(
      "MasterProblemBlock::set_l_bar: kUpperLower stabilization is not "
      "yet implemented" ) );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_linear_part( const std::vector< double > & b )
{
 if( int( b.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_linear_part: b must have NumVars entries" ) );

 f_linear_part = b;

 // Under the primal MP the linear part belongs directly to the Objective.
 if( IsPrimal ) {
  set_b( b );
  return;
  }

 // Dual MP: the linear part b of the original sum-function enters the
 // j-th coupling row *scaled by the master multiplier lambda*
 //     z_j - lambda * b_j + s^+_j - s^-_j - sum_{k, i} theta^k_i a^k_{i,j} = 0
 // and is therefore NOT a constant RHS offset. The reason is the master
 // normalization: each hard component's simplex sum_i theta^k_i + gamma^k
 // equals lambda (with K * lambda + r - omega = 1), so the theta-side
 // terms are carried at "mass lambda". The linear part is the gradient of
 // the super-easy 0-th component and must ride at the *same* mass, i.e.
 // lambda * b_j, otherwise it would be over-weighted by a factor K = 1/lambda
 // relative to every f_k and the aggregate z* (= the search direction)
 // would be a wrong convex combination -- which, when b is comparable to
 // the component subgradients, can even flip the sign of z* and send the
 // bundle uphill. The lambda * b_j term lives as the coefficient -b_j on
 // Var_lambda, kept at position 1 of every coupling LinearFunction by
 // CreateDualMP; the RHS stays 0
 int j = 0;
 for( auto & cns : CouplingCns ) {
  auto * lf = static_cast< LinearFunction * >( cns.get_function() );
  lf->modify_coefficient( 1 , - b[ j ] , eModBlck );
  ++j;
  }
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::add_vars( int n )
{
 if( n <= 0 )
  return;

 // Structural change: resizing NumVars while the master problem is live
 // requires a coordinated rewrite of:
 //   - the static-variable groups Var_d / Var_v_hard / Var_z
 //     (the Block internally keeps pointers into the underlying vector,
 //     and a plain resize() would invalidate them);
 //   - the dynamic constraint group CouplingCns (one new row z_j = 0
 //     per added coordinate);
 //   - every PolyhedralFunctionBlock sub-Block (its f_polyf
 //     active-variables list, plus the set_conjugate_constraint
 //     bookkeeping, must mirror the new NumVars);
 //   - the DQuadFunction triples laid out as [z..r..omega] in the
 //     master Objective (the new z_j entries must be inserted before
 //     r_obj_idx and omega_obj_idx shifted accordingly).
 //
 // Until the full plumbing is wired in this is a no-op, which is safe
 // whenever NumVars is stable across the algorithm (the typical case
 // for the bundle pipelines exercised by the in-tree tests). A one-shot
 // warning on std::cerr makes the mismatch visible the first time the
 // method is called with a non-zero count, so silent miscomputation can
 // be diagnosed without crashing the surrounding solver.
 static bool warned = false;
 if( ! warned ) {
  std::cerr << "WARNING: MasterProblemBlock::add_vars( " << n
            << " ): structural resize not implemented yet; the master "
               "will not track new coordinates of the original "
               "sum-function. (This warning is shown once.)"
            << std::endl;
  warned = true;
  }
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::remove_vars( const int * , int sz )
{
 if( sz <= 0 )
  return;
 // see add_vars: structural shrink mirrors the grow path and is left
 // as a no-op with a one-shot warning until the same plumbing lands
 static bool warned = false;
 if( ! warned ) {
  std::cerr << "WARNING: MasterProblemBlock::remove_vars( ..., " << sz
            << " ): structural shrink not implemented yet; the master "
               "will not drop removed coordinates of the original "
               "sum-function. (This warning is shown once.)"
            << std::endl;
  warned = true;
  }
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::forward_log( std::ostream * log_stream )
{
 const auto & solvers = this->get_registered_solvers();
 if( solvers.empty() )
  return;
 solvers.front()->set_log( log_stream );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_max_time( double t )
{
 const auto & solvers = this->get_registered_solvers();
 if( solvers.empty() )
  return;
 auto * slv = solvers.front();
 const auto idx = slv->dbl_par_str2idx( "dblMaxTime" );
 if( idx < slv->get_num_dbl_par() )
  slv->set_par( idx , t );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::refresh_primal_objective( void )
{
 if( ! IsPrimal || ! f_primal_objective_dirty )
  return;

 if( NumVars == 0 ) {
  f_primal_objective_dirty = false;
  return;
  }

 auto * obj = dynamic_cast< FRealObjective * >( get_objective() );
 auto * dqf = obj ? dynamic_cast< DQuadFunction * >( obj->get_function() )
                  : nullptr;
 if( ! dqf )
  throw( std::logic_error(
       "MasterProblemBlock::refresh_primal_objective: expected "
       "DQuadFunction" ) );

 const bool has_quad = ( StblType == kProximal ||
                         StblType == kDoublyStabilized );
 const double quad = has_quad ? 1.0 / ( 2.0 * t_stab ) : 0.0;

 std::vector< double > linear( NumVars );
 std::vector< double > quadratic( NumVars , quad );
 bool changed = false;

 for( int j = 0 ; j < NumVars ; ++j ) {
  linear[ j ] = f_linear_part[ j ];
  if( f_v2_form && has_quad )
   linear[ j ] -= f_x_bar[ j ] / t_stab;

  changed = changed ||
   ( linear[ j ] !=
     dqf->get_linear_coefficient( DQuadFunction::Index( j ) ) ) ||
   ( quadratic[ j ] !=
     dqf->get_quadratic_coefficient( DQuadFunction::Index( j ) ) );
  }

 if( changed )
  dqf->modify_terms( quadratic.cbegin() , linear.cbegin() ,
                     Range( 0 , DQuadFunction::Index( NumVars ) ) );

 f_primal_objective_dirty = false;
}

/*--------------------------------------------------------------------------*/

int MasterProblemBlock::solve_master( void )
{
 const auto & solvers = this->get_registered_solvers();
 if( solvers.empty() )
  throw( std::logic_error(
       "MasterProblemBlock::solve_master: no Solver registered" ) );
 auto * slv = solvers.front();

 // Empty-bundle short-circuit, symmetric for primal and dual.
 //
 // Primal MP with an empty bundle (no cuts pushed yet by add_cut) is
 // intrinsically unbounded below in the v_k variables: the master has no
 // v_k >= a_i.d + b_i row constraining v_k from below, so min sum v_k +
 // (1/(2t))||d||^2 -> -INF.
 //
 // Dual MP with an empty bundle is structurally infeasible under the
 // K * lambda + r - omega = 1 normalization once r and omega are
 // pinned to 0 (their default state when no global LB / level is
 // installed): the per-PFB rows  sum_i theta^k_i + gamma^k = lambda
 // collapse to gamma^k = lambda (theta empty), and any hard component
 // whose PFB has no global bound has gamma^k fixed to 0, forcing
 // lambda = 0 and contradicting K * lambda = 1.
 //
 // The classic Bundle algorithm convention is that at the very first
 // call (bundle empty, no Fi(.) value known yet) no master needs to
 // be solved at all: the surrounding driver will compute
 // Fi(Lambda1 = Lambda = 0) and push the first round of subgradients,
 // then re-enter solve_master with a non-empty bundle. We therefore
 // short-circuit here, returning d = 0 (no movement) and signaling
 // kOK to the caller
 if( is_bundle_empty() ) {
  if( IsPrimal )
   for( int i = 0 ; i < int( Var_d.size() ) ; ++i )
    Var_d[ i ].set_value( f_v2_form ? f_x_bar[ i ] : 0.0 );
  else
   for( auto & zi : Var_z ) zi.set_value( 0.0 );
  return( Solver::kOK );
  }

 refresh_primal_objective();

 const int rc = slv->compute();

 // SMS++ pattern: compute() only writes the solution to the Solver's internal
 // buffers; the ColVariable on the Block stay at their stale values until
 // get_var_solution() is called. Without this push, the driver
 // would read d* / z* / theta as zeros after every master solve
 if( rc == Solver::kOK || rc == Solver::kLowPrecision ) {
  slv->get_var_solution( nullptr );
  // In the primal linearized PFB representation the bundle multipliers are
  // the dual values of the cut constraints, rather than explicit theta
  // variables. Bundle management and aggregation therefore need both sides
  // of the QP solution.
  if( IsPrimal )
   if( auto * cda = dynamic_cast< CDASolver * >( slv ) )
    cda->get_dual_solution( nullptr );
  }

 return( rc );
 }  // end( MasterProblemBlock::solve_master )

/*--------------------------------------------------------------------------*/
/*------------------------- STABILIZATION PARAMETER ------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_t( double t )
{
 if( t <= 0.0 )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_t: t must be strictly positive" ) );

 t_stab = t;

 if( IsPrimal ) {
  f_primal_objective_dirty = true;
  return;
  }

 // From here on only the dual MP is handled. Pure #kNone has no
 // stabilization and pure #kLevel only uses the level row.
 if( StblType == kNone )
  return;
 if( StblType == kLevel )
  return;

 // The z_j quadratic terms occupy the first NumVars entries.
 if( Var_z.empty() )
  return;

 auto obj = dynamic_cast< FRealObjective * >( get_objective() );
 if( ! obj )
  return;
 auto dqf = dynamic_cast< DQuadFunction * >( obj->get_function() );
 if( ! dqf )
  return;

 // Refresh only the quadratic coefficient; the linear coefficient carries
 // the centre-dependent z term installed by set_x_bar().
 const double sgn = IsConvex ? -1.0 : 1.0;
 const double quad_coeff = - sgn * t_stab / 2.0;
 for( int i = 0 ; i < NumVars ; ++i ) {
  const double lin_coeff = dqf->get_linear_coefficient(
                                       DQuadFunction::Index( i ) );
  dqf->modify_term( DQuadFunction::Index( i ) , lin_coeff , quad_coeff );
  }

 }  // end( MasterProblemBlock::set_t )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_f_lev( double f )
{
 f_lev = f;

 if( IsPrimal ) {
  // primal MP: f_lev is the RHS of the level constraint sum_k v^k <= f_lev,
  // present only under #kLevel / #kDoublyStabilized.
  if( StblType == kLevel || StblType == kDoublyStabilized )
   LevelCns.set_rhs( f_lev );
  return;
  }

 // dual MP: f_lev is the linear coefficient on omega in the master Objective;
 // nothing to update under #kProximal (omega is fixed to 0 and is not part
 // of the Objective triples).
 if( omega_obj_idx < 0 )
  return;

 auto obj = dynamic_cast< FRealObjective * >( get_objective() );
 if( ! obj )
  return;
 auto dqf = dynamic_cast< DQuadFunction * >( obj->get_function() );
 if( ! dqf )
  return;

 // : the omega-side contribution to the dual objective is
 //     -omega * Lvl_xbar  ,    Lvl_xbar = Lvl + f_C
 // in the textbook eMax form ; the convex
 // case flips the whole row sign. A non-finite f_lev (= no level set
 // yet) collapses the term to zero rather than leaking Inf into the
 // Objective.
 const bool finite = std::isfinite( f_lev );
 const double sgn = IsConvex ? -1.0 : 1.0;
 const double coeff = finite ? - sgn * ( f_lev + f_C ) : 0.0;
 dqf->modify_term( DQuadFunction::Index( omega_obj_idx ) , coeff , 0.0 );

 // mirror set_global_LB(): Var_omega is meaningful only when f_lev is
 // finite. Unfix it now or re-pin it to 0 otherwise, so the master
 // normalization K * lambda + r - omega = 1 stays bounded
 if( finite ) {
  if( Var_omega.is_fixed() )
   Var_omega.is_fixed( false , eNoMod );
  }
 else {
  Var_omega.set_value( 0 );
  if( ! Var_omega.is_fixed() )
   Var_omega.is_fixed( true , eNoMod );
  }

 }  // end( MasterProblemBlock::set_f_lev )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_b( const std::vector< double > & b )
{
 if( int( b.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_b: b must have NumVars entries" ) );

 f_linear_part = b;

 // Only the primal MP exposes a linear b*d / b*x term: in the dual MP the
 // contribution x_bar*b lives in the driver-managed x_bar*z linear
 // part, which is updated through a different API.
 if( ! IsPrimal || Var_d.empty() )
  return;

 f_primal_objective_dirty = true;

 }  // end( MasterProblemBlock::set_b )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_LB( int k , double LB )
{
 if( k < 0 || k >= NoHardCmps )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_LB: hard-component index out of range" ) );

 if( k >= int( HardCmps.size() ) || ! HardCmps[ k ] )
  return;
 auto pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return;

 if( k < int( f_LB_raw.size() ) )
  f_LB_raw[ k ] = LB;

 auto & poly = pfb->get_PolyhedralFunction();
 const double no_bound = poly.is_convex()
                         ? - Inf< Function::FunctionValue >()
                         :   Inf< Function::FunctionValue >();

 // Primal raw MP: the PFB epigraph variable is the physical value v^k, so its
 // native global bound is exactly the row v^k >= LB^k. The PFB dispatcher
 // propagates modify_bound() to its f_bcv BoxConstraint.
 if( IsPrimal ) {
  poly.modify_bound( std::isfinite( LB ) ? LB : no_bound );
  return;
  }

 // Dual MP: LB^k is the *raw* native lower bound of F_k handed over by
 // the driver. MPB caches it (so subsequent set_reference() calls can
 // refresh the translated form without the driver having to re-call
 // set_LB), then propagates the linearization-error form
 //   LB^k - F_k( x_bar )
 // to the per-PFB gamma * (.) contribution in the dual Objective via
 // PolyhedralFunction::modify_bound(). When LB^k is at -INFshift the
 // bound is genuinely absent and gamma stays fixed to 0
 const double F_xb = ( k < int( f_F_at_x_bar.size() ) )
                     ? f_F_at_x_bar[ k ] : 0.0;

 // PolyhedralFunction::modify_bound() accepts a finite value or the
 // *correct* infinity sentinel (-Inf for convex / max representation,
 // +Inf for concave / min representation): any other infinity is
 // rejected with "wrong INF value to global bound". Route a non-finite
 // LB through the side-appropriate "no bound" sentinel: this collapses
 // the gamma^k * LB^k term in the dual Objective to zero, the same
 // effect as set_global_LB() under !isfinite(LB)
 if( std::isfinite( LB ) )
  // the lower bound LB^k is a *horizontal* cut ( g = 0, constant = LB ): its
  // master-side stored value is therefore the same linearization-error form
  // add_cut() uses for the diagonal cuts ( incoming_b = F_k(x_bar) - alpha +
  // g . x_bar, here with g = 0 and alpha = LB ), i.e. F_k(x_bar) - LB >= 0.
  // The per-PFB dual Objective carries the error with the sign of the master
  // sense: +error in the convex/min representation, -error in the concave/max
  // representation. Thus a slack bound is unattractive in either sense and
  // gamma^k can grow only when the bound is tight ( F_k = LB ), as required by
  // complementary slackness.
  poly.modify_bound( ( IsConvex ? 1.0 : -1.0 ) * ( F_xb - LB ) );
 else {
  poly.modify_bound( no_bound );
  }

 }  // end( MasterProblemBlock::set_LB )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_fictitious_LB( int k , bool on )
{
 if( k < 0 || k >= NoHardCmps )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_fictitious_LB: index out of range" ) );

 // A component k whose bundle is empty models F_k as the max over an
 // empty set of affine pieces, i.e. the improper constant -infinity:
 // the primal master is then unbounded below and the dual master
 // infeasible (the per-PFB simplex row sum_i theta^k_i + gamma^k =
 // lambda collapses to gamma^k = lambda with gamma^k fixed to 0).
 // Following the historical QP-bundle device, the driver installs a
 // *fictitious* model lower bound  v_k >= 0  which makes the empty
 // component "disappear": v_k is pinned to 0, gamma^k becomes a free
 // multiplier with a zero objective coefficient that absorbs the
 // simplex mass lambda, and the master regains feasibility. The
 // bound is *fictitious* in the sense that 0 is expressed directly
 // in the d-space / linearization-error frame (no F_k(x_bar) shift,
 // unlike a genuine native bound routed through set_LB): the model
 // value v_k lives in that frame already. The driver is responsible
 // for removing it (on == false) as soon as the component receives
 // a real cut, restoring the gamma^k-fixed "no bound" state
 // the fictitious model lower bound value: 0 when no genuine global LB
 // exists, otherwise strictly below the translated global LB
 // (global_LB_xbar - 1), so the per-component fictitious constraint
 //   v^k >= fict
 // never becomes more stringent than the genuine aggregate
 //   v   >= global_LB_xbar
 // which would corrupt the master direction (legacy OSIMPSolver device:
 // gamma_i cost = global_LB - 1, cf. OSIMPSolver::SetLowerBound)
 const double fict = std::isfinite( f_global_LB_xbar )
                     ? ( f_global_LB_xbar - 1.0 ) : 0.0;

 if( k >= int( HardCmps.size() ) || ! HardCmps[ k ] )
  return;
 auto pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return;
 auto & poly = pfb->get_PolyhedralFunction();
 const double no_bound = poly.is_convex()
                         ? - Inf< Function::FunctionValue >()
                         :   Inf< Function::FunctionValue >();

 if( IsPrimal ) {
  // In the translated primal frame v_k is already relative to F_k(x_bar), so
  // the empty component disappears with v_k >= 0. In the raw frame v_k is
  // the physical function value instead: use F_k(x_bar), which again gives
  // zero after get_FiBLambda() subtracts the reference.
  const double primal_fict =
   ( f_v2_form && k < int( f_F_at_x_bar.size() ) )
   ? f_F_at_x_bar[ k ] : 0.0;
  poly.modify_bound( on ? primal_fict : no_bound );
  }
 else
  poly.modify_bound( on ? ( IsConvex ? 1.0 : -1.0 ) * fict : no_bound );
 }  // end( MasterProblemBlock::set_fictitious_LB )

/*--------------------------------------------------------------------------*/
/*------------------------- LOAD INTO THE SOLVER ---------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::load_problem( void )
{
 if( this->get_registered_solvers().empty() )
  throw( std::logic_error(
       "MasterProblemBlock::load_problem: no Solver has been registered "
       "to the MasterProblemBlock" ) );

 // Ordinary [MILP]Solver-s do not expose a public load_problem(): they
 // ingest the abstract representation lazily on the first compute(), so
 // here we only verify that at least one Solver has been registered.
 // Should an implementation-specific bridge between the abstract
 // representation and the attached Solver ever be required (e.g. an
 // ad-hoc "compact" embedding of the easy components into the master
 // LP/QP), this is the natural hook to trigger it.

 }  // end( MasterProblemBlock::load_problem )

/*--------------------------------------------------------------------------*/
/*----------------------- End File MasterProblemBlock.cpp ------------------*/
/*--------------------------------------------------------------------------*/
