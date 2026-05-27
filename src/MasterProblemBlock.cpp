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
                          bool convex )
{
 // - - - sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( max_bundle_size < 0 || num_vars < 0 || num_hard_cmps < 0 )
  throw( std::invalid_argument(
       "MasterProblemBlock::configure: negative size" ) );

 const int n_easy  = int( easy_components.size() );
 const int n_total = num_hard_cmps + n_easy;

 // - - - sizes / form / stabilisation - - - - - - - - - - - - - - - - - - -
 // SetDim() starts from a clean slate (it calls clear()) and re-initialises
 // MaxBSize / NumVars / NoTotCmps / NoEasyCmps / NoHardCmps
 SetDim( max_bundle_size , num_vars , n_total , n_easy );

 IsPrimal = primal;
 IsConvex = convex;
 StblType = reg;

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
 for( int k = 0 ; k < n_easy ; ++k ) {
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
   absorb_LBF_into_dual_MP( lbf );
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
 // allocated here and the surrounding BundleSolver then
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
 Var_d.resize( NumVars );          // d is free by default
 if( NumVars > 0 )
  add_static_variable( Var_d , "MPB_d" );

 // ---- one PolyhedralFunctionBlock sub-Block per "hard" component ---------
 //
 // The PFB is wired in its *linearized primal* representation (rep == 1).
 // Each PFB's "x" variables are bound to Var_d (the master step), so each
 // cut v_k >= a_i . d + b_i pushed by add_row hits the right master-side
 // variables. We do NOT call set_lambda() here: that is dual-form specific
 // and would introduce an extra slack in the bundle's normalization; the
 // linearized-primal rep does not have a lambda at all (no dual
 // normalization row).
 //
 // The PolyhedralFunction interior bundle is empty here: the
 // BundleSolver feeds rows (g, alpha) into each f_polyf via
 // the Modification interface as new linearizations are produced

 HardCmps.clear();
 HardCmps.reserve( NoHardCmps );
 const SimpleConfiguration< int > rep_lin_primal( 1 );  // bit 0 = 1, bit 1 = 0

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
 //   min  b * d + sum_k f_v[ k ] + (1/(2t)) || d ||^2_2
 //
 // The proximal quadratic term (1/(2t))||d||^2_2 is present only under
 // kProximal / kDoublyStabilized; under pure kLevel the master is an LP
 // and the d_i contribute only through the linear b*d part.
 //
 // The linear coefficient on d (the "b" of the paper, i.e. the constant
 // gradient of the linear part of the original sum-function) is left at 0
 // here; the BundleSolver is expected to install it through the dedicated
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
 for( int i = 0 ; i < NumVars ; ++i )
  triples.emplace_back( & Var_d[ i ] , 0.0 , quad_coeff );

 FRealObjective * obj;
 if( triples.empty() ) {
  obj = new FRealObjective( this , new LinearFunction() );
  }
 else {
  obj = new FRealObjective( this , new DQuadFunction( std::move( triples ) ) );
  }
 obj->set_sense( Objective::eMin , eNoMod );
 set_objective( obj , eNoMod );

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
 Var_r.is_positive( true , eNoMod );
 if( Stbl == kLevel || Stbl == kDoublyStabilized ) {
  Var_omega.is_positive( true , eNoMod );
  }
 else {
  Var_omega.set_value( 0 );
  Var_omega.is_fixed( true , eNoMod );
  }
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

 // ---- global normalization row: lambda + r - omega = 1 -------------------
 {
  LinearFunction::v_coeff_pair norm_terms;
  norm_terms.reserve( 3 );
  norm_terms.emplace_back( & Var_lambda ,  1.0 );
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
 // BundleSolver knows the linear part b of the original sum-function;
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
 // BundleSolver feeds rows (g, alpha) into each f_polyf via
 // the Modification interface as new linearizations are produced.

 HardCmps.clear();
 HardCmps.reserve( NoHardCmps );
 const SimpleConfiguration< int > rep_dual( 3 );  // bit 0 = 1, bit 1 = 1

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
 // The other linear part x_bar * z is left at 0 here: the BundleSolver
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
 // by the surrounding BundleSolver via register_easy_component(),
 // which adds the easy-cmp sub-Block and augments every CouplingCns[j]
 // with the +A^k_{i,j} u^k_i terms produced by that component.

 // Two coefficients are intentionally left unset by CreateDualMP and must
 // be filled in by the surrounding BundleSolver as soon as
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

void MasterProblemBlock::absorb_LBF_into_dual_MP( LagBFunction * lbf )
{
 if( ! lbf )
  throw( std::invalid_argument(
       "MasterProblemBlock::absorb_LBF_into_dual_MP: null LagBFunction" ) );

 // walk every RowConstraint of the inner Block of the LagBFunction
 // (already transferred under *this* by configure()) and re-install it
 // on the master as the stationarity row of the dual MP at pi^k
 // ():
 //
 //     E^k_i u^k + lambda * e^k_i = 0    for every row i.
 //
 // The original RowConstraint of the inner Block is relaxed in place
 // (LHS = -INF, RHS = +INF) so that the inner :MILPSolver loaded by the
 // dual MP back-end does not enforce it twice.

 auto * inner = lbf->get_inner_block();
 if( ! inner )
  throw( std::invalid_argument(
       "MasterProblemBlock::absorb_LBF_into_dual_MP: null inner Block" ) );

 auto process_one = [ this ]( FRowConstraint & ci ) {
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
  // pairs plus the +e^k_i coupling on Var_lambda
  LinearFunction::v_coeff_pair pairs;
  pairs.reserve( lf->get_num_active_var() + 1 );
  for( Function::Index j = 0 ; j < lf->get_num_active_var() ; ++j )
   pairs.emplace_back(
        static_cast< ColVariable * >( lf->get_active_var( j ) ) ,
        lf->get_coefficient( j ) );
  pairs.emplace_back( & Var_lambda , e_i );

  auto * new_fun = new LinearFunction( std::move( pairs ) );
  auto * new_cns = new FRowConstraint();
  new_cns->set_function( new_fun , eNoMod );
  new_cns->set_lhs( 0.0 , eNoMod );
  new_cns->set_rhs( 0.0 , eNoMod );
  add_static_constraint( *new_cns , "MPB_LBF_easy" );

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
 // or a std::list of them (dynamic only).
 auto walk_any = [ & process_one ]( boost::any & any ) {
  if( auto * one = boost::any_cast< FRowConstraint >( & any ) ) {
   process_one( *one );
   }
  else if( auto * vec = boost::any_cast< std::vector< FRowConstraint > >( & any ) ) {
   for( auto & ci : *vec ) process_one( ci );
   }
  else if( auto * lst = boost::any_cast< std::list< FRowConstraint > >( & any ) ) {
   for( auto & ci : *lst ) process_one( ci );
   }
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
 if( IsPrimal )
  return( 0.0 );

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

 // duplicate detection: walk the bundle of HardCmps[k] and look for an
 // entry that matches (g, alpha) within a relative tolerance. A
 // duplicate brings no new information to the master problem, so the
 // insertion is suppressed and the matching slot is reported back; the
 // surrounding bundle solver consults this return value to refrain
 // from declaring the master as "changed", which in turn lets the
 // noise-reduction machinery kick in when the oracle keeps returning
 // the same subgradient at a frozen Lambda
 const auto & A = pfb->get_PolyhedralFunction().get_A();
 const auto & b = pfb->get_PolyhedralFunction().get_b();
 constexpr double rel_tol = 1e-12;
 const auto cuts_match = [ & ]( std::size_t i ) -> bool {
  if( i >= b.size() || i >= A.size() )
   return( false );
  const double a_tol = rel_tol * std::max( { std::abs( alpha ) ,
                                              std::abs( b[ i ] ) ,
                                              double( 1 ) } );
  if( std::abs( alpha - b[ i ] ) > a_tol )
   return( false );
  const auto & gi = A[ i ];
  const std::size_t n = std::min( g.size() , gi.size() );
  for( std::size_t j = 0 ; j < n ; ++j ) {
   const double g_tol = rel_tol * std::max( { std::abs( g[ j ] ) ,
                                               std::abs( gi[ j ] ) ,
                                               double( 1 ) } );
   if( std::abs( g[ j ] - gi[ j ] ) > g_tol )
    return( false );
   }
  return( true );
  };

 for( int s = 0 ; s < int( slot_to_local[ k ].size() ) ; ++s ) {
  const int local = slot_to_local[ k ][ s ];
  if( local < 0 )
   continue;
  if( cuts_match( std::size_t( local ) ) )
   return( s );  // duplicate: do not touch the bundle, report the slot
  }

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
 const int new_local = int( pfb->get_PolyhedralFunction().get_nrows() );
 pfb->get_PolyhedralFunction().add_row( std::move( g ) , alpha ,
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

namespace {

PolyhedralFunctionBlock * pfb_at( const std::vector< Block * > & HardCmps ,
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

}  // namespace ( anonymous )

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

 pfb->get_PolyhedralFunction().modify_row(
       PolyhedralFunction::Index( slot_to_local[ k ][ slot ] ) ,
       std::move( g ) , alpha );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::modify_alpha( int k , int slot , double alpha )
{
 auto pfb = pfb_at( HardCmps , k , "modify_alpha" );
 if( slot < 0 || slot >= int( slot_to_local[ k ].size() ) ||
     slot_to_local[ k ][ slot ] < 0 )
  throw( std::invalid_argument(
       "MasterProblemBlock::modify_alpha: slot empty or out of range" ) );

 pfb->get_PolyhedralFunction().modify_constant(
       PolyhedralFunction::Index( slot_to_local[ k ][ slot ] ) , alpha );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_theta( int k , int slot ) const
{
 if( IsPrimal )
  return( 0.0 );
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
 const auto thetas = const_cast< PolyhedralFunctionBlock * >( pfb )
                       ->get_dynamic_variable< ColVariable >( "PolyF_theta" );
 if( ! thetas || loc >= int( thetas->size() ) )
  return( 0.0 );

 // f_theta is a std::list, so linear traversal is needed; the bundle
 // sizes are bounded by intBPar2 (typically a few hundred per component)
 auto it = thetas->cbegin();
 std::advance( it , loc );
 return( it->get_value() );
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
  auto v_of = [ this ]( int kk ) -> double {
   if( kk < 0 || kk >= int( HardCmps.size() ) )
    return( 0.0 );
   auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ kk ] );
   if( ! pfb )
    return( 0.0 );
   const auto * v = pfb->get_v();
   return( v ? v->get_value() : 0.0 );
   };
  if( k >= 0 )
   return( v_of( k ) );
  double sum = 0.0;
  for( int kk = 0 ; kk < int( HardCmps.size() ) ; ++kk )
   sum += v_of( kk );
  return( sum );
  }

 // dual MP: by LP duality the per-cmp Objective contribution is
 //  sum_i theta^k_i b^k_i + gamma^k * LB^k
 // which (modulo sign convention) is exactly v*[k]; the gamma^k LB^k
 // piece is currently not separately read out -- if it becomes
 // load-bearing the impl can be extended to add it explicitly.
 return( get_aggregated_alpha( k ) );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > MasterProblemBlock::get_z_vector( void ) const
{
 std::vector< double > out;
 if( IsPrimal ) {
  // primal MP: z* = sum_k sum_i theta^k_i g^k_i across every hard cmp
  out.assign( NumVars , 0.0 );
  for( int k = 0 ; k < int( HardCmps.size() ) ; ++k ) {
   const auto zk = get_aggregated_subgradient( k );
   for( int j = 0 ; j < NumVars && j < int( zk.size() ) ; ++j )
    out[ j ] += zk[ j ];
   }
  return( out );
  }
 out.reserve( Var_z.size() );
 for( const auto & zj : Var_z )
  out.push_back( zj.get_value() );
 return( out );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > MasterProblemBlock::get_aggregated_subgradient( int k ) const
{
 std::vector< double > out( NumVars , 0.0 );

 if( k < 0 || k >= int( HardCmps.size() ) )
  return( out );
 const auto pfb = dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( out );

 // The aggregated subgradient is the convex combination
 //
 //   z^k = sum_i theta^k_i * g^k_i
 //
 // where g^k_i is row i of f_polyf.A and theta^k_i is the multiplier
 // attached to the i-th bundle row. The dual MP and the linearized
 // primal MP store theta in different places:
 //
 //  - in the dual MP, the PolyhedralFunctionBlock allocates a
 //    "PolyF_theta" dynamic ColVariable group (one entry per bundle row)
 //    and the master Solver writes the optimal multipliers there;
 //  - in the linearized primal MP, theta is the dual multiplier of the
 //    bundle row v_k >= g^k_i . d + b^k_i carried by the f_const
 //    dynamic FRowConstraint list, and is recovered via
 //    RowConstraint::get_dual().
 //
 // Both branches share the same g^k_i source (f_polyf.A)
 const auto & A = const_cast< PolyhedralFunctionBlock * >( pfb )
                       ->get_PolyhedralFunction().get_A();

 auto add_row = [ & ]( std::size_t i , double th ) {
  if( th == 0.0 )
   return;
  if( i >= A.size() )
   return;
  const auto & gi = A[ i ];
  const std::size_t n = std::min( gi.size() , out.size() );
  for( std::size_t j = 0 ; j < n ; ++j )
   out[ j ] += th * gi[ j ];
  };

 if( IsPrimal ) {
  const auto cgroup = const_cast< PolyhedralFunctionBlock * >( pfb )
                        ->get_dynamic_constraint< FRowConstraint >( "" );
  if( ! cgroup )
   return( out );
  std::size_t i = 0;
  for( const auto & ci : *cgroup ) {
   add_row( i , ci.get_dual() );
   ++i;
   }
  return( out );
  }

 const auto thetas = const_cast< PolyhedralFunctionBlock * >( pfb )
                       ->get_dynamic_variable< ColVariable >( "PolyF_theta" );
 if( ! thetas )
  return( out );
 if( thetas->size() != A.size() )
  return( out );
 auto it = thetas->cbegin();
 for( std::size_t i = 0 ; i < A.size() ; ++i , ++it )
  add_row( i , it->get_value() );
 return( out );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > MasterProblemBlock::get_d_vector( void ) const
{
 std::vector< double > out;

 if( IsPrimal ) {
  out.reserve( Var_d.size() );
  for( const auto & di : Var_d )
   out.push_back( di.get_value() );
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
  if( z.empty() || Var_d.empty() )
   return( 0.0 );
  const std::size_t n = std::min( z.size() , Var_d.size() );
  double s = 0.0;
  for( std::size_t j = 0 ; j < n ; ++j )
   s += z[ j ] * Var_d[ j ].get_value();
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
 return( std::inner_product( g.begin() , g.end() , d.begin() , 0.0 ) );
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
 if( IsPrimal )
  return( 0.0 );

 auto contrib = [ this ]( int kk ) -> double {
  if( kk < 0 || kk >= int( HardCmps.size() ) )
   return( 0.0 );
  const auto pfb = dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ kk ] );
  if( ! pfb )
   return( 0.0 );
  auto * th = const_cast< PolyhedralFunctionBlock * >( pfb )
                ->get_dynamic_variable< ColVariable >( "PolyF_theta" );
  if( ! th )
   return( 0.0 );
  const auto & b = const_cast< PolyhedralFunctionBlock * >( pfb )
                       ->get_PolyhedralFunction().get_b();
  if( th->size() != b.size() )
   return( 0.0 );
  double s = 0.0;
  std::size_t i = 0;
  for( const auto & v : *th ) {
   s += v.get_value() * b[ i ];
   ++i;
   }
  return( s );
  };

 if( k >= 0 )
  return( contrib( k ) );

 double total = 0.0;
 for( int kk = 0 ; kk < int( HardCmps.size() ) ; ++kk )
  total += contrib( kk );
 return( total );
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
 if( IsPrimal )
  return( out );
 if( k < 0 || k >= int( HardCmps.size() ) )
  return( out );
 const auto pfb = dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return( out );

 const auto thetas = const_cast< PolyhedralFunctionBlock * >( pfb )
                       ->get_dynamic_variable< ColVariable >( "PolyF_theta" );
 if( ! thetas )
  return( out );

 out.reserve( thetas->size() );
 for( const auto & v : *thetas )
  out.push_back( v.get_value() );
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
 if( IsPrimal )
  return;
 if( ! L.empty() && int( L.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: L must be empty or of size NumVars" ) );
 if( ! U.empty() && int( U.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: U must be empty or of size NumVars" ) );

 f_L = L;
 f_U = U;

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
 const double sgn = IsConvex ? -1.0 : 1.0;
 const double coeff = std::isfinite( LB ) ? sgn * ( LB + f_C ) : 0.0;
 dqf->modify_term( DQuadFunction::Index( r_obj_idx ) , coeff , 0.0 );
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_x_bar( const std::vector< double > & x_bar )
{
 if( int( x_bar.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_x_bar: x_bar must have NumVars entries" ) );

 // refresh the absorbed BendersBFunction RowConstraints (primal MP):
 // their right-hand side(s) are  A_i . x_bar + b_i  for the side(s)
 // marked dynamic by the original ConstraintSide; the static side(s)
 // keep their snapshot from absorb_BBF_into_primal_MP()
 if( IsPrimal ) {
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

 // refresh only the linear coefficient on every z_j; the quadratic
 // coefficient ±t/2 (or 0 under kLevel) is left consistent with the
 // sense chosen at CreateDualMP time. The +x̄·z textbook contribution
 // is negated under IsConvex (whole Objective negated).
 const double sgn = IsConvex ? -1.0 : 1.0;
 const double quad_coeff = ( StblType == kProximal ||
                             StblType == kDoublyStabilized )
                           ? - sgn * t_stab / 2.0 : 0.0;
 for( int j = 0 ; j < NumVars ; ++j )
  dqf->modify_term( DQuadFunction::Index( z_obj_idx + j ) ,
                    sgn * x_bar[ j ] , quad_coeff );
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

 // under the primal MP the linear part is folded directly into the master
 // Objective by set_b and there is nothing to install in the coupling rows
 if( IsPrimal )
  return;

 // dual MP: CreateDualMP wires every CouplingCns[ j ] with Var_lambda
 // at position 1 of its LinearFunction, initial coefficient 0. The
 // linear part of the original sum-function enters the equation as the
 // term  -lambda * b_j, so we refresh that coefficient (and only that
 // one) here. The change is
 // broadcast to every attached Solver via eModBlck.
 constexpr int LAMBDA_POS = 1;
 int j = 0;
 for( auto & cns : CouplingCns ) {
  auto * lf = static_cast< LinearFunction * >( cns.get_function() );
  if( lf )
   lf->modify_coefficient( LAMBDA_POS , - b[ j ] , eModBlck );
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

int MasterProblemBlock::solve_master( void )
{
 const auto & solvers = this->get_registered_solvers();
 if( solvers.empty() )
  throw( std::logic_error(
       "MasterProblemBlock::solve_master: no Solver registered" ) );
 auto * slv = solvers.front();

 // Primal MP with an empty bundle (no cuts pushed yet by add_cut) is
 // intrinsically unbounded below in the v_k variables: the master has no
 // v_k >= a_i.d + b_i row constraining v_k from below, so min sum v_k +
 // (1/(2t))||d||^2 -> -INF. The classic Bundle algorithm convention is
 // that at the very first call (bundle empty, no Fi(.) value known yet)
 // no master needs to be solved at all: the surrounding BundleSolver
 // will compute Fi(Lambda1 = Lambda = 0) and push the first
 // round of subgradients, then re-enter solve_master with a non-empty
 // bundle. We therefore short-circuit here, returning d = 0 (no movement)
 // and signaling kOK to the caller.
 //
 // The dual MP does not have this issue: an empty bundle is feasible
 // (theta empty, gamma fixed, lambda absorbs the unit) and the solver
 // returns a well-defined zero solution
 if( IsPrimal && is_bundle_empty() ) {
  for( auto & di : Var_d ) di.set_value( 0.0 );
  return( Solver::kOK );
  }

 const int rc = slv->compute();

 // SMS++ pattern: compute() only writes the solution to the Solver's internal
 // buffers; the ColVariable on the Block stay at their stale values until
 // get_var_solution() is called. Without this push, the BundleSolver
 // would read d* / z* / theta as zeros after every master solve
 if( rc == Solver::kOK || rc == Solver::kLowPrecision )
  slv->get_var_solution( nullptr );

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

 // Pure #kNone has no stabilization at all and pure #kLevel (in the
 // dual MP) only uses the level row; in both cases the Objective is a
 // plain LinearFunction (or a DQuadFunction with no quadratic d/z
 // entries), so there is nothing to refresh.
 if( StblType == kNone )
  return;
 if( ! IsPrimal && StblType == kLevel )
  return;

 // The diagonal quadratic coefficient lives in the DQuadFunction wrapped
 // by the FRealObjective set by CreatePrimal/DualMP. Both layouts place
 // the quadratic terms on the first NumVars triple entries:
 //   - primal:  d_i with coefficient  +1/(2t)
 //   - dual:    z_j with coefficient  -t/2
 // so the refresh is a straight loop over [0, NumVars).
 if( ( IsPrimal && Var_d.empty() ) || ( ! IsPrimal && Var_z.empty() ) )
  return;

 auto obj = dynamic_cast< FRealObjective * >( get_objective() );
 if( ! obj )
  return;
 auto dqf = dynamic_cast< DQuadFunction * >( obj->get_function() );
 if( ! dqf )
  return;

 // refresh only the quadratic coefficient on every d_i / z_j; the linear
 // coefficient carries the b*d term (primal, via set_b) or the ±x_bar*z
 // term (dual, via set_x_bar) and must be preserved as-is. The primal is
 // in min form (quad +1/(2t) on d); the dual carries quad -t/2 on z in
 // textbook max form, negated under IsConvex (whole Objective in min).
 const double sgn = IsConvex ? -1.0 : 1.0;
 const double quad_coeff = IsPrimal ? 1.0 / ( 2.0 * t_stab )
                                    : - sgn * t_stab / 2.0;
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
 const double sgn = IsConvex ? -1.0 : 1.0;
 const double coeff = std::isfinite( f_lev )
                       ? - sgn * ( f_lev + f_C ) : 0.0;
 dqf->modify_term( DQuadFunction::Index( omega_obj_idx ) , coeff , 0.0 );

 }  // end( MasterProblemBlock::set_f_lev )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_b( const std::vector< double > & b )
{
 if( int( b.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_b: b must have NumVars entries" ) );

 // Only the primal MP exposes a linear b*d term: in the dual MP the
 // contribution x_bar*b lives in the BundleSolver-managed x_bar*z linear
 // part, which is updated through a different API.
 if( ! IsPrimal || Var_d.empty() )
  return;

 auto obj = dynamic_cast< FRealObjective * >( get_objective() );
 if( ! obj )
  return;
 auto dqf = dynamic_cast< DQuadFunction * >( obj->get_function() );
 if( ! dqf )
  return;

 // d_i sit at the first NumVars triple entries of the DQuadFunction
 const double quad_coeff = ( StblType == kProximal ||
                             StblType == kDoublyStabilized )
                           ? 1.0 / ( 2.0 * t_stab ) : 0.0;
 for( int i = 0 ; i < NumVars ; ++i )
  dqf->modify_term( DQuadFunction::Index( i ) , b[ i ] , quad_coeff );

 }  // end( MasterProblemBlock::set_b )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_LB( int k , double LB )
{
 if( k < 0 || k >= NoHardCmps )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_LB: hard-component index out of range" ) );

 // Primal MP: v^k >= LB^k via BoxConstraint LHS.
 if( IsPrimal ) {
  if( int( Bounds_v_hard.size() ) != NoHardCmps )
   throw( std::logic_error(
        "MasterProblemBlock::set_LB: primal MP not built or bounds "
        "group has the wrong size" ) );
  Bounds_v_hard[ k ].set_lhs( LB );
  return;
  }

 // Dual MP: LB^k is the global lower bound of the underlying f_polyf of
 // HardCmps[k]; it propagates to the per-PFB gamma * LB^k contribution in
 // the dual Objective via PolyhedralFunction::modify_bound().
 if( k >= int( HardCmps.size() ) || ! HardCmps[ k ] )
  return;
 auto pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
 if( ! pfb )
  return;
 pfb->get_PolyhedralFunction().modify_bound( LB );

 }  // end( MasterProblemBlock::set_LB )

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
