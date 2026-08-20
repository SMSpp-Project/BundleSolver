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

#include <cstdlib>
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
/*-------------------- ABSTRACT REPRESENTATION FLAGS -----------------------*/
/*--------------------------------------------------------------------------*/

// Same convention used by PolyhedralFunctionBlock: the low bits are left free
// for representation choices, while these bits track which generated slices
// are already available.
static constexpr char k_mpb_built_var  = 0x04;
static constexpr char k_mpb_built_cnst = 0x08;
static constexpr char k_mpb_built_obj  = 0x10;
static constexpr char k_mpb_built_mask =
 k_mpb_built_var | k_mpb_built_cnst | k_mpb_built_obj;

static char mpb_built_stage( Configuration * cfg )
{
 auto * scfg = dynamic_cast< SimpleConfiguration< int > * >( cfg );
 return( scfg ? char( scfg->f_value & k_mpb_built_mask ) : 0 );
}

/*--------------------------------------------------------------------------*/
/*----------------------- CLEAR / REINITIALIZE -----------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::clear()
{
 // forget any per-component sub-Block (the sub-Block objects themselves
 // are owned by the base Block, which will dispose of them in due time)
 EasyCmps.clear();
 EasyCmps_SB.clear();
 EasyObjVars.clear();
 EasyObjCoeffs.clear();
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
 EasyLocal2Global.clear();

 // back to the "default" MP type
 IsPrimal      = false;
 StblType      = kDoublyStabilized;
 HardCmpScaling = 0;
 t_stab        = 1.0;
 f_lev         = 0.0;
 z_obj_idx     = -1;
 r_obj_idx     = -1;
 omega_obj_idx = -1;
 easy_obj_idx  = -1;
 level_model_obj_idx = -1;
 f_dual_level_probe_active = false;
 f_abs_rep = 0;

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
 EasyCmps_SB.reserve( NoEasyCmps );
 EasyLocal2Global.reserve( NoEasyCmps );
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
                          int hard_cmp_scaling ,
                          const std::vector< std::vector< Index > > & easy_local2global )
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
 if( ! is_easy.empty() && int( is_easy.size() ) == n_total ){
  IsEasyCmp = is_easy;
  EasyLocal2Global = easy_local2global;
 }

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
 // inside register_Solver) TODO
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

   /* IMPORTANT NOTE: In the dual version, the per-row stationarity
    * constraints of each easy components would read:
    *       
    *       E^k_i u^k + lambda * e^k_i = 0
    * 
    * However, in the current implementation we simply register the inner
    * block of the LagBFunction to *this, directly importing the constraints
    * 
    *       E^k_i u^k + e^k_i = 0
    * 
    * Adding the contribution of \lambda would require "hacking" the internal
    * representation of LagBFunction, hence contradicting the general idea 
    * of SMS++. For this reason, if easy components are available in the dual
    * master problem, we simply force \lambda = 1, making the two formulations
    * equivalent. Hopefully, this will be addressed in the future with some
    * copy or scaling mechanism. */

   inner->transfer_ownership_to( this );
   EasyCmps.push_back( lbf );
   EasyCmps_SB.push_back( inner );
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
   EasyCmps_SB.push_back( inner );
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
/*------------------ Abstract representation methods------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::generate_abstract_variables( Configuration * stvv )
{
 Configuration * cfg = stvv;
 if( ( ! cfg ) && f_BlockConfig )
  cfg = f_BlockConfig->f_static_variables_Configuration;
 f_abs_rep |= mpb_built_stage( cfg );

 if( f_abs_rep & k_mpb_built_var )
  return;

 if( IsPrimal )
  generate_primal_abstract_variables();
 else
  generate_dual_abstract_variables();

 f_abs_rep |= k_mpb_built_var;
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::generate_abstract_constraints( Configuration * stcc )
{
 Configuration * cfg = stcc;
 if( ( ! cfg ) && f_BlockConfig )
  cfg = f_BlockConfig->f_static_constraints_Configuration;
 f_abs_rep |= mpb_built_stage( cfg );

 if( f_abs_rep & k_mpb_built_cnst )
  return;

 if( ! ( f_abs_rep & k_mpb_built_var ) )
  throw( std::logic_error(
       "MasterProblemBlock::generate_abstract_constraints: variables "
       "must be generated first" ) );

 if( IsPrimal )
  generate_primal_abstract_constraints();
 else
  generate_dual_abstract_constraints();

 f_abs_rep |= k_mpb_built_cnst;
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::generate_objective( Configuration * objc )
{
 Configuration * cfg = objc;
 if( ( ! cfg ) && f_BlockConfig )
  cfg = f_BlockConfig->f_objective_Configuration;
 f_abs_rep |= mpb_built_stage( cfg );

 if( f_abs_rep & k_mpb_built_obj )
  return;

 if( ! ( f_abs_rep & k_mpb_built_var ) )
  throw( std::logic_error(
       "MasterProblemBlock::generate_objective: variables must be "
       "generated first" ) );

 if( IsPrimal )
  generate_primal_objective();
 else
  generate_dual_objective();

 f_abs_rep |= k_mpb_built_obj;
}

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIMAL MP ---------------------------------*/
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
 f_dual_level_probe_active = false;

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

 // PFB stvv bit 4 disables the objective function. This is needed in the
 // pure level case because the model terms are managed by the root
 // objective during the one-shot proximal probe.
 const int built_obj_cfg = ( Stbl == kLevel ) ? 0x10 : 0;

 const SimpleConfiguration< int > rep_lin_primal( 1 | scaling_cfg |
                                   built_obj_cfg );
 for( int k = 0 ; k < NoHardCmps ; ++k ) {
  auto * pfb = new PolyhedralFunctionBlock( this );

  /* Even tough in this phase we are creating the physical representation 
   * of the MP, we already call the method for generating the abstract variables 
   * of each PFB. This is needed because currently SMS++ writes solution 
   * vectors inside the Block variables, which therefore should be initialized.
   * Hopefully, this will go away someday. */
  pfb->generate_abstract_variables(
               const_cast< SimpleConfiguration< int > * >( & rep_lin_primal ) );

  HardCmps.push_back( pfb );
  add_nested_Block( pfb );
  }

 // bookkeeping fields: Bounds_v_hard / Var_v_hard are *unused* in this
 // refactor; we keep them clean so the dual MP code path (which still
 // owns them) stays unaffected
 Bounds_v_hard.clear();
 Var_v_hard.clear();

 }  // end( MasterProblemBlock::CreatePrimalMP )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::generate_primal_abstract_variables( void )
{
 // ---- static Variable: the step d ---------------------------------------
 //
 // NOTE: The epigraph variable v^k of each hard component is NOT a static
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

 Var_d.clear();
 Var_d.resize( NumVars );          // x (raw form) or d (translated form)
 if( NumVars > 0 )
  add_static_variable( Var_d , f_v2_form ? "MPB_x" : "MPB_d" );


 // Now handle each hard component
 for( int k = 0 ; k < NoHardCmps ; ++k ) {
  auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );

  // wire the "x" variables of f_polyf to the master-side Var_d, so every
  // cut v_k >= a . d + b pushed via add_row lands on the right columns
  PolyhedralFunction::VarVector vv;
  vv.reserve( NumVars );
  for( auto & di : Var_d )
   vv.push_back( & di );
  pfb->get_PolyhedralFunction().set_variables( std::move( vv ) );

  // abstract variables have already been generated in the CreatePrimal method
 }

}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::generate_primal_abstract_constraints( void )
{
 // Initialize the possible bounds on the step d
 Bounds_d.clear();
 Bounds_d.resize( NumVars );
 for( int j = 0 ; j < NumVars ; ++j ) {
  Bounds_d[ j ].set_variable( & Var_d[ j ] , eNoMod );
  Bounds_d[ j ].set_lhs( - Inf< double >() , eNoMod );
  Bounds_d[ j ].set_rhs( Inf< double >() , eNoMod );
  }
 if( NumVars > 0 )
  add_static_constraint( Bounds_d , "MPB_box" );

 // Now handle each hard component
 for( int k = 0 ; k < NoHardCmps ; ++k ) {
  auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
  
  // ensure f_polyf has NO global lower bound (for the default convex
  // case, f_bound defaults to +INF which would translate the linearized
  // primal box-constraint f_bcv to "+INF <= f_v <= +INF" -> infeasible
  // master; we want f_v to be free). modify_bound is called with eNoMod
  // because we are still in setup phase and no Solver is attached yet
  pfb->get_PolyhedralFunction().modify_bound(
                  - Inf< Function::FunctionValue >() , eNoMod );

  // Now generate the abstract constraints of each PFB
  pfb->generate_abstract_constraints();
 }

 // ---- level constraint b*d + sum_k v^k <= f_lev -------------------------
 //
 // The first NumVars coefficients are reserved for the linear component b;
 // set_linear_part() refreshes them after construction. The remaining terms
 // reference each PFB's f_v directly via get_v(). LevelCns is a master-side
 // static FRowConstraint so the inner Solver picks it up.

 if( ( StblType == kLevel || StblType == kDoublyStabilized ) 
          && NoHardCmps > 0 ) {
  LinearFunction::v_coeff_pair lvl_terms;
  lvl_terms.reserve( NumVars + NoHardCmps );
  for( int j = 0 ; j < NumVars ; ++j ) {
   const double coeff = ( j < int( f_linear_part.size() ) )
                        ? f_linear_part[ j ] : 0.0;
   lvl_terms.emplace_back( & Var_d[ j ] , coeff );
   }

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

 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::generate_primal_objective( void )
{
 // Generate the objective of each hard component
 for( int k = 0 ; k < NoHardCmps ; ++k ) {
  auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );

  pfb->generate_objective();
 }

 // ---- Objective ---------------------------------------------------------
 //
 // kProximal / kDoublyStabilized:
 //   translated: min b*d + sum_k f_v[k] + (1/(2t)) ||d||^2
 //   raw:        min b*x + sum_k f_v[k] + (1/(2t)) ||x-x_bar||^2
 //
 // kLevel:
 //   first non-empty solve:
 //     translated: min b*d + sum_k f_v[k] + (1/(2t)) ||d||^2
 //     raw:        min b*x + sum_k f_v[k] + (1/(2t)) ||x-x_bar||^2
 //   later solves:
 //     translated: min (1/2) ||d||^2
 //     raw:        min (1/2) ||x-x_bar||^2
 //   subject to the level row b*d + sum_k f_v[k] <= f_lev
 //
 // Each hard PolyhedralFunctionBlock contributes +f_v[k] through its own
 // Objective in the SMS++ objective scan, except in pure kLevel where that
 // nested Objective is suppressed. The one-shot proximal probe therefore
 // places the f_v[k] terms in the root Objective, then removes them.

 const bool pure_level = ( StblType == kLevel );
 const bool level_probe = pure_level;
 const bool has_quad =
  pure_level || ( StblType == kProximal ) || ( StblType == kDoublyStabilized );

 DQuadFunction::v_coeff_triple triples;
 triples.reserve( NumVars + ( level_probe ? NoHardCmps : 0 ) );

 const double quad_coeff = level_probe ? 1.0 / ( 2.0 * t_stab )
                                      : pure_level ? 0.5
                                      : ( has_quad ? 1.0 / ( 2.0 * t_stab )
                                                   : 0.0 );
 for( int i = 0 ; i < NumVars ; ++i ) {
  const double lin_coeff = level_probe
                           ? ( f_linear_part[ i ] +
                               ( f_v2_form ? - f_x_bar[ i ] / t_stab : 0.0 ) )
                           : pure_level ? ( f_v2_form ? - f_x_bar[ i ] : 0.0 )
                           : ( ( f_v2_form && has_quad )
                               ? - f_x_bar[ i ] / t_stab : 0.0 );
  triples.emplace_back( & Var_d[ i ] , lin_coeff , quad_coeff );
  }

 level_model_obj_idx = -1;
 if( level_probe ) {
  level_model_obj_idx = int( triples.size() );
  for( int k = 0 ; k < NoHardCmps ; ++k ) {
   auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
   if( pfb )
    triples.emplace_back( pfb->get_v() , 1.0 , 0.0 );
   }
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
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DUAL MP -----------------------------------*/
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
 // Pure level starts with the same proximal seed used by the primal form. The
 // flag is cleared by remove_initial_level_objective() once a level is known.
 f_dual_level_probe_active = ( Stbl == kLevel );

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

  /* Even tough in this phase we are creating the physical representation 
   * of the MP, we already call the method for generating the abstract variables 
   * of each PFB. This is needed because currently SMS++ writes solution 
   * vectors inside the Block variables, which therefore should be initialized.
   * Hopefully, this will go away someday. */
  pfb->generate_abstract_variables(
               const_cast< SimpleConfiguration< int > * >( & rep_dual ) );

  HardCmps.push_back( pfb );
  add_nested_Block( pfb );
  }

 }  // end( MasterProblemBlock::CreateDualMP )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::generate_dual_abstract_variables( void )
{
 // ---- static dual multipliers: lambda (global), r, omega -----------------
 //
 // lambda is the *single* global non-negative multiplier paired with the
 // model-value equation
 //
 //     v = d b + sum_E pi^k e^k + sum_H v_x^k
 //
 // of the lower model and the stationarity condition (ii) at 
 //     
 //     v_x^k: lambda = sum_i theta^k_i + gamma^k for every k in H. 
 //
 // Hence the same lambda enters the simplex (= normalization) row of *every* 
 // hard component PFB sub-Block with coefficient +1 
 // (cf. PolyhedralFunctionBlock::set_lambda, called below with a single 
 // shared pointer), so each per-PFB row reads
 //
 //     sum_i theta^k_i + gamma^k = lambda.
 //
 // r is the non-negative multiplier of the global lower bound v >= LB;
 // omega is the non-negative multiplier of the level row Lvl >= v,
 // live only under #kLevel / #kDoublyStabilized. The master-side
 // stationarity at v then reads
 //
 //     lambda + r - omega = 1      (proximal / doubly stabilized)
 //     lambda + r - omega = 0      (pure level)
 //
 // The corresponding row becomes the single static "normalization" constraint
 // installed below; the per-PFB rows are owned by the PFB sub-Blocks themselves.
 //
 // NOTE: when easy components are considered, \lambda is fixed to 1 (see 
 // MasterProblemBlock.353 for further details). This means that the above
 // equation reads
 //
 //      r = omega.
 //
 // Therefore the global lower-bound multiplier r and the level multiplier 
 // omega can only appear in a perfectly balanced way. In particular, if 
 // omega is absent or fixed to zero, then r is forced to zero as well, 
 // so the global lower bound cannot contribute through its dual multiplier.

 if( NoEasyCmps > 0 ){
  Var_lambda.set_value( 1 );
  Var_lambda.is_fixed( true , eNoMod );
 }
 else
  Var_lambda.is_positive( true , eNoMod );

 // Var_r is the dual multiplier of the global LB row. It is structurally
 // present in every stabilization type, but only carries an objective
 // contribution when set_global_LB() is called with a finite LB; until
 // then we pin it to 0 so the proximal/doubly stabilized master normalization
 //     lambda + r - omega = 1
 // collapses to  lambda = 1 , and the per-PFB
 //     sum_i theta^k_i + gamma^k = lambda
 // gives every hard component its full convexity mass. Leaving r free
 // with a 0 objective coefficient would allow it to pick any value in
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
 // pin, the lambda + r - omega normalization would let omega drift to
 // +infinity, pulling lambda along the same direction, with the same
 // gamma^k unboundedness consequence described above for Var_r
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
 // The box belongs to the physical representation and may already have been
 // stored in f_L / f_U before abstract variables are generated, so initialize
 // the fixed status from that state. Missing sides stay fixed to 0.
 Var_s_plus.clear();
 Var_s_minus.clear();
 Var_s_plus.resize( NumVars );
 Var_s_minus.resize( NumVars );
 for( int j = 0 ; j < NumVars ; ++j ) {
  const bool has_L = ! f_L.empty() && std::isfinite( f_L[ j ] );
  const bool has_U = ! f_U.empty() && std::isfinite( f_U[ j ] );

  Var_s_plus[ j ].is_positive( true , eNoMod );
  Var_s_plus[ j ].set_value( 0.0 );
  Var_s_plus[ j ].is_fixed( ! has_L , eNoMod );
  Var_s_minus[ j ].is_positive( true , eNoMod );
  Var_s_minus[ j ].set_value( 0.0 );
  Var_s_minus[ j ].is_fixed( ! has_U , eNoMod );
  }
 if( NumVars > 0 ) {
  add_static_variable( Var_s_plus  , "MPB_s_plus"  );
  add_static_variable( Var_s_minus , "MPB_s_minus" );
  }

 // Now handle each hard component
 for( int k = 0 ; k < NoHardCmps ; ++k ) {
  auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );

  // wire the "x" variables of f_polyf to the master-side Var_d, so every
  // cut v_k >= a . d + b pushed via add_row lands on the right columns
  PolyhedralFunction::VarVector vv;
  vv.reserve( NumVars );
  for( auto & zj : Var_z )
   vv.push_back( & zj );
  pfb->get_PolyhedralFunction().set_variables( std::move( vv ) );

  // abstract variables have already been generated in the CreatePrimal method
 }

}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::generate_dual_abstract_constraints( void )
{
 // ---- global normalization row: lambda + r - omega = rhs ---------------
 // The disaggregated proximal dual stationarity in v_k (the per-component
 // model decrease) reads
 //     dL/dv_k = 1 - sum_i theta^k_i - gamma^k - r + omega = 0
 // hence per component sum_i theta^k_i + gamma^k = 1 - r + omega =: lambda.
 // In pure level there is no +v_k term in the primal objective: the level
 // multiplier omega is the only multiplier of the model-value row, so the
 // same stationarity loses the constant 1 and gives
 //     sum_i theta^k_i + gamma^k = omega - r =: lambda.
 // The per-PFB simplex row enforces sum_i theta^k_i + gamma^k = lambda for
 // every hard k, so this global row only has to pin lambda + r - omega = rhs.
 // The one-shot level probe is proximal and keeps rhs = 1; true pure level
 // switches to rhs = 0 when remove_initial_level_objective() is called
 // (coefficient 1 on lambda, NOT NoHardCmps): the dual aggregates the
 // *sum* of the K components f = sum_k f_k, and each one carries unit
 // convexity mass independently. A K factor here would instead force
 // lambda = 1 / K when r = omega = 0, i.e. it would average the components
 // rather than sum them, shrinking the aggregate subgradient z = sum theta g
 // by 1 / K and the proximal step d = - t z with it.
 {
  const bool pure_level =
   ( StblType == kLevel ) && ( ! has_initial_level_objective() );
  const double norm_rhs = pure_level ? 0.0 : 1.0;
  LinearFunction::v_coeff_pair norm_terms;
  norm_terms.reserve( 3 );
  norm_terms.emplace_back( & Var_lambda , 1.0 );
  norm_terms.emplace_back( & Var_r      ,  1.0 );
  norm_terms.emplace_back( & Var_omega  , -1.0 );

  NormalizationCns.set_lhs( norm_rhs , eNoMod );
  NormalizationCns.set_rhs( norm_rhs , eNoMod );
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
 // u^k A^k terms (easy components in the dual MP) must be set explicity
 // using the mapping provided by the LagBFunction class.
 // CouplingCns is exposed as a *dynamic* group because
 // PolyhedralFunctionBlock::set_conjugate_constraint takes
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
  // Prepare a vector that will contain the pairs of (var,coeff) for each coupling
  // constraint.
  std::vector< LinearFunction::v_coeff_pair > vp_Cns( NumVars );
  for( int j = 0 ; j < NumVars ; ++j ) {
   LinearFunction::v_coeff_pair vp;
   vp.reserve( 4 );
   vp.emplace_back( & Var_z[ j ]        ,  1.0 );  // pos 0
   vp.emplace_back( & Var_lambda        ,  0.0 );  // pos 1 (set by set_linear_part)
   vp.emplace_back( & Var_s_plus[ j ]   ,  1.0 );  // pos 2
   vp.emplace_back( & Var_s_minus[ j ]  , -1.0 );  // pos 3

   vp_Cns[ j ] = vp;
  }

  // Now append easy-component terms to the vectors
  if( NoEasyCmps > 0 )
   add_LBF_to_coupling_rows( vp_Cns );

  // Finally, initialize the Linear function with each prepared list
  auto it = CouplingCns.begin();
  for( int j = 0 ; j < NumVars ; ++j , ++it ) {
   it->set_function( new LinearFunction( std::move( vp_Cns[ j ] ) , 0.0 ) , eNoMod );
   it->set_lhs( 0.0 , eNoMod );
   it->set_rhs( 0.0 , eNoMod );
  }
 }

 if( NumVars > 0 )
  add_dynamic_constraint( CouplingCns , "MPB_coupling" );

 // Now handle each hard component
 for( int k = 0 ; k < NoHardCmps ; ++k ) {
  auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );

  // Now generate the abstract constraints of each PFB
  pfb->generate_abstract_constraints();

  // inject the lambda multiplier into the normalization constraint
  pfb->set_lambda( & Var_lambda );

  // add the PFB coupling terms to the MPB coupling rows
  pfb->set_conjugate_constraint( CouplingCns );
 }

 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::generate_dual_objective( void )
{
 // Generate the objective of each hard component
 for( int k = 0 ; k < NoHardCmps ; ++k ) {
  auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );

  pfb->generate_objective();
 }
 
 // ---- master-side FRealObjective: quadratic z term + omega * f_lev --------
 //
 // The bundle-summing terms theta^k_i b^k_i + gamma^k LB^k of every hard
 // component already live in the sub-PFB Objectives and are accumulated by
 // the SMS++ engine when the master is solved. Here we only need to add
 // the master-side stabilization terms:
 //
 //  - the quadratic stabilization is -(t/2) || z ||^2_2 under #kProximal and
 //    #kDoublyStabilized, and -1/2 || z ||^2_2 under pure #kLevel because
 //    the primal level master minimizes 1/2 || d ||^2;
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
 //                                   for proximal/doubly and for the one-shot
 //                                   level probe, -1/2 for true level, else 0);
 //                                   set_x_bar moves the linear coefficient to
 //                                   x_bar[j]
 //   triples[ NumVars ]            :  r with (linear = 0, quad = 0);
 //                                   set_global_LB moves the linear
 //                                   coefficient to LB
 //   triples[ NumVars + 1 ]        :  omega with (linear = f_lev, quad = 0)
 //                                   IF has_omega_lin; set_f_lev refreshes
 //                                   the linear coefficient
 const bool pure_level = ( StblType == kLevel );
 const bool level_probe = has_initial_level_objective();
 const bool true_level = pure_level && ( ! level_probe );
 const bool has_quad =
  pure_level || ( StblType == kProximal ) || ( StblType == kDoublyStabilized );
 const bool has_omega_lin =
  ( StblType == kLevel ) || ( StblType == kDoublyStabilized );

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
 const double quad_coeff = has_quad ? - sgn * ( true_level ? 0.5
                                                               : t_stab / 2.0 )
                                      : 0.0;
 for( int j = 0 ; j < NumVars ; ++j )
  triples.emplace_back( & Var_z[ j ] , 0.0 , quad_coeff );

 r_obj_idx = int( triples.size() );
 triples.emplace_back( & Var_r , 0.0 , 0.0 );

 omega_obj_idx = -1;
 if( has_omega_lin ) {
  omega_obj_idx = int( triples.size() );
  // The one-shot level probe is proximal and has no active level multiplier;
  // omega becomes meaningful only after remove_initial_level_objective().
  const double omega_lin =
   ( ( ! level_probe ) && std::isfinite( f_lev ) )
   ? - sgn * ( f_lev + f_C ) : 0.0;
  triples.emplace_back( & Var_omega , omega_lin , 0.0 );
  }

 // : the box slacks s^+ / s^- contribute
 //     + s^+ ( L - x_bar ) - s^- ( U - x_bar )
 // in the textbook eMax form. Coefficients are initialized from the physical
 // box state f_L / f_U if already available, and set_box() / set_x_bar()
 // refresh them later. The s^+ / s^- triples are laid out contiguously after
 // omega so that the refresh logic can address them by a fixed base offset.
 s_plus_obj_idx  = -1;
 s_minus_obj_idx = -1;
 if( NumVars > 0 ) {
  const bool iterate = ( f_v2_form != 0 );
  s_plus_obj_idx = int( triples.size() );
  for( int j = 0 ; j < NumVars ; ++j ) {
   const double xj = ( ( ! iterate ) && j < int( f_x_bar.size() ) )
                     ? f_x_bar[ j ] : 0.0;
   const bool has_L = ! f_L.empty() && std::isfinite( f_L[ j ] );
   triples.emplace_back( & Var_s_plus[ j ] ,
                         has_L ? sgn * ( f_L[ j ] - xj ) : 0.0 , 0.0 );
   }
  s_minus_obj_idx = int( triples.size() );
  for( int j = 0 ; j < NumVars ; ++j ) {
   const double xj = ( ( ! iterate ) && j < int( f_x_bar.size() ) )
                     ? f_x_bar[ j ] : 0.0;
   const bool has_U = ! f_U.empty() && std::isfinite( f_U[ j ] );
   triples.emplace_back( & Var_s_minus[ j ] ,
                         has_U ? - sgn * ( f_U[ j ] - xj ) : 0.0 , 0.0 );
   }
  }

 // In iterate form, x_bar * z already contains x_bar * g^k(u^k) for every
 // exact easy component through the coupling equations. Displacement form
 // has no x_bar * z term, so carry the equivalent contribution explicitly:
 //
 //     sum_k sum_j x_bar_j g^k_j(u^k).
 //
 // Register each involved inner variable once. set_x_bar() will update only
 // these linear coefficients when the stability centre changes.
 EasyObjVars.clear();
 EasyObjCoeffs.clear();
 easy_obj_idx = -1;
 if( ! f_v2_form ) {
  for( Index easy_id = 0 ; easy_id < EasyCmps.size() ; ++easy_id ) {
   auto * lbf = EasyCmps[ easy_id ];
   for( Index i = 0 ; i < lbf->get_num_active_var() ; ++i ) {
    const Index j = easy_local_to_global( easy_id , i );
    auto * gi = dynamic_cast< LinearFunction * >(
                                      lbf->get_Lagrangian_term( i ) );
    if( ! gi )
     throw( std::logic_error(
          "MasterProblemBlock::CreateDualMP: LagBFunction term is not "
          "a LinearFunction" ) );

    for( Function::Index h = 0 ; h < gi->get_num_active_var() ; ++h ) {
     auto * u = static_cast< ColVariable * >( gi->get_active_var( h ) );
     auto it = std::find( EasyObjVars.begin() , EasyObjVars.end() , u );
     std::size_t pos;
     if( it == EasyObjVars.end() ) {
      pos = EasyObjVars.size();
      EasyObjVars.push_back( u );
      EasyObjCoeffs.emplace_back();
      }
     else
      pos = std::size_t( std::distance( EasyObjVars.begin() , it ) );

     EasyObjCoeffs[ pos ].emplace_back( j , gi->get_coefficient( h ) );
     }
    }
   }

  if( ! EasyObjVars.empty() ) {
   easy_obj_idx = int( triples.size() );
   for( auto * u : EasyObjVars )
    triples.emplace_back( u , 0.0 , 0.0 );
   }
  }

 auto obj = new FRealObjective( this ,
                                new DQuadFunction( std::move( triples ) ) );
 obj->set_sense( IsConvex ? Objective::eMin : Objective::eMax , eNoMod );
 set_objective( obj , eNoMod );

 // Easy components are not allocated here: they are registered one by one
 // by the surrounding driver via register_easy_component(),
 // which adds the easy-cmp sub-Block and augments every CouplingCns[j]
 // with the +A^k_{i,j} u^k_i terms produced by that component.

 // Two coefficients are intentionally left unset by generate_duaò_objective 
 // and must be filled in by the surrounding driver as soon as
 // the corresponding pieces of the sum-function become known:
 //  - the omega * h coefficient on the level row, whenever X is a
 //    polyhedron with explicit right-hand side h;
 //  - the constant term b_j on every CouplingCns[j] rhs, coming from the
 //    linear part of the original sum-function.
}

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

void MasterProblemBlock::add_LBF_to_coupling_rows(
     std::vector< LinearFunction::v_coeff_pair > & vp_Cns )
{
 if( EasyCmps.empty() )
  throw( std::invalid_argument(
       "MasterProblemBlock::add_LBF_to_coupling_rows: there are no "
       "easy components" ) );

 for( Index easy_id = 0 ; easy_id < EasyCmps.size() ; ++easy_id ) {
  // Retrieve specific easy component
  auto * lbf = EasyCmps[ easy_id ];

  for( Index i = 0 ; i < lbf->get_num_active_var() ; ++i ){
   // Retrieve global index of i-th local variable
   const Index j = easy_local_to_global( easy_id , i );

   if( j >= Index( NumVars ) )
    throw( std::logic_error(
        "MasterProblemBlock::add_LBF_to_coupling_rows: "
        "global index outside master dimension" ) );

   // Collect lagrangian terms associated with the variable
   auto * gi = dynamic_cast< LinearFunction * >(
     lbf->get_Lagrangian_term( i ) );

   if( ! gi )
    throw( std::logic_error(
        "MasterProblemBlock::add_LBF_to_coupling_rows: "
        "LagBFunction term is not a LinearFunction" ) );

   // Coupling convention:
   //
   //   z_j - lambda b_j + s+_j - s-_j
   //       - sign(F_internal) g_i(u) - hard_terms = 0.
   //
   // For a convex minimisation F_internal = F, hence append -g_i(u).
   // For a concave maximisation BundleSolver minimises -F internally,
   // hence the easy-component subgradient is -g_i(u) and we append +g_i(u).
   const double easy_sign = IsConvex ? 1.0 : -1.0;
   for( Function::Index h = 0 ; h < gi->get_num_active_var() ; ++h ) {
    auto * u = static_cast< ColVariable * >( gi->get_active_var( h ) );
    const double a = gi->get_coefficient( h );

    // Append the new term to the coupling constraint terms associated
    // with the j-th global variable
    vp_Cns[ j ].emplace_back( u , easy_sign * a );
   }

   /*
   const double cst = gi->get_constant_term();
   if( cst != 0.0 ) {
    // Since the row receives -g_i(u), it also receives -cst.
    row_lf->modify_constant( row_lf->get_constant_term() - cst , eNoMod );
   }*/
  }
 }
}

/*--------------------------------------------------------------------------*/

Index MasterProblemBlock::easy_local_to_global(
    Index easy_id,
    Index local_i ) const
{
 // Dense case: no sparse map provided.
 if( EasyLocal2Global.empty() )
  return local_i;

 if( easy_id >= EasyLocal2Global.size() )
  throw( std::logic_error(
     "easy_local_to_global: invalid easy component index" ) );

 // Retrieve easy components specific map
 const auto & map = EasyLocal2Global[ easy_id ];

 if( local_i >= map.size() )
  throw( std::logic_error( "easy_local_to_global: invalid local variable "
     "index" ) );
 // return global index
 return map[ local_i ];
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
 if( k < 0 || k >= int( EasyCmps_SB.size() ) )
  return( nullptr );
 return( EasyCmps_SB[ k ] );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_dual_norm_squared( void ) const
{
 const auto z = get_z_vector();
 return( std::transform_reduce(
            z.cbegin() , z.cend() , 0.0 , std::plus<>() ,
            []( double zj ) {
             return( zj * zj );
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

  double rhs_mass = 1.0;
  if( StblType == kLevel )
   rhs_mass = has_initial_level_objective() ? 1.0 : 0.0;
  if( StblType == kLevel || StblType == kDoublyStabilized )
   rhs_mass += get_level_multiplier();

  const double gamma = rhs_mass - diagonal_mass;
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
 // Easy components are represented exactly in the dual master. Therefore,
 // for a global easy-component index k, return its true value at the
 // tentative point x_bar + d*, rather than a translated bundle-model value.
 if( ( ! IsPrimal ) && k >= 0 && k < int( IsEasyCmp.size() ) &&
     IsEasyCmp[ k ] ) {
  const int easy_k = int( std::count( IsEasyCmp.cbegin() ,
                                     IsEasyCmp.cbegin() + k , true ) );
  if( easy_k >= int( EasyCmps.size() ) ||
      easy_k >= int( EasyCmps_SB.size() ) )
   return( Inf< double >() );

  double value = 0.0;
  std::function< void( Block * ) > add_objectives =
   [ & value , & add_objectives ]( Block * block ) {
    if( ! block )
     return;

    if( auto * obj = dynamic_cast< RealObjective * >(
                                              block->get_objective() ) ) {
     obj->compute();
     value += obj->value();
     }

    for( auto * sub_block : block->get_nested_Blocks() )
     add_objectives( sub_block );
    };

  add_objectives( EasyCmps_SB[ easy_k ] );

  const auto d = get_d_vector();
  auto * lbf = EasyCmps[ easy_k ];
  for( Index i = 0 ; i < lbf->get_num_active_var() ; ++i ) {
   const Index j = easy_local_to_global( Index( easy_k ) , i );
   if( j >= f_x_bar.size() || j >= d.size() )
    return( Inf< double >() );

   auto * gi = lbf->get_Lagrangian_term( i );
   if( ! gi )
    return( Inf< double >() );

   gi->compute();
   value += ( f_x_bar[ j ] + d[ j ] ) * gi->get_value();
   }

  // BundleSolver represents a concave maximisation as the minimisation of
  // -F. Return the easy value in those same internal units.
  return( IsConvex ? - value : value );
  }

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

 // dual MP: v*[k] is the cutting-plane model's *predicted decrease*,
 // not the aggregated linearization error Sigma_k. In the proximal case
 // d* = -t z*, while in the pure-level case Var_z stores eta z* and hence
 // d* = -eta z*. The scalar below is therefore the displacement mass that
 // converts the normalized aggregate z* into the actual step.
 const double step_scale = uses_pure_level_aggregation()
                           ? get_level_multiplier() : t_stab;
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
  return( - sigma_k - step_scale * dot );
  }

 for( int kk = 0 ; kk < int( HardCmps.size() ) ; ++kk )
  if( ! has_model_row( kk ) )
   return( Inf< double >() );

 // total v* = -( Sigma + step_scale ||z*||^2 ), with step_scale equal to
 // t in proximal mode and eta in pure-level mode.
 return( - ( get_aggregated_alpha( -1 ) +
             step_scale * get_dual_norm_squared() ) );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > MasterProblemBlock::get_z_vector( void ) const
{
 std::vector< double > out;
 if( IsPrimal ) {
  // In primal form, stationarity of the stabilized master gives
  // z* + d*/t = 0. This recovers the complete essential subgradient,
  // including the linear part and active-domain/box multipliers.
  out = get_d_vector();

  if( StblType == kLevel && ! has_initial_level_objective() ) {
   // For the pure level MP, min 1/2 ||d||^2 s.t. m(d) <= L, stationarity is
   // d* + eta z* = 0, where eta is the multiplier of the level row.
   const double eta = get_level_multiplier();
   if( eta > 0.0 )
    for( auto & zj : out )
     zj = - zj / eta;
   else
    std::fill( out.begin() , out.end() , 0.0 );
   return( out );
   }

  if( t_stab > 0.0 &&
      ( StblType == kProximal || StblType == kDoublyStabilized ||
        StblType == kLevel ) )
   for( auto & zj : out )
    zj = - zj / t_stab;
  else
   std::fill( out.begin() , out.end() , 0.0 );
  return( out );
  }
 const bool normalize_level_z = uses_pure_level_aggregation();
 const double eta = normalize_level_z ? get_level_multiplier() : 1.0;
 out.reserve( Var_z.size() );

 if( normalize_level_z && eta <= 0.0 ) {
  out.assign( Var_z.size() , 0.0 );
  return( out );
  }

 for( const auto & zj : Var_z )
  // In pure level the dual stationarity vector is eta z*, so expose the
  // normalized aggregate expected by BundleSolver stopping and cut logic.
  out.push_back( normalize_level_z ? zj.get_value() / eta : zj.get_value() );
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

 const bool normalize_level_theta = uses_pure_level_aggregation();
 const double level_eta = normalize_level_theta ? get_level_multiplier() : 1.0;
 if( normalize_level_theta && level_eta <= 0.0 )
  return( out );

 for( std::size_t i = 0 ; i < A.size() ; ++i ) {
  double theta =
   pfb->get_row_multiplier( PolyhedralFunctionBlock::Index( i ) );

  if( normalize_level_theta ) {
   // Pure-level theta values are level-multiplier masses, so divide by eta
   // before forming the aggregate cut in both primal and dual master forms.
   theta /= level_eta;
   }

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

 out.reserve( Var_z.size() );

 if( uses_pure_level_aggregation() ) {
  // In pure level Var_z already stores eta z*, so the physical step is
  // d* = -eta z* = -Var_z rather than the proximal -t z* identity.
  for( const auto & zj : Var_z )
   out.push_back( - zj.get_value() );
  return( out );
  }

 // dual MP, proximal stabilization: d* = -t * z*
 for( const auto & zj : Var_z )
  out.push_back( - t_stab * zj.get_value() );
 return( out );
 }

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_level_multiplier( void ) const
{
 if( StblType != kLevel && StblType != kDoublyStabilized )
  return( 0.0 );

 if( ! IsPrimal ) {
  // In the dual master the level-row multiplier is represented explicitly by
  // omega, while in the primal master it is the dual value of LevelCns below.
  return( std::max( 0.0 , Var_omega.get_value() ) );
  }

 if( ! ( f_abs_rep & k_mpb_built_cnst ) )
  return( 0.0 );

 return( std::max( 0.0 , LevelCns.get_dual() ) );
}

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_Gid_aggregate( void ) const
{
 // Once get_z_vector() and get_d_vector() have translated the concrete master
 // representation into physical BundleSolver quantities, Gid is the same
 // scalar product in primal, dual, proximal, and level forms.
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

  if( uses_pure_level_aggregation() ) {
   // Pure-level theta/gamma masses are scaled by eta in the dual rows; divide
   // here so Sigma remains the normalized bundle-method linearization error.
   const double eta = get_level_multiplier();
   s = ( eta > 0.0 ) ? s / eta : 0.0;
   }

  return( s );
  };

 if( k >= 0 )
  return( contrib( k ) );

 const double mp_obj = get_master_objective_value();
 if( ( ! uses_pure_level_aggregation() ) && std::isfinite( mp_obj ) ) {
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
 // With exact easy components the master is not empty even before the first
 // hard cut: their inner variables, objective and coupling terms already
 // define a nontrivial optimization problem.
 if( is_bundle_empty() && ( NoEasyCmps == 0 ) ){
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

void MasterProblemBlock::refresh_box_coordinate( Index j ,
                                                 DQuadFunction * dqf )
{
 if( IsPrimal ) {
  if( int( Bounds_d.size() ) != NumVars )
   return;
  const double lower = f_L.empty() ? - Inf< double >() : f_L[ j ];
  const double upper = f_U.empty() ? Inf< double >() : f_U[ j ];
  double lhs = std::isfinite( lower ) ? lower : - Inf< double >();
  double rhs = std::isfinite( upper ) ? upper : Inf< double >();
  if( ! f_v2_form ) {
   const double xj = j < f_x_bar.size() ? f_x_bar[ j ] : 0.0;
   if( std::isfinite( lhs ) )
    lhs -= xj;
   if( std::isfinite( rhs ) )
    rhs -= xj;
   }
  Bounds_d[ j ].set_lhs( lhs );
  Bounds_d[ j ].set_rhs( rhs );
  return;
  }

 if( ! dqf || s_plus_obj_idx < 0 || s_minus_obj_idx < 0 )
  return;

 const double sgn = IsConvex ? -1.0 : 1.0;
 const bool iterate = ( f_v2_form != 0 );
 const double xj = ( ( ! iterate ) && j < f_x_bar.size() )
                   ? f_x_bar[ j ] : 0.0;

 const double lower = f_L.empty() ? - Inf< double >() : f_L[ j ];
 const bool has_L = std::isfinite( lower );
 Var_s_plus[ j ].is_fixed( ! has_L , eNoMod );
 if( ! has_L )
  Var_s_plus[ j ].set_value( 0.0 );
 dqf->modify_term( DQuadFunction::Index( s_plus_obj_idx + int( j ) ) ,
                   has_L ? sgn * ( lower - xj ) : 0.0 , 0.0 );

 const double upper = f_U.empty() ? Inf< double >() : f_U[ j ];
 const bool has_U = std::isfinite( upper );
 Var_s_minus[ j ].is_fixed( ! has_U , eNoMod );
 if( ! has_U )
  Var_s_minus[ j ].set_value( 0.0 );
 dqf->modify_term( DQuadFunction::Index( s_minus_obj_idx + int( j ) ) ,
                   has_U ? - sgn * ( upper - xj ) : 0.0 , 0.0 );
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

 const std::vector< double > lower = L.empty()
  ? std::vector< double >( NumVars , - Inf< double >() ) : L;
 const std::vector< double > upper = U.empty()
  ? std::vector< double >( NumVars , Inf< double >() ) : U;

 // BundleSolver may resubmit the complete box before every master solve.
 // Compare semantically so that an empty side and an explicit vector of
 // infinities are treated alike after a preceding partial update.
 bool lower_changed = false;
 bool upper_changed = false;
 for( Index j = 0 ; j < Index( NumVars ) ; ++j ) {
  const double old_L = f_L.empty() ? - Inf< double >() : f_L[ j ];
  const double old_U = f_U.empty() ? Inf< double >() : f_U[ j ];
  lower_changed = lower_changed || ( old_L != lower[ j ] );
  upper_changed = upper_changed || ( old_U != upper[ j ] );
  if( lower_changed && upper_changed )
   break;
  }
 if( ! lower_changed && ! upper_changed )
  return;

 f_L = L;
 f_U = U;

 auto issue_box_modifications = [ this , & lower , & upper ,
                                  lower_changed , upper_changed ]( void ) {
  if( ! anyone_there() )
   return;
  const Range range( 0 , Index( NumVars ) );
  if( lower_changed )
   add_Modification( std::make_shared< MasterProblemRngdMod >(
                     this , MasterProblemMod::BoxLowerChanged , range ,
                     MasterProblemRngdMod::Values( lower ) ) );
  if( upper_changed )
   add_Modification( std::make_shared< MasterProblemRngdMod >(
                     this , MasterProblemMod::BoxUpperChanged , range ,
                     MasterProblemRngdMod::Values( upper ) ) );
  };

 if( IsPrimal ) {
  for( Index j = 0 ; j < Index( NumVars ) ; ++j )
   refresh_box_coordinate( j , nullptr );
  issue_box_modifications();
  return;
  }

 auto obj = dynamic_cast< FRealObjective * >( get_objective() );
 auto dqf = obj ? dynamic_cast< DQuadFunction * >( obj->get_function() )
                : nullptr;
 for( Index j = 0 ; j < Index( NumVars ) ; ++j )
  refresh_box_coordinate( j , dqf );

 issue_box_modifications();
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_box( std::vector< double > L ,
                                  std::vector< double > U , Range range )
{
 if( range.second < range.first || range.second > Index( NumVars ) )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: invalid range" ) );
 const auto size = range.second - range.first;
 if( ! L.empty() && L.size() != size )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: L size does not match range" ) );
 if( ! U.empty() && U.size() != size )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: U size does not match range" ) );
 if( size == 0 )
  return;

 if( L.empty() )
  L.assign( size , - Inf< double >() );
 if( U.empty() )
  U.assign( size , Inf< double >() );

 bool lower_changed = false;
 bool upper_changed = false;
 for( Index i = 0 ; i < size ; ++i ) {
  const Index j = range.first + i;
  const double old_L = f_L.empty() ? - Inf< double >() : f_L[ j ];
  const double old_U = f_U.empty() ? Inf< double >() : f_U[ j ];
  lower_changed = lower_changed || ( old_L != L[ i ] );
  upper_changed = upper_changed || ( old_U != U[ i ] );
  if( lower_changed && upper_changed )
   break;
  }
 if( ! lower_changed && ! upper_changed )
  return;

 if( f_L.empty() )
  f_L.assign( NumVars , - Inf< double >() );
 if( f_U.empty() )
  f_U.assign( NumVars , Inf< double >() );
 std::copy( L.begin() , L.end() , f_L.begin() + range.first );
 std::copy( U.begin() , U.end() , f_U.begin() + range.first );

 DQuadFunction * dqf = nullptr;
 if( ! IsPrimal ) {
  auto obj = dynamic_cast< FRealObjective * >( get_objective() );
  dqf = obj ? dynamic_cast< DQuadFunction * >( obj->get_function() )
            : nullptr;
  }
 for( Index j = range.first ; j < range.second ; ++j )
  refresh_box_coordinate( j , dqf );

 if( anyone_there() ) {
  if( lower_changed )
   add_Modification( std::make_shared< MasterProblemRngdMod >(
                     this , MasterProblemMod::BoxLowerChanged , range ,
                     std::move( L ) ) );
  if( upper_changed )
   add_Modification( std::make_shared< MasterProblemRngdMod >(
                     this , MasterProblemMod::BoxUpperChanged , range ,
                     std::move( U ) ) );
  }
 }

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_box( std::vector< double > L ,
                                  std::vector< double > U , Subset subset ,
                                  bool ordered )
{
 const auto size = subset.size();
 if( ! L.empty() && L.size() != size )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: L size does not match subset" ) );
 if( ! U.empty() && U.size() != size )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: U size does not match subset" ) );
 if( size == 0 )
  return;

 if( L.empty() )
  L.assign( size , - Inf< double >() );
 if( U.empty() )
  U.assign( size , Inf< double >() );

 if( ( ! ordered ) && size > 1 ) {
  std::vector< std::size_t > order( size );
  std::iota( order.begin() , order.end() , std::size_t( 0 ) );
  std::sort( order.begin() , order.end() ,
             [ & subset ]( auto i , auto j )
             { return( subset[ i ] < subset[ j ] ); } );
  Subset sorted_subset( size );
  std::vector< double > sorted_L( size );
  std::vector< double > sorted_U( size );
  for( std::size_t i = 0 ; i < size ; ++i ) {
   sorted_subset[ i ] = subset[ order[ i ] ];
   sorted_L[ i ] = L[ order[ i ] ];
   sorted_U[ i ] = U[ order[ i ] ];
   }
  subset = std::move( sorted_subset );
  L = std::move( sorted_L );
  U = std::move( sorted_U );
  }

 if( subset.back() >= Index( NumVars ) )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_box: subset index out of range" ) );
 for( std::size_t i = 1 ; i < size ; ++i )
  if( subset[ i - 1 ] >= subset[ i ] )
   throw( std::invalid_argument(
        "MasterProblemBlock::set_box: unordered or repeated subset" ) );

 bool lower_changed = false;
 bool upper_changed = false;
 for( std::size_t i = 0 ; i < size ; ++i ) {
  const auto j = subset[ i ];
  const double old_L = f_L.empty() ? - Inf< double >() : f_L[ j ];
  const double old_U = f_U.empty() ? Inf< double >() : f_U[ j ];
  lower_changed = lower_changed || ( old_L != L[ i ] );
  upper_changed = upper_changed || ( old_U != U[ i ] );
  if( lower_changed && upper_changed )
   break;
  }
 if( ! lower_changed && ! upper_changed )
  return;

 if( f_L.empty() )
  f_L.assign( NumVars , - Inf< double >() );
 if( f_U.empty() )
  f_U.assign( NumVars , Inf< double >() );
 for( std::size_t i = 0 ; i < size ; ++i ) {
  f_L[ subset[ i ] ] = L[ i ];
  f_U[ subset[ i ] ] = U[ i ];
  }

 DQuadFunction * dqf = nullptr;
 if( ! IsPrimal ) {
  auto obj = dynamic_cast< FRealObjective * >( get_objective() );
  dqf = obj ? dynamic_cast< DQuadFunction * >( obj->get_function() )
            : nullptr;
  }
 for( auto j : subset )
  refresh_box_coordinate( j , dqf );

 if( anyone_there() ) {
  if( lower_changed ) {
   auto lower_subset = subset;
   add_Modification( std::make_shared< MasterProblemSbstMod >(
                     this , MasterProblemMod::BoxLowerChanged ,
                     std::move( lower_subset ) , std::move( L ) , true ) );
   }
  if( upper_changed )
   add_Modification( std::make_shared< MasterProblemSbstMod >(
                     this , MasterProblemMod::BoxUpperChanged ,
                     std::move( subset ) , std::move( U ) , true ) );
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

 if( anyone_there() )
  add_Modification( std::make_shared< MasterProblemLowerBoundMod >(
                    this , MasterProblemLowerBoundMod::GlobalLowerBound ,
                    LB ) );
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

 auto issue_reference_mod = [ this ]( void ) {
  if( anyone_there() )
   add_Modification( std::make_shared< MasterProblemMod >(
                     this , MasterProblemMod::ReferenceChanged ) );
  };

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
 if( IsPrimal ) {
  for( int k = 0 ; k < int( HardCmps.size() ) ; ++k ) {
   auto * pfb = dynamic_cast< PolyhedralFunctionBlock * >( HardCmps[ k ] );
   if( ! pfb )
    continue;
   auto & poly = pfb->get_PolyhedralFunction();

   if( ! f_v2_form && old_x_bar.size() == f_x_bar.size() ) {
    const double dF = ( k < int( old_F.size() ) )
                      ? ( F_at_x_bar[ k ] - old_F[ k ] )
                      : F_at_x_bar[ k ];
    const auto & A = poly.get_A();
    const auto b = poly.get_b();
    for( PolyhedralFunction::Index i = 0 ; i < b.size() ; ++i ) {
     if( poly.is_row_vertical( i ) )
      continue;
     double dot = 0.0;
     if( i < A.size() ) {
      const auto & Ai = A[ i ];
      const std::size_t n = std::min( Ai.size() , f_x_bar.size() );
      for( std::size_t j = 0 ; j < n ; ++j )
       dot += Ai[ j ] * ( f_x_bar[ j ] - old_x_bar[ j ] );
      }
     poly.modify_constant( i , b[ i ] + dot - dF );
     }
    }

   // Refresh finite native lower bounds from their cached physical value.
   // Absent bounds may be temporary fictitious bounds managed by BundleSolver.
   if( k < int( f_LB_raw.size() ) && std::isfinite( f_LB_raw[ k ] ) )
    poly.modify_bound( get_stored_constant( k , {} , f_LB_raw[ k ] , false ) );
   }
  issue_reference_mod();
  return;
  }

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

 issue_reference_mod();
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

 // Refresh the z quadratic consistently with generate_dual_objective().
 // The level probe is a proximal master and must keep the t-dependent
 // curvature; true pure level switches to the fixed 1/2 ||z||^2 curvature.
 // Dropping this term during set_x_bar() leaves the dual level master nearly
 // linear and can make the first non-empty QP numerically pathological.
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
 const bool level_probe = has_initial_level_objective();
 const bool true_level = ( StblType == kLevel ) && ( ! level_probe );
 const bool has_quad = ( StblType == kProximal ) ||
                       ( StblType == kDoublyStabilized ) ||
                       ( StblType == kLevel );
 const double quad_coeff = has_quad ? - sgn * ( true_level ? 0.5
                                                            : t_stab / 2.0 )
                                    : 0.0;
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

 // Displacement form has no explicit x_bar * z term. Restore its exact easy
 // component part directly as x_bar * g^k(u^k), using the cached sparse map
 // built in CreateDualMP(). Iterate form needs no such correction.
 if( ( ! iterate ) && easy_obj_idx >= 0 )
  for( std::size_t h = 0 ; h < EasyObjCoeffs.size() ; ++h ) {
   double coeff = 0.0;
   for( const auto & [ j , a ] : EasyObjCoeffs[ h ] )
    if( j < f_x_bar.size() )
     coeff += f_x_bar[ j ] * a;
   dqf->modify_term( DQuadFunction::Index( easy_obj_idx + int( h ) ) ,
                     coeff , 0.0 );
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

 auto issue_linear_part_mod = [ this ]( void ) {
  if( anyone_there() )
   add_Modification( std::make_shared< MasterProblemMod >(
                     this , MasterProblemMod::LinearPartChanged ) );
  };

 // Under the primal MP the linear part belongs directly to the Objective.
 if( IsPrimal ) {
  set_b( b );
  issue_linear_part_mod();
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

 issue_linear_part_mod();
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

bool MasterProblemBlock::has_initial_level_objective( void ) const
{
 if( StblType != kLevel )
  return( false );

 return( IsPrimal ? level_model_obj_idx >= 0 : f_dual_level_probe_active );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::remove_initial_level_objective( void )
{
 if( ! has_initial_level_objective() )
  return;

 auto * obj = dynamic_cast< FRealObjective * >( get_objective() );
 auto * dqf = obj ? dynamic_cast< DQuadFunction * >( obj->get_function() )
                  : nullptr;

 if( ! IsPrimal ) {
  // The dual probe was a proximal master: lambda + r - omega = 1, omega fixed,
  // and z quadratic -t/2 ||z||^2. Switch it to true level: lambda + r - omega
  // = 0, omega active when the level is finite, and z quadratic -1/2 ||z||^2.
  f_dual_level_probe_active = false;

  if( f_abs_rep & k_mpb_built_cnst ) {
   NormalizationCns.set_lhs( 0.0 , eNoBlck );
   NormalizationCns.set_rhs( 0.0 , eNoBlck );
   }

  if( f_abs_rep & k_mpb_built_var ) {
   const bool finite = std::isfinite( f_lev );
   if( finite ) {
    if( Var_omega.is_fixed() )
     Var_omega.is_fixed( false , eNoBlck );
    }
   else {
    Var_omega.set_value( 0 );
    if( ! Var_omega.is_fixed() )
     Var_omega.is_fixed( true , eNoBlck );
    }
   }

  if( dqf ) {
   const double sgn = IsConvex ? -1.0 : 1.0;
   for( int j = 0 ; j < NumVars ; ++j ) {
    const double lin_coeff = dqf->get_linear_coefficient(
                                      DQuadFunction::Index( z_obj_idx + j ) );
    dqf->modify_term( DQuadFunction::Index( z_obj_idx + j ) ,
                      lin_coeff , - sgn * 0.5 , eNoBlck );
    }

   if( omega_obj_idx >= 0 ) {
    const double coeff = std::isfinite( f_lev )
                         ? - sgn * ( f_lev + f_C ) : 0.0;
    dqf->modify_term( DQuadFunction::Index( omega_obj_idx ) ,
                      coeff , 0.0 , eNoBlck );
    }
   }

  return;
  }

 if( ! dqf )
  return;

 const auto n_terms = dqf->get_num_active_var();
 if( DQuadFunction::Index( level_model_obj_idx ) < n_terms ) {
  const auto n_v_terms = n_terms - DQuadFunction::Index( level_model_obj_idx );
  std::vector< double > linear( n_v_terms , 0.0 );
  std::vector< double > quadratic( n_v_terms , 0.0 );
  dqf->modify_terms( quadratic.cbegin() , linear.cbegin() ,
                     Range( DQuadFunction::Index( level_model_obj_idx ) ,
                            n_terms ) );
  }

 std::vector< double > linear( NumVars );
 std::vector< double > quadratic( NumVars , 0.5 );
 for( int j = 0 ; j < NumVars ; ++j )
  linear[ j ] = f_v2_form ? - f_x_bar[ j ] : 0.0;

 if( NumVars > 0 )
  dqf->modify_terms( quadratic.cbegin() , linear.cbegin() ,
                     Range( 0 , DQuadFunction::Index( NumVars ) ) );

 f_primal_objective_dirty = false;
 level_model_obj_idx = -1;
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

 const bool pure_level = ( StblType == kLevel );
 const bool level_probe = has_initial_level_objective();
 const bool has_quad = pure_level || ( StblType == kProximal ||
                                      StblType == kDoublyStabilized );
 const double quad = level_probe ? 1.0 / ( 2.0 * t_stab )
                                 : pure_level ? 0.5
                                : ( has_quad ? 1.0 / ( 2.0 * t_stab )
                                             : 0.0 );

 std::vector< double > linear( NumVars );
 std::vector< double > quadratic( NumVars , quad );
 bool changed = false;

 for( int j = 0 ; j < NumVars ; ++j ) {
  linear[ j ] = ( pure_level && ! level_probe ) ? 0.0 : f_linear_part[ j ];
  if( f_v2_form && has_quad ) {
   if( pure_level && ! level_probe )
    linear[ j ] -= f_x_bar[ j ];
   else
    linear[ j ] -= f_x_bar[ j ] / t_stab;
   }

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

void MasterProblemBlock::refresh_primal_level_linear_part( void )
{
 if( ! IsPrimal || ( StblType != kLevel && StblType != kDoublyStabilized ) )
  return;
 if( NumVars <= 0 )
  return;

 auto * lf = dynamic_cast< LinearFunction * >( LevelCns.get_function() );
 if( ! lf )
  return;

 LinearFunction::Vec_FunctionValue coeff( NumVars , 0.0 );
 const auto n = std::min( std::size_t( NumVars ) , f_linear_part.size() );
 for( std::size_t j = 0 ; j < n ; ++j )
  coeff[ j ] = f_linear_part[ j ];

 lf->modify_coefficients( std::move( coeff ) ,
                          Range( 0 , Index( NumVars ) ) , eModBlck );
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
 // then re-enter solve_master with a non-empty bundle. This does not apply
 // when exact easy components are present: their objective and coupling
 // terms already make the master non-empty and must be optimized.
 if( is_bundle_empty() && ( NoEasyCmps == 0 ) ) {
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

  if( std::getenv( "BS_PRINT_MULTIPLIERS" ) ) {
   std::cerr << "MPB_MULT lambda=" << get_lambda()
             << " r=" << get_r()
             << " omega=" << get_omega()
             << " eta=" << get_level_multiplier();

   for( int k = 0 ; k < int( HardCmps.size() ) ; ++k ) {
    const auto * pfb =
     dynamic_cast< const PolyhedralFunctionBlock * >( HardCmps[ k ] );
    if( ! pfb )
     continue;

    auto & poly = const_cast< PolyhedralFunctionBlock * >( pfb )
                   ->get_PolyhedralFunction();
    double theta_sum = 0.0;
    for( PolyhedralFunction::Index i = 0 ; i < poly.get_nrows() ; ++i )
     theta_sum += pfb->get_row_multiplier( i );

    const double gamma = get_gamma( k );
    std::cerr << " c" << k
              << "_theta=" << theta_sum
              << " c" << k
              << "_gamma=" << gamma
              << " c" << k
              << "_mass=" << theta_sum + gamma;
    }
   std::cerr << std::endl;
   }
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

 auto issue_t_mod = [ this ]( void ) {
  if( anyone_there() )
   add_Modification( std::make_shared< MasterProblemParamMod >(
                    this , MasterProblemMod::TChanged , t_stab ) );
  };

 if( ! ( f_abs_rep & k_mpb_built_obj ) ) {
  if( IsPrimal )
   f_primal_objective_dirty = true;
  issue_t_mod();
  return;
  }

 if( IsPrimal ) {
  f_primal_objective_dirty = true;
  issue_t_mod();
  return;
  }

 // From here on only the dual MP is handled. Pure #kNone has no
 // stabilization; true pure #kLevel uses no t, but the one-shot level probe is
 // proximal and must keep the z quadratic synchronized with t_stab.
 if( StblType == kNone ) {
  issue_t_mod();
  return;
  }
 if( StblType == kLevel && ! has_initial_level_objective() ) {
  issue_t_mod();
  return;
  }

 // The z_j quadratic terms occupy the first NumVars entries.
 if( Var_z.empty() ) {
  issue_t_mod();
  return;
  }

 auto obj = dynamic_cast< FRealObjective * >( get_objective() );
 if( ! obj ) {
  issue_t_mod();
  return;
  }
 auto dqf = dynamic_cast< DQuadFunction * >( obj->get_function() );
 if( ! dqf ) {
  issue_t_mod();
  return;
  }

 // Refresh only the quadratic coefficient; the linear coefficient carries
 // the centre-dependent z term installed by set_x_bar().
 const double sgn = IsConvex ? -1.0 : 1.0;
 const double quad_coeff = - sgn * t_stab / 2.0;
 for( int i = 0 ; i < NumVars ; ++i ) {
  const double lin_coeff = dqf->get_linear_coefficient(
                                       DQuadFunction::Index( i ) );
  dqf->modify_term( DQuadFunction::Index( i ) , lin_coeff , quad_coeff ,
                    eNoBlck );
  }

 issue_t_mod();

 }  // end( MasterProblemBlock::set_t )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_f_lev( double f )
{
 f_lev = f;

 const bool uses_level =
  ( StblType == kLevel || StblType == kDoublyStabilized );

 if( IsPrimal ) {
  // primal MP: f_lev is the RHS of the level constraint sum_k v^k <= f_lev,
  // present only under #kLevel / #kDoublyStabilized and only after the
  // abstract constraints have been generated.
  if( ( f_abs_rep & k_mpb_built_cnst ) && uses_level )
   LevelCns.set_rhs( f_lev , eNoBlck );
  }
 else if( uses_level ) {
  // dual MP: f_lev controls the omega side. The variable state belongs to
  // the abstract variables stage, while the omega coefficient belongs to the
  // abstract objective stage.
  const bool finite = std::isfinite( f_lev );

  if( f_abs_rep & k_mpb_built_var ) {
   // During the one-shot level probe the dual master is proximal, so omega is
   // not part of the active normalization even if BundleSolver installs a
   // temporary level value before the probe is removed.
   if( has_initial_level_objective() ) {
    Var_omega.set_value( 0 );
    if( ! Var_omega.is_fixed() )
     Var_omega.is_fixed( true , eNoBlck );
    }
   else if( finite ) {
    if( Var_omega.is_fixed() )
     Var_omega.is_fixed( false , eNoBlck );
    }
   else {
    Var_omega.set_value( 0 );
    if( ! Var_omega.is_fixed() )
     Var_omega.is_fixed( true , eNoBlck );
    }
   }

  if( ( f_abs_rep & k_mpb_built_obj ) && omega_obj_idx >= 0 ) {
   auto obj = dynamic_cast< FRealObjective * >( get_objective() );
   auto dqf = obj ? dynamic_cast< DQuadFunction * >( obj->get_function() )
                  : nullptr;
   if( dqf ) {
    // The omega-side contribution to the dual objective is
    //     -omega * Lvl_xbar  ,    Lvl_xbar = Lvl + f_C
    // in the textbook eMax form; the convex case flips the whole row sign.
    // A non-finite f_lev (= no level set yet) collapses the term to zero
    // rather than leaking Inf into the Objective.
    const double sgn = IsConvex ? -1.0 : 1.0;
    const double coeff = ( finite && ( ! has_initial_level_objective() ) )
                         ? - sgn * ( f_lev + f_C ) : 0.0;
    dqf->modify_term( DQuadFunction::Index( omega_obj_idx ) , coeff , 0.0 ,
                      eNoBlck );
    }
   }
  }

 if( anyone_there() )
  add_Modification( std::make_shared< MasterProblemParamMod >(
                    this , MasterProblemMod::LevelChanged , f_lev ) );

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
 refresh_primal_level_linear_part();

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

 auto issue_lower_bound_mod = [ this , k , LB ]( void ) {
  if( anyone_there() )
   add_Modification( std::make_shared< MasterProblemLowerBoundMod >(
                     this , k , LB ) );
  };

 auto & poly = pfb->get_PolyhedralFunction();
 const double no_bound = poly.is_convex()
                         ? - Inf< Function::FunctionValue >()
                         :   Inf< Function::FunctionValue >();

 // Primal MP: route the native lower bound through the same storage-frame
 // conversion used by diagonal cuts. In iterate form v^k is the physical value
 // and the bound stays LB^k; in displacement form v^k is relative to
 // F_k(x_bar), hence the bound is LB^k - F_k(x_bar).
 if( IsPrimal ) {
  const double primal_bound = std::isfinite( LB )
                              ? get_stored_constant( k , {} , LB , false )
                              : no_bound;
  poly.modify_bound( primal_bound );
  issue_lower_bound_mod();
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

 issue_lower_bound_mod();

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
