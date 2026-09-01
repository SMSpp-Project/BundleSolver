/*--------------------------------------------------------------------------*/
/*----------------------- File MILPMPSolver.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MILPMPSolver class.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/

#include "MILPMPSolver.h"

#include <AbstractBlock.h>
#include <BlockSolverConfig.h>
#include <ColVariable.h>
#include <DQuadFunction.h>
#include <FRealObjective.h>
#include <FRowConstraint.h>
#include <LinearFunction.h>
#include <MILPSolver.h>
#include <OneVarConstraint.h>

#include <stdexcept>

// note: deliberately no `using namespace SMSpp_di_unipi_it` at file scope —
// both NDO_di_unipi_it and SMSpp_di_unipi_it define a templated Inf<T>()
// helper, which would create ambiguity. SMSpp types are qualified inline.

/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- HELPER MACRO -------------------------------*/
/*--------------------------------------------------------------------------*/

#define NYI( fname )                                                          \
 throw( std::logic_error( "MILPMPSolver::" fname ": not yet implemented" ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTRUCTOR -------------------------------*/
/*--------------------------------------------------------------------------*/

MILPMPSolver::MILPMPSolver( std::istream * iStrm )
 : MPSolver()
{
 // intentionally minimal: the actual master Block and MILPSolver are
 // built lazily in SetDim() once the FiOracle dimensions are known.
 // Defaults: OptEps = FsbEps = 1e-10.
 }

/*--------------------------------------------------------------------------*/

MILPMPSolver::~MILPMPSolver()
{
 // tear down in the reverse of SetDim's order: solver_config knows how
 // to clear the Solver(s) it registered, then the master Block goes away
 if( solver_config && master_block )
  master_block->unregister_Solvers( true );
 delete solver_config;
 solver_config = nullptr;
 milp_solver = nullptr;  // owned by the Block via the registration list
 master_block.reset();
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

// builds the master AbstractBlock skeleton (y, v, empty cuts list,
// proximal DQuadFunction objective) and applies the user-supplied
// BlockSolverConfig to register the concrete MILPSolver.
void MILPMPSolver::SetDim( cIndex MxBSz , FiOracle * Oracle ,
                           const bool UsAvSt )
{
 if( ! MxBSz ) {
  // zero bundle size: deallocate everything and quietly wait
  if( solver_config && master_block )
   master_block->unregister_Solvers( true );
  milp_solver = nullptr;
  master_block.reset();
  MaxBSize = 0;
  NrFi = 0;
  MaxSGLen = 0;
  FIO = nullptr;
  comp_lb.clear();
  return;
  }

 if( ! Oracle )
  // a previous SetDim with a nullptr Oracle keeps the structure but
  // does no setup
  return;

 // record dimensions from the oracle
 MaxBSize = MxBSz;
 FIO = Oracle;
 MaxSGLen = Oracle->GetMaxNumVar();
 NrFi = Oracle->GetNrFi();
 useactiveset = UsAvSt;

 // size the lower-bounds vector (slot 0 reserved for the global LB)
 comp_lb.assign( NrFi + 1 , - Inf< HpNum >() );

 // per-slot bookkeeping for the bundle (Nm in [0, MaxBSize))
 items.assign( MaxBSize , ItemSlot{} );
 item_maxname = 0;

 // ----- build the master AbstractBlock skeleton -------------------------
 // teardown any previous master (e.g. SetDim called a second time with
 // different dimensions)
 if( solver_config && master_block )
  master_block->unregister_Solvers( true );
 milp_solver = nullptr;
 master_block = std::make_unique< SMSpp_di_unipi_it::AbstractBlock >();

 // y[NumVar]: Lambda dual multipliers, free continuous
 const Index num_y = Oracle->GetNumVar();
 if( num_y ) {
  auto y_vec = new std::vector< SMSpp_di_unipi_it::ColVariable >( num_y );
  // y are free continuous: ColVariable default-constructs to kContinuous,
  // so no set_type call is required
  master_block->add_static_variable( *y_vec , "y" );
  }

 // v[NrFi]: per-component epigraph values, continuous, lower-bounded by
 // a generous finite default so that the *initial* solve (before any
 // cut has been added) is bounded. Tighter per-component bounds arrive
 // through SetLowerBound (a follow-up).
 if( NrFi ) {
  auto v_vec = new std::vector< SMSpp_di_unipi_it::ColVariable >( NrFi );
  master_block->add_static_variable( *v_vec , "v" );

  auto v_lb = new std::vector< SMSpp_di_unipi_it::BoxConstraint >( NrFi );
  for( Index j = 0 ; j < NrFi ; ++j ) {
   ( *v_lb )[ j ].set_variable( &( ( *v_vec )[ j ] ) );
   // -1e+9 is enough to keep the master bounded below until cuts arrive,
   // without forcing CPX's QP barrier into the "large objective shift"
   // numerical regime that -1e+20 triggers (since the proximal coefficient
   // 1/(2t) can be > 1e+9 when t is small, mixing 1e+20 and 1e+9 in the
   // objective wrecks the optimality criterion). SetLowerBound() tightens
   // this per-component once the bundle provides real component LBs.
   ( *v_lb )[ j ].set_lhs( -1e+9 );
   ( *v_lb )[ j ].set_rhs(  Inf< double >() );
   }
  master_block->add_static_constraint( *v_lb , "v_lb" );
  }

 // cuts: empty dynamic list of FRowConstraints, populated by SetItem
 auto cuts = new std::list< SMSpp_di_unipi_it::FRowConstraint >();
 master_block->add_dynamic_constraint( *cuts , "cuts" );

 // objective: sum_wFi v[wFi] + (1/(2t)) sum_i (y_i - yc_i)^2
 // expanded as a DQuadFunction with triples per variable:
 //   y[i]: linear = -yc_i / t , quadratic = 1/(2t)
 //   v[j]: linear = 1         , quadratic = 0
 // plus a constant term (1/(2t)) * sum_i yc_i^2 . initially y_center is
 // zero, so the y linear and the constant are zero. Sett() updates the
 // quadratic coefficients on t changes; ChangeCurrPoint() updates the y
 // linear coefficients (and the constant) on center moves.
 SMSpp_di_unipi_it::DQuadFunction::v_coeff_triple obj_triples;
 obj_triples.reserve( num_y + NrFi );

 // ensure y_center has the right size and is initialised to zero
 y_center.assign( num_y , LMNum( 0 ) );

 const double inv_2t = ( t > 0 ) ? ( 1.0 / ( 2.0 * t ) ) : 0.0;

 if( num_y ) {
  auto y_vec = master_block->get_static_variable_v<
                                 SMSpp_di_unipi_it::ColVariable >( "y" );
  for( Index i = 0 ; i < num_y ; ++i )
   obj_triples.emplace_back( &( ( *y_vec )[ i ] ) ,
                             /* linear   */ 0.0 ,
                             /* quadratic*/ inv_2t );
  }

 if( NrFi ) {
  auto v_vec = master_block->get_static_variable_v<
                                 SMSpp_di_unipi_it::ColVariable >( "v" );
  for( Index j = 0 ; j < NrFi ; ++j )
   obj_triples.emplace_back( &( ( *v_vec )[ j ] ) ,
                             /* linear   */ 1.0 ,
                             /* quadratic*/ 0.0 );
  }

 auto obj = new SMSpp_di_unipi_it::FRealObjective();
 obj->set_function( new SMSpp_di_unipi_it::DQuadFunction(
                                              std::move( obj_triples ) ) );
 obj->set_sense( SMSpp_di_unipi_it::Objective::eMin ,
                 SMSpp_di_unipi_it::eNoMod );
 master_block->set_objective( obj );

 // ----- apply the BlockSolverConfig to register the chosen MILPSolver ---
 if( solver_config ) {
  solver_config->apply( master_block.get() );
  const auto & slvrs = master_block->get_registered_solvers();
  if( slvrs.empty() )
   throw( NDOException(
    "MILPMPSolver::SetDim: BlockSolverConfig registered no Solver" ) );
  milp_solver = dynamic_cast< SMSpp_di_unipi_it::MILPSolver * >(
                                                       slvrs.front() );
  if( ! milp_solver )
   throw( NDOException(
    "MILPMPSolver::SetDim: registered Solver is not a MILPSolver" ) );
  }
 // if no config was provided, milp_solver stays nullptr — SolveMP will
 // throw later. the LegacyBundleSolver instantiation path is expected to call
 // SetSolverConfig() before the first SolveMP() (a follow-up wiring).
 }

/*--------------------------------------------------------------------------*/

void MILPMPSolver::Sett( cHpNum tt )
{
 if( t == tt )
  return;

 t = tt;
 // refresh the proximal coefficients on the y entries of the objective
 // DQuadFunction. v entries stay (linear 1, quad 0) so they are skipped.
 if( ! master_block )
  return;

 auto * obj = dynamic_cast< SMSpp_di_unipi_it::FRealObjective * >(
                                          master_block->get_objective() );
 if( ! obj )
  return;
 auto * dq = dynamic_cast< SMSpp_di_unipi_it::DQuadFunction * >(
                                                     obj->get_function() );
 if( ! dq )
  return;

 const double inv_t = ( t > 0 ) ? ( 1.0 / t ) : 0.0;
 const double inv_2t = inv_t * 0.5;

 // y entries occupy the first MaxSGLen positions in the v_triples vector
 // (see SetDim); for each i update both the linear (-yc[i]/t) and the
 // quadratic (1/(2t)) coefficient
 const Index num_y = static_cast< Index >( y_center.size() );
 for( Index i = 0 ; i < num_y ; ++i )
  dq->modify_term( i ,
                   /* linear   */ - y_center[ i ] * inv_t ,
                   /* quadratic*/ inv_2t );
 }

/*--------------------------------------------------------------------------*/

void MILPMPSolver::SetPar( const int wp , cHpNum value )
{
 switch( wp ) {
  case kMaxTme: MaxTime = value; break;

  case kOptEps:
   // because the dual formulation is used, the role of primal and dual
   // tolerances is reversed
   OptEps = value;
   // propagate to milp_solver via ComputeConfig
   break;

  case kFsbEps:
   FsbEps = value;
   // propagate to milp_solver via ComputeConfig
   break;

  default:
   throw( NDOException( "MILPMPSolver::SetPar( HpNum ): unknown parameter" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void MILPMPSolver::SetThreads( int nthreads )
{
 f_nthreads = nthreads;
 // propagate to milp_solver via intThreads ComputeConfig
 }

/*--------------------------------------------------------------------------*/

// needs the slack/\rho/\gamma_i column layout in the master
// Block; for a follow-up just record the value
void MILPMPSolver::SetLowerBound( cHpNum LwBnd , cIndex wFi )
{
 const Index h = wFi > NrFi ? 0 : wFi;
 if( comp_lb.size() <= h )
  comp_lb.resize( h + 1 , - Inf< HpNum >() );
 comp_lb[ h ] = LwBnd;
 }

/*--------------------------------------------------------------------------*/

void MILPMPSolver::SetMPLog( std::ostream * outs , const char lvl )
{
 MPSolver::SetMPLog( outs , lvl );
 // milp_solver verbosity is pushed via SetMPLog overrides on the SMS++
 // Solver in a follow-up (intLogVerb ComputeConfig)
 }

/*--------------------------------------------------------------------------*/
/*--------------------- MILPMPSolver-SPECIFIC METHODS ----------------------*/
/*--------------------------------------------------------------------------*/

void MILPMPSolver::SetSolverConfig(
                          SMSpp_di_unipi_it::BlockSolverConfig * bsc )
{
 // if a Block already exists, unregister the currently-attached Solvers
 // (the previous config's apply effect) before switching configs
 if( master_block )
  master_block->unregister_Solvers( true );
 delete solver_config;
 solver_config = bsc;
 milp_solver = nullptr;

 // if SetDim has already built the master, apply the new config now;
 // otherwise it will be applied at the end of SetDim
 if( solver_config && master_block ) {
  solver_config->apply( master_block.get() );
  // pick the first registered Solver as our handle for convenience
  const auto & slvrs = master_block->get_registered_solvers();
  milp_solver = slvrs.empty() ? nullptr
                            : dynamic_cast< SMSpp_di_unipi_it::MILPSolver * >(
                                                            slvrs.front() );
  }
 }

/*--------------------------------------------------------------------------*/

void MILPMPSolver::SetStabType( const StabFun sf )
{
 stab_type = sf;
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

// invoke the configured MILPSolver and translate its
// SMS++-style return code (Solver::sol_type) to the MPSolver-side
// MPStatus enum.
//
// Note: at this stage the master Block only contains the y / v variables,
// an empty cut list and a placeholder linear objective; the proximal
// QuadFunction and the convexity / lower-bound rows will be wired by
// a follow-up. SolveMP will therefore not produce useful results until then —
// but the call dispatch path is correct from now on.
MPSolver::MPStatus MILPMPSolver::SolveMP( void )
{
 if( ! milp_solver )
  throw( NDOException(
   "MILPMPSolver::SolveMP: no Solver configured "
   "(SetSolverConfig must be called before the first SolveMP)" ) );

 // empty bundle: the v[wFi] epigraph variables are bounded below only by
 // a placeholder -1e+20, so the LP would be unbounded (CPX/GRB will
 // either return kUnbounded or, more commonly, flag numerical failure).
 // QPPenaltyMP / OSIMPSolver paper over this by formulating the dual
 // master, where the empty case is trivially feasible with multipliers
 // equal to zero. Replicate that here: when no item is in the bundle,
 // short-circuit to a "no-movement" solution (d* = 0, Fi = 0).
 bool any_used = false;
 for( const auto & s : items )
  if( s.used ) { any_used = true; break; }
 if( ! any_used ) {
  d_scratch.assign( MaxSGLen , LMNum( 0 ) );
  return( kOK );
  }

 using S = SMSpp_di_unipi_it::Solver;
 const int sc = milp_solver->compute();

 // pull the primal solution into the ColVariables and (when meaningful)
 // the dual values into the FRowConstraints — readers downstream
 // (Readd, ReadFiBLambda, ReadMult, ReadGid, CheckSubG's ScPri, ...)
 // read straight from the Block, so they need this to have happened
 const bool solved_ok = ( sc == S::kOK )
                     || ( sc == S::kLowPrecision )
                     || ( sc == S::kStopTime )
                     || ( sc == S::kStopIter );
 if( solved_ok ) {
  if( milp_solver->has_var_solution() )
   milp_solver->get_var_solution( nullptr );
  if( milp_solver->has_dual_solution() )
   milp_solver->get_dual_solution( nullptr );
  }

 // also refresh d_scratch so the readers (and CheckSubG's ScPri) have
 // an immediately-usable d* = y* - y_center
 if( solved_ok ) {
  const Index num_y = MaxSGLen;
  d_scratch.assign( num_y , 0 );
  auto y_vec = master_block->get_static_variable_v<
                                    SMSpp_di_unipi_it::ColVariable >( "y" );
  if( y_vec ) {
   const Index yc_size = static_cast< Index >( y_center.size() );
   for( Index i = 0 ; i < num_y && i < y_vec->size() ; ++i ) {
    const LMNum yc_i = ( i < yc_size ) ? y_center[ i ] : LMNum( 0 );
    d_scratch[ i ] = ( *y_vec )[ i ].get_value() - yc_i;
    }
   }
  }

 if( sc == S::kOK || sc == S::kLowPrecision ) return( kOK );
 if( sc == S::kInfeasible )                   return( kUnfsbl );
 if( sc == S::kUnbounded )                    return( kUnbndd );
 if( sc == S::kStopTime || sc == S::kStopIter ) return( kStppd );
 return( kError );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

// solution readers. Roadmap:
//   immediately implementable (no d* needed): ReadLowerBound, ReadLinErr,
//     MakeLambda1 (closed-form), SensitAnals (no-op)
//   needs SolveMP + master state: ReadFiBLambda (master obj), Readd (y* - yc),
//     ReadZ (d*/t), ReadMult (cut duals), ReadDt / ReadSigma / ReadDStart
//     (proximal-term values, a follow-up once the QuadFunction is wired)
//   easy-component readers: ReadDualEasy / ReadReducedCostsEasy are a
//   follow-up
//

cHpRow MILPMPSolver::refresh_lin_err_scratch()
{
 lin_err_scratch.assign( item_maxname , 0 );
 for( Index k = 0 ; k < item_maxname ; ++k )
  if( items[ k ].used )
   lin_err_scratch[ k ] = items[ k ].alpha;
 return( lin_err_scratch.data() );
 }

// Cutting-plane model value Fi[wFi]_B(Lambda). For wFi in [1, NrFi]
// returns v[wFi-1].get_value() (the epigraph variable's optimal value).
// For wFi == Inf returns the total cutting-plane value = sum over wFi.
// Caller must have triggered a SolveMP first.
HpNum MILPMPSolver::ReadFiBLambda( cIndex wFi )
{
 if( ! master_block )
  return( 0 );
 auto v_vec = master_block->get_static_variable_v<
                                   SMSpp_di_unipi_it::ColVariable >( "v" );
 if( ! v_vec )
  return( 0 );

 if( wFi == NDO_di_unipi_it::Inf< Index >() ) {
  HpNum tot = 0;
  for( const auto & v : *v_vec )
   tot += v.get_value();
  return( tot );
  }

 if( ( wFi >= 1 ) && ( wFi <= NrFi ) && ( ( wFi - 1 ) < v_vec->size() ) )
  return( ( *v_vec )[ wFi - 1 ].get_value() );

 return( 0 );
 }

// Value of the proximal stabilisation term (1/(2tt)) * ||d*||^2 at the
// current d_scratch.
HpNum MILPMPSolver::ReadDt( cHpNum tt )
{
 if( tt <= 0 )
  return( 0 );
 HpNum sq = 0;
 for( const auto v : d_scratch )
  sq += v * v;
 return( sq / ( 2.0 * tt ) );
 }

// Predicted decrease at d* for component wFi. The cutting-plane lower
// bound Fi[wFi]_B(Lambda1) at the next iterate is captured by the
// epigraph value v[wFi], so "predicted decrease" is the gap between
// the model value at the center (which equals -alpha_min over the
// active cuts) and at Lambda1; reporting v[wFi].get_value() is the
// usual convention since the bundle compares against the current Fi.
// For wFi == Inf returns the master's objective at d* (sum v +
// proximal), which is the canonical Sigma in the bundle's eU test.
HpNum MILPMPSolver::ReadSigma( cIndex wFi )
{
 if( wFi != NDO_di_unipi_it::Inf< Index >() )
  return( ReadFiBLambda( wFi ) );

 // global Sigma = -(master objective at d*). The master minimises
 //   sum_wFi v[wFi]  +  (1/(2t)) || y* - y_center ||^2
 // We compute it directly from v* and d_scratch instead of relying on
 // FRealObjective::get_value(), which caches the function value and may
 // not have been recomputed at y* / v* by the MILPSolver pipeline (the
 // QPSolver computes its own objective internally and the Function's
 // cache stays at its last explicit compute()).
 if( ! master_block )
  return( 0 );
 HpNum total = 0;
 auto v_vec = master_block->get_static_variable_v<
                                  SMSpp_di_unipi_it::ColVariable >( "v" );
 if( v_vec )
  for( const auto & vi : *v_vec )
   total += vi.get_value();
 total += ReadDt( t );  // proximal (1/(2t)) ||d*||^2
 return( - total );
 }

// ReadDStart( t ) is documented by LegacyBundleSolver as
//   ReadDStart( t ) == t * || z* ||_2^2 / 2
// which, given d = -t z*, is exactly the proximal stabilisation term
// value at the master optimum (i.e. ||d*||^2 / (2t)). LegacyBundleSolver uses
// it to back-derive ||d*||_2 (and from there ||z*||_2) in FormD, so it
// must NOT default to zero — that's what was making ||z*|| stick at 0
// and the outer bundle spin forever in NR steps.
HpNum MILPMPSolver::ReadDStart( cHpNum tt )       { return( ReadDt( tt ) ); }

// d* = y* - y_center. Read it on demand from the master Block. If the
// master has not been solved yet (no y values populated) the buffer is
// returned as all-zeros. Fulld is honoured only structurally for now.
cLMRow MILPMPSolver::Readd( bool Fulld )
{
 const Index num_y = MaxSGLen;  // dense representation length
 d_scratch.assign( num_y , 0 );
 if( ! master_block )
  return( d_scratch.data() );

 auto y_vec = master_block->get_static_variable_v<
                                   SMSpp_di_unipi_it::ColVariable >( "y" );
 if( ! y_vec )
  return( d_scratch.data() );

 // align: yc is empty (a follow-up fills it). treat missing entries as 0.
 const Index yc_size = static_cast< Index >( y_center.size() );
 for( Index i = 0 ; i < num_y && i < y_vec->size() ; ++i ) {
  const LMNum yc_i = ( i < yc_size ) ? y_center[ i ] : LMNum( 0 );
  d_scratch[ i ] = ( *y_vec )[ i ].get_value() - yc_i;
  }
 return( d_scratch.data() );
 }

// ReadZ: z* = d* / t in the primal-quadratic master (sign matches the
// bundle's convention of pointing towards the minimiser).
void MILPMPSolver::ReadZ( LMRow tz , cIndex_Set & I , Index & D , cIndex )
{
 D = MaxSGLen;
 I = nullptr;
 const Index ds_size = static_cast< Index >( d_scratch.size() );
 const HpNum inv_t = ( t > 0 ) ? ( 1.0 / t ) : 0.0;
 for( Index i = 0 ; i < D ; ++i )
  tz[ i ] = ( i < ds_size ) ? ( d_scratch[ i ] * inv_t ) : LMNum( 0 );
 }

// Cut multipliers: read the dual values of the FRowConstraints stored
// in items[] (or only those of the requested component). The returned
// buffer is allocated in mult_scratch; mult_idx_scratch carries the
// corresponding slot indices when sparse.
//
// IncldCnst toggles whether constraint-type items are included alongside
// subgradient-type ones.
cHpRow MILPMPSolver::ReadMult( cIndex_Set & I , Index & D ,
                               cIndex wFi , const bool IncldCnst )
{
 mult_scratch.clear();
 mult_idx_scratch.clear();

 const bool per_comp = ( wFi != NDO_di_unipi_it::Inf< Index >() );
 for( Index k = 0 ; k < items.size() ; ++k ) {
  const auto & slot = items[ k ];
  if( ! slot.used ) continue;
  if( per_comp && slot.wFi != wFi ) continue;
  if( ! slot.is_subg && ! IncldCnst ) continue;

  mult_scratch.push_back( slot.row_it->get_dual() );
  mult_idx_scratch.push_back( k );
  }

 D = static_cast< Index >( mult_scratch.size() );

 // LegacyBundleSolver's sparse-MBse consumers (e.g. UpdtCntrs) iterate the
 // index vector with `for( ; *(MBse++) < InINF ; )` — i.e. they expect
 // it to be Inf< Index >()-terminated, ignoring D. Append the sentinel
 // (plus a dummy multiplier entry to keep the two arrays in lockstep
 // for any consumer that does pair-walk them).
 if( ! mult_idx_scratch.empty() ) {
  mult_idx_scratch.push_back( NDO_di_unipi_it::Inf< Index >() );
  mult_scratch.push_back( HpNum( 0 ) );
  }

 I = mult_idx_scratch.empty() ? nullptr : mult_idx_scratch.data();
 return( mult_scratch.empty() ? nullptr : mult_scratch.data() );
 }

// Component lower-bound multiplier. a follow-up wires the LB row.
HpNum MILPMPSolver::ReadLBMult( cIndex wFi )      { return( 0 ); }

// Easy-component duals — a follow-up.
cHpRow MILPMPSolver::ReadDualEasy( cIndex wFi )   { return( nullptr ); }
cHpRow MILPMPSolver::ReadReducedCostsEasy( cIndex wFi )
                                                  { return( nullptr ); }

// Scalar product < g_Nm , d* > for item Nm. Uses the stored g_terms of
// the slot (Index-coefficient pairs) and the current d_scratch (last
// d* populated by SolveMP).
HpNum MILPMPSolver::ReadGid( cIndex Nm )
{
 if( Nm == NDO_di_unipi_it::Inf< Index >() ) {
  // aggregated request: sum over all used slots
  HpNum tot = 0;
  for( const auto & slot : items )
   if( slot.used )
    for( const auto & p : slot.g_terms )
     if( p.first < d_scratch.size() )
      tot += p.second * d_scratch[ p.first ];
  return( tot );
  }
 if( Nm >= items.size() || ! items[ Nm ].used )
  return( 0 );
 HpNum gd = 0;
 for( const auto & p : items[ Nm ].g_terms )
  if( p.first < d_scratch.size() )
   gd += p.second * d_scratch[ p.first ];
 return( gd );
 }

// Lambda1 = Lambda + Tau * d* — closed-form, no master state needed.
// Uses the d* currently in scratch (Readd should have populated it for
// this iteration; if not, Lambda1 == Lambda).
void MILPMPSolver::MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau )
{
 const Index num_y = MaxSGLen;
 const Index ds_size = static_cast< Index >( d_scratch.size() );
 for( Index i = 0 ; i < num_y ; ++i ) {
  const LMNum di = ( i < ds_size ) ? d_scratch[ i ] : LMNum( 0 );
  Lmbd1[ i ] = Lmbd[ i ] + Tau * di;
  }
 }

// Sensitivity analysis on t. Without a quadratic stab there is nothing
// to be sensitive to; return both endpoints as the current t (i.e. no
// allowed perturbation range). a follow-up supplies real values.
void MILPMPSolver::SensitAnals( HpNum & lp , HpNum & cp )
{
 lp = t;
 cp = t;
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

// read-only queries over the bundle slot table
//
// BSize(wFi)    — number of items in component wFi (sub + constraints).
//                 wFi in [1,NrFi] is per-component; wFi == Inf is the total.
// BCSize(wFi)   — number of constraint items in wFi (or total).
// MaxName(wFi)  — highest slot index +1 occupied by items of wFi.
//                 wFi == Inf returns item_maxname directly.
// WComponent(i) — items[i].wFi, or Inf if slot is free.
// IsSubG(i)     — true if items[i] is a subgradient (not a constraint).
//
// Variable-bound queries return zero / false: the master here has y free
// continuous and v free continuous, so no non-negative or bounded vars.

Index MILPMPSolver::BSize( cIndex wFi )
{
 Index n = 0;
 const bool per_comp = ( wFi > 0 ) && ( wFi <= NrFi );
 for( const auto & slot : items ) {
  if( ! slot.used ) continue;
  if( per_comp && slot.wFi != wFi ) continue;
  ++n;
  }
 return( n );
 }

Index MILPMPSolver::BCSize( cIndex wFi )
{
 Index n = 0;
 const bool per_comp = ( wFi > 0 ) && ( wFi <= NrFi );
 for( const auto & slot : items ) {
  if( ! slot.used || slot.is_subg ) continue;
  if( per_comp && slot.wFi != wFi ) continue;
  ++n;
  }
 return( n );
 }

Index MILPMPSolver::MaxName( cIndex wFi )
{
 if( wFi == NDO_di_unipi_it::Inf< Index >() )
  return( item_maxname );
 // scan backwards from item_maxname for the highest slot belonging to wFi
 for( Index k = item_maxname ; k > 0 ; --k )
  if( items[ k - 1 ].used && items[ k - 1 ].wFi == wFi )
   return( k );
 return( 0 );
 }

Index MILPMPSolver::WComponent( cIndex i )
{
 if( i >= items.size() || ! items[ i ].used )
  return( NDO_di_unipi_it::Inf< Index >() );
 return( items[ i ].wFi );
 }

bool MILPMPSolver::IsSubG( cIndex i )
{
 if( i >= items.size() || ! items[ i ].used )
  return( false );
 return( items[ i ].is_subg );
 }

// y is free continuous in our master, so no non-negative or bounded vars
Index MILPMPSolver::NumNNVars( void )     { return( 0 ); }
Index MILPMPSolver::NumBxdVars( void )    { return( 0 ); }
bool  MILPMPSolver::IsNN( cIndex )        { return( false ); }

// duplicate-detection flag consumed by CheckSubG / CheckCnst
void MILPMPSolver::CheckIdentical( const bool Chk ) { check_id = Chk; }

cHpRow MILPMPSolver::ReadLinErr( void )
{
 // build a contiguous vector of HpNum from the per-item alpha (lazy)
 return( refresh_lin_err_scratch() );
 }

HpNum MILPMPSolver::ReadLowerBound( cIndex wFi )
{
 // mirror SetLowerBound's index convention: wFi > NrFi (i.e. Inf) means
 // the global LB stored at slot 0
 const Index h = wFi > NrFi ? 0 : wFi;
 if( h >= comp_lb.size() )
  return( - Inf< HpNum >() );
 return( comp_lb[ h ] );
 }
HpNum  MILPMPSolver::EpsilonD( void )                  { return FsbEps; }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

// GetItem returns a writable buffer of length MaxSGLen the
// caller (LegacyBundleSolver::GetGi) fills with the next subgradient. The
// buffer is reused across calls; SetItem consumes it to produce a new
// FRowConstraint on the master Block. wFi (1-based component index in
// the MPSolver world, or NrFi+1 for the linear part) is recorded so
// SetItem knows which v[wFi] the cut applies to.
SgRow MILPMPSolver::GetItem( cIndex wFi )
{
 if( g1k_buffer.size() < MaxSGLen )
  g1k_buffer.resize( MaxSGLen );
 next_item_wFi = wFi;
 // default to subgradient mode; CheckCnst will flip the flag if the
 // caller is actually adding a constraint
 next_item_is_subg = true;
 return( g1k_buffer.data() );
 }

// SetItemBse stores the sparsity base for the next SetItem.
// Passing nullptr means dense (the whole NumVar range).
void MILPMPSolver::SetItemBse( cIndex_Set SGBse , cIndex SGBDm )
{
 g1k_base.clear();
 if( SGBse && SGBDm )
  g1k_base.assign( SGBse , SGBse + SGBDm );
 // SetItem reads g1k_base to build the LinearFunction terms (dense if empty)
 }

// CheckSubG translates the linearization error Ai of the
// subgradient stored in g1k_buffer from "at Lambda1" to "at the current
// point Lambda", and computes the scalar product ScPri = <g, d*> where
// d* is the current primal optimal of the master.
//
// At the current scaffold level d* is not yet available (a follow-up
// implements Readd), so ScPri is set to 0 — this is mathematically wrong
// but keeps the dispatch path runnable until a follow-up fills in d*. The
// duplicate-check return is also stubbed to "no duplicate" (Inf).
//
// Both values are cached for SetItem to consume.
// CheckSubG: ScPri = <g, d*> where g is the new subgradient currently
// in g1k_buffer (sparse over g1k_base if non-empty) and d* is the most
// recent direction stored in d_scratch by SolveMP / Readd. Ai is then
// translated by the standard formula Ai_out = Ai_in - DFi + (Tau/t)*ScPri.
Index MILPMPSolver::CheckSubG( cHpNum DFi , cHpNum Tau ,
                               HpNum & Ai , HpNum & ScPri )
{
 // compute ScPri = <g, d*>
 ScPri = 0;
 const Index ds_size = static_cast< Index >( d_scratch.size() );
 if( g1k_base.empty() ) {
  for( Index i = 0 ; i < MaxSGLen && i < ds_size ; ++i )
   ScPri += g1k_buffer[ i ] * d_scratch[ i ];
  }
 else {
  for( Index k = 0 ; k < g1k_base.size() ; ++k ) {
   const Index i = g1k_base[ k ];
   if( i < ds_size )
    ScPri += g1k_buffer[ k ] * d_scratch[ i ];
   }
  }

 // Ai_out = Ai_in - DFi + (Tau/t)*ScPri
 if( ( Tau > 0 ) && ( t > 0 ) )
  Ai = Ai - DFi + ( Tau / t ) * ScPri;

 next_item_Ai      = Ai;
 next_item_ScPri   = ScPri;
 next_item_is_subg = true;

 return( Inf< Index >() );  // duplicate detection deferred
 }

// CheckCnst handles a *constraint* item (vertical/horizontal
// linearisation) instead of a subgradient. The semantics differ from
// CheckSubG: Ai is interpreted as the rhs of the constraint, not as a
// linearization error to translate. CrrPnt is the current point (Lambda
// or Lambda1) and is unused at the scaffold level.
Index MILPMPSolver::CheckCnst( HpNum & Ai , HpNum & ScPri , cHpRow )
{
 // same ScPri computation as CheckSubG — the bundle's d* and the
 // constraint's coefficient row produce the same scalar-product term
 ScPri = 0;
 const Index ds_size = static_cast< Index >( d_scratch.size() );
 if( g1k_base.empty() ) {
  for( Index i = 0 ; i < MaxSGLen && i < ds_size ; ++i )
   ScPri += g1k_buffer[ i ] * d_scratch[ i ];
  }
 else {
  for( Index k = 0 ; k < g1k_base.size() ; ++k ) {
   const Index i = g1k_base[ k ];
   if( i < ds_size )
    ScPri += g1k_buffer[ k ] * d_scratch[ i ];
   }
  }

 next_item_Ai      = Ai;
 next_item_ScPri   = ScPri;
 next_item_is_subg = false;

 return( Inf< Index >() );  // duplicate detection deferred
 }
// SetItem above will trigger Modifications via add_dynamic_
// constraints, so the master MILPSolver re-solves from a perturbed state.
// Sett / ChangeCurrPoint similarly modify the DQuadFunction. Hence the
// MP solution always changes between solves except for trivial no-ops.
bool MILPMPSolver::ChangesMPSol( void )             { return( true ); }

// commit the subgradient/constraint currently held in
// g1k_buffer (filled by GetItem + the caller, with sparsity in g1k_base
// from SetItemBse) as a new FRowConstraint on the master Block, stored
// under bundle name Nm.
//
// Cut form for a subgradient on component wFi:
//   v[wFi] >= alpha + g^T y       i.e.   alpha <= v[wFi] - g^T y <= +INF
// where alpha is the *post-CheckSubG* linearization error at Lambda
// (which CheckSubG cached in next_item_Ai). LinearFunction terms are
//   (+1, v[wFi-1])
//   (-g[i], y[i])  for each non-zero coefficient
//
// Constraint variant (next_item_is_subg == false) has the same shape
// minus the v term:
//   alpha <= - g^T y <= +INF   (a follow-up will refine the sign/rhs handling)
void MILPMPSolver::SetItem( cIndex Nm )
{
 if( ! master_block )
  throw( NDOException(
   "MILPMPSolver::SetItem: master Block not yet built (call SetDim)" ) );

 if( Nm < NDO_di_unipi_it::Inf< Index >() && Nm >= items.size() )
  throw( NDOException( "MILPMPSolver::SetItem: Nm out of range" ) );

 // ----- assemble the LinearFunction coefficient pairs -------------------
 const Index wFi_local = next_item_wFi;  // 1-based in MPSolver world
 SMSpp_di_unipi_it::LinearFunction::v_coeff_pair coeffs;

 // (+1, v[wFi_local - 1]) only for subgradient items; constraints have
 // no v-term in this convention (they bound y directly)
 if( next_item_is_subg && wFi_local && ( wFi_local <= NrFi ) ) {
  auto v_vec = master_block->get_static_variable_v<
                                   SMSpp_di_unipi_it::ColVariable >( "v" );
  if( v_vec )
   coeffs.emplace_back( &( (*v_vec)[ wFi_local - 1 ] ) , 1.0 );
  }

 // (-g[i], y[i]) for each non-zero coefficient (sparse if g1k_base
 // populated by SetItemBse, otherwise dense over the whole y range).
 // we also collect the (y_idx, g[i]) pairs into a temporary g_pairs so
 // SetItem can later store them in the slot for ChangeCurrPoint /
 // ChgSubG / ReadGid use.
 auto y_vec = master_block->get_static_variable_v<
                                  SMSpp_di_unipi_it::ColVariable >( "y" );
 if( ! y_vec )
  throw( NDOException(
   "MILPMPSolver::SetItem: master Block missing 'y' static variables" ) );

 std::vector< std::pair< Index , double > > g_pairs;

 if( g1k_base.empty() ) {
  // dense: g1k_buffer[i] is the coefficient on y[i] for i in [0, MaxSGLen)
  for( Index i = 0 ; i < MaxSGLen && i < y_vec->size() ; ++i ) {
   const double gi = g1k_buffer[ i ];
   if( gi != 0.0 ) {
    coeffs.emplace_back( &( (*y_vec)[ i ] ) , - gi );
    g_pairs.emplace_back( i , gi );
    }
   }
  }
 else {
  // sparse: g1k_buffer[k] is the coefficient on y[ g1k_base[k] ]
  for( Index k = 0 ; k < g1k_base.size() ; ++k ) {
   const Index i = g1k_base[ k ];
   const double gi = g1k_buffer[ k ];
   if( gi != 0.0 && i < y_vec->size() ) {
    coeffs.emplace_back( &( (*y_vec)[ i ] ) , - gi );
    g_pairs.emplace_back( i , gi );
    }
   }
  }

 // ----- build the FRowConstraint ----------------------------------------
 // it'll be placed in a one-element temporary list and add_dynamic_constraints
 // splices it into the master's "cuts" list, so the iterator stays valid
 auto tmp = new std::list< SMSpp_di_unipi_it::FRowConstraint >( 1 );
 auto & cut = tmp->front();
 cut.set_function( new SMSpp_di_unipi_it::LinearFunction(
                                                   std::move( coeffs ) ) );
 cut.set_lhs( next_item_Ai );
 cut.set_rhs(   Inf< double >() );

 auto cuts_list = master_block->get_dynamic_constraint<
                              SMSpp_di_unipi_it::FRowConstraint >( "cuts" );
 if( ! cuts_list )
  throw( NDOException(
   "MILPMPSolver::SetItem: master Block missing 'cuts' dynamic list" ) );

 auto it_before_splice = tmp->begin();  // iterator stays valid after splice
 master_block->add_dynamic_constraints( *cuts_list , *tmp );
 delete tmp;  // tmp is now empty, just free the container

 // ----- bookkeeping -----------------------------------------------------
 if( Nm < NDO_di_unipi_it::Inf< Index >() ) {
  auto & slot = items[ Nm ];
  slot.used    = true;
  slot.is_subg = next_item_is_subg;
  slot.wFi     = wFi_local;
  slot.alpha   = next_item_Ai;
  slot.g_terms = std::move( g_pairs );
  slot.row_it  = it_before_splice;
  if( Nm + 1 > item_maxname )
   item_maxname = Nm + 1;
  }
 // Nm == InINF: anonymous insertion (not addressed via items[]); used by
 // aggregation, which is a follow-up
 }
// SubstItem replaces the cut at slot Nm with the one currently
// pending in g1k_buffer / g1k_base / next_item_*. The simplest reliable
// implementation is "remove then add" via the dynamic-constraint plumbing
// in SetItem; the master MILPSolver sees the equivalent BlockModRmv +
// BlockModAdd Modifications. A direct in-place ChgSubG-style refresh
// would be cheaper but needs the per-row LinearFunction term ordering
// known — deferred with ChgSubG above.
void MILPMPSolver::SubstItem( cIndex Nm )
{
 if( Nm >= items.size() || ! items[ Nm ].used )
  return;
 RmvItem( Nm );
 SetItem( Nm );
 }

// remove the cut at slot Nm from the master's dynamic cuts list
// via Block::remove_dynamic_constraint (issuing a BlockModRmv so the
// attached MILPSolver re-syncs).
void MILPMPSolver::RmvItem( cIndex Nm )
{
 if( Nm >= items.size() || ! items[ Nm ].used )
  return;
 auto cuts_list = master_block->get_dynamic_constraint<
                              SMSpp_di_unipi_it::FRowConstraint >( "cuts" );
 if( ! cuts_list )
  return;
 // remove the single FRowConstraint pointed to by row_it. SMS++ Block
 // exposes remove_dynamic_constraint( list , iterator , issueMod ).
 master_block->remove_dynamic_constraint( *cuts_list , items[ Nm ].row_it );
 items[ Nm ] = ItemSlot{};
 }

// bulk teardown of every cut currently in the dynamic list.
void MILPMPSolver::RmvItems( void )
{
 if( ! master_block )
  return;
 auto cuts_list = master_block->get_dynamic_constraint<
                              SMSpp_di_unipi_it::FRowConstraint >( "cuts" );
 if( cuts_list && ! cuts_list->empty() ) {
  // empty-subset overload removes the entire list and issues a single
  // BlockModRmvSbst with subset().empty() (see Block.h:3994-4009).
  // SMSpp_di_unipi_it::Block::Subset is a std::vector< Index > alias.
  std::vector< Index > all;
  master_block->remove_dynamic_constraints<
                              SMSpp_di_unipi_it::FRowConstraint >(
                          *cuts_list , std::move( all ) );
  }
 for( auto & slot : items )
  slot = ItemSlot{};
 item_maxname = 0;
 }

// dimension-changing operations
//
// AddVars / RmvVars manipulate the Lambda dimensionality (the number of y
// variables in the master). In the current MILPMPSolver layout y is
// allocated as a *static* std::vector<ColVariable> by SetDim, so adding
// or removing entries at runtime would invalidate all the LinearFunction
// pointers held by existing cuts. Supporting this cleanly requires
// promoting y to a dynamic std::list<ColVariable> (matching the pattern
// in MILPSolver/test_dynamic.cpp). Deferred to a follow-up: the bundle
// usage in TSSB EC / LDS+LDS+MILP, which is the driver for this work,
// keeps the Lambda dimension constant for the whole solve.
void MILPMPSolver::AddVars( cIndex NNwVrs )
{
 if( NNwVrs == 0 ) return;
 // promote y to a dynamic variable group and
 // append NNwVrs new free ColVariables here, then extend y_center with
 // zeros (this becomes the new center for the freshly added dimensions).
 // For now silently no-op — call sites in LegacyBundleSolver only fire AddVars
 // when the FiOracle dimension grows mid-solve, which the bundle paths
 // exercised by TSSB EC do not do.
 }

void MILPMPSolver::RmvVars( cIndex_Set whch , Index hwmny )
{
 if( ! hwmny ) return;
 // same reason as AddVars (dynamic-list y).
 // Removing y[i] would also need to walk every cut and drop the term
 // associated with that variable, plus shrink y_center.
 }

// Sparse-Lambda (active-set) management. In sparse-Lambda mode only a
// subset of the y variables is "active" at any one time; cuts reference
// only the active ones. The classic implementation is via a bitmask
// (Aset[]) and reconfigures the OSI columns on every change.
//
// For our SMS++-backed master the natural translation is to:
//   - hold a sparsity vector (already in g1k_base for the current item)
//   - keep y_center / d_scratch full-length, but treat inactive entries
//     as effectively zero so the proximal term doesn't penalise them
//   - mask out cut contributions on inactive indices (handled at SetItem
//     time, which already respects g1k_base)
//
// The structural change to y (static-vec) is the same blocker as AddVars,
// so for now SetActvSt / AddActvSt / RmvActvSt just record the active
// indices without touching the master — sufficient for the LegacyBundleSolver
// dense-Lambda path which is the only one tested.
void MILPMPSolver::SetActvSt( cIndex_Set AVrs , cIndex AVDm )
{
 // when sparse Lambda is finally supported,
 //   - track active indices here (probably a std::vector<Index>)
 //   - update the master's y bounds: inactive y_i fixed to 0, active y_i
 //     free
 //   - re-run ChangeCurrPoint to refresh proximal centres
 useactiveset = ( AVrs != nullptr );
 }

void MILPMPSolver::AddActvSt( cIndex_Set Addd , cIndex AdDm , cIndex_Set )
{
 // same blocker as SetActvSt
 }

void MILPMPSolver::RmvActvSt( cIndex_Set Rmvd , cIndex RmDm , cIndex_Set )
{
 // same blocker as SetActvSt
 }

// ChgAlfa overloads
//   (DeltaAlfa)              — shift every cut's alpha by DeltaAlfa[wFi+1]
//                              (or DeltaAlfa[0] for the linear-part slot)
//   (NewAlfa, wFi)           — assign NewAlfa[k] to every cut k of wFi
//   (i, Ai)                  — set alpha of cut at slot i to Ai
// All variants must push the new value to the FRowConstraint's lhs so the
// master MILPSolver re-syncs.
void MILPMPSolver::ChgAlfa( cHpRow DeltaAlfa )
{
 if( ! DeltaAlfa )
  return;
 for( Index k = 0 ; k < items.size() ; ++k ) {
  auto & slot = items[ k ];
  if( ! slot.used )
   continue;
  slot.alpha += DeltaAlfa[ slot.wFi ];
  slot.row_it->set_lhs( slot.alpha );
  }
 }

void MILPMPSolver::ChgAlfa( cHpRow NewAlfa , cIndex wFi )
{
 if( ! NewAlfa )
  return;
 Index k_in_comp = 0;
 for( Index k = 0 ; k < items.size() ; ++k ) {
  auto & slot = items[ k ];
  if( ! slot.used || slot.wFi != wFi )
   continue;
  slot.alpha = NewAlfa[ k_in_comp++ ];
  slot.row_it->set_lhs( slot.alpha );
  }
 }

void MILPMPSolver::ChgAlfa( cIndex i , cHpNum Ai )
{
 if( i >= items.size() || ! items[ i ].used )
  return;
 items[ i ].alpha = Ai;
 items[ i ].row_it->set_lhs( Ai );
 }

// proximal-center move helpers
//
// The master objective is
//     sum_wFi v[wFi] + (1/(2t)) sum_i (y_i - yc_i)^2
// stored as a DQuadFunction in expanded form. When the center moves
//     yc_new = yc_old + DLambda
// the y linear coefficients become -yc_new[i]/t (the quadratic stays at
// 1/(2t) unless Sett() runs). Cut linearization errors also shift:
//     alpha_new = alpha_old + DFi[wFi+1] - <g_k, DLambda>
// since the cut v[wFi] - g^T y >= alpha has a center-dependent alpha.
//
// DFi format (mirror of LegacyBundleSolver.cpp:GotoLambda1 / GotoLambda):
//   DFi[0]   = total Fi-value delta across all components
//   DFi[i+1] = per-component delta for component i (i in [0,NrFi))

namespace {

// pull the y-linear from a DQuadFunction obj, update y_center[i] += dy,
// push the new linear back. constant term is also kept in sync (a follow-up
// extension once we wire constant updates via DQuadFunction).
static void refresh_y_linear( SMSpp_di_unipi_it::DQuadFunction * dq ,
                              const std::vector< LMNum > & y_center ,
                              double inv_t )
{
 const Index num_y = static_cast< Index >( y_center.size() );
 for( Index i = 0 ; i < num_y ; ++i )
  dq->modify_linear_coefficient( i , - y_center[ i ] * inv_t );
 }

}  // anonymous namespace

void MILPMPSolver::ChangeCurrPoint( cLMRow DLambda , cHpRow DFi )
{
 // update y_center: yc_new = yc_old + DLambda
 const Index num_y = static_cast< Index >( y_center.size() );
 for( Index i = 0 ; i < num_y ; ++i )
  y_center[ i ] += DLambda[ i ];

 // propagate to the master proximal y-linears
 if( master_block && t > 0 ) {
  auto * obj = dynamic_cast< SMSpp_di_unipi_it::FRealObjective * >(
                                          master_block->get_objective() );
  if( obj ) {
   auto * dq = dynamic_cast< SMSpp_di_unipi_it::DQuadFunction * >(
                                                      obj->get_function() );
   if( dq )
    refresh_y_linear( dq , y_center , 1.0 / t );
   }
  }

 // shift each stored cut's alpha by the standard bundle update rule
 //   alpha_new = alpha_old + DFi[wFi] - <g_k, DLambda>
 // where g_k is the cut's stored subgradient (from g_terms) and DLambda
 // is the center move. for constraint items the DFi correction does
 // not apply (their alpha represents a rhs, not a linearisation error).
 for( Index k = 0 ; k < items.size() ; ++k ) {
  auto & slot = items[ k ];
  if( ! slot.used )
   continue;

  HpNum delta = 0;
  if( DFi && slot.is_subg )
   delta += DFi[ slot.wFi ];

  // <g_k, DLambda> = sum over stored (idx, g[i]) pairs
  double gd = 0;
  for( const auto & p : slot.g_terms )
   if( p.first < num_y )
    gd += p.second * DLambda[ p.first ];
  delta -= gd;

  if( delta != 0 ) {
   slot.alpha += delta;
   slot.row_it->set_lhs( slot.alpha );
   }
  }
 }

void MILPMPSolver::ChangeCurrPoint( cHpNum Tau , cHpRow DFi )
{
 // DLambda = Tau * d* where d* was populated by the latest Readd() call.
 // d_scratch is the natural source for d* (see a follow-up).
 const Index num_y = static_cast< Index >( y_center.size() );
 const Index ds_size = static_cast< Index >( d_scratch.size() );

 // build DLambda inline so we can reuse the DLambda-form path
 std::vector< LMNum > dlambda( num_y , LMNum( 0 ) );
 for( Index i = 0 ; i < num_y && i < ds_size ; ++i )
  dlambda[ i ] = Tau * d_scratch[ i ];

 ChangeCurrPoint( dlambda.data() , DFi );
 }

// change a single coefficient of a stored subgradient. The
// FRowConstraint's LinearFunction holds (-g[i]) at the y[i] term; the
// Modification is pushed by LinearFunction::modify_coefficient.
//
// Currently a no-op: walking back from
// (strt, stp, wFi) to the right v_coeff_pair indices requires either
// caching the LinearFunction term order per slot or scanning by variable
// pointer. Deferred to the follow-up that also stores g per slot.
void MILPMPSolver::ChgSubG( cIndex strt , Index stp , cIndex wFi )
{
  // no-op for now: changes already flow to the MP via add_Modification
 // outside this method
 }

/*--------------------------------------------------------------------------*/

#undef NYI

}  // end( namespace NDO_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*---------------------- End File MILPMPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
