/*--------------------------------------------------------------------------*/
/*------------------- File MasterProblemBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MasterProblemBlock class.
 *
 * The current state is a *work-in-progress* skeleton: the structural part
 * (Block plumbing, dimension setup, Solver registration, abstract /
 * factory hooks) is in place, while the actual generation of the static
 * Variable, Constraint and Objective of the primal and dual MP is still
 * TODO, see CreatePrimalMP() and CreateDualMP().
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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MasterProblemBlock.h"

#include "BlockSolverConfig.h"
#include "Configuration.h"
#include "DQuadFunction.h"
#include "FRealObjective.h"
#include "LinearFunction.h"
#include "PolyhedralFunctionBlock.h"

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
 IsPrimal = false;
 StblType = kDoublyStabilized;
 t_stab   = 1.0;

 // drop the dynamically-sized variable groups (primal d and v^k, dual z
 // and CouplingCns); the scalar members (Var_lambda, Var_r, Var_omega,
 // NormalizationCns) keep their default state and are released by their
 // own destructors when *this is destroyed
 Var_d.clear();
 Var_v_hard.clear();
 Var_z.clear();
 CouplingCns.clear();

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
       "MasterProblemBlock::register_Solver: empty configuration filename; "
       "a BlockSolverConfig is required to attach a Solver to the Master "
       "Problem Block" ) );

 auto cfg = Configuration::deserialize( solv_cfg_filename );
 auto MPBSC = dynamic_cast< BlockSolverConfig * >( cfg );
 if( ! MPBSC ) {
  delete cfg;
  throw( std::invalid_argument(
       "MasterProblemBlock::register_Solver: the provided Configuration "
       "file is not a BlockSolverConfig" ) );
  }

 MPBSC->apply( this );
 MPBSC->clear();
 delete MPBSC;

 }  // end( MasterProblemBlock::register_Solver )

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
 // we only support proximal stabilization in the primal form at the moment;
 // level / doubly-stabilized variants need either a separate level row or a
 // dedicated sub-block and are not implemented yet
 if( Stbl != kProximal )
  throw( std::logic_error(
       "MasterProblemBlock::CreatePrimalMP: only kProximal is supported "
       "in the primal form for now" ) );

 StblType = Stbl;
 IsPrimal = true;

 // ---- static Variable: the step d and the per-hard-component epigraph v^k --

 Var_d.clear();
 Var_d.resize( NumVars );          // d is free by default
 if( NumVars > 0 )
  add_static_variable( Var_d , "MPB_d" );

 Var_v_hard.clear();
 Var_v_hard.resize( NoHardCmps );  // v^k is free; the v >= g*d + alpha cuts
                                   // are added later as dynamic constraints
 if( NoHardCmps > 0 )
  add_static_variable( Var_v_hard , "MPB_v" );

 // ---- Objective: min  b*d + sum_k v^k  +  (1/(2t)) || d ||^2_2 ------------
 //
 // the linear coefficient on d (the "b" of the paper, i.e. the constant
 // gradient of the linear part of the original sum-function) is left at 0
 // here; the BundleSolver is expected to install it through the standard
 // ModBlck interface (the LinearFunction member of the DQuadFunction tracks
 // its variables and reacts to coefficient changes via Modification).

 const double quad_coeff = 1.0 / ( 2.0 * t_stab );

 DQuadFunction::v_coeff_triple triples;
 triples.reserve( NumVars + NoHardCmps );

 for( int i = 0 ; i < NumVars ; ++i )
  triples.emplace_back( & Var_d[ i ] , 0.0 , quad_coeff );

 for( int k = 0 ; k < NoHardCmps ; ++k )
  triples.emplace_back( & Var_v_hard[ k ] , 1.0 , 0.0 );

 auto obj = new FRealObjective( this , new DQuadFunction( std::move( triples ) ) );
 obj->set_sense( Objective::eMin , eNoMod );
 set_objective( obj , eNoMod );

 }  // end( MasterProblemBlock::CreatePrimalMP )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::CreateDualMP( stabilization_type Stbl )
{
 StblType = Stbl;
 IsPrimal = false;

 // ---- static dual multipliers: lambda >= 0, r >= 0, omega ----------------
 //
 // omega is the multiplier of the level constraint v <= f_lev and is
 // non-negative under #kLevel / #kDoublyStabilized; with pure proximal
 // stabilization the level row does not exist, so omega is fixed to 0.

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
 add_static_variable( Var_r , "MPB_r" );
 add_static_variable( Var_omega , "MPB_omega" );

 // ---- z auxiliary variables (free, one per coordinate) -------------------

 Var_z.clear();
 Var_z.resize( NumVars );
 if( NumVars > 0 )
  add_static_variable( Var_z , "MPB_z" );

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

 // ---- coupling rows z_j = b_j (b = 0 until BundleSolver sets the linear --
 // ---- part of the original sum-function)                              ----
 //
 // CouplingCns is exposed as a *dynamic* group because that is the only
 // overload of Block::add_*_constraint() that accepts std::list, and the
 // PolyhedralFunctionBlock::set_conjugate_constraint API takes a
 // std::list< FRowConstraint > & by reference. The list size itself is in
 // fact constant (NumVars) and never grows during the algorithm.

 CouplingCns.clear();
 CouplingCns.resize( NumVars );
 {
  auto it = CouplingCns.begin();
  for( int j = 0 ; j < NumVars ; ++j , ++it ) {
   LinearFunction::v_coeff_pair vp{ { & Var_z[ j ] , 1.0 } };
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
 // (Generalized)BundleSolver feeds rows (g, alpha) into each f_polyf via
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

  pfb->generate_abstract_variables(
                              const_cast< SimpleConfiguration< int > * >( & rep_dual ) );
  pfb->generate_abstract_constraints();
  pfb->generate_objective();

  pfb->set_lambda( & Var_lambda );
  pfb->set_conjugate_constraint( CouplingCns );

  HardCmps.push_back( pfb );
  add_nested_Block( pfb );
  }

 // ---- master-side FRealObjective: -(t/2) || z ||^2_2 ---------------------
 //
 // The bundle-summing terms theta^k_i b^k_i + gamma^k LB^k of every hard
 // component already live in the sub-PFB Objectives and are accumulated by
 // the SMS++ engine when the master is solved. Here we only need to add
 // the master-side quadratic stabilization on z and the (currently zero)
 // linear part x_bar * z, which the BundleSolver will later modify via
 // LinearFunction::modify_coefficient as the stability centre changes.

 const double quad_coeff = - t_stab / 2.0;

 DQuadFunction::v_coeff_triple triples;
 triples.reserve( NumVars );
 for( int j = 0 ; j < NumVars ; ++j )
  triples.emplace_back( & Var_z[ j ] , 0.0 , quad_coeff );

 auto obj = new FRealObjective( this , new DQuadFunction( std::move( triples ) ) );
 obj->set_sense( Objective::eMax , eNoMod );
 set_objective( obj , eNoMod );

 // Easy components are not allocated here: they are registered one by one
 // by the surrounding (Generalized)BundleSolver via register_easy_component(),
 // which adds the easy-cmp sub-Block and augments every CouplingCns[j]
 // with the +A^k_{i,j} u^k_i terms produced by that component.

 // TODO (post-MVP):
 //  - the omega * h coefficient on the level row when X is a polyhedron;
 //  - inject the constant term b_j on every CouplingCns[j] rhs as soon as
 //    the linear part of the original sum-function is known.

 }  // end( MasterProblemBlock::CreateDualMP )

/*--------------------------------------------------------------------------*/
/*--------------------------- EASY COMPONENTS ------------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::register_easy_component(
                  Block * easy_blk ,
                  std::vector< LinearFunction::v_coeff_pair > && A_rows )
{
 if( IsPrimal )
  throw( std::logic_error(
       "MasterProblemBlock::register_easy_component: easy components are "
       "only supported in the dual MP" ) );

 if( ! easy_blk )
  throw( std::invalid_argument(
       "MasterProblemBlock::register_easy_component: null sub-Block" ) );

 if( int( EasyCmps.size() ) >= NoEasyCmps )
  throw( std::logic_error(
       "MasterProblemBlock::register_easy_component: all NoEasyCmps slots "
       "have been registered already" ) );

 if( int( A_rows.size() ) != NumVars )
  throw( std::invalid_argument(
       "MasterProblemBlock::register_easy_component: A_rows must have "
       "NumVars entries" ) );

 // 1. attach the sub-Block (and transfer ownership)
 add_nested_Block( easy_blk );
 EasyCmps.push_back( easy_blk );

 // 2. augment every CouplingCns[j] with the contributed +A^k_{i,j} u^k_i
 //    terms, leaving the rhs untouched
 auto it = CouplingCns.begin();
 for( int j = 0 ; j < NumVars ; ++j , ++it ) {
  if( A_rows[ j ].empty() )
   continue;
  auto lf = dynamic_cast< LinearFunction * >( it->get_function() );
  if( ! lf )
   throw( std::logic_error(
        "MasterProblemBlock::register_easy_component: CouplingCns row "
        "does not carry a LinearFunction" ) );
  lf->add_variables( std::move( A_rows[ j ] ) , eNoMod );
  }

 }  // end( MasterProblemBlock::register_easy_component )

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
/*------------------------- STABILIZATION PARAMETER ------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_t( double t )
{
 if( t <= 0.0 )
  throw( std::invalid_argument(
       "MasterProblemBlock::set_t: t must be strictly positive" ) );

 t_stab = t;

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

 const double quad_coeff = IsPrimal ?   1.0 / ( 2.0 * t_stab )
                                    : - t_stab / 2.0;
 for( int i = 0 ; i < NumVars ; ++i )
  dqf->modify_term( DQuadFunction::Index( i ) , 0.0 , quad_coeff );

 }  // end( MasterProblemBlock::set_t )

/*--------------------------------------------------------------------------*/
/*------------------------- LOAD INTO THE SOLVER ---------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::load_problem( void )
{
 const auto & solvers_list = this->get_registered_solvers();
 if( solvers_list.empty() )
  throw( std::logic_error(
       "MasterProblemBlock::load_problem: no Solver has been registered "
       "to the Master Problem Block" ) );

 // TODO: once the abstract representation is fully generated by
 //       CreatePrimalMP() / CreateDualMP(), trigger here any
 //       implementation-specific bridge between the abstract
 //       representation and the registered Solver (e.g. the per-easy
 //       component "compact" embedding into the master LP/QP).
 //
 // For the moment, just verify that at least one Solver is attached:
 // ordinary [MILP]Solver-s do not expose a public load_problem(), they
 // ingest the abstract representation lazily on the first compute().
 (void) solvers_list;

 }  // end( MasterProblemBlock::load_problem )

/*--------------------------------------------------------------------------*/
/*----------------------- End File MasterProblemBlock.cpp ------------------*/
/*--------------------------------------------------------------------------*/
