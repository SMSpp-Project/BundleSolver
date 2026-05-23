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

 // drop any group of dynamic / static gamma multipliers that might have
 // been previously allocated; the rest of the abstract representation
 // (Var_lambda, Var_r, Var_omega, NormalizationCns) is owned by the
 // class as plain members and is automatically released by the
 // corresponding destructors when *this is destroyed
 Var_gamma.clear();

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

 // reserve the slots for the per-component sub-Blocks
 EasyCmps.reserve( NoEasyCmps );
 HardCmps.reserve( NoHardCmps );

 // pre-size the per-hard-component LB multipliers
 Var_gamma.resize( NoHardCmps );

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
 StblType = Stbl;

 // TODO: build the primal MP described in (P), i.e.
 //
 //  - one ColVariable per component of d (NumVars of them);
 //  - one ColVariable v^k per "hard" component, plus the global v;
 //  - one PolyhedralFunctionBlock sub-Block per "hard" component
 //    (in its *primal* linearised representation), feeding the
 //    sub-gradient cuts v^k >= g^k_i * d + alpha^k_i;
 //  - the global epigraph constraint v >= b * d + sum_k v^k;
 //  - the FRealObjective wrapping a DQuadFunction with the proximal
 //    quadratic term (1/2t) ||d||^2_2 plus the linear part (v).

 throw( std::logic_error(
      "MasterProblemBlock::CreatePrimalMP: not implemented yet" ) );

 }  // end( MasterProblemBlock::CreatePrimalMP )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::CreateDualMP( stabilization_type Stbl )
{
 StblType = Stbl;

 // TODO: build the dual MP described in (D), i.e.
 //
 //  - the static non-negative ColVariable lambda / r / omega (the latter
 //    fixed to 0 unless Stbl == kLevel or Stbl == kDoublyStabilized);
 //  - one non-negative ColVariable gamma^k per "hard" component (already
 //    pre-sized in Var_gamma by SetDim());
 //  - one PolyhedralFunctionBlock sub-Block per "hard" component (in its
 //    *dual* linearised representation), feeding the theta^k_i multipliers
 //    and the normalization sum_i theta^k_i + gamma^k = lambda;
 //  - one sub-Block per "easy" component (e.g. the inner Block of the
 //    corresponding LagBFunction) for the compact u^k variables and the
 //    A^k / G^k constraints;
 //  - the global normalization row lambda + r - omega = 1;
 //  - the FRealObjective wrapping a DQuadFunction with the dual proximal
 //    term -(t/2) ||z||^2_2 plus the linear part.

 throw( std::logic_error(
      "MasterProblemBlock::CreateDualMP: not implemented yet" ) );

 }  // end( MasterProblemBlock::CreateDualMP )

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
