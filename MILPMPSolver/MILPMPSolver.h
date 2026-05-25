/*--------------------------------------------------------------------------*/
/*------------------------ File MILPMPSolver.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * MILPMPSolver is an MPSolver-derived class that solves the (stabilized)
 * Master Problem of the Generalized Bundle algorithm using an SMS++
 * MILPSolver registered to an AbstractBlock holding the master.
 *
 * The master problem is represented as an AbstractBlock with:
 *   - a vector of free ColVariables y of size NumVar (the dual
 *     multipliers driven by the bundle);
 *   - a vector of free ColVariables v of size NrFi (per-component
 *     epigraph values);
 *   - one FRowConstraint per bundle item, of the form
 *         alpha <= v[wFi] - g^T y <= +INF   (subgradient cuts)
 *   - an FRealObjective minimising  sum_i v[i]  +  proximal stabilisation
 *     term (a 2-norm quadratic penalty on y centred at the current
 *     stability centre).
 *
 * Master modifications (add cut, remove cut, change LB, change current
 * point, etc.) translate to SMS++ Modifications on the AbstractBlock that
 * the underlying MILPSolver consumes via its standard Modification
 * pipeline. SolveMP() calls solver->compute(); primal/dual readers go
 * through the Block's Variable/Constraint accessors. The concrete
 * MILPSolver backend is selected at runtime through a BlockSolverConfig
 * passed via SetSolverConfig.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MILPMPSolver
 #define __MILPMPSolver

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MPSolver.h"

#include <FRowConstraint.h>

#include <list>
#include <memory>
#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------- FORWARD DECLARATIONS ---------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{
 class AbstractBlock;
 class BlockSolverConfig;
 class MILPSolver;
 class ColVariable;
 class FRealObjective;
 class LinearFunction;
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS MILPMPSolver -----------------------------*/
/*--------------------------------------------------------------------------*/

class MILPMPSolver : public MPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// stabilisation function applied to the proximal term
 enum StabFun { unset = 0 , none , boxstep , quadratic };

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTRUCTOR -------------------------------*/
/*--------------------------------------------------------------------------*/

 MILPMPSolver( std::istream * iStrm = nullptr );

 ~MILPMPSolver() override;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 void SetDim( cIndex MxBSz = 0 , FiOracle * Oracle = nullptr ,
              const bool UsAvSt = false ) override;

 void Sett( cHpNum tt = 1 ) override;

 void SetPar( const int wp , cHpNum value ) override;

 void SetThreads( int nthreads ) override;

 void SetLowerBound( cHpNum LwBnd = - Inf< HpNum >() ,
                     cIndex wFi = Inf< Index >() ) override;

 void SetMPLog( std::ostream * outs , const char lvl ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- MILPMPSolver-SPECIFIC METHODS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// configure the solver attached to the master Block
 /** Hands the master a fully-formed `BlockSolverConfig` whose `apply()`
  * will instantiate and register the concrete solver (CPX / GRB / SCIP /
  * HiGHS / …). `MILPMPSolver` takes ownership of the passed config —
  * it stores it and invokes `apply()` from `SetDim()` once the master
  * `AbstractBlock` is built. Passing `nullptr` clears the previously
  * configured solver.
  *
  * This is the single point of solver-choice configuration; the class
  * never instantiates a concrete `MILPSolver` subclass directly. */

 void SetSolverConfig( SMSpp_di_unipi_it::BlockSolverConfig * bsc );

 /// pick the stabilisation type (none / boxstep / quadratic)
 void SetStabType( const StabFun sf = none );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

 MPStatus SolveMP( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 HpNum ReadFiBLambda( cIndex wFi = Inf< Index >() ) override;

 HpNum ReadDt( cHpNum tt = 1 ) override;

 HpNum ReadSigma( cIndex wFi = Inf< Index >() ) override;

 HpNum ReadDStart( cHpNum tt = 1 ) override;

 cLMRow Readd( bool Fulld = false ) override;

 void ReadZ( LMRow tz , cIndex_Set & I , Index & D ,
             cIndex wFi = Inf< Index >() ) override;

 cHpRow ReadMult( cIndex_Set & I , Index & D ,
                  cIndex wFi = Inf< Index >() ,
                  const bool IncldCnst = true ) override;

 HpNum ReadLBMult( cIndex wFi = Inf< Index >() ) override;

 cHpRow ReadDualEasy( cIndex wFi ) override;

 cHpRow ReadReducedCostsEasy( cIndex wFi ) override;

 HpNum ReadGid( cIndex Nm = Inf< Index >() ) override;

 void MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau ) override;

 void SensitAnals( HpNum & lp , HpNum & cp ) override;

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

 Index BSize( cIndex wFi = Inf< Index >() ) override;

 Index BCSize( cIndex wFi = Inf< Index >() ) override;

 Index MaxName( cIndex wFi = Inf< Index >() ) override;

 Index WComponent( cIndex i ) override;

 bool IsSubG( cIndex i ) override;

 Index NumNNVars( void ) override;

 Index NumBxdVars( void ) override;

 bool IsNN( cIndex i ) override;

 void CheckIdentical( const bool Chk ) override;

 cHpRow ReadLinErr( void ) override;

 HpNum ReadLowerBound( cIndex wFi = Inf< Index >() ) override;

 HpNum EpsilonD( void ) override;

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 SgRow GetItem( cIndex wFi = Inf< Index >() ) override;

 void SetItemBse( cIndex_Set SGBse , cIndex SGBDm ) override;

 Index CheckSubG( cHpNum DFi , cHpNum Tau , HpNum & Ai ,
                  HpNum & ScPri ) override;

 Index CheckCnst( HpNum & Ai , HpNum & ScPri , cHpRow CrrPnt ) override;

 bool ChangesMPSol( void ) override;

 void SetItem( cIndex Nm = Inf< Index >() ) override;

 void SubstItem( cIndex Nm ) override;

 void RmvItem( cIndex i ) override;

 void RmvItems( void ) override;

 void SetActvSt( cIndex_Set AVrs , cIndex AVDm ) override;

 void AddActvSt( cIndex_Set Addd , cIndex AdDm ,
                 cIndex_Set AVrs ) override;

 void RmvActvSt( cIndex_Set Rmvd , cIndex RmDm ,
                 cIndex_Set AVrs ) override;

 void AddVars( cIndex NNwVrs ) override;

 void RmvVars( cIndex_Set whch , Index hwmny ) override;

 void ChgAlfa( cHpRow DeltaAlfa ) override;

 void ChgAlfa( cHpRow NewAlfa , cIndex wFi ) override;

 void ChgAlfa( cIndex i , cHpNum Ai ) override;

 void ChangeCurrPoint( cLMRow DLambda , cHpRow DFi ) override;

 void ChangeCurrPoint( cHpNum Tau , cHpRow DFi ) override;

 void ChgSubG( cIndex strt , Index stp , cIndex wFi ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 // master Block + attached MILPSolver
 std::unique_ptr< SMSpp_di_unipi_it::AbstractBlock > master_block;
 SMSpp_di_unipi_it::MILPSolver * milp_solver = nullptr;

 // solver configuration (owned; applied to master_block in SetDim) and
 // stabilisation
 SMSpp_di_unipi_it::BlockSolverConfig * solver_config = nullptr;
 StabFun stab_type = unset;

 // bundle parameters
 HpNum    t = 1;            ///< proximity parameter
 HpNum    MaxTime = 0;      ///< per-call time limit (0 = none)
 HpNum    OptEps = 1e-10;   ///< optimality tolerance
 HpNum    FsbEps = 1e-10;   ///< feasibility tolerance
 int      f_nthreads = 0;   ///< thread count for the MILPSolver backend
 Index    MaxBSize = 0;     ///< maximum bundle dimension
 Index    NrFi = 0;         ///< number of components
 Index    MaxSGLen = 0;     ///< maximum subgradient length
 FiOracle * FIO = nullptr;  ///< oracle handle
 bool     useactiveset = false;
 bool     check_id = false; ///< CheckIdentical flag

 // per-component lower bounds (stored verbatim; tightening of the
 // corresponding row constraints happens in SetLowerBound)
 std::vector< HpNum > comp_lb;

 // working data populated by the GetItem / SetItemBse / CheckSubG /
 // CheckCnst → SetItem pipeline
 std::vector< double > g1k_buffer;   // GetItem scratch
 std::vector< Index >  g1k_base;     // SetItemBse copy (empty == dense)
 Index   next_item_wFi   = 0;        // recorded by GetItem
 HpNum   next_item_Ai    = 0;        // computed by CheckSubG / CheckCnst
 HpNum   next_item_ScPri = 0;        // ditto
 bool    next_item_is_subg = true;   // subgradient vs constraint
 Index   item_maxname     = 0;       // highest Nm ever passed to SetItem

 // per-slot bookkeeping: for each item name Nm in [0, MaxBSize), where
 // it lives in the master's dynamic cuts list (used == false means the
 // slot is free). g_terms stores the (y_index, subgradient_coefficient)
 // pairs of the cut — required by ChangeCurrPoint(DLambda) to compute
 // the - <g_k, DLambda> shift on alpha, by ReadGid to compute
 // <g_Nm, d*>, and by ChgSubG to update LinearFunction coefficients in
 // place (instead of remove+add).
 struct ItemSlot {
  bool  used = false;
  bool  is_subg = true;             // cut from subgradient vs constraint
  Index wFi = 0;                    // component index (1..NrFi)
  HpNum alpha = 0;                  // current linearization error
  std::vector< std::pair< Index , double > > g_terms; // (y_idx, g[i])
  std::list< SMSpp_di_unipi_it::FRowConstraint >::iterator row_it{};
  };
 std::vector< ItemSlot > items;

 // current proximal center y_c, populated by ChangeCurrPoint. The
 // default-zero value is the correct starting centre for the very first
 // SolveMP call. Scratch buffers below are sized lazily by their
 // respective Read* methods.
 std::vector< LMNum > y_center;
 std::vector< LMNum > d_scratch;     // Readd buffer (y* - y_center)
 std::vector< LMNum > z_scratch;     // ReadZ buffer (= d_scratch / t)
 std::vector< HpNum > lin_err_scratch;  // ReadLinErr buffer
 std::vector< HpNum > mult_scratch;     // ReadMult buffer (cut duals)
 std::vector< Index > mult_idx_scratch; // ReadMult sparse index map

 // linearization-error vector accessor: lazily rebuilds the contiguous
 // HpNum array from items[].alpha for ReadLinErr()
 cHpRow refresh_lin_err_scratch();

 };  // end( class MILPMPSolver )

/*--------------------------------------------------------------------------*/

}  // end( namespace NDO_di_unipi_it )

#endif  /* __MILPMPSolver */

/*--------------------------------------------------------------------------*/
/*------------------------ End File MILPMPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
