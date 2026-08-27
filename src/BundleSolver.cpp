/*--------------------------------------------------------------------------*/
/*---------------------- File BundleSolver.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BundleSolver class.
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BundleSolver.h"

#include "LagBFunction.h"

#include "OneVarConstraint.h"

#include <iomanip>
#include <unordered_set>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define BLOG( l , x ) if( f_log && ( LogVerb > l ) ) *f_log << x

#define BLOG2( l , c , x ) if( f_log && ( LogVerb > l ) && c ) *f_log << x

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Perform long and costly checks on the data structures, coded bit-wise:
  *
  * - CHECK_DS & 1 == checks the data structures representing the bundle and
  *                   the global pools against the MPSolver and the
  *                   C05Function(s)
  *
  * - CHECK_DS & 2 == checks that the aggregated linearization produced by
  *                   the C05Function agrees with that produced by the
  *                   MPSolver
  *
  * - CHECK_DS & 4 == checks that the aggregated linearization errors
  *                   directly computed with fresh data (linearization +
  *                   constant) out of the C05Function agree with these
  *                   stored in the MPSolver
  *
  * - CHECK_DS & 8 == checks that the lower bounds out of the C05Function
  *                   agree with these stored in the MPSolver */

#else
 #define CHECK_DS 0
 // never change this
#endif

#define CHECK_BAD_F 0
/* Bundle methods are supposed to work on convex functions. Technically,
 * this boils down to the fact that each (eps-)subgradient produced by
 * each orcale must be a linear lower approximation of the corresponding
 * function on all the space. This is immediately tested right away for the
 * current stability centre Lambda by computing the linearization error of
 * the subgradient (for the corresponding component) w.r.t. that point. If
 *
 * - the function is convex to start with
 *
 * - the function values and (eps-)subgradients (and the value of eps) are
 *   correctly computed
 *
 * - the oracle is "faithful", i.e., it correctly reports the subgradients
 *   as eps-ones (with the correct value of eps) rather than pretending that
 *   they are exact, i.e., it correctly computes and returns the upper and
 *   lower bounds on the function value
 *
 * then no negative linearization error should ever appear. Sometimes this is
 * not the case. In Lagrangian optimization, for instance, some oracles may
 * not solve the Lagrangian subproblem exactly and they may not be capable
 * (or willing) to compute correct upper/lower bounds on the objective value
 * so as to correctly declare the subgradient as an eps-one and provide a
 * correct estimate of the eps; rather, these "cheating" oracles may just
 * pretend the subgradient to be an exact one and leave the poor Bundle
 * method to fend off with the consequences. These may be particularly
 * nefarious in that the linearization errors are then used to compute the
 * aggregate linearization error (Sigma) which enters in the crucial
 * stopping criterion of the algorithm. Negative linearization errors may
 * lead to a negative Sigma, which breaks the stopping criterion.
 *
 * BundleSolver is engineered to be "resistant" to Negative linearization
 * errors by detecting negative Sigma and performing "noise reduction steps"
 * to try to make them go away. However, in general one may expect that, for
 * some applications, this should never happen as the functions are convex
 * linearization and especially negative Sigma, would be a sign that the
 * oracles are not behaving as expected. This macro, coded bitwose, causes
 * checks on negative linearization errors and/or negative Sigma to be
 * performed and warnings to be printed on std::cerr if "negative enough"
 * values are found. The exact coding is:
 *
 * - CHECK_BAD_F & 1 == checks the sign of the aggregate linearization error
 *                      (Sigma) as soon as produced by the Master Problem
 *
 * - CHECK_BAD_F & 2 == checks the sign of any linearization error of any
 *                      new subgradient w.r.t. the current stability centre
 *                      Lambda as soon as the subgradient is extracted from
 *                      the corresponding oracle; the check is disable at
 *                      the first iteration and in general whenever the
 *                      reference value of the corresponding component is
 *                      undefined, as in this case linearization errors are
 *                      computed w.r.t. an arbitrary reference value (say, 0)
 *                      and them being negative does not mean anything
 *                      untowards having happened.
 *
 * Note that the warnings are printed whatever the value of the log verbosity
 * and on std::cerr rather than on the BundleSolver log stream (if any), hence
 * they are "rather invasive". */

#if CHECK_BAD_F
 #define CHECK_BAD_F_EPS 1e-6
 // the relative threshold for checking negative linearization errors
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Index = BundleSolver::Index;
using SIndex = int;
using VarValue = BundleSolver::VarValue;
using Subset = BundleSolver::Subset;
using c_Subset = BundleSolver::c_Subset;
using Range = BundleSolver::Range;
using LinearCombination = BundleSolver::LinearCombination;
using Vec_VarValue = BundleSolver::Vec_VarValue;

std::pair< double , double >
BundleSolver::effective_bounds( const ColVariable * var )
{
 double lb = var->is_positive()
             ? 0.0 : - std::numeric_limits< double >::infinity();
 double ub = var->get_ub();

 for( Index i = 0 ; i < var->get_num_active() ; ++i ) {
  const auto c = var->get_active( i );
  const auto b = dynamic_cast< const OneVarConstraint * >( c );
  if( ! b )
   continue;

  // Only one-variable rows can tighten a column bound. Read their generic
  // interval once instead of trying each concrete bound subclass in turn.
  lb = std::max( lb , double( b->get_lhs() ) );
  ub = std::min( ub , double( b->get_rhs() ) );
  }

 return( std::pair< double , double >( lb , ub ) );
}

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr double Nearly  = 1.01;
static constexpr double Nearly2 = 1.02;

static constexpr char LogBnd = 16;        // log Bundle changes
static constexpr char LogVar = 32;        // log variables changes

static constexpr Index tSP1Msk = 12;  // mask for tSPar1: the long-term t-s
static constexpr Index kSLTTS =  4;   // "soft" long-term t-strategy
static constexpr Index kHLTTS =  8;   // "hard" long-term t-strategy
static constexpr Index kBLTTS = 12;   // "balancing" long-term t-strategy
static constexpr Index kEGTTS = 16;   // "endgame" long-term t-strategy
static constexpr Index tSPHMsk1 = 192;  // mask for heuristics: bits 6 and 7
static constexpr Index tSPHMsk2 = 768;  // mask for heuristics: bits 7 and 8

static constexpr unsigned char RstAlg = 1;  // don't reset algorithmic params
static constexpr unsigned char RstCrr = 2;  // don't reset current point to
                                            // all-0, use Variable value()

static constexpr auto InINF = SMSpp_di_unipi_it::Inf< Index >();

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/
// set precision for long floats (10 digits) in scientific notation

static inline std::ostream & def( std::ostream & os ) {
 os.setf( std::ios::scientific , std::ios::floatfield );
 os << std::setprecision( 10 );
 return( os );
 }

/*--------------------------------------------------------------------------*/
// set precision for short floats (2 digits) in scientific notation

static inline std::ostream & shrt( std::ostream & os ) {
 os.setf( std::ios::scientific , std::ios::floatfield );
 os << std::setprecision( 2 );
 return( os );
 }

/*--------------------------------------------------------------------------*/
// set precision for short floats (4 digits) in fixed notation

static inline std::ostream & fixd( std::ostream & os ) {
 os.setf( std::ios::fixed , std::ios::floatfield );
 os << std::setprecision( 4 );
 return( os );
 }

/*--------------------------------------------------------------------------*/
// set precision for short floats (4 digits) in scientific notation

static inline std::ostream & shrt4( std::ostream & os ) {
 os.setf( std::ios::scientific , std::ios::floatfield );
 os << std::setprecision( 4 );
 return( os );
 }

/*--------------------------------------------------------------------------*/
// cleanly print +/-INF

static inline void pval( std::ostream & os , double val ) {
 if( val == BundleSolver::INFshift )
  os << "INF";
 else
  if( val == -BundleSolver::INFshift )
   os << "-INF";
  else
   os << val;
 }

/*--------------------------------------------------------------------------*/

static void Compact( BundleSolver::Vec_VarValue & g ,
                     BundleSolver::c_Subset & B )
{
 // takes a "dense" n-vector g and "compacts" it deleting the elements whose
 // indices are in B; all elements of B must be in the range 0 .. n, B must
 // be ordered in increasing sense
 // the remaining entries in g are shifted left of the minimum possible
 // amount in order to fill the holes left by the deleted ones
 // g is *not* resized in here

 auto Bit = B.begin();
 auto i = *(Bit++);
 auto git = g.begin() + (i++);

 for( ; Bit != B.end() ; ++i ) {
  auto h = *(Bit++);
  while( i < h )
   *(git++) = g[ i++ ];
  }

 std::copy( g.begin() + i , g.end() , git );

 }  // end( Compact )

/*--------------------------------------------------------------------------*/

static void set_difference_in_place( BundleSolver::Subset & S1 ,
                                     BundleSolver::c_Subset & S2 )
{
 // removes from S1 all elements in S2, resizing it accordingly
 // both S1 and S2 are assumed to be ordered and with unique elements

 if( S1.empty() )  // nothing to delete from
  return;          // nothing to do

 auto S1it = S1.begin();
 auto S2it = S2.begin();

 // first phase: find the first element present in both S1 and S2

 for( ; ; ) {
  while( ( S1it != S1.end() ) && ( *S1it < *S2it ) )
   ++S1it;
  if( S1it == S1.end() )
   break;
  while( ( S2it != S2.end() ) && ( *S1it > *S2it ) )
   ++S2it;
  if( S2it == S2.end() )
   break;
  if( *S1it == *S2it )
   break;
  }

 if( ( S1it == S1.end() ) || ( S2it == S2.end() ) ) // if there are none
  return;                                           // all done

 // now S1it points to the first element in S1 == than the first in S2
 // elements in S1 after the common one(s) will have to be moved
 auto S1wit = S1it++;  // skip the first equal element
 S2it++;

 for( ; ( S1it != S1.end() ) && ( S2it != S2.end() ) ; ) {
  while( ( S1it != S1.end() ) && ( *S1it < *S2it ) )
   *(S1wit++) = *(S1it++);
  if( S1it == S1.end() )
   break;
  while( ( S2it != S2.end() ) && ( *S1it > *S2it ) )
   ++S2it;
  if( S2it == S2.end() )
   break;
  if( *S1it == *S2it ) { ++S1it; ++S2it; }
  }

 while( S1it != S1.end() )  // copy the part remaining after the end of S2
  *(S1wit++) = *(S1it++);

 S1.resize( std::distance( S1.begin() , S1wit ) );

 }  // end( set_difference_in_place )

/*--------------------------------------------------------------------------*/

static void set_union_in_place( BundleSolver::Subset & S1 ,
                                BundleSolver::c_Subset & S2 )
{
 // make S1 to be the union of S1 and S2
 if( S2.empty() )
  return;

 if( S1.empty() )
  S1 = S2;
 else {
  BundleSolver::Subset tmp;
  std::set_union( S1.begin() , S1.end() , S2.begin() , S2.end() ,
                  std::back_inserter( tmp ) );
  S1 = std::move( tmp );
  }
 }  // end( set_union_in_place )

/*--------------------------------------------------------------------------*/

static void set_union_in_place( BundleSolver::Subset & S1 ,
                                BundleSolver::Subset && S2 )
{
 // make S1 to be the union of S1 and S2, if useful destroy S2 in the process
 if( S2.empty() )
  return;

 if( S1.empty() )
  S1 = std::move( S2 );
 else {
  BundleSolver::Subset tmp;
  std::set_union( S1.begin() , S1.end() , S2.begin() , S2.end() ,
                  std::back_inserter( tmp ) );
  S1 = std::move( tmp );
  }
 }  // end( set_union_in_place )

/*--------------------------------------------------------------------------*/

static double norm( const BundleSolver::Vec_VarValue & v , char t )
{
 double res = 0;
 if( t == 0 ) {    // INF-norm
  for( auto el : v )
   if( std::abs( el ) > res )
    res = std::abs( el );
  }
 else
  if( t == 1 )    // 1-norm
   for( auto el : v )
    res += std::abs( el );
  else {          // 2-norm
   for( auto el : v )
    res += el * el;

   res = sqrt(  res );
   }

 return( res );
 }

/*--------------------------------------------------------------------------*/

static void vect_sum( BundleSolver::Vec_VarValue & v1 , double * v2 )
{
 for( auto & el : v1 )
  el += *(v2++);
 }

/*--------------------------------------------------------------------------*/

static void chgsign( double * v , Index n )
{
 for( const auto ev = v + n ; v < ev ; ++v )
  *v = - *v;
 }

/*--------------------------------------------------------------------------*/

static std::string ps_insert( const std::string & name ,
                              const std::string & insert )
{
 auto pos = name.rfind('.');
 if( pos != std::string::npos )
  return( name.substr( 0 , pos ) + insert + name.substr( pos ) );
 else
  return( name + insert );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BundleSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0( BundleSolver );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register BundleSolverState to the State factory

SMSpp_insert_in_factory_cpp_0( BundleSolverState );

/*--------------------------------------------------------------------------*/
/*------------- METHODS OF BundleSolver -------------------------*/
/*--------------------------------------------------------------------------*/

BundleSolver::VarValue
BundleSolver::read_alpha_global( Index name ) const
{
 if( MasterPB && name < ItemVcblr.size() ) {
  const auto & loc = ItemVcblr[ name ];
  // ItemVcblr[ name ].second == InINF means "empty slot"; nothing to read
  if( loc.second < Inf< Index >() )
   return( MasterPB->get_alpha( int( loc.first ) , int( name ) ) );
  }
 return( 0 );
 }

/*--------------------------------------------------------------------------*/

BundleSolver::VarValue
BundleSolver::read_theta_global( Index name ) const
{
 if( MasterPB && name < ItemVcblr.size() ) {
  const auto & loc = ItemVcblr[ name ];
  if( loc.second < Inf< Index >() )
   return( MasterPB->get_theta( int( loc.first ) , int( name ) ) );
  }
 return( 0 );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::remove_cut_global( Index name )
{
 if( MasterPB && name < ItemVcblr.size() ) {
  const auto & loc = ItemVcblr[ name ];
  if( loc.second < Inf< Index >() )
   MasterPB->remove_cut( int( loc.first ) , int( name ) );
  }
 }

/*--------------------------------------------------------------------------*/

int BundleSolver::compute( bool changedvars )
{
 // ensure no concurrent accesses - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 lock();                        // ... either from other threads

 if( Result == kStillRunning )  // ... or from the same
  throw( std::logic_error( "BundleSolver::compute() called within itself" )
         );

 // basic sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! f_Block ) {    // no Block is there to compute
  Result = kBlockLocked;

  BundleSolver_error_return:
  f_mutex.unlock();  // unlock the mutex
  return( Result );
  }

 if( MaxIter == 0 ) {  // no iteration must be performed
  Result = kStopIter;
  goto BundleSolver_error_return;
  }

 Result = kStillRunning;    // still working

 // start timer now (so that processing Modification is included) - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 c_start = std::chrono::system_clock::now();

 double tot_time = 0;  // independently keep function evaluation time
 long tot_NrEvls = 0;  // total number of C05Function evaluations

 // first, process any outstanding Modification - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // v_mod is atomically copied in a temporary data structure to be processed,
 // but while the latter happens new Modification may come in; hence,
 // process_outstanding_Modification() may be called more than once

 while( num_outstanding_Modification() ) {
  bool owned = f_Block->is_owned_by( f_id );       // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) ) {  // if not read_lock now
   Result = kBlockLocked;                           // return error on failure
   goto BundleSolver_error_return;
   }

  // process any other Modification
  process_outstanding_Modification();

  if( ! owned )             // if the Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it
  }

 // initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there is a 0-th component and the value had not been computed before
 // (maybe because the 0-th component has changed), do it now. Fi0Lmb is
 // reset to INFshift both at the very first call and whenever a
 // C05FunctionModLin* mutates the linear part's coefficients, so it
 // doubles as the "the master's linear part is out of date" signal below
 const bool lin_recompute = f_lf && ( Fi0Lmb == INFshift );
 if( lin_recompute ) {
  f_lf->compute( true );
  Fi0Lmb = rs( f_lf->get_upper_estimate() );
  if( UpFiLmbdef == NrFi ) {  // ready to compute the total upper bound
   ++UpFiLmbdef;              // do so
   UpFiLmb.back() = std::accumulate( UpFiLmb.begin() , --(UpFiLmb.end()) ,
                                     Fi0Lmb );
   }
  if( LwFiLmbdef == NrFi ) {  // ready to compute the total lower bound
   ++LwFiLmbdef;              // do so
   LwFiLmb.back() = std::accumulate( LwFiLmb.begin() , --(LwFiLmb.end()) ,
                                     Fi0Lmb );
   }
  }

 // hand the constant gradient of the linear 0-th component f_lf to the
 // master problem so that the dual coupling rows pick up the affine drift
 //   z_j = b_j - sum_k sum_i theta^k_i * A^k_{i,j}
 // This must be (re)done at the first call and again every time the
 // LinearFunction's coefficients change (e.g. a "change linear objective"
 // oracle Modification): lin_recompute, captured above, is exactly that
 // signal. A stale linear part silently corrupts the dual coupling RHS and
 // hence the search direction. For problems treated as concave maxima
 // (f_convex == false) the gradient is flipped in sign to match the
 // convex-min convention used inside MasterProblemBlock
 if( MasterPB && f_lf && ( ( ! f_linear_part_set ) || lin_recompute ) ) {
  std::vector< double > b( NumVar , 0.0 );
  const auto & cf = f_lf->get_v_var();
  for( Index i = 0 ; i < cf.size() && i < NumVar ; ++i )
   b[ i ] = f_convex ? cf[ i ].second : - cf[ i ].second;
  MasterPB->set_linear_part( b );
  f_linear_part_set = true;
  }

 f_wFi = NrFi - 1;  // since not all components are necessarily evaluated
                    // at all iterations, the order in which they are seen
 // may be important; keep track of the last evaluated component so as to
 // proceed round-robin-like across multiple iterations
 double lastETT = 0;  // last "time" eEveryTTime events have been called
 ParIter = 0;         // number of iterations in this call
 ++SCalls;            // one more call
 RifeqFi = ( UpRifFi == UpFiLmb );  // true if the reference values are right

 if( NeedsG1() )
  G1.resize( NrFi );
 else
  G1.clear();

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( ; ; ) {
  // check if time is over- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( MaxTime < INFshift ) && ( get_elapsed_time() > MaxTime ) ) {
   BLOG( 1 , " ~ stop due to max time" << std::endl );
   Result = kStopTime;
   break;
   }

  // run time-periodic events - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( EveryTTm && ( get_elapsed_time() >= lastETT + EveryTTm ) ) {
   for( auto & ev : v_events[ eEveryTTime ] ) {
    auto res = ev();
    if( res == eStopOK ) { Result = kOK; break; }
    if( res == eStopError ) { Result = kError; break; }
    }

   lastETT = get_elapsed_time();  // reset counter
   }

  if( Result != kStillRunning ) {
   BLOG( 1 , " ~ stop due to time-periodic event" << std::endl );
   break;
   }

  // construct the direction d- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!! PrintBundle();

  FormD();

  // some log - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log1();

  // another iteration (master problem solution)- - - - - - - - - - - - - - -

  ++ParIter;

  // check for "bad" termination- - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( Result == kStopTime )  // time is up
   continue;                 // return at the start and stop

  if( Result == kInfeasible ) {  // the Master Problem is infeasible
   const bool level_empty = UsesPrimalMaster();
   if( UsesLevelStabilization() && f_level_initialized &&
       ( f_level_value < INFshift ) && ( UpFiLmb.back() < INFshift ) &&
       ( ! get_bc_size() ) && level_empty ) {
    BLOG( 1 , " ~ level empty: LB = " << def << f_level_value
              << std::endl );
    record_level_lower_bound( f_level_value );
    if( ! refresh_level_after_master( true ) ) {
     BLOG( 1 , " ~ stop (empty level refresh made no progress)"
               << std::endl );
     Result = kLowPrecision;
     break;
     }
    continue;
    }
   BLOG( 1 , " ~ stop (infeasible)" << std::endl );
   break;
   }

  if( Result == kUnbounded ) {  // the Master Problem is unbounded
   const bool level_empty = ! UsesPrimalMaster();
   if( UsesLevelStabilization() && f_level_initialized &&
       ( f_level_value < INFshift ) && ( UpFiLmb.back() < INFshift ) &&
       level_empty ) {
    BLOG( 1 , " ~ level empty: LB = " << def << f_level_value
              << std::endl );
    record_level_lower_bound( f_level_value );
    if( ! refresh_level_after_master( true ) ) {
     BLOG( 1 , " ~ stop (empty level refresh made no progress)"
               << std::endl );
     Result = kLowPrecision;
     break;
     }
    continue;
    }
   BLOG( 1 , " ~ stop (MP unbounded)" << std::endl );
   break;
   }

  if( Result >= kError ) {  // problems in the Master Problem solver
   BLOG( 1 , " ~ error in the MPSolver" << std::endl );
   break;
   }

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( IsOptimal() ) {  // if optimality is detected
   // run optimality events - - - - - - - - - - - - - - - - - - - - - - - - -
   int res = eContinue;
   for( auto & ev : v_events[ eBeforeTermination ] )
    if( ( res = ev() ) != eContinue )
     break;

   if( res == eForceContinue ) {
    BLOG( 1 , " ~ optimal stop aborted by optimality event" << std::endl );
    continue;  // go back to master problem solution
    }

   if( res == eStopError ) {
    BLOG( 1 , " ~ stop (error) by optimality event" << std::endl );
    Result = kError;
    break;
    }

   BLOG( 1 , " ~ stop (optimal)" << std::endl );
   Result = kOK;
   break;
   }

  // run iteration-periodic events- - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( EverykIt && ( ! ( ParIter % EverykIt ) ) )
   for( auto & ev : v_events[ eEverykIteration ] ) {
    auto res = ev();
    if( res == eStopOK ) { Result = kOK; break; }
    if( res == eStopError ) { Result = kError; break; }
    }

  if( Result != kStillRunning ) {
   BLOG( 1 , " ~ stop due to time-periodic event" << std::endl );
   break;
   }

  // check if "ex-ante" Noise Reduction is needed - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ensure that the Sigma* is "not too negative", if it is increase t (if
  // possible) and re-solve the MP; note that this kind of NR only happens if
  // the oracle is "unfaithful", i.e., it pretends to provide information with
  // the required accuracy but in fact it does not
  //
  // however, avoid doing any of this if the linearization errors are not
  // computed w.r.t. the "true" value of UpFiLmb but w.r.t. a "random"
  // reference value, since then the fact that linearization errors are
  // negative is not meaningful

  if( ( ! UsesPureLevelStabilization() ) &&
      RifeqFi && ( vStar.back() < INFshift ) &&
      ( Sigma < - max_error( UpRifFi.back() , RelAcc ) ) &&
      ( Sigma <= - m3 * DST ) ) {
   if( t >= tMaior ) {
    BLOG( 1 , " ~ stop: NR required but t maximum" << std::endl );
    Result = kLowPrecision;
    break;
    }

   t = std::min( t * mxIncr , tMaior );
   BLOG( 2 , " ~ NR: t increased to " << shrt << t << std::endl );
   tHasChgd = true;
   continue;
   }

  // update out-of-base counters- - - - - - - - - - - - - - - - - - - - - - -

  UpdtCntrs();

  // Hard Long-Term t-strategy for quadratic stabilization- - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the hard long-term t-strategy requires t to increase if the step is too
  // small, and therefore has to be checked before the others
  // however, it is only viable under a quadratic stabilization
  //
  // however, avoid doing any of this if the linearization errors are not
  // computed w.r.t. the "true" value of UpFiLmb (the GBS master is always
  // a quadratically stabilized MasterProblemBlock, so the test is
  // unconditional)

  if( ( ! UsesPureLevelStabilization() ) &&
      ( tStar > 0 ) && ( ( tSPar1 & tSP1Msk ) == kHLTTS ) && RifeqFi ) {

   double AFL = std::abs( UpFiLmb.back() );
   if( AFL < 1 )
    AFL = 1;

   if( abs( vStar.back() ) <= tSPar2 * EpsU * AFL ) {
    BLOG( 1 , "small v => increase t" << std::endl << "           " );

    // collect two numbers vc and vl such that v( tNew ) >= vc + tNew * vl
    // we require that v( tNew ) >= vc + tNew * vl = tSPar2 * EpsU * AFL
    // ==> tNew = ( tSPar2 * EpsU * AFL - vc ) / vl

    double vl , vc;
    if( MasterPB )
     MasterPB->sensitivity_analysis( vl , vc );
    else
     { vl = 0; vc = 0; }

    double tt;
    if( - vl < 1e-15 )  // v( t ) is [~] constant ==> D*_t [~]= 0
     tt = tStar;                 // ==> the CP model is [~]bounded
    else
     tt = std::min( tStar , ( tSPar2 * EpsU * AFL * Nearly + vc ) /
                            ( - vl ) );

    if( ( tHasChgd = ( tt != t ) ) ) {
     t = tt;
     continue;         // loop only if t changes
     }
    }
   }  // end if( Hard t-strategy )  - - - - - - - - - - - - - - - - - - - - -
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // compute Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormLambda1( t );

  // update the number of items to be fetched from the oracle - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UpdtaBP3();

  // eliminate outdated info- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // This is done *after* the call to Master->SensitAnals() in the Hard
  // Long-Term t-strategy and to FormLambda1(), because elimination of items
  // from the bundle may make the current solution of the master problem
  // invalid, and therefore all solution information may be lost. In theory
  // this should not happen, since only items "out of base" are eliminated,
  // and therefore the solution remains optimal; however, not all MPSolvers
  // may behave in this respect.

  SimpleBStrat();

  // run the inner loop - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // first some initializations - - - - - - - - - - - - - - - - - - - - - - -
  // all stuff that must be computed/changed inside InnerLoop()

  Alfa1 = 0;
  ScPr1 = NeedsScPr1() ? read_Gid_aggregate() : 0;
  if( NeedsG1() ) {
   G1Norm = INFshift;
   G1.assign( NrFi , double( 0 ) );
   }

  CurrNrEvls.assign( NrFi , Index( 0 ) );
  MPchgs = 0;  // != 0 if the MP is guaranteed to change enough after the
               // insertion of new information to ensure convergence

  auto start = std::chrono::system_clock::now();

  auto cnt = InnerLoop();

  auto end = std::chrono::system_clock::now();
  std::chrono::duration< double > elapsed = end - start;

  tot_time += elapsed.count();
  tot_NrEvls += std::accumulate( CurrNrEvls.begin() , CurrNrEvls.end() , 0 );

  CmptdinL = false;

  // compute DeltaFi- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1.back() == INFshift )
   DeltaFi = INFshift;
  else
   if( UpFiLmb1.back() == -INFshift )
    DeltaFi = -INFshift;
   else
    DeltaFi = UpFiLmb1.back() - UpRifFi.back();

  // update FiBest- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1.back() < UpFiBest ) {
   UpFiBest = UpFiLmb1.back();
   if( MaxSol > 1 )
    LmbdBst = Lambda1;
   }

  // update the "aggregated" Alfa1 and ScPr1- - - - - - - - - - - - - - - - -

  UpdateHeuristicInfo();

  // some log about the newly obtained information- - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log2( tot_time );

  // check whether either any error has occurred or time has expired- - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1.back() == -INFshift ) {
   BLOG2( 1 , f_convex , " ~ stop (Fi = -INF)" << std::endl );
   BLOG2( 1 , ! f_convex , " ~ stop (Fi = INF)" << std::endl );
   Result = kUnbounded;
   break;
   }

  if( Result == kError ) {
   BLOG( 1 , " ~ stop (error)" << std::endl );
   break;
   }

  if( Result == kStopTime )  // time has ran up inside InnerLoop()
   continue;                 // go back at the beginning to stop

  // check for the conditional lower bound- - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( ! TrueLB ) &&
      ( UpFiBest <= LowerBound.back() *
                    ( 1 - ( LowerBound.back() > 0 ? RelAcc : - RelAcc ) ) )
      ) {
   BLOG( 1 , "            FiBest " );
   BLOG2( 1 , f_convex , "< conditional LB" );
   BLOG2( 1 , ! f_convex , "> conditional UB" );
   BLOG( 1 , ": unbounded " << std::endl );
   if( UpFiLmb1.back() < UpFiLmb.back() )  // Lambda1 is better than Lambda
    GotoLambda1();                         // go to Lambda1
   else                                    // if not
    if( ! RifeqFi )    // and the alfas are not computed w.r.t. UpFiLmb
     GotoLambda();     // ensure they are so

   Result = kUnbounded;
   break;
   }

  // avoid the t-changing phase if a vertical linearization has been found- -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is because the vertical linearization making Lambda1 unfeasible,
  // surely "change enough the master problem already"
  // yet, one possible t-strategy would be to set t to the largest value
  // that would have produced a feasible point: t := Alfa1 / ( - ScPr1 )
  // (with Alfa1 and ScPr1 of that particular constraint, though, not the
  // "global" ones)

  if( MPchgs > 1 ) {
   if( ParIter >= MaxIter ) {  // if we have done too many iterations
    BLOG( 1 , " ~ stop due to max iter" << std::endl );
    Result = kStopIter;        // stop already
    break;
    }
   else                        // otherwise
    continue;                  // go to the next one
   }

  // avoid the t-changing phase if the linearization errors are not reliable-
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is because we have firmly established one feasible (finite) upper
  // estimate in Lambda1, which ends the "phase 0" in which the linearization
  // errors were computed against an arbitrary value and starts the "phase 1"
  // in which the real optimization takes place

  if( ( ! RifeqFi ) && ( UpFiLmb1.back() < INFshift ) ) {
   // if we are still in "phase 0", and we just found a point where the
   // function value is finite, end the "phase 0" by immediately jumping
   // there. note that one may expect the thing on the function value to be
   // redundant since any component evaluating to +INF should generate a
   // vertical linearization and therefore set MPchgs = 2, which is acted
   // upon right above, but this may not happen. which is a problem if
   // MPchgs == 0 (but this is acted upon right below) but not otherwise,
   // since a "normal" NS will be done which is the right thing to do
   BLOG( 1 , "            Fi1 defined ==> SS " << std::endl );
   GotoLambda1();              // go to the feasible point
   if( ParIter >= MaxIter ) {  // if this was the last possible iteration
    BLOG( 1 , " ~ stop due to max iter" << std::endl );
    Result = kStopIter;
    break;                     // main loop ends here
    }
   else
    continue;                  // go start the actual minimization of Fi()
   }

  // check if noise reduction has to be done- - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( ! MPchgs ) && ( ! UsesPureLevelStabilization() ) ) {
   if( t >= tMaior ) {
    BLOG( 1 , "            stop: NR required but t maximum" << std::endl );
    Result = kLowPrecision;
    break;
    }
   t = std::min( t * mxIncr , tMaior );
   BLOG( 1 , "            NR: t increased to " << shrt << t << std::endl );
   tHasChgd = true;
   continue;
   }

  // Check if we exceeded the maximum noise reduction steps for the level
  if( UsesPureLevelStabilization() && LevelNRCntr >= MaxLevelNR ) {
   BLOG( 1 , "            stop: NR required but maximum nummber of "
              "level NR has been reached" << std::endl );
   Result = kLowPrecision;
   break;
  }

  // the NS / SS decision - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // note again the "<" in the SS condition below (which means this is ever
  // so slightly stronger than it should), which is there to avoid the
  // condition to fire when UpFiLmb1.back() == INF == UpTrgt

  SSDone = ( UpFiLmb1.back() < UpTrgt ) ? true : false;

  VarValue tt = t , tm = t , tp = t;  // setup for the heuristic t

  if( SSDone ) {  // SS - - - - - - - - - - - - - - - - - - - - - - - - - - -

   BLOG( 1 , std::endl << " SS[" << CSSCntr << "]: DFi = " << shrt );
   if( f_convex ) {
    BLOG( 1 , DeltaFi << def << " ~ Up1(" << UpFiLmb1.back()
              << ") <= UpTrgt(" << UpTrgt << ")" );
    }
   else
    BLOG( 1 , - DeltaFi << def << " ~ Lw1(" << - UpFiLmb1.back()
              << ") >= LwTrgt(" << - UpTrgt << ")" );

   if( ( ! UsesPureLevelStabilization() ) && ( tSPar1 & 1 ) ) {
    tt = Heuristic( tSPar1 >> 6 );
    BLOG( 1 , " ~ Ht = " << shrt << tt );
    }

   if( ( ! UsesPureLevelStabilization() ) && tSPar3 ) {
    tp *= std::abs( tSPar3 );
    if( tSPar3 > 0 )
     tm /= tSPar3;
    }

   const bool gated_level_update = CSSCntr + 1 > MnSSC;
   if( ( ++CSSCntr > MnSSC ) &&
       ( ! UsesPureLevelStabilization() ) ) {
    // due to the fact that the counter has just been increased
    if( ( ( tSPar1 & tSP1Msk ) == kBLTTS )  &&
        ( DSTS <= tSPar2 * Sigma ) && ( CSSCntr < 10 ) ) {  //!! 10!
     // if the "balancing" long-term t-strategy is active and D*_t( 1 )
     // is small already, inhibit t increases (but not small heuristic
     // decreases, if active) unless "too many SS happened"
     BLOG( 1 , " ~ small D*_t( 1 )" );
     tp = t;
     }
    else {
     tm = t * mnIncr;  // minimum significant increase
     tp = t * mxIncr;  // maximum significant increase
     CSSCntr = 0;      // a significant increase happened, reset counter
     }
    }

   BLOG( 1 , std::endl );

   GotoLambda1();
   update_level_after_step( true , gated_level_update );
   CNSCntr = 0;
   CmptdinL = ( cnt == NrFi - NrEasy );
   }
  else {        // NS - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   BLOG( 1 , std::endl << " NS[" << CNSCntr << "]: " );
   BLOG2( 1 , DeltaFi < INFshift , "DFi = " << shrt << rs( DeltaFi )
              <<  " ~ " << def );
   if( f_convex ) {
    BLOG( 1 , "Lw1(" << def << LwFiLmb1.back() << ") >= LwTrgt(" << LwTrgt
              << ")" );
    }
   else
    BLOG( 1 , "Up1(" << - LwFiLmb1.back() << ") <= UpTrgt(" << - LwTrgt
              << ")" );

   if( ( ! UsesPureLevelStabilization() ) && ( tSPar1 & 2 ) ) {
    tt = Heuristic( tSPar1 >> 8 );
    BLOG( 1 , " ~ Ht = " << shrt << tt );
    }

   if( ( ! UsesPureLevelStabilization() ) && tSPar3 ) {
    tm /= std::abs( tSPar3 );
    if( tSPar3 > 0 )
     tp *= tSPar3;
    }

   const bool gated_level_update = CNSCntr + 1 > MnNSC;
   if( ( ++CNSCntr > MnNSC ) &&
       ( ! UsesPureLevelStabilization() ) ) {
    // due to the fact that the counter has just been increased
    if( ( ( ( tSPar1 & tSP1Msk ) == kSLTTS ) ||
          ( ( tSPar1 & tSP1Msk ) == kHLTTS ) ) &&
        ( abs( vStar.back() ) <= tSPar2 * EpsU * max_error() ) ) {
     // if either the "hard" or the "soft" long-term t-strategy is active
     // and v* is small already, inhibit t decreases (but not small
     // heuristic increases, if active)
     BLOG( 1 , " small v" );
     tm = t;
     }
    else
     if( ( ( tSPar1 & tSP1Msk ) == kBLTTS ) &&
         ( tSPar2 * DSTS >= abs( Sigma ) ) ) {
      // if the "balancing" long-term t-strategy is active and D*_t( 1 )
      // is large already, inhibit t decreases (but not small heuristic
      // increases, if active); note that one may add the clause "unless
      // too many NS happened", i.e., "&& ( CNSCntr < 20 )": this version
      // avoids problems which may occur with ill-set tStar or tSPar2, but
      // it may give worse performances with "difficult" problems
      // also note the "abs( Sigma )": Sigma should be positive, but in
      // case it is not the control would always be true irrespectively of
      // the magnitude of tSPar2 and tStar just because of the sign
      BLOG( 1 , " ~ large D*_t( 1 )" );
      tm = t;
      }
     else {
      tm = t * mxDecr;  // maximum significant decrease
      tp = t * mnDecr;  // minimum significant decrease
      CNSCntr = 0;      // a significant decrease happened, reset counter
      }
    }


   BLOG( 1 , std::endl );
   update_level_after_step( false , gated_level_update );
   CSSCntr = 0;

   }   // end else( NS )- - - - - - - - - - - - - - - - - - - - - - - - - - -

  // actually update t- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // if the endgame t-strategy fires (note the "/ 10"!!), the regular
  // t-updating mechanism is superseeded
  if( ( ! UsesPureLevelStabilization() ) &&
      ( tSPar1 & kEGTTS ) &&
      ( UpFiLmb.back() < INFshift ) &&
      ( DSTS < max_error() / 10 ) ) {
    tt = std::max( t * ( mxDecr + mnDecr ) / 2 , tMinor );
    BLOG( 1 , " ~ endgame, t = " << shrt << tt );
    //!! the reverse should also be done: if sigma is small and D*( t* ) is
    //!! large, t should be increased --> but this would happen surely at
    //!! the beginning, it should be done only near the end
    }
  else             // regular update mechanism
   if( tm != tp ){ // if t can change, select it in [ tm , tp ]
    tt = std::min( std::min( tMaior , tp ) ,
                   std::max( std::max( tMinor , tm ) , tt ) );
    tt = std::max( std::min( tp , tt ) , tm );
    tt = std::max( std::min( tMaior , tt ) , tMinor );
    }
   else            // else
    tt = t;        // keep it as it is

  if( ( tHasChgd = ( t != tt ) ) )
   t = tt;

  // check max number of iterations - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ParIter >= MaxIter ) {
   BLOG( 1 , " ~ stop due to max iter" << std::endl );
   Result = kStopIter;
   break;
   }
  }  // end( main loop )

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if necessary, force one last SS to the stability center - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( FrcLstSS & 1 ) && ( ! CmptdinL ) &&
     ( ( Result == kOK ) || ( Result == kStopIter ) ||
       ( Result == kLowPrecision ) ) ) {
  BLOG( 1 , "            Recomputing the current point" << std::endl );

  UpFiLmb1 = UpFiLmb;
  LwFiLmb1 = LwFiLmb;
  Fi0Lmb1 = Fi0Lmb;
  UpTrgt = UpFiLmb1.back();
  LwTrgt = LwFiLmb1.back();

  FiStatus.assign( NrFi , kUnEval );
  for( Index i = 0 ; i < NumVar ; i++ )
   LamVcblr[ i ]->set_value( Lambda[ i ] );

  // note that Alfa1, ScPr1, G1 are computed inside GetGi() that is not
  // called inside this call to InnerLoop(), so they are not initialised
  CurrNrEvls.assign( NrFi , Index( 0 ) );
  MPchgs = 0;  // != 0 if the MP is guaranteed to change enough after the
               // insertion of new information to ensure convergence

  auto start = std::chrono::system_clock::now();

  auto cnt = InnerLoop( true );

  auto end = std::chrono::system_clock::now();
  std::chrono::duration< double > elapsed = end - start;

  tot_time += elapsed.count();

  CmptdinL = ( cnt == NrFi - NrEasy );

  // not being able to compute all non-easy components is an error
  if( ( ! CmptdinL ) && ( Result != kStopTime ) )
   Result = kError;
  }

 // final printouts - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_log && ( LogVerb >= 1 ) ) {
  *f_log << std::endl << "Call " << SCalls << ": "  << fixd << ParIter
         << " ~ " << tot_NrEvls  << " ~ " << get_elapsed_time() << " ~ "
         << tot_time << " -> ";
  switch( Result ) {
   case( kOK ):           *f_log << "optimal"; break;
   case( kStopTime ):     *f_log << "max time"; break;
   case( kStopIter ):     *f_log << "max iter"; break;
   case( kInfeasible ):   *f_log << "infeasible"; break;
   case( kUnbounded ):    *f_log << "unbounded"; break;
   case( kLowPrecision ): *f_log << "inexact oracle"; break;
   default:               *f_log << "error";
   }
  if( ( Result != kInfeasible ) && ( Result != kUnbounded ) ) {
   *f_log << " ~ Fi* = " << def;
   pval( *f_log , rs( UpRifFi.back() ) );
   }
  *f_log << std::endl;
  }

 unlock();  // unlock the mutex

 return( Result );

 }  // end( BundleSolver::compute )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::set_Block( Block * block )
{
 if( f_Block == block )  // registering to the same Block
  return;                // cowardly and silently return

 if( f_Block ) {  // changing from a previous oracle - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  guts_of_destructor();   // deallocate memory
  }

 Solver::set_Block( block );  // attach to the new Block

 if( ! f_Block )  // that was actually clearing the Block
  return;         // all done

 /* Immediately create the MasterProblemBlock. This is needed because
  * in this phase two things will happen involving MPB:
  *   - it will populate the structure needed for the master problem
  *     using all the information coming from the Block;
  *   - it will silently register itself as father of the easy components.
  *     This is needed because MPB will have to respond to all the
  *     Modifications affecting such components, and hence must be informed.
 .*/
 CreateMPB();

 // lock the Block - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
             "BundleSolver::set_Block: unable to lock the Block" ) );

 // generate the abstract representation
 f_Block->generate_abstract_variables();
 f_Block->generate_abstract_constraints();
 f_Block->generate_objective();

 /* Two types of block can be handled by the BundleSolver:

     1. Only one single non-smooth function
     2. a sum of some non-smooth functions

    The algorithm here developed aims at solving non-constrained non-smooth
    optimization. The block can have box constraints at the most.
    It is expected the block to have in the first case a FRealObjective
    whose the function is a C05Function one and having no children, while in
    the second case a FRealObjective one whose the function is a
    LinearFunction and having as many sub-blocks as the number of components.
    In the latter case, each sub-block must not contain any Variable or
    Constraint. cVariable may have a lower and upper bound. If the lower
    bound  has a finite value, it must be 0. */

 const auto & sb = f_Block->get_nested_Blocks();

 if( sb.empty() ) {  // no sub-Block
  // the objective function of the Block must be a C05Function  - - - - - - -

  auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
  if( ! obj )
   throw( std::invalid_argument(
              "BundleSolver::set_Block: "
              "objective is not a FRealObjective" ) );

  auto c05f = dynamic_cast< C05Function * >( obj->get_function() );
  if( ! c05f )
   throw( std::invalid_argument(
              "BundleSolver::set_Block: "
              "the objective is not a C05Function" ) );

  f_convex = c05f->is_convex();
  if( ( ! f_convex ) && ( ! c05f->is_concave() ) )
   throw( std::invalid_argument(
              "BundleSolver::set_Block: "
              "only convex or concave objectives allowed" ) );

  if( ( f_convex && ( obj->get_sense() == Objective::eMax ) ) ||
      ( ( ! f_convex ) && ( obj->get_sense() == Objective::eMin ) ) )
   throw( std::invalid_argument(
              "BundleSolver::set_Block: "
              "can only minimize convex / maximize concave" ) );
  v_c05f.push_back( c05f );
  f_lf = nullptr;
  }
 else {  // there are sub-Block
  // the objective function of the block must be a LinearFunction- - - - - - -

  if( ! f_Block->get_objective() )  // there is no Objective
   f_lf = nullptr;
  else {
   auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
   if( ! obj )
    throw( std::logic_error(
               "BundleSolver::set_Block: "
               "the objective is not a real function" ) );

   if( ! obj->get_function() ) { // the FRealObjective has no Function
    f_lf = nullptr;
    }
   else {
    f_lf = dynamic_cast< LinearFunction * >( obj->get_function() );
    if( ! f_lf )
     throw( std::logic_error(
                "BundleSolver::set_Block: "
                "the objective is not a LinearFunction" ) );

    if( ! f_lf->get_num_active_var() )  // the LinearFunction has no Variable
     f_lf = nullptr;

    }
   }

  v_c05f.resize( sb.size() );

  for( Index i = 0 ; i < sb.size() ; ++i ) {  // for each sub-block
   // the objective function of each sub-block must be a C05Function - - - - -

   auto obj = dynamic_cast< FRealObjective * >( sb[ i ]->get_objective() );
   if( ! obj )
    throw( std::logic_error(
               "BundleSolver::set_Block: sub-Block "
               "objective is not a real function" ) );

   auto c05f = dynamic_cast< C05Function * >( obj->get_function() );
   if( ! c05f )
    throw( std::logic_error(
               "BundleSolver::set_Block: sub-Block "
               "objective is not a C05Function" ) );

   // all have to be the same convexity - - - - - - - - - - - - - - - - - - -
   v_c05f[ i ] = c05f;
   if( ! i )
    f_convex = c05f->is_convex();
   else
    if( c05f->is_convex() != f_convex )
     throw( std::invalid_argument(
                "BundleSolver::set_Block: "
                "objectives must be all convex or all concave" ) );

   // all have to be max/min in the right way- - - - - - - - - - - - - - -
   if( ( ! f_convex ) && ( ! c05f->is_concave() ) )
    throw( std::invalid_argument(
               "BundleSolver::set_Block: "
               "only convex or concave objectives allowed" ) );

   if( ( f_convex && ( obj->get_sense() == Objective::eMax ) ) ||
       ( ( ! f_convex ) && ( obj->get_sense() == Objective::eMin ) ) )
    throw( std::invalid_argument(
               "BundleSolver::set_Block: "
               "can only minimize convex / maximize concave" ) );

   // nephews are not allowed- - - - - - - - - - - - - - - - - - - - - - - - -
   if( sb[ i ]->get_nested_Blocks().size() )
    throw( std::logic_error(
               "BundleSolver::set_Block: "
               "nephew sub-Blocks are not allowed" ) );

   // Variable not allowed - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( sb[ i ]->get_static_variables().size() )
    throw( std::logic_error(
               "BundleSolver::set_Block: "
               "static Variable in sub-Block are not allowed" ) );

   if( sb[ i ]->get_dynamic_variables().size() )
    throw( std::logic_error(
               "BundleSolver::set_Block: "
               "dynamic Variable in sub-Block are not allowed" ) );

   // neither are Constraint - - - - - - - - - - - - - - - - - - - - - - - - -
   if( sb[ i ]->get_static_constraints().size() )
    throw( std::logic_error(
               "BundleSolver::set_Block: "
               "static Constraint in sub-Block are not allowed" ) );

   if( sb[ i ]->get_dynamic_constraints().size() )
    throw( std::logic_error(
               "BundleSolver::set_Block: "
               "dynamic Constraint in sub-Block are not allowed" ) );
   }  // end( for each sub-Block )
  }  // end( there are sub-Block )

 // build LamVcblr as the union of "active" Variables across all v_c05f[ h ]
 // (and f_lf, if any), in first-encounter order. Each v_c05f[ h ] is
 // allowed to expose either the full union (dense path) or a strict
 // subset (sparse path) of LamVcblr. v_local2global[ h ] records,
 // for each h, the index in LamVcblr of h's i-th active Variable in the
 // order get_linearization_coefficients writes them, and is left empty
 // when the sparse path is not needed.- - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Lambda2Idx.clear();
 LamVcblr.clear();
 v_ref_count.clear();
 v_local2global.assign( v_c05f.size() , {} );
 f_sparse_lambda = false;

 auto register_active = [ & ]( C05Function * f ,
                               std::vector< Index > * map_out ) {
  const Index loc_NV = f->get_num_active_var();
  if( map_out )
   map_out->reserve( loc_NV + 1 );
  auto it_end = std::as_const( f )->end();
  for( auto it_v = std::as_const( f )->begin() ; it_v != it_end ; ++it_v ) {
   auto p = static_cast< ColVariable * >( & ( *it_v ) );
   auto [ jt , inserted ] = Lambda2Idx.try_emplace( p , LamVcblr.size() );
   if( inserted ) {
    LamVcblr.push_back( p );
    v_ref_count.push_back( 1 );
    }
   else
    ++v_ref_count[ jt->second ];
   if( map_out )
    map_out->push_back( jt->second );
   }
  };

 if( f_lf )
  register_active( f_lf , nullptr );

 for( Index i = 0 ; i < v_c05f.size() ; ++i )
  register_active( v_c05f[ i ] , & v_local2global[ i ] );

 NumVar = LamVcblr.size();

 // detect whether any v_c05f[ h ] (or f_lf) exposes a proper subset of
 // LamVcblr, or the full set in a non-identity order. In both cases
 // switch to the sparse Lambda path.
 if( f_lf && f_lf->get_num_active_var() != NumVar )
  f_sparse_lambda = true;
 for( Index h = 0 ; ! f_sparse_lambda && h < v_local2global.size() ; ++h ) {
  const auto & m = v_local2global[ h ];
  if( m.size() != NumVar ) {
   f_sparse_lambda = true;
   break;
   }
  for( Index i = 0 ; i < NumVar ; ++i )
   if( m[ i ] != i ) {
    f_sparse_lambda = true;
    break;
    }
  }

 if( f_sparse_lambda ) {
  // sanity: f_lf must cover the full LamVcblr in identity order when
  // sparse is engaged, otherwise the f_lf gather paths would also need
  // translation. Typical sparse-producing callers leave f_lf == nullptr
  // (the linear term is empty), so this check is mostly defensive.
  if( f_lf ) {
   if( f_lf->get_num_active_var() != NumVar )
    throw( std::logic_error(
               "BundleSolver::set_Block: "
               "sparse Lambda mode requires f_lf to cover "
               "the full union of active variables" ) );
   auto v = f_lf->begin();
   for( auto vi = LamVcblr.begin() ; vi != LamVcblr.end() ; ++v , ++vi )
    if( static_cast< ColVariable * >( & ( *v ) ) != *vi )
     throw( std::logic_error(
                "BundleSolver::set_Block: "
                "sparse Lambda mode requires f_lf to follow "
                "the LamVcblr order" ) );
   }

  // each v_c05f[ h ]'s active Variables must be presented in an order
  // that is strictly increasing with respect to their position in
  // LamVcblr (= first-encounter order in the union across all v_c05f
  // and f_lf). If the caller broke this invariant, throw rather than
  // silently sort the dual pairs. The terminator Inf< Index >() at the
  // end of each map makes its .data() usable as an Inf-terminated index
  // array by the gather paths.
  for( Index h = 0 ; h < v_local2global.size() ; ++h ) {
   auto & m = v_local2global[ h ];
   for( Index li = 1 ; li < m.size() ; ++li )
    if( m[ li ] <= m[ li - 1 ] )
     throw( std::logic_error(
                "BundleSolver::set_Block: sparse Lambda mode: "
                "v_c05f[" + std::to_string( h ) + "] active Variables are "
                "not in strictly increasing LamVcblr order; the caller "
                "must present dual pairs sorted by global Variable "
                "position" ) );
   m.push_back( Inf< Index >() );
   }
  }
 else {
  // dense path: drop the per-component maps, the global lookup, and the
  // refcount. This is just defensive — the maps would all be the
  // identity, but we avoid keeping ~ f_nsb * NumVar of redundant Index
  // data live.
  v_local2global.clear();
  v_local2global.shrink_to_fit();
  Lambda2Idx.clear();
  Lambda2Idx.rehash( 0 );  // shrink the bucket array to 0
  v_ref_count.clear();
  v_ref_count.shrink_to_fit();
  }

 // if some Variable are present, they are of the ColVariable type - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index NumBVar = 0;  // count the number of Variable in the Block
 auto v_s_Variable = f_Block->get_static_variables();
 for( auto & el : v_s_Variable ) {
  auto sz = un_any_thing_count_static( ColVariable , el );
  if( sz == Inf< std::size_t >() )
   throw( std::logic_error(
              "BundleSolver::set_Block: "
              "some static Variable is not a ColVariable" ) );
  NumBVar += sz;
  }

 auto v_d_Variable = f_Block->get_dynamic_variables();
 for( auto & el : v_d_Variable ) {
  auto sz = un_any_thing_count_dynamic( ColVariable , el );
  if( sz == Inf< std::size_t >() )
   throw( std::logic_error(
              "BundleSolver::set_Block: "
              "some dynamic Variable is not a ColVariable" ) );
  NumBVar += sz;
  }

 if( NumBVar < NumVar )
  throw( std::logic_error(
             "BundleSolver::set_Block: "
             "too few ColVariable in the Block" ) );

 // check that the Variable in the Block agree with that in the C05Function- -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector< ColVariable * > LamBVcblr( NumBVar );

 Index cnt = 0;
 for( auto & el : v_s_Variable )
  un_any_static( el , [ & ]( ColVariable & sv ) { LamBVcblr[ cnt++ ] = & sv;
                  } , un_any_type< ColVariable >() );

 for( auto & el : v_d_Variable )
  un_any_dynamic( el , [ & ]( ColVariable & sv ) { LamBVcblr[ cnt++ ] = & sv;
                  } , un_any_type< ColVariable >() );

 std::sort( LamBVcblr.begin() , LamBVcblr.end() ); // These are the block variables

 std::vector< ColVariable * > LamVcblrO( LamVcblr ); // These are the Objective variables (only of the first subblock)
 std::sort( LamVcblrO.begin() , LamVcblrO.end() );

 if( ! std::includes( LamBVcblr.begin() , LamBVcblr.end() ,
                      LamVcblrO.begin() , LamVcblrO.end() ) )
  throw( std::logic_error(
             "BundleSolver::set_Block: "
             "some ColVariable in C05Function are not in the Block" ) );

 LamVcblrO.clear();
 LamBVcblr.clear();

 // if some Constraints are present, their can only be either BoxConstraints
 // (with LHS == 0), LB0Constraints or NNConstraints - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // one day general linear constraints will be allowed
 //
 // note that un_any_thing_*() only serves to verify that the stuff is of the
 // right type, and therefore it has to do nothing; this is obtained by
 // passing it as, the "function" argument, a void --> void lambda doing
 // nothing immediately applied to nothing, cue the curios list of
 // parentheses "[](){}()"

 for( auto & el : f_Block->get_static_constraints() ) {
  if( un_any_thing_static( BoxConstraint , el , [](){}() ) )
   continue;
  if( un_any_thing_static( LB0Constraint , el , [](){}() ) )
   continue;
  if( un_any_thing_static( NNConstraint , el , [](){}() ) )
   continue;
  //!! this should never have been needed in the first place
  //!!if( un_any_const_static( el , []( BoxConstraint & b ){} ,
  //!!                         un_any_type< BoxConstraint >() ) )
  //!! continue;
  throw( std::logic_error(
             "BundleSolver::set_Block: "
             "unsupported type of static Constraint" ) );
  }

 for( auto & el : f_Block->get_dynamic_constraints() ) {
  if( un_any_thing_dynamic( BoxConstraint , el , [](){}() ) )
   continue;
  if( un_any_thing_dynamic( LB0Constraint , el , [](){}() ) )
   continue;
  if( un_any_thing_dynamic( NNConstraint , el , [](){}() ) )
   continue;
  throw( std::logic_error(
             "BundleSolver::set_Block: "
             "unsupported type of dynamic Constraint" ) );
  }

 // read information about the C05Function - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NrFi = v_c05f.size();

 // check if there are "easy" components and deal with them- - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NrEasy = 0;
 if( DoEasy ) {
 // retrieve the ComputeConfig for the "easy" components, if any
 ComputeConfig * eCC = nullptr;
 if( ! EasyCfg.empty() ) {
  auto cfg = Configuration::deserialize( EasyCfg );
  if( ! ( eCC = dynamic_cast< ComputeConfig * >( cfg ) ) )
   delete cfg;
  }

 IsEasy.resize( NrFi , 0 );
 auto NEit = NoEasy.begin();
 for( Index k = 0 ; k < NrFi ; ++k ) {
  // check if component k is marked as being forbidden to be easy (this
  // assumes NoEasy ordered in increasing sense and without duplications)
  if( ( NEit != NoEasy.end() ) && ( Index( *NEit ) == k ) ) {
   ++NEit;
   continue;
  }

   auto LagB = dynamic_cast< LagBFunction * >( v_c05f[ k ] );
   if( LagB ) {
    auto MILPs = new MILPSolver();
    try {  // check if the inner Block of the LagBFunction is all-linear
     // do this by trying to register the MILPSolver to the inner Block; if
     // the operation succeeds than the component may be easy (provided that
     // also all variables are continuous), otherwise it surely is not,
     // which is captured by the fact that exception is thrown; note that
     // [MILP]Solver::set_Block() does *not* call Block::register_Solver(),
     // which therefore may have to be done later
     MILPs->set_Block( LagB->get_inner_block() );
     // the component is easy only if it is a real LP: all variables are
     // continuous and there is no quadratic constraint (the master problem
     // cannot absorb quadratic rows; also, with them the coefficient matrix
     // is stored row-wise, while GetBDesc() hands the master the column-wise
     // description, see GetBNC())
     if( ( ! MILPs->get_num_integer_vars() ) &&
	 ( ! MILPs->get_numquadrows() ) ) {
      IsEasy[ k ] = MILPs;
      ++NrEasy;

     // the master MP cannot carry the constant term of an easy
     // component on its own, so it is folded into constant_value and
     // added back when reporting Fi values
     constant_value += LagB->get_constant_term();
     }
    }
   catch( ... ) {  // exception means that something nonlinear is there
    }
   delete MILPs;
   }
  }

 if( ! NrEasy )
  IsEasy.clear();
 else {
  if( NrEasy == NrFi )
   throw( std::logic_error(
          "BundleSolver: all components are easy, this is no supported" ) );

 // ComputeConfig-ure the easy components: clone eCC once per component
 // that needs the shared template, except for the last such component
 // which receives eCC by ownership transfer (saving one clone + delete)
 Index last_using_eCC = NrFi;
 if( eCC )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( IsEasy[ k ] &&
       ( k >= Index( CmpCfg.size() ) || CmpCfg[ k ].empty() ) )
    last_using_eCC = k;

 if( eCC || ( ! CmpCfg.empty() ) )
  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( ! IsEasy[ k ] )
    continue;

   ComputeConfig * cfg = nullptr;
   if( ( k < Index( CmpCfg.size() ) ) && ( ! CmpCfg[ k ].empty() ) ) {
    auto tcfg = Configuration::deserialize( CmpCfg[ k ] );
    if( ! ( cfg = dynamic_cast< ComputeConfig * >( tcfg ) ) )
      delete tcfg;
   }

   if( ( ! cfg ) && eCC ) {
    if( k == last_using_eCC ) {
     cfg = eCC;
     eCC = nullptr;
     }
    else
     cfg = eCC->clone();
    }

   if( cfg )
    v_c05f[ k ]->set_ComputeConfig( cfg );
  }

 }

 delete eCC;  // no-op when ownership has been transferred above

 }  // end( if( DoEasy ) )

 // configure all non-easy components- - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if the non-easy ComputeConfig is provided, apply it.
 //
 // in all cases, set the global pool size, *after* having configured them.
 // this is necessary in that with BPar2 == 0, BundleSolver just takes
 // whatever size of the global pool it finds in the C05Function (that may
 // have been just set by the non-easy ComputeConfig). with BPar2 > 0,
 // instead, BundleSolver ensures that the size of the global pool is *at
 // least* BPar2 by increasing it if it is below. this means that:
 // - BundleSolver never *decreases* the size of the global pool
 // - BundleSolver only uses the first BPar2 linearizations in each global
 //   pool; if there are more, the other ones are ignored
 //
 // meanwhile, also set the accuracy of multipliers: the GBS master is a
 // generic [MILP]Solver attached to a MasterProblemBlock, whose effective
 // accuracy is bounded by the back-end (e.g. ~1e-8 for Gurobi barrier);
 // we therefore give the C05Function some slack with 1e-7
 const double eps = 1e-7;

 vBPar2.resize( NrFi + 1, 0 );
 InvItemVcblr.resize( NrFi );

 // retrieve the ComputeConfig for the non-easy components, if any
 ComputeConfig * hCC = nullptr;
 if( ! HardCfg.empty() ) {
  auto cfg = Configuration::deserialize( HardCfg );
  if( ! ( hCC = dynamic_cast< ComputeConfig * >( cfg ) ) )
   delete cfg;
  }

 // mirror of the trick used above for eCC: the last hard component that
 // would borrow hCC takes it by ownership instead, sparing one clone
 Index last_using_hCC = NrFi;
 if( hCC )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ( ! NrEasy || ! IsEasy[ k ] ) &&
       ( k >= Index( CmpCfg.size() ) || CmpCfg[ k ].empty() ) )
    last_using_hCC = k;

 for( Index k = 0 ; k < NrFi ; ++k ) {
  if( NrEasy && IsEasy[ k ] )
   continue;

  // ComputeConfig-ure the non-easy component
  ComputeConfig * cfg = nullptr;
  if( ( k < Index( CmpCfg.size() ) ) && ( ! CmpCfg[ k ].empty() ) ) {
   auto tcfg = Configuration::deserialize( CmpCfg[ k ] );
   if( ! ( cfg = dynamic_cast< ComputeConfig * >( tcfg ) ) )
    delete tcfg;
   }

  if( ( ! cfg ) && hCC ) {
   if( k == last_using_hCC ) {
    cfg = hCC;
    hCC = nullptr;
    }
   else
    cfg = hCC->clone();
   }

  if( cfg )
   v_c05f[ k ]->set_ComputeConfig( cfg );

  // ensure that the accuracy of multipliers is at least eps
  if( v_c05f[ k ]->get_dbl_par( C05Function::dblAAccMlt ) > eps )
   v_c05f[ k ]->set_par( C05Function::dblAAccMlt , eps );

  // manage the global pool size
  Index gps = v_c05f[ k ]->get_int_par( C05Function::intGPMaxSz );
  if( BPar2 == 0 ) {  // use the current global pool size
   if( gps < 2 )
    throw( std::logic_error(
               "BundleSolver::set_Block: "
               "BPar2 == 0 but global pool is too small" ) );
   vBPar2.back() += gps;
   vBPar2[ k ] = gps;
   }
  else {              // force the global pool size to be *at least* BPar2
   if( gps < BPar2 )
    v_c05f[ k ]->set_par( C05Function::intGPMaxSz , int( BPar2 ) );
   vBPar2.back() += BPar2;
   vBPar2[ k ] = BPar2;
   }

  InvItemVcblr[ k ].resize( vBPar2[ k ] , InINF );
  }

 delete hCC;  // no-op when ownership has been transferred above

 // set the component-specific string parameters, if any - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_C05_SPAR_Names.empty() ) {
  if( v_C05_SPAR_Names.size() > v_C05_SPAR_Vals.size() )
   throw( std::logic_error(
              "BundleSolver::set_Block: "
              "vstr_C05_SPAR_Names.size() > vstr_C05_SPAR_Vals.size()" ) );

  for( Index k = 0 ; k < NrFi ; ++k ) {
   ComputeConfig Ck;
   Ck.set_diff( true );
   Ck.set_relax( true );
   auto Vit = v_C05_SPAR_Vals.begin();
   for( const auto & name : v_C05_SPAR_Names ) {
    auto par = ps_insert( *(Vit++) , "_" + std::to_string( k ) );
    if( ( name.size() > 4 ) && ( name.substr( 0 , 4 ) == "vstr" ) )
     Ck.set_par( std::string( name ) ,
                 std::vector< std::string >( { par } ) );
    else
     Ck.set_par( std::string( name ) , std::move( par ) );
    }

   v_c05f[ k ]->set_ComputeConfig( & Ck );
   }
  }

 // allocate memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 t = tInit;
 Prevt = INFshift;

 Lambda.resize( NumVar );    // the default starting point
 Lambda1.resize( NumVar );   // the tentative point

 if( MaxSol > 1 )  // best point found so far
  LmbdBst.resize( NumVar );

 OOBase.resize( vBPar2.back() , Inf< SIndex >() );
 // counter for eliminating outdated items: Inf< SIndex >() means empty

 ItemVcblr.resize( vBPar2.back() , std::make_pair( InINF , InINF ) );

 NrItems.resize( NrFi + 1 , 0 );
 FrFItem.resize( NrFi , 0 );
 MaxItem.resize( NrFi , 0 );
 FictLB.assign( NrFi , false );  // no fictitious LB installed yet

 FreList = {};
 whisZ.resize( NrFi , InINF );
 Zvalid.resize( NrFi , false );

 CurrNrEvls.resize( NrFi , 0 );

 FiStatus.resize( NrFi , kUnEval );
 TrueLB = false;

 UpFiBest = INFshift;      // best, ...
 UpRifFi.resize( NrFi + 1 , 0 );  // and reference Fi() values
 RifeqFi = false;                 // reference values != UpFiLmb
 UpFiLmb1.resize( NrFi + 1 );     // upper and lower function value
 LwFiLmb1.resize( NrFi + 1 );     // ... at the tentative point
 UpFiLmb.resize( NrFi + 1 ,  INFshift );  // upper
 LwFiLmb.resize( NrFi + 1 , -INFshift );  // ... and lower Fi-value
 UpFiLmbdef = LwFiLmbdef = 0;             // ... at the current point
 LowerBound.resize( NrFi + 1 , -INFshift );  // global lower bounds
 f_global_LB = -INFshift;         // algorithmic global LB

 vStar.resize( NrFi + 1 , 0 );
 whisG1.resize( NrFi , InINF );  // no representative yet

 Result = kError;
 SSDone = false;

 // initialize the MasterProblemBlock - - - - - - - - - - - - - -  - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 InitMPB();
 reset_level_stabilization();

 // install the per-coordinate box  L <= Lambda <= U  on the dual master, so
 // that the proximal direction d* = Lambda1 - Lambda is projected onto the
 // feasible Lambda region (exactly the box FormLambda1 clamps Lambda1
 // against). Without it the dual MasterProblemBlock returns the
 // *unconstrained* d* = -t z*, which inflates || z* || (and hence v*) and
 // stalls the bundle on finite optima: the legacy NDOFiOracle MPSolvers
 // receive the same box from the NDOSolver interface and project d*. A
 // effective bounds combine the ColVariable's own bounds with any supported
 // active bound constraint; set_box() treats any non-finite entry as "no bound"
 // (the matching slack stays fixed to 0)
 if( MasterPB ) {
  std::vector< double > Lbox( NumVar ) , Ubox( NumVar );
  for( Index i = 0 ; i < NumVar ; ++i ) {
   const auto bounds = effective_bounds( LamVcblr[ i ] );
   Lbox[ i ] = bounds.first;
   Ubox[ i ] = bounds.second;
   }
  MasterPB->set_box( Lbox , Ubox );
  }

 // After the physical representation have been cretaed into the MPB, register
 // the solver. This will allow also solver using abstract representation to
 // correctly generate it.
 MasterPB->register_Solver( std::string( MPBSolverCfg ) );

 // cleanup easy-component caches - - - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now that the MPSolver has read all the information it needs out of the
 // MILPSolver, cleanup all un-necessary data; note that clear_problem() tells
 // which data to delete with exactly the same bit mapping as DoEasy tells
 // which to keep, hence an xor works; the first bit gets zeroed (hence
 // ultimately set to 1) so that the coefficient matrix is always discarded,
 // as changing it is not managed by the MPSolver (and neither really is by
 // MILPSolver, currently)
 //
 // note that if "easy" components are fully static one may want to get rid
 // of the MILPSolver entirely; however, if new variables are added then
 // GetADesc() can be called, which relies on the MILPSolver. getting rid of
 // the MILPSolver would require knowing that the variables set is also
 // static, which would require one parameter; doable, but not now

 if( NrEasy && MasterPB ) {
  // The easy sub-Block is solved as part of the master MP, so the bulk
  // of the cleanup happens through the master Solver. No dedicated
  // per-easy clear_problem is needed here: the master Solver will
  // re-ingest the easy sub-Block on the next compute() through the
  // Modification it has already received.
  }

 // reset algorithm  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this has to be done before the next step since it sets Lambda, and
 //       if linearizations are added to the bundle their linearization error
 //       depends on Lambda

 ReSetAlg( RstAlgPrm );  // Lambda is reset inside

 // deal with existing linearizations- - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if BPar7 & 8, read from all the C05Function and immediately add to the
 // master problem each and every linearization found in their global pools
 //
 // otherwise read from all the C05Function and mark into InvItemVcblr each
 // and every linearization found in their global pools
 //
 // note that if ( BPar7 & 3 ) >= 2, BundleSolver will happily delete from
 // the global pool any linearization it deletes from the bundle; yet we do
 // not immediately delete existing linearizations from the global pools here.
 // these will likely be overwritten during the optimization, and if memory
 // is a problem they can be cleaned up by the user before set_Block() is
 // called. besides, in many scenarios there will be no linearizations anyway
 //
 // however, if ( BPar7 & 3 ) == 3 then BundleSolver does not care at all
 // about what linearizations are there in the global pool because it will
 // treat any position in the global pool as available for it regardless to
 // if there is anything there; hence, in this case we do not bother to even
 // look

 if( BPar7 & 8 ) {
  for( Index k = 0 ; k < NrFi ; ++k )
   for( Index i = 0 ; i < vBPar2[ k ] ; ++i )
    if( v_c05f[ k ]->is_linearization_there( i ) )
     add_to_bundle( k , i );
  }
 else
  if( ( BPar7 & 3 ) < 3 ) {
   for( Index k = 0 ; k < NrFi ; ++k )
    for( Index i = 0 ; i < vBPar2[ k ] ; ++i )
     if( v_c05f[ k ]->is_linearization_there( i ) )
      add_to_global_pool( k , i );
   }

 //!! PrintBundle();
 #if CHECK_DS & 1
  CheckBundle();
 #endif
 #if CHECK_DS & 4
  CheckAlpha();
 #endif

 // finally, release the Block - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! owned )
  f_Block->unlock( f_id );

 }  // end( BundleSolver::set_Block )

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ):
   if( value < 0 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: MaxIter must be >= 0" ) );
   MaxIter = value;
   break;
  case( intMaxSol ):
   if( value < 1 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: MaxSol must be >= 1" ) );
   MaxSol = value;
   break;
  case( intEverykIt ):
   EverykIt = value;
   break;
  case( intLogVerb ): LogVerb = value; break;
  case( intBPar1 ):
   if( value < 1 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: BPar1 must be >= 1" ) );
   BPar1 = value;
   break;
  case( intBPar2 ):
   if( value < 2 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: BPar2 must be >= 2" ) );
   if( BPar2 == Index( value ) )
    break;
   if( f_Block )
    throw( std::invalid_argument(
               "BundleSolver::set_par: "
               "changing BPar2 not supported yet" ) );
   BPar2 = value;
   break;
  case( intBPar3 ):
   if( Index( value ) < BPar4 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: BPar3 must be >= BPar4" ) );
   BPar3 = value;
   break;
  case( intBPar4 ):
   if( value < 1 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: BPar4 must be >= 1" ) );
   BPar4 = value;
   if( BPar4 > BPar3 )
    BPar3 = BPar4;
   break;
  case( intBPar6 ): BPar6 = value; break;
  case( intBPar7 ): BPar7 = value; break;
  case( intMnSSC ): MnSSC = value; break;
  case( intMnNSC ): MnNSC = value; break;
  case( inttSPar1 ): tSPar1 = value; break;
  case( intMaxNrEvls ): MaxNrEvls = value; break;
  case( intDoEasy ): DoEasy = value; break;
  case( intWZNorm ):
   if( WZNorm != char( value ) ) {
    WZNorm = char( value );
    if( ! ( WZNorm & ~3 ) )  // the easy case, constant
     NrmZFctr = 1;           // factor is known
    else
     NrmZFctr = INFshift;    // factor to be computed
    }
   break;
  case( intFrcLstSS ): FrcLstSS = value; break;
  case( intTrgtMng ):  TrgtMng = Index( value ); break;
  case( intMPStbl ):
   MPStbl = static_cast< MasterProblemBlock::stabilization_type >( value );
   break;
  case( intMPPrimal ): IsMPPrimal = bool( value ); break;
  case( intRstAlg ): RstAlgPrm = value; break;
  case( intMPV2Form ):
   if( value < 0 || value > 1 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: MPV2Form must be either 0 or 1" ) );
   MPV2Form = value;
   break;
  case( intMPHScaling ):
   if( value < 0 || value > 3 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: MPHScaling must be between 0 and 3" ) );
   MPHScaling = value;
   break;
  case( intMaxLevelNR ): MaxLevelNR = value; break;
  default: CDASolver::set_par( par , value );
  }
 }  // end( BundleSolver::set_par( int ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ):
   if( value <= 0 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: "
               "dblMaxTime must be > 0" ) );
   MaxTime = value;
   break;
  case( dblRelAcc ):
   if( ( value <= 0 ) || ( value >= INFshift ) )
    throw( std::invalid_argument(
               "BundleSolver::set_par: "
               "RelAcc must be > 0 and finite" ) );
   RelAcc = value;
   break;
  case( dblAbsAcc ):
   if( value <= 0 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: AbsAcc must be > 0" ) );
   AbsAcc = value;
   break;
  case( dblEveryTTm ):
   if( value < 0 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: EveryTTm must be >= 0" ) );
   EveryTTm = value;
   break;
  case( dblNZEps ):
   if( value < 0 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: NZEps must be >= 0" ) );
   NZEps = value;
   break;
  case( dbltStar ):
   tStar = value;
   break;
  case( dblMinNrEvls ):
   MinNrEvls = std::max( double( -1 ) , value );
   break;
  case( dblBPar5 ):
   BPar5 = value;
   break;
  case( dblm1 ):
   if( std::abs( value ) >= 1 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: "
               "| m1 | must be in (0, 1)" ) );
   m1 = value;
   break;
  case( dblm2 ):
   if( ( value < std::abs( m1 ) ) || ( value >= 1 ) )
    throw( std::invalid_argument(
               "BundleSolver::set_par: "
               "m2 must be in [ | m1 |, 1)" ) );
   m2 = value;
   break;
  case( dblm3 ):
  if( ( value <= 0 ) || ( value >= 1 ) )
    throw( std::invalid_argument(
               "BundleSolver::set_par: m3 must be in (0, 1)" ) );
   m3 = value;
   break;
  case( dblmxIncr ):
  if( value <= 1 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: mxIncr must be > 1" ) );
   mxIncr = value;
   break;
  case( dblmnIncr ):
   if( value <= 1 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: mnIncr must be > 1" ) );
   mnIncr = std::min( value , mxIncr );
   break;
  case( dblmxDecr ):
   if( ( value <= 0 ) || ( value > 1 ) )
    throw( std::invalid_argument(
               "BundleSolver::set_par: "
               "mxDecr must be in (0, 1)" ) );
   mxDecr = value;
   break;
  case( dblmnDecr ):
   if( ( value <= 0 ) || ( value > 1 ) )
    throw( std::invalid_argument(
               "BundleSolver::set_par: "
               "mnDecr must be in (0, 1)" ) );
   mnDecr = std::max( value , mxDecr );
   break;
  case( dbltMaior ):
   if( value <= 0 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: tMaior must be > 0" ) );
   tMaior = value;
   break;
  case( dbltMinor ):
   if( ( value <= 0 ) || ( value > tMaior ) )
    throw( std::invalid_argument(
               "BundleSolver::set_par: "
               "tMinor must be in (0, tMaior]" ) );
   tMinor = value;
   break;
  case( dbltInit ):
   if( ( value < tMinor ) || ( value > tMaior ) )
    throw( std::invalid_argument(
               "BundleSolver::set_par: "
               "tInit must be in [tMinor, tMaior]" ) );
   tInit = value;
   break;
  case( dbltSPar2 ):
   if( value <= 0 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: tSPar2 must be > 0" ) );
   tSPar2 = value;
   break;
  case( dbltSPar3 ):
   tSPar3 = std::abs( value ) > 1 ? value : 0;
   break;
  case( dblLStabM ):
   if( ( value <= 0 ) || ( value >= 1 ) )
    throw( std::invalid_argument(
               "BundleSolver::set_par: LStabM must be in (0, 1)" ) );
   LStabM = value;
   break;
  case( dblLStabDlt ):
   if( value <= 0 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: LStabDlt must be > 0" ) );
   LStabDlt = value;
   break;
  case( dblLStabIncr ):
   if( value <= 1 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: LStabIncr must be > 1" ) );
   LStabIncr = value;
   break;
  case( dblLStabSmall ):
   if( value <= 0 )
    throw( std::invalid_argument(
               "BundleSolver::set_par: LStabSmall must be > 0" ) );
   LStabSmall = value;
   break;
  default:
   CDASolver::set_par( par , value );
  }
 }  // end( BundleSolver::set_par( double ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par , std::string && value )
{
 switch( par ) {
  case( strEasyCfg ):
   EasyCfg = std::move( value );
   break;
  case( strHardCfg ):
   HardCfg = std::move( value );
   break;
  case( strMPBSolverCfg ):
   MPBSolverCfg = std::move( value );
   break;
  default:
   CDASolver::set_par( par , std::move( value ) );
  }
 }  // end( BundleSolver::set_par( std::string && ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par , std::vector< int > && value )
{
 if( par == vintNoEasy )
  NoEasy = std::move( value );
 else
  CDASolver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par ,
                            std::vector< std::string > && value )
{
 switch( par ) {
  case( vstrCmpCfg ): CmpCfg = std::move( value ); break;
  case( vstr_C05_SPAR_Names ):    v_C05_SPAR_Names = std::move( value );
                                  break;
  case( vstr_C05_SPAR_Vals ):     v_C05_SPAR_Vals = std::move( value );
                                  break;
  case( vstr_C05_EI_SPAR_Names ): v_C05_EI_SPAR_Names = std::move( value );
                                  break;
  case( vstr_C05_EI_SPAR_Vals ):  v_C05_EI_SPAR_Vals = std::move( value );
                                  break;
  default: CDASolver::set_par( par , std::move( value ) );
  }
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::set_log( std::ostream * log_stream )
{
 f_log = log_stream;
 if( MasterPB )
  MasterPB->forward_log( f_log );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::get_var_solution( Configuration *solc )
{
 // first take the values of the ColVariable in the Block

 if( ( MaxSol > 1 ) && ( UpFiBest < UpRifFi.back() ) ) {
  for( Index i = 0 ; i < NumVar ; i++ )
   LamVcblr[ i ]->set_value( LmbdBst[ i ] );
  }
 else
  for( Index i = 0 ; i < NumVar ; i++ )
   LamVcblr[ i ]->set_value( Lambda[ i ] );

 // now, if so instructed, also take the dual optimal solutions of the
 // (chosen) easy components

 // case where a subset of (easy) components is specified
 if( auto c = dynamic_cast< SimpleConfiguration< std::vector<
                                   std::pair< int , int > > > * >( solc ) ) {
  for( auto p : c->f_value ) {
   if( ( p.first < 0 ) || ( p.first >= int( NrFi ) ) )
    throw( std::invalid_argument(
               "BundleSolver::get_var_solution: invalid index " +
               std::to_string( p.first ) ) );
   if( ! IsEasy[ p.first ] )
    throw( std::invalid_argument(
               "BundleSolver::get_var_solution: component " +
               std::to_string( p.first ) + " is not easy" ) );
   if( p.second & 1 )
    get_var_solution_easy_pi( p.first );

   if( p.second & 2 )
    get_var_solution_easy_rc( p.first );
   }

  return;
  }

 // case "all the easy components"
 if( auto c = dynamic_cast< SimpleConfiguration< int > * >( solc ) ) {
  auto h = c->f_value;
  if( ! ( h & 3 ) )  // it'd be funny, but ...
   return;

  for( Index i = 0 ; i < NrFi ; ++i )
   if( IsEasy[ i ] ) {
    if( h & 1 )
     get_var_solution_easy_pi( i );

    if( h & 2 )
     get_var_solution_easy_rc( i );
    }

  return;
  }
 }  // end( BundleSolver::get_var_solution )

/*--------------------------------------------------------------------------*/

void BundleSolver::get_var_solution_easy_pi( Index k )
{
 if( ! ( ( DoEasy & 12 ) == 12 ) )
  throw( std::logic_error(
       "BundleSolver::get_var_solution_easy_pi: "
       "intDoEasy & 12 == 12 required to get easy components pi" ) );

 // The k-th easy sub-Block registered into MasterPB is the very same
 // inner Block exposed by v_c05f[ k ]->build_easy_master_block(); when
 // the master Solver computes(), the dual values of its RowConstraint
 // are written back in place, so there is nothing to copy here.

 if( ! MasterPB || ! MasterPB->get_easy_component( int( k ) ) )
  throw( std::logic_error(
       "BundleSolver::get_var_solution_easy_pi: "
       "easy component " + std::to_string( k ) + " not registered "
       "in MasterProblemBlock" ) );

 }  // end( BundleSolver::get_var_solution_easy_pi() )

/*--------------------------------------------------------------------------*/

void BundleSolver::get_var_solution_easy_rc( Index k )
{
 if( ! ( DoEasy & 8 ) )
  throw( std::logic_error(
       "BundleSolver::get_var_solution_easy_rc: "
       "intDoEasy & 8 required to get easy components rc" ) );

 // Same observation as in get_var_solution_easy_pi: the reduced costs of
 // the easy sub-Block ColVariable are produced in place by the master
 // Solver, no explicit propagation is needed.

 if( ! MasterPB || ! MasterPB->get_easy_component( int( k ) ) )
  throw( std::logic_error(
       "BundleSolver::get_var_solution_easy_rc: "
       "easy component " + std::to_string( k ) + " not registered "
       "in MasterProblemBlock" ) );

 }  // end( BundleSolver::get_var_solution_easy_rc() )

/*--------------------------------------------------------------------------*/

void BundleSolver::get_dual_solution( Configuration * solc )
{
 if( auto c = dynamic_cast< SimpleConfiguration< std::vector< int > > * >(
                                                                  solc ) ) {
  for( auto k : c->f_value ) {
   if( ( k < 0 ) || ( k >= int( NrFi ) ) )
    throw( std::invalid_argument(
               "BundleSolver::get_dual_solution: invalid index " +
               std::to_string( k ) ) );

   if( NrEasy && ( ! IsEasy[ k ] ) )
    get_dual_solution_hard( k );
   else
    get_dual_solution_easy( k );
   }

  return;
  }

 if( NrEasy ) {
  for( Index k = 0 ; k < NrFi ; ++k )
   if( IsEasy[ k ] )
    get_dual_solution_easy( k );
   else
    get_dual_solution_hard( k );
  }
 else
  for( Index k = 0 ; k < NrFi ; ++k )
   get_dual_solution_hard( k );

 }  // end( BundleSolver::get_dual_solution() )

/*--------------------------------------------------------------------------*/

void BundleSolver::get_dual_solution_easy( Index k )
{
 // The optimal primal u^k of the k-th easy component is the primal
 // solution of the easy sub-Block registered into MasterPB; since that
 // sub-Block is the very same inner Block exposed by
 // v_c05f[ k ]->build_easy_master_block(), the ColVariable values are
 // already in place after the master Solver has solve()d.

 if( ! MasterPB || ! MasterPB->get_easy_component( int( k ) ) )
  throw( std::logic_error(
       "BundleSolver::get_dual_solution_easy: "
       "easy component " + std::to_string( k ) + " not registered "
       "in MasterProblemBlock" ) );

 }  // end( BundleSolver::get_dual_solution_easy() )

/*--------------------------------------------------------------------------*/

void BundleSolver::get_dual_solution_hard( Index k )
{
 // construct the important linearization for the non-easy component (unless
 // it is already there, and signal to the C05Functions which one it is

 C05Function::LinearCombination lc;

 if( Zvalid[ k ] ) {
  // the optimal aggregated linearization for component k is in the
  // bundle (and, therefore, global pool) already: the optimal
  // coefficients are very simple, it's just that one
  lc.resize( 1 );
  lc[ 0 ].first = ItemVcblr[ whisZ[ k ] ].second;
  lc[ 0 ].second = 1;
  }
 else if( MasterPB ) {
  // retrieve optimal multipliers from MasterPB and copy them into the
  // LinearCombination
  const auto thetas = MasterPB->get_thetas( int( k ) );
  lc.reserve( thetas.size() );
  for( int slot = 0 ; slot < int( thetas.size() ) ; ++slot ) {
   if( thetas[ slot ] == 0 )
    continue;
   // find the global name corresponding to ( k , slot )
   for( Index name = 0 ; name < Index( ItemVcblr.size() ) ; ++name )
    if( ItemVcblr[ name ].first == k &&
        ItemVcblr[ name ].second == Index( slot ) ) {
     lc.emplace_back( ItemVcblr[ name ].second , thetas[ slot ] );
     break;
     }
   }
  }

 v_c05f[ k ]->set_important_linearization( std::move( lc ) );

 }  // end( BundleSolver::get_dual_solution_hard() )

/*--------------------------------------------------------------------------*/

int BundleSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):   return( MaxIter );
  case( intMaxSol ):    return( MaxSol );
  case( intEverykIt ):  return( EverykIt );
  case( intLogVerb ):   return( LogVerb );
  case( intBPar1 ):     return( BPar1 );
  case( intBPar2 ):     return( BPar2 );
  case( intBPar3 ):     return( BPar3 );
  case( intBPar4 ):     return( BPar4 );
  case( intBPar6 ):     return( BPar6 );
  case( intBPar7 ):     return( BPar7 );
  case( intMnSSC ):     return( MnSSC );
  case( intMnNSC ):     return( MnNSC );
  case( inttSPar1 ):    return( tSPar1 );
  case( intMaxNrEvls ): return( MaxNrEvls );
  case( intDoEasy ):    return( DoEasy );
  case( intWZNorm ):    return( WZNorm );
  case( intFrcLstSS ):  return( FrcLstSS );
  case( intTrgtMng ):   return( TrgtMng );
  case( intMPStbl ):    return( MPStbl );
  case( intMPPrimal ):  return( int( IsMPPrimal ) );
  case( intRstAlg ):    return( RstAlgPrm  );
  case( intMPV2Form ):  return( MPV2Form );
  case( intMPHScaling ): return( MPHScaling );
  case( intMaxLevelNR ): return( MaxLevelNR );
  default:              return( CDASolver::get_int_par( par ) );
  }
 }  // end( BundleSolver::get_int_par )

/*--------------------------------------------------------------------------*/

double BundleSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ):   return( MaxTime );
  case( dblRelAcc ):    return( RelAcc );
  case( dblAbsAcc ):    return( AbsAcc );
  case( dblEveryTTm ):  return( EveryTTm );
  case( dblNZEps ):     return( NZEps );
  case( dbltStar ):     return( tStar );
  case( dblMinNrEvls ): return( MinNrEvls );
  case( dblBPar5 ):     return( BPar5 );
  case( dblm1 ):        return( m1 );
  case( dblm2 ):        return( m2 );
  case( dblm3 ):        return( m3 );
  case( dblmxIncr ):    return( mxIncr );
  case( dblmnIncr ):    return( mnIncr );
  case( dblmxDecr ):    return( mxDecr );
  case( dblmnDecr ):    return( mnDecr );
  case( dbltMaior ):    return( tMaior );
  case( dbltMinor ):    return( tMinor );
  case( dbltInit ):     return( tInit );
  case( dbltSPar2 ):    return( tSPar2 );
  case( dbltSPar3 ):    return( tSPar3 );
  case( dblLStabM ):     return( LStabM );
  case( dblLStabDlt ):   return( LStabDlt );
  case( dblLStabIncr ):  return( LStabIncr );
  case( dblLStabSmall ): return( LStabSmall );
  default:               return( CDASolver::get_dbl_par( par ) );
  }
 }  // end( BundleSolver::get_dbl_par )

/*--------------------------------------------------------------------------*/

const std::string & BundleSolver::get_str_par( idx_type par ) const
{
 switch( par ) {
  case( strEasyCfg ):       return( EasyCfg );
  case( strHardCfg ):       return( HardCfg );
  case( strMPBSolverCfg ):  return( MPBSolverCfg );
  default:                  return( CDASolver::get_str_par( par ) );
  }
 }  // end( BundleSolver::get_str_par )

/*--------------------------------------------------------------------------*/

const std::vector< int > & BundleSolver::get_vint_par( idx_type par ) const
{
 if( par == vintNoEasy )
  return( NoEasy );

 return( CDASolver::get_vint_par( par ) );

 }  // end( BundleSolver::get_vint_par )

/*--------------------------------------------------------------------------*/

const std::vector< std::string > & BundleSolver::get_vstr_par( idx_type par )
 const
{
 switch( par ) {
  case( vstrCmpCfg ):             return( CmpCfg );
  case( vstr_C05_SPAR_Names ):    return( v_C05_SPAR_Names );
  case( vstr_C05_SPAR_Vals ):     return( v_C05_SPAR_Vals );
  case( vstr_C05_EI_SPAR_Names ): return( v_C05_EI_SPAR_Names );
  case( vstr_C05_EI_SPAR_Vals ):  return( v_C05_EI_SPAR_Vals );
  }

 return( CDASolver::get_vstr_par( par ) );

 }  // end( BundleSolver::get_vstr_par )

/*--------------------------------------------------------------------------*/
/*----------- METHODS FOR HANDLING THE State OF THE BundleSolver -----------*/
/*--------------------------------------------------------------------------*/

State * BundleSolver::get_State( void ) const {
  return( new BundleSolverState( this ) );
  }  // end( BundleSolver::get_State )

/*--------------------------------------------------------------------------*/

void BundleSolver::put_State( const State & state )
{
 // if state is not a BundleSolverState &, exception will be thrown
 const auto & s = dynamic_cast< const BundleSolverState & >( state );

 guts_of_put_State( s );

 Lambda = s.Lambda;

 for( Index i = 0 ; i < NrFi ; ++i )
  if( s.v_comp_State[ i ] )
   v_c05f[ i ]->put_State( *(s.v_comp_State[ i ]) );

 }  // end( BundleSolver::put_State( const & ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::put_State( State && state )
{
 // if state is not a BundleSolverState &&, exception will be thrown
 auto && s = dynamic_cast< BundleSolverState && >( state );

 guts_of_put_State( s );

 Lambda = std::move( s.Lambda );

 for( Index i = 0 ; i < NrFi ; ++i )
  if( s.v_comp_State[ i ] ) {
   v_c05f[ i ]->put_State( std::move( *(s.v_comp_State[ i ]) ) );
   delete s.v_comp_State[ i ];
   }

 s.v_comp_State.clear();

 }  // end( BundleSolver::put_State( && ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::serialize_State( netCDF::NcGroup & group ,
                                    const std::string & sub_group_name ) const
{
 if( ! sub_group_name.empty() ) {
  auto gr = group.addGroup( sub_group_name );
  serialize_State( gr );
  return;
  }

 // do it "by hand" since there is no BundleSolverState available to call
 // State::serialize() from
 group.putAtt( "type", "BundleSolverState" );

 auto nv = group.addDim( "BundleSolver_NumVar" , NumVar );

 ( group.addVar( "BundleSolver_Lambda" , netCDF::NcDouble() , nv ) ).putVar(
                                       { 0 } , {  NumVar } , Lambda.data() );

 ( group.addVar( "BundleSolver_t" , netCDF::NcDouble() ) ).putVar( &t );

 auto nfi = group.addDim( "BundleSolver_NrFi" , NrFi + 1 );

 if( UpFiLmbdef ) {
  group.addDim( "BundleSolver_UpFiLmbdef" , UpFiLmbdef );
  ( group.addVar( "BundleSolver_UpFiLmb" , netCDF::NcDouble() , nfi )
    ).putVar( { 0 } , { NrFi + 1 } , UpFiLmb.data() );
  }

 if( LwFiLmbdef ) {
  group.addDim( "BundleSolver_LwFiLmbdef" , LwFiLmbdef );
  ( group.addVar( "BundleSolver_LwFiLmb" , netCDF::NcDouble() , nfi )
    ).putVar( { 0 } , { NrFi + 1 } , LwFiLmb.data() );
  }

 if( Fi0Lmb != 0 )
  ( group.addVar( "BundleSolver_Fi0Lmb" , netCDF::NcDouble() )
    ).putVar( & Fi0Lmb );

 if( f_global_LB > - INFshift )
  ( group.addVar( "BundleSolver_global_LB" , netCDF::NcDouble() )
    ).putVar( & f_global_LB );

 for( Index i = 0 ; i < NrFi ; ++i ) {
  if( NrEasy && IsEasy[ i ] )
   continue;

  v_c05f[ i ]->serialize_State( group ,
                                "Component_State_" + std::to_string( i ) );
  }
 }  // end( BundleSolver::serialize_State )

/*--------------------------------------------------------------------------*/
/*----------------------- OTHER PROTECTED METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::guts_of_put_State( const BundleSolverState & state )
{
 if( Result == kStillRunning )
  throw( std::logic_error(
        "BundleSolver::put_State() called within BundleSolver::compute()" ) );

 if( ( NrFi != state.NrFi ) || ( NumVar != state.NumVar ) )
  throw( std::invalid_argument(
                          "BundleSolver::put_State(): inconsistent State" ) );

 if( t != state.t ) {
  t = state.t;
  tHasChgd = true;
  }

 if( state.UpFiLmbdef ) {
  UpFiLmbdef = state.UpFiLmbdef;
  UpFiLmb = state.UpFiLmb;
  }
 else {
  UpFiLmbdef = 0;
  std::fill( UpFiLmb.begin() , UpFiLmb.end() ,  INFshift );
  }

 if( state.LwFiLmbdef ) {
  LwFiLmbdef = state.LwFiLmbdef;
  LwFiLmb = state.LwFiLmb;
  }
 else {
  LwFiLmbdef = 0;
  std::fill( LwFiLmb.begin() , LwFiLmb.end() ,  -INFshift );
  }

 // if Lambda has changed, the master problem need be informed
 if( Lambda != state.Lambda ) {
  Vec_VarValue foo( NrFi + 1 , 0 );

  if( UpFiLmbdef == NrFi + 1 ) {
   // the function value is known, so it has to be passed to the master:
   // see the comments inside GotoLambda1() for the rationale of the
   // curios definition
   foo.front() = UpFiLmb.back() - UpRifFi.back();
   std::transform( UpFiLmb.begin() , --(UpFiLmb.end()) , UpRifFi.begin() ,
                   ++(foo.begin()) , std::minus< double >() );
   UpRifFi = UpFiLmb;
   RifeqFi = true;
   }
  // else no change in the (unknown) f-values, so all-0 is OK

  if( MasterPB ) {
   std::vector< double > F_hard;
   F_hard.reserve( NrFi );
   for( Index k = 0 ; k < NrFi ; ++k ) {
    if( NrEasy && IsEasy[ k ] )
     continue;
    F_hard.push_back( UpRifFi[ k ] );
    }
   MasterPB->set_reference( state.Lambda , F_hard );
   }
  }

 if( FrcLstSS & 2 ) {
  // something may have changed in the Block in the meantime, force the
  // recomputation of all components before trusting the current state
  UpFiLmbdef = 0;
  std::fill( UpFiLmb.begin() , UpFiLmb.end() ,  INFshift );
  RifeqFi = false;
  }

 Fi0Lmb = state.Fi0Lmb;
 f_global_LB = state.global_LB;

 }  // end( BundleSolver::guts_of_put_State )

/*--------------------------------------------------------------------------*/

BundleSolver::VarValue BundleSolver::reliable_level_LB( void ) const
{
 VarValue lb = -INFshift;
 if( f_global_LB > lb )
  lb = f_global_LB;
 if( TrueLB && ( LowerBound.back() > lb ) )
  lb = LowerBound.back();
 return( lb );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::reset_level_stabilization( void )
{
 f_level_Delta = 0;
 f_level_value = INFshift;
 f_level_LB = -INFshift;
 f_level_reliable_LB = false;
 f_level_initialized = false;
 if( MasterPB && UsesLevelStabilization() )
  MasterPB->set_f_lev( INFshift );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::install_level_stabilization( void )
{
 if( ! ( MasterPB && UsesLevelStabilization() ) )
  return;
 if( ( UpFiLmb.back() >= INFshift ) || ( f_level_value >= INFshift ) ) {
  MasterPB->set_f_lev( INFshift );
  return;
  }

 auto lev = f_level_value;
 if( ! ( UsesPrimalMaster() && MPV2Form ) ) {
  VarValue rf = 0;
  if( ( ! UsesPrimalMaster() ) && MPV2Form ) {
   // In the dual iterate frame the explicit x_bar . z objective term already
   // translates the linear component (and every exact easy component). The
   // PFB constants still carry F_k(x_bar), so only hard component references
   // must be removed from the absolute level.
   for( Index k = 0 ; k < NrFi ; ++k )
    if( ( ! NrEasy ) || ( ! IsEasy[ k ] ) )
     rf += UpRifFi[ k ];
   }
  else {
   rf = UpRifFi.back();
   if( NrEasy )
    for( Index k = 0 ; k < NrFi ; ++k )
     if( IsEasy[ k ] )
      rf -= UpRifFi[ k ];
   }
  lev -= rf;
  }
 MasterPB->set_f_lev( lev );
 }

/*--------------------------------------------------------------------------*/

bool BundleSolver::refresh_level_after_master( bool force )
{
 const auto old_level = f_level_value;

 if( ! UsesLevelStabilization() )
  return( false );
 if( UpFiLmb.back() >= INFshift ) {
  if( MasterPB )
   MasterPB->set_f_lev( INFshift );
  return( false );
  }

 const auto lb = reliable_level_LB();
 if( lb > -INFshift ) {
  if( force || ( ! f_level_initialized ) || ( ! f_level_reliable_LB ) ||
      ( lb > f_level_LB ) ) {
   const auto gap = UpFiLmb.back() - lb;
   f_level_Delta = gap > 0 ? ( 1.0 - LStabM ) * gap : 0.0;
   f_level_LB = lb;
   }
  f_level_reliable_LB = true;
  }
 else if( ( ! f_level_initialized ) || f_level_reliable_LB ||
          ( f_level_Delta <= 0 ) ||
          ( UpFiLmb.back() - f_level_Delta >= UpFiLmb.back() ) ) {
  f_level_Delta = LStabDlt * std::max( std::abs( UpFiLmb.back() ) , 1.0 );
  f_level_reliable_LB = false;
  }

 f_level_value = UpFiLmb.back() - f_level_Delta;
 f_level_initialized = true;
 install_level_stabilization();

 if( old_level == f_level_value )
  return( false );

 if( ( ! std::isfinite( old_level ) ) ||
     ( ! std::isfinite( f_level_value ) ) )
  return( true );

 const auto scale = std::max( { std::abs( old_level ) ,
                                std::abs( f_level_value ) , VarValue( 1 ) } );
 const auto tol = 16 * std::numeric_limits< VarValue >::epsilon() * scale;
 return( std::abs( f_level_value - old_level ) > tol );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::update_level_after_step( bool serious_step ,
                                           bool gated_update )
{
 if( ! ( UsesLevelStabilization() && f_level_initialized ) )
  return;

 const auto lb = reliable_level_LB();
 if( lb > -INFshift ) {
  const auto gap = UpFiLmb.back() - lb;
  const auto from_lb = gap > 0 ? ( 1.0 - LStabM ) * gap : 0.0;
  if( ( ! f_level_reliable_LB ) || f_level_Delta <= 0 )
   f_level_Delta = from_lb;
  else if( serious_step ) {
   if( gated_update )
    f_level_Delta = std::min( f_level_Delta , from_lb );
   }

  f_level_LB = lb;
  f_level_reliable_LB = true;
  }

 // Pure-level safeguard:
 // If this is a null step but the master problem has not actually changed,
 // do not apply the ordinary NS rule Delta <- m_l Delta.
 // Instead, override it and enlarge Delta.
 //
 // Rationale: shrinking Delta raises the level and relaxes the target. This
 // is meaningful after an informative NS, but not after a non-informative
 // one, where no useful model information was generated.
 // The safeguard is only needed when the level multiplier is small: this is
 // the ill-conditioned regime where pure-level aggregation divides by eta.
 // With eta already sizeable, increasing Delta can over-tighten the level and
 // cause an increase/shrink cycle; let the ordinary NS rule relax it instead.
 const auto level_eta = MasterPB ? MasterPB->get_level_multiplier() : 0.0;
 const bool small_level_eta = level_eta <= 1.0;
 if( ( ! serious_step ) && UsesPureLevelStabilization() && 
        ( ! MPchgs ) && small_level_eta ) {
  f_level_Delta *= LStabIncr;

  if( f_level_reliable_LB ) {
   const auto rlb = reliable_level_LB();
   if( rlb > -INFshift ) {
    const auto gap = UpFiLmb.back() - rlb;
    if( gap > 0 )
     f_level_Delta = std::min( f_level_Delta , gap );
    }
  }

  // Count this as a special level-noise / non-informative NS event.
  ++LevelNRCntr;

  BLOG( 1 , " ~ level NR: MP unchanged, Delta increased to "
          << shrt << f_level_Delta << std::endl );

  if( f_level_Delta <= 0 &&
      ( ! f_level_reliable_LB || reliable_level_LB() <= -INFshift ) )
   f_level_Delta = LStabDlt * std::max( std::abs( UpFiLmb.back() ) , 1.0 );

  f_level_value = UpFiLmb.back() - f_level_Delta;
  install_level_stabilization();
  return;
  }

 if( serious_step && gated_update && UsesPureLevelStabilization() &&
     ( lb <= -INFshift ) && ( ! f_level_reliable_LB ) &&
     ( UpFiLmb.back() < INFshift ) && ( UpRifFi.back() < INFshift ) &&
     ( vStar.back() < INFshift ) ) {
  const auto level_model_value = UpRifFi.back() + vStar.back();
  const auto level_model_gap =
   std::max( UpFiLmb.back() - level_model_value , VarValue( 0 ) );
  const auto level_model_ratio =
   level_model_gap / std::max( std::abs( level_model_value ) , VarValue( 1 ) );

  if( level_model_ratio <= LStabSmall ) {
   f_level_Delta *= LStabIncr;
   CSSCntr = 0;
   }
  }
  
 // reset the NR counter if needed
 if( serious_step || MPchgs )
  LevelNRCntr = 0;

 if( ( ! serious_step ) && gated_update )
  f_level_Delta *= LStabM;

 if( f_level_Delta <= 0 &&
     ( ! f_level_reliable_LB || reliable_level_LB() <= -INFshift ) )
  f_level_Delta = LStabDlt * std::max( std::abs( UpFiLmb.back() ) , 1.0 );

 f_level_value = UpFiLmb.back() - f_level_Delta;
 install_level_stabilization();
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::record_level_lower_bound( VarValue lb )
{
 if( lb <= -INFshift )
  return;
 if( f_global_LB < lb )
  f_global_LB = lb;
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FormD( void )
{
 // initialize the Master Problem Solver- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // change/set t as required- - - - - - - - - - - - - - - - - - - - - - - - -

 // Special treatment of the "empty Master Problem" case: there are no
 // subgradients, so this is only a feasibility problem, which should give
 // a feasible point as close as possible to the starting one. This is
 // "free" with some stabilizing terms (e.g. the quadratic one), but not
 // necessarily so with others (e.g. the trust region). In order to "stay
 // as close as possible", t is temporarily decreased to its minimum value.
 // As soon as there is something in the bundle, the current value of t is
 // restored (Prevt is used to hold it).

 // bundle-is-empty check: "no subgradient cut has been added yet to any
 // hard component". Reads is_bundle_empty() from MasterPB.
 const bool empty_bundle = MasterPB ? MasterPB->is_bundle_empty() : true;
 if( empty_bundle ) {
  if( ( t > tMinor ) && ( Prevt == INFshift ) ) {
   Prevt = t;
   t = tMinor;
   tHasChgd = true;
   }
  }
 else
  if( Prevt < INFshift ) {
   if( t != Prevt ) {
    t = Prevt;
    tHasChgd = true;
    }
   Prevt = INFshift;
   }

 if( tHasChgd ) {
  if( MasterPB )
   MasterPB->set_t( t );
  tHasChgd = false;
  }

 // collect and set individual lower bounds, if any - - - - - - - - - - - - -
 // if the MPSolver accepts them, collect and if necessary set the individual
 // lower bounds. note that if all of them are finite and the 0-th component
 // is not there their sum would give an alternative valid global lower bound.
 // however, the same information is already encoded in the individual bounds,
 // hence it's of no use. The MasterProblemBlock accepts individual lower
 // bounds on every component, so they are forwarded unconditionally.

 for( Index k = 0 ; k < NrFi ; ++k ) {
  if( NrEasy && IsEasy[ k ] )  // skip easy components
   continue;

  // get the lower bound out of the C05Function
  auto LwrBndk = f_convex ?   v_c05f[ k ]->get_global_lower_bound()
                          : - v_c05f[ k ]->get_global_upper_bound();

  // use it to update the lower estimate in Lambda
  update_LwFiLambd( k , LwrBndk );

  if( LwrBndk != LowerBound[ k ] ) {  // if it has changed
   LowerBound[ k ] = LwrBndk;         // record the new value

   // pass the raw lower bound to the master: MPB stores it as the
   // native bound of its PFB and handles any required translation
   // (gamma * LB term coefficient) internally as the reference
   // F_k( x_bar ) evolves via set_reference. The legacy pre-shift
   // " LwrBndk -= UpRifFi[ k ] " is gone, paralleling the removal of
   // the per-cut alpha promotion in GetGi
   if( MasterPB )
    MasterPB->set_LB( int( k ) , LwrBndk );
   }
  }

 // collect and set the global lower bound, if any- - - - - - - - - - - - - -
 // note: this used to be done elsewhere, in particular after each function
 //       evaluation. the rationale was that in the Lagrangian case one would
 //       do heuristics as a part of the computation, and these could produce
 //       a better lower bound that one may immediately want to check. but
 //       the truth is that one does not: the real way in which a lower bound
 //       is useful is when it enters the master problem and helps in getting
 //       the optimality conditions. thus, the right place to check and
 //       update the lower bound(s) is right before the master problem is
 //       solved. the conditional lower bound ( TrueLB == false ) is indeed
 //       useful right after the function computation to prove unboundedness,
 //       but it is available then, and anyway it typically does not change
 //       when the function is computed, unlike the "hard" one

 auto LwrBnd = f_convex ?   f_Block->get_valid_lower_bound( false )
                        : - f_Block->get_valid_upper_bound( false );

 if( LwrBnd != LowerBound.back() ) {
  if( TrueLB ) {  // if the global lower bound was a "true" one, its value
   // is "baked in" the total lower and upper estimate of the function value
   // in Lambda: ensure it is recomputed. this works both if the bound was
   // there and it is changed and if it is reset
   if( UpFiLmbdef > NrFi ) {
    --UpFiLmbdef;
    UpFiLmb.back() = INFshift;
    }
   if( LwFiLmbdef > NrFi ) {
    --LwFiLmbdef;
    LwFiLmb.back() = -INFshift;
    }
   }

  TrueLB = ( LwrBnd > -INFshift );  // if it's finite it's a true LB

  // if a true global lower bound was not there and it is set, then ensure
  // it will be properly "baked in" the total lower and upper estimate
  if( TrueLB ) {
   if( UpFiLmbdef > NrFi ) {
    --UpFiLmbdef;
    UpFiLmb.back() = INFshift;
   }
   if( LwFiLmbdef > NrFi ) {
    --LwFiLmbdef;
    LwFiLmb.back() = -INFshift;
    }
   }

  // if the total upper estimate in Lambda needs be recomputed, do it now
  // the only role of the total upper estimate in Lambda in the master
  // problem is to translate the global lower bound (if any); this will
  // be done during the final call to SetLowerBound()
  if( UpFiLmbdef == NrFi ) {
   ++UpFiLmbdef;  // all components + the sum computed
   UpFiLmb.back() = std::accumulate( UpFiLmb.begin() , --(UpFiLmb.end()) ,
                                     Fi0Lmb );
   // note that here the global lower bound (if any) is "baked in" the total
   // upper function estimate
   if( TrueLB && ( UpFiLmb.back() < LwrBnd ) )
    UpFiLmb.back() = LwrBnd;
   }

  // if the total lower estimate in Lambda needs be recomputed, do it now
  if( LwFiLmbdef == NrFi ) {
   ++LwFiLmbdef;  // all components + the sum computed
   LwFiLmb.back() = std::accumulate( LwFiLmb.begin() , --(LwFiLmb.end()) ,
                                     Fi0Lmb );
   // note that here the global lower bound (if any) is "baked in" the total
   // lower function estimate
   if( TrueLB && ( LwFiLmb.back() < LwrBnd ) )
    LwFiLmb.back() = LwrBnd;
   }

  LowerBound.back() = LwrBnd;        // in all cases, record it
  if( TrueLB ) {   // if the bound value is finite
   // Translate it using the active storage frame. The displacement frame
   // removes the linear and hard-component reference values (but not exact
   // easy components). In the dual iterate frame x_bar . z already translates
   // the linear and exact-easy terms, so only hard-component references are
   // removed here.
   VarValue rf = 0;
   if( ( ! UsesPrimalMaster() ) && MPV2Form ) {
    // As for the level row, x_bar . z already accounts for the linear and
    // exact-easy reference terms in the dual iterate frame.
    for( Index k = 0 ; k < NrFi ; ++k )
     if( ( ! NrEasy ) || ( ! IsEasy[ k ] ) )
      rf += UpRifFi[ k ];
    }
   else {
    rf = UpRifFi.back();
    if( NrEasy )
     for( Index k = 0 ; k < NrFi ; ++k )
      if( IsEasy[ k ] )
       rf -= UpRifFi[ k ];
    }

   LwrBnd -= rf;
   }

  // set the global lower bound in the master problem (translated if it
  // is finite); the bound has to be set even if it is -INF, because
  // before it was not so, hence it has to be reset
  if( MasterPB )
   MasterPB->set_global_LB( LwrBnd );

  }  // end( if( the global lower bound has changed ) )

 if( ! TrueLB )  // if no true LB, see if "conditional" one is there
  LowerBound.back() = f_convex
                      ?   f_Block->get_valid_lower_bound( true )
                      : - f_Block->get_valid_upper_bound( true );

 if( MasterPB && MasterPB->has_initial_level_objective() &&
     f_level_initialized ) {
  MasterPB->remove_initial_level_objective();
  install_level_stabilization();
  }

 const bool initial_level_probe =
  MasterPB && UsesLevelStabilization() &&
  MasterPB->has_initial_level_objective() &&
  ( ParIter > 0 ) &&
  RifeqFi &&
  ( UpFiLmb.back() < INFshift );

 if( initial_level_probe )
  MasterPB->set_f_lev( INFshift );
 else
  refresh_level_after_master();

 // fictitious-LB sync for empty components - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // A hard component whose bundle is empty models its F_k as the max over
 // an empty set of affine pieces, i.e. the improper constant -infinity:
 // the dual master is then infeasible (the per-component simplex row
 // sum_i theta^k_i + gamma^k = lambda collapses to gamma^k = lambda with
 // gamma^k structurally fixed to 0). Following the historical QP-bundle /
 // OSIMPSolver device, install a *fictitious* model lower bound on every
 // such component, which makes it "disappear" (gamma^k becomes a free
 // multiplier absorbing the simplex mass) and restores master
 // feasibility. The fictitious bound value is handled inside
 // MasterProblemBlock::set_fictitious_LB: 0 when no global LB exists, or
 // strictly below the genuine global LB otherwise, so it never tightens
 // the master beyond the aggregate bound (the subtle point OSIMPSolver
 // makes via "gamma_i cost = global_LB - 1"). The subsequent Fi(.)
 // evaluation then either refills the bundle (proper function -> cut
 // added, fictitious removed next time) or returns F_k = -INF (improper
 // function -> the UpFiLmb1.back() == -INFshift path right after
 // InnerLoop propagates kUnbounded, *not* masked by the fictitious
 // bound). Components owning a genuine individual LB never need this.
 // With no easy components, the all-empty case is left to the
 // MasterProblemBlock::solve_master short-circuit. With exact easy
 // components, however, that master must be solved immediately; therefore
 // every empty hard component also needs its fictitious bound from the first
 // iteration, otherwise gamma^k = 0 forces lambda = 0 in its normalization
 // row while the easy master fixes lambda = 1.
 if( MasterPB && ( NrEasy || ( ! MasterPB->is_bundle_empty() ) ) )
  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( NrEasy && IsEasy[ k ] )
    continue;
   const bool want = ( NrItems[ k ] == 0 ) &&
                      ( LowerBound[ k ] <= - INFshift );
   // Refresh an active fictitious bound after every possible reference-point
   // change. set_reference() rebuilds the PFB bound from the physical lower
   // bound and can therefore overwrite the fictitious value even though the
   // BundleSolver-side FictLB flag is still true.
   if( want ) {
    MasterPB->set_fictitious_LB( int( k ) , want );
    FictLB[ k ] = true;
   }
   else
    if( FictLB[ k ] ) {
     // A finite physical bound collected by FormD() has already replaced
     // the fictitious one through set_LB(). In that case only clear our
     // bookkeeping flag: calling set_fictitious_LB(false) would erase the
     // genuine bound that has just been installed.
     if( LowerBound[ k ] <= - INFshift )
      MasterPB->set_fictitious_LB( int( k ) , false );
     FictLB[ k ] = false;
    }
  }

 for( ; ; )  // error-handling loop - - - - - - - - - - - - - - - - - - - - - -
 {           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // ensure the master Solver will not take too much time
  if( MaxTime < INFshift ) {
   if( MasterPB )
    MasterPB->set_max_time( MaxTime - get_elapsed_time() );
   }

  int mps = Solver::kOK;
  if( MasterPB ) {
   std::vector< double > Lbox( NumVar ) , Ubox( NumVar );
   for( Index i = 0 ; i < NumVar ; ++i ) {
    const auto bounds = effective_bounds( LamVcblr[ i ] );
    Lbox[ i ] = bounds.first;
    Ubox[ i ] = bounds.second;
    }
   MasterPB->set_box( Lbox , Ubox );

   const auto rc = MasterPB->solve_master();
   // an OK or a low-precision OK from the inner Solver counts as kOK;
   // anything else is forwarded so the surrounding error-handling kicks in
   mps = ( rc == Solver::kOK || rc == Solver::kLowPrecision )
         ? Solver::kOK : rc;
   }

  if( mps == Solver::kOK )           // everything's alright
   break;

  const bool level_empty =
   UsesPrimalMaster()
   ? ( ( mps == Solver::kInfeasible ) && ( ! get_bc_size() ) )
   : ( mps == Solver::kUnbounded );

  if( UsesLevelStabilization() && f_level_initialized &&
      ( f_level_value < INFshift ) && level_empty ) {
   record_level_lower_bound( f_level_value );
   if( ! refresh_level_after_master( true ) ) {
    BLOG( 1 , std::endl
              << "Bundle::FormD: empty level refresh made no progress" );
    Result = kError;
    return;
    }
   continue;
   }

  /* If it's not OK, three things can happen: unfeasible, unbounded, or a
   * numerical error.
   *
   * The primal Master Problem is
   *
   *    P_{B,Lambda,t}:   inf{ Fi_{B,Lambda}( d ) + D_t( d ) }
   *
   * where Fi_{B,Lambda}() is the cutting-plane model and D_t() is the
   * stabilisation term, although with easy components the most natural
   * implementation is actually in terms of the dual
   *
   *    D_{B,Lambda,t}:   inf{ Fi_{B,Lambda}*( z ) + D_t*( -z ) }
   *
   * The primal Master Problem can therefore be unfeasible only if there
   * are vertical linearizations. With Quadratic or BoxStep stabilisation
   * the primal Master Problem can never be "naturally" unbounded. Yet,
   * it can still be unbounded if the dual Master Problem is empty,
   * which can happen if there are easy components. If there are not,
   * the inner Solver returning kInfeasible/unbounded can only be a
   * numerical error in disguise. */

  if( mps == Solver::kInfeasible ) {  // the MP is (primal) empty
   if( ! get_bc_size() )              // there are no vertical linearizations
    mps = Solver::kError;             // it must be a numerical error
   else {                             // there are vertical linearizations
    Result = kInfeasible;             // the MP can really be infeasible
    return;                           // nothing else to do
    }
   }

  if( mps == Solver::kUnbounded ) {   // the MP is (primal) unbounded
   if( ! NrEasy )                     // there are no easy components
    mps = Solver::kError;             // it must be a numerical error
   else {                             // there are easy components
    Result = kUnbounded;              // the MP can really be unbounded
    return;                           // nothing else to do
    }
   }

  if( mps == Solver::kStopTime ) {  // stopped by time limit
   Result = kStopTime;
   break;
   }

  // mps == Solver::kError, i.e., there has been a numerical problem in
  // the inner Solver; it's not yet time to despair, as by eliminating
  // items it may be possible to solve the problem.

  BLOG( 2 , std::endl << "Bundle::FormD: error in MP, emergency delete" );

  Index i = InINF;

  // the last *removable* item in the bundle is eliminated; walk all
  // global names from the top, pick the first removable one with
  // positive theta ("in base, non-zero multiplier")
  if( MasterPB ) {
   for( Index j = Index( ItemVcblr.size() ) ; j-- ; )
    if( ( OOBase[ j ] >= 0 ) &&
        ( read_theta_global( j ) >= 1e-15 ) ) {
     i = j;
     break;
     }
   }

  if( i == InINF )  // there are no *removable* items in Base - - - - - - - -
  for( Index j = get_max_name() ; j-- ; )  // pick any removable item
    if( ( OOBase[ j ] >= 0 ) && ( OOBase[ j ] < Inf< SIndex >() ) ) {
     i = j;
     break;
     }

  if( i == InINF ) {  // there are no removable items at all - - - - - - - -
   BLOG( 1 , std::endl << "Bundle::FormD: unrecoverable MP failure." );
   Result = kError;
   return;
   }

  // i can be deleted, do it and try again hoping this will solve the
  // numerical issue in the master problem

  if( ( BPar7 & 3 ) == 3 ) {
   // if BPar7 tells that the removal must be "hard", actually remove the
   // linearization from the global pool, as Delete() does not do that
   inhibit_Modification( true );
   v_c05f[ ItemVcblr[ i ].first ]->delete_linearization(
                                                   ItemVcblr[ i ].second );
   inhibit_Modification( false );
   }
  Delete( i );

  }  // end ( error-handling loop )- - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // read Sigma* (aggregated linearization error)
 if( MasterPB )
  Sigma = MasterPB->get_aggregated_alpha();
 // read the total v* (predicted decrease at Lambda1)
 if( MasterPB )
  vStar.back() = MasterPB->get_FiBLambda();

 // v* is the predicted decrease in  Lambda1 w.r.t. the value in Lambda;
 // however, if any (non-easy) component does not have any subgradient
 // in the bundle this value is not well-defined (the master problem is
 // added an artificial constraint to make v[ k ] bounded) and +INF is
 // returned for that component, and therefore for the total v* (the sum)

 // now retrieve vStar[ k ] for each component: for easy ones, the master
 // problem produces the *exact* Fi-value (up = lw) at Lambda1, which we
 // store unmodified in vStar[ k ] (that is not used anyway) for it to be
 // retrieved later by FormLambda1()
 for( Index k = 0 ; k < NrFi ; ++k )
  if( MasterPB )
   vStar[ k ] = MasterPB->get_FiBLambda( int( k ) );

 if( NrEasy ) {  // there are easy components
  // adjust the contribution of the easy components to v* and Sigma*.
  // get_FiBLambda() / get_aggregated_alpha() already return v* and Sigma*
  // *translated* by the reference Fi-value for the HARD components (their
  // cuts are stored in linearization-error form, refreshed against the
  // current reference by set_reference). The EASY components, instead, have
  // their value entered into the model un-translated, so only their
  // reference Fi-value has to be corrected here (note that v* and Sigma*
  // carry opposite signs, hence the "-=" and the "+=")
  VarValue EasyRifFi = 0;
  for( Index k = 0 ; k < NrFi ; ++k )
   if( IsEasy[ k ] )
    EasyRifFi += UpRifFi[ k ];

  if( vStar.back() < INFshift )
   vStar.back() -= EasyRifFi;
  Sigma += EasyRifFi;
  }

 #if CHECK_BAD_F & 1
  // if so required, check for "very negative Sigma" and print a warning on
  // std::cerr if it is found; however, do that only if (an upper bound on
  // the correct value of) Fi is known, for otherwise the linearization
  // errors have been computed against an arbitrary reference value and
  // they easily be negative "by chance", so the same holds for Sigma,
  // without this indicating a problem; this may happen either at the very
  // first call to compute(), or at the first iteration of a subsequent
  // call to compute() where Modification have made the Fi-value invalid
  if( UpFiLmb.back() < INFshift ) {
   auto rel_error = ( Sigma / std::max( std::abs( UpRifFi.back() ) , 1.0 ) );
   if( rel_error < - CHECK_BAD_F_EPS )
    std::cerr << std::endl << "Warning[ " << SCalls << ", " << ParIter
              << " ]: Sigma = " << rel_error << std::endl;
   }
 #endif

 DSTS = read_DStart( std::abs( tStar ) );  // D_{t*}( z* )

 // Sigma* + D*_{t*}( -z* ) is the "maximum expected increase" used in
 // the stopping criterion, EpsU is that relative to Fi( Lambda )
 if( UpFiLmb.back() < INFshift )
  EpsU = ( DSTS + Sigma ) / std::max( std::abs( UpFiLmb.back() ) ,
                                      double( 1 ) );
 else
  EpsU = 1;  // ensure EpsU is initialized somehow

 Zvalid.assign( NrFi , false );    // the z[ i ] are no longer valid

 DST = read_DStart( t );  // D_t( z* )

 // Delta* = D_t( z* ) + Sigma* is <= - v*, and a weaker requirement about
 // how much the (total) function must increase for a NS to be declared;
 // however, this only holds if v* is "true", which means that all the
 // components have some diagonal linearization in their bundle; otherwise
 // v* is "fake" and so is Delta* (but anyway, this means that any finite
 // lower bound is much better than what we currently have)

 // compute || d* ||_2
 // Under quadratic stabilisation, ReadDStart( t ) == t || z* ||_2^2 / 2
 // and d = - t z*, hence || d* ||_2 == t || z* ||_2 ==
 //                                     t * sqrt( 2 * DST / t )
 // For non-quadratic stabilisations (e.g. boxstep) we materialise d*.
 NrmD = 0;
 if( MasterPB ) {
  const auto d = MasterPB->get_d_vector();
  for( const auto v : d )
   NrmD += v * v;
  }
 NrmD = sqrt( NrmD );

 // in the easy case NrmD also gives the (2-)norm of z*; otherwise we
 // compute it explicitly. The relationship d* = - t z* is NOT VALID
 // WHEN THERE ARE CONSTRAINTS if z* is to be interpreted as the
 // aggregate subgradient of the objective; it IS VALID IF z* IS TO BE
 // INTERPRETED AS THE AGGREGATE SUBGRADIENT OF THE ESSENTIAL OBJECTIVE
 // ( f + i_X ), which is exactly what is needed here.
 if( ( WZNorm & 3 ) == 2 )
  NrmZ = NrmD / t;
 else {  // otherwise compute the norm directly from the z* vector
  std::vector< double > tZ;
  if( MasterPB )
   tZ = MasterPB->get_z_vector();
  NrmZ = ::norm( tZ , WZNorm & 3 );
  }

 // if still needed, compute the scaling factor for z*
 if( NrmZFctr == INFshift )
  compute_NrmZFctr();

 // if the scaling factor could be computed one can check if z* == 0
 // has happened and declare a globally valid LB
 // note: the update of f_global_LB used to be under guard
 // if( f_global_LB < UpFiLmb.back() + vStar.back() )
 // i.e., one would always report the largest f_global_LB ever found.
 // however, declaring a global LB is slippery, as it requires to set NZEps
 // "small enough" and no-one really knows how to do that. As a consequence,
 // one may end up with the final LB being higher than the final UB, which
 // is not something any Solver should ever report. We rather take the
 // conservative stance where the final reported LB is the one of the
 // stopping iteration: since the UB is "that one + v^*" and v^* is negative,
 // this ensures that UB >= LB
 if( ( UpFiLmb.back() < INFshift ) && ( vStar.back() < INFshift ) &&
     ( NrmZFctr < INFshift ) && ( NrmZ <= NrmZFctr * NZEps ) ) {
  f_global_LB = UpFiLmb.back() + vStar.back();
  refresh_level_after_master();
 }

 if( initial_level_probe && MasterPB ) {
  if( ( UpFiLmb.back() < INFshift ) && ( vStar.back() < INFshift ) ) {
   f_level_Delta = std::max( - vStar.back() , VarValue( 0 ) );
   if( f_level_Delta <= 0 )
    f_level_Delta = LStabDlt * std::max( std::abs( UpFiLmb.back() ) ,
                                         VarValue( 1 ) );
   f_level_LB = -INFshift;
   f_level_reliable_LB = false;
   f_level_value = UpFiLmb.back() - f_level_Delta;
   f_level_initialized = true;
   }
  else
   refresh_level_after_master();
  }

 }  // end( BundleSolver::FormD )

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdtCntrs( void )
{
 // increase all the OOBase[] counters but those == +/-Inf< SIndex >() - - - - -
 // items whose OOBase[] becomes 0 (e.g. the newly entered items, which have
 // OOBase[] == -1) are set to +1, in such a way that only the items in the
 // optimal base have OOBase[] == 0; note that the converse is not true, as
 // items in the optimal base may have OOBase[] < 0 instead

 for( auto OOit = OOBase.begin() ;
      OOit != OOBase.begin() + get_max_name() ; ++OOit )
  if( ( *OOit < Inf< SIndex >() ) && ( *OOit > -Inf< SIndex >() ) ) {
   ++(*OOit);
   if( ! *OOit )
    ++(*OOit);
   }

 // set to 0 the OOBase[] counter for items in base (if not < 0)
 // checking if the multiplier is strictly positive should be redundant
 // if one is trusting the inner Solver
 constexpr double eps = 1e-12;

 if( MasterPB ) {
  // walk every global name; OOBase[i] > 0 means "out of base, removable";
  // a positive theta means "in base", so reset OOBase[i] to 0 (in base)
  for( Index i = 0 ; i < Index( ItemVcblr.size() ) ; ++i )
   if( ( OOBase[ i ] > 0 ) && ( read_theta_global( i ) >= eps ) )
    OOBase[ i ] = 0;
  }

 /*!!
 // note that there is a case in which a component wFi has Z[ wFi ] "for free"
 // in the bundle: this is when wFi only has *one* subgradient in base (or, in
 // practice, a subgradient with multiplier very close to one). This is
 // checked here (it is basically for free), and in case whisZ[] is properly
 // set so as to avoid pointless aggregations and OOBase[] is set to -1,
 // because under no circumstances such a subgradient can ever be removed
 // from the bundle
 //
 // it now seems to me that this is stupid, since if there is only one
 // subgradient in base no aggregation is ever performed; the only issue
 // is if the base is all (but one) taken by constraints, but then even
 // aggregating does not help

 if( MBse ) {
  for( Index i ; ( i = *(MBse++) ) < InINF ; Mlt++ )
   if( *Mlt >= 1e-15 ) {
    if( ( *Mlt >= 1 - RMPAccSol ) && is_subgradient_global( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ wcomponent_global( i ) - 1 ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }
  }
 else
  for( Index i = 0 ; i < MBDim ; i++ , Mlt++ )
   if( *Mlt >= 1e-15 ) {
    if( ( *Mlt >= 1 - RAccSol ) && is_subgradient_global( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ wcomponent_global( i ) - 1 ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }
    !!*/

 }  // end( UpdtCntrs ) - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FormLambda1( double Tau )
{
 // Lambda1 = Lambda + ( Tau / t_master ) * d_master, with
 //   d_master = MasterPB->get_d_vector() : the *displacement* already
 //              scaled by the master's current t (i.e. d_master = t * d_true
 //              for the primal MP, and d_master = -t * z* for the dual MP)
 //   t_master = MasterPB->get_t() : the proximal parameter t that the
 //              master used when computing d_master
 // the ( Tau / t ) normalization supports t-strategies that pass a
 // Tau != t_master; in the common case Tau == t_master the scale is 1
 // and Lambda1 = Lambda + d_master, i.e. the master-prescribed step
 if( MasterPB ) {
  const auto d = MasterPB->get_d_vector();
  const double t_master = MasterPB->get_t();
  const double scale = ( t_master != 0.0 ) ? ( Tau / t_master ) : 1.0;
  for( Index i = 0 ; i < NumVar ; ++i )
   Lambda1[ i ] = Lambda[ i ] + scale * ( i < d.size() ? d[ i ] : 0.0 );
  }

 // walk LamVcblr to clamp Lambda1[ i ] into the effective [ lb , ub ]:
 // ColVariable bounds combined with any supported active bound constraint
 if( MasterPB )
  for( Index i = 0 ; i < NumVar ; ++i ) {
   const auto bounds = effective_bounds( LamVcblr[ i ] );
   if( Lambda1[ i ] < bounds.first )
    Lambda1[ i ] = bounds.first;
   if( Lambda1[ i ] > bounds.second )
    Lambda1[ i ] = bounds.second;
   }

 // move the value from Lambda1 to the ColVariable - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiStatus.assign( NrFi , kUnEval );
 whisG1.assign( NrFi , InINF );

 for( Index i = 0 ; i < NumVar ; i++ )
  LamVcblr[ i ]->set_value( Lambda1[ i ] );

 // compute the upper and lower model at the tentative point - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 UpFiLmb1def = LwFiLmb1def = 0;
 for( Index k = 0 ; k < NrFi ; ++k )
  if( NrEasy && IsEasy[ k ] ) {  // k is an easy component
   UpFiLmb1[ k ] = LwFiLmb1[ k ] = vStar[ k ];  // we know the exact value
   ++UpFiLmb1def;                               // thus both UpFi1 and LwFi1
   ++LwFiLmb1def;                               // are known (and equal)
   }
  else {                         // k is a hard component
   // compute upper and lower bound for k (possibly +/- INF)

   UpFiLmb1[ k ] = INFshift;

   // computing the upper bound requires the Lipschitz constant and is
   // only done if bit 4 of TrgtMng == 1
   if( ( TrgtMng & 16 ) && ( UpFiLmb[ k ] < INFshift ) ) {
    c_VarValue Lk = v_c05f[ k ]->get_Lipschitz_constant();
    if( Lk < INFshift ) {
     UpFiLmb1[ k ] = UpFiLmb[ k ] + Lk * NrmD;
     ++UpFiLmb1def;
     }
    }

   if( vStar[ k ] < INFshift ) {
    LwFiLmb1[ k ] = std::max( UpRifFi[ k ] + vStar[ k ] , LowerBound[ k ] );
    ++LwFiLmb1def;
    }
   else
    LwFiLmb1[ k ] = -INFshift;
   }

 // now compute total upper and lower bound (possibly +/- INF)
 // this requires the value of the linear function, so ensure it is computed
 if( f_lf ) {
  f_lf->compute( true );
  Fi0Lmb1 = rs( f_lf->get_upper_estimate() );
  }
 else
  Fi0Lmb1 = 0;

 if( UpFiLmb1def == NrFi ) {
  ++UpFiLmb1def;  // all components + the sum computed
  // can now compute the total value:
  UpFiLmb1.back() = std::accumulate( UpFiLmb1.begin() , --(UpFiLmb1.end()) ,
                                     Fi0Lmb1 );
  // note that this is the point where the lower bound (if any) is taken
  // into account when computing the total upper function estimate
  if( TrueLB && ( UpFiLmb1.back() < LowerBound.back() ) )
   UpFiLmb1.back() = LowerBound.back();
  }
 else
  UpFiLmb1.back() = INFshift;

 // note that this computation does not make explicit use of the global
 // lower bound; this is not necessary since if it is a "true" lower bound
 // then it is inserted in the master problem, and therefore it determines
 // the values of the vStar[ k ] in such a way that their sum (with all the
 // required corrections) is never < than the global bound
 if( LwFiLmb1def == NrFi ) {
  ++LwFiLmb1def;  // all components + the sum computed
  LwFiLmb1.back() = std::accumulate( LwFiLmb1.begin() , --(LwFiLmb1.end()) ,
                                     Fi0Lmb1 );
  }
 else
  LwFiLmb1.back() = -INFshift;

 // update the upper and lower targets - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that vStar.back() == INFshift if the bundle of any of the hard
 // components is empty, in which case there are no targets: any "finite"
 // information about them is better than what we currently have

 if( vStar.back() < INFshift ) {
  UpTrgt = UpRifFi.back() + ( 1.0 - m2 ) * vStar.back();

  if( m1 > 0 )
   LwTrgt = UpRifFi.back() + vStar.back() + m1 * ( DST + Sigma );
  else
   LwTrgt = UpRifFi.back() + ( 1.0 + m1 ) * vStar.back();
  }
 else {
  UpTrgt = INFshift;
  LwTrgt = -INFshift;
  }

 // print the new information - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_log && ( LogVerb > 4 ) ) {
  *f_log << std::endl << "    Lambda1 = [ ";
  for( auto el : Lambda1 )
   *f_log << el << " ";
  *f_log << "]";

  if( LogVerb > 5 )
   for( Index k = 0 ; k < NrFi ; ++k ) {
    *f_log << std::endl << "    UB[ " << k << " ] = " << def;
    if( f_convex )
     pval( *f_log , UpFiLmb1[ k ] );
    else
     pval( *f_log , - LwFiLmb1[ k ] );
    *f_log << ", LB[ " << k << " ] = ";
    if( f_convex )
     pval( *f_log , LwFiLmb1[ k ] );
    else
     pval( *f_log , - UpFiLmb1[ k ] );
    }
  }
 }  // end( BundleSolver::FormLambda1 )

/*--------------------------------------------------------------------------*/

BundleSolver::Index BundleSolver::InnerLoop( bool extrastep )
{
 // here one might change the value of wFi, corresponding to the first
 // component to be evaluated, if a non-strictly-round-robin order is
 // sought for

 // compute the minimum number of components to evaluate
 Index minceval = ( MinNrEvls >= 0 ? Index( MinNrEvls )
                                   : ( NrFi - NrEasy ) * ( - MinNrEvls ) );
 Index ceval = 0;  // how many components have been evaluated so far

 for( bool insrtd = false ; ; ) {
  // round-robin-like loop between the different components
  if( ! FindNext() ) {  // find next component
   if( ! ceval )        // no component found to evaluate
    Result = kError;    // this is an error (or is it?)
   break;               // anyway, nothing else to do but stop
   }

  if( FiAndGi( f_wFi , ! extrastep ) )
   insrtd = true;

  // return if an unrecoverable error happens
  if( ( FiStatus[ f_wFi ] <= kUnEval ) || ( FiStatus[ f_wFi ] >= kError ) ) {
   Result = kError;
   break;
   }

  // if any component evaluates to -INF, then the whole problem evaluates to
  // -INF and it is therefore unbounded below; due to convexity, it is "very
  // seriously unbounded" in that a convex function being -INF anywhere is
  // -INF everywhere, hence if this ever happens it will do it "very soon"
  // (the very first time the offending component is evaluated)
  if( UpFiLmb1[ f_wFi ] == -INFshift )  {
   UpFiLmb1.back() = -INFshift;
   break;
   }

  if( ! CurrNrEvls[ f_wFi ] )  // not evaluated before
   ++ceval;                    // one more evaluated
  ++CurrNrEvls[ f_wFi ];       // evaluated once more

  // a SS can be performed: note the "<" in the SS condition below (which
  // means it is ever so slightly stronger than it should), which is there
  // to avoid the condition to work when UpFiLmb1.back() == INF == UpTrgt
  if( ( ! MPchgs ) && ( UpFiLmb1.back() < UpTrgt ) )
   MPchgs = 1;

  if( ( ! MPchgs ) && insrtd && RifeqFi && ( LwFiLmb1.back() > LwTrgt ) )
   // doing a NS without possibly evaluating all the components is inhibited
   // if the linearization errors are not computed w.r.t. the "true" value
   // of (every component of); this corresponds to the assumption in the
   // theory that a finite upper bound is known for every component. this
   // implies that eventually all components will be evaluated, which will
   // typically yield a SS
   //
   // for a NS to be performed, LwFiLmb1 must be > than the lower target;
   // again, note the ">" instead of the ">=" (which means this is ever so
   // slightly stronger than it should), which is there to avoid the
   // condition to work when LwFiLmb1.back() == -INF == LwTrgt
   //
   // however, for a NS to guarantee no cycling, at least something must
   // have been inserted (on top of all the other conditions)
   MPchgs = 1;

  if( ( MaxTime < INFshift ) && ( get_elapsed_time() > MaxTime ) ) {
   Result = kStopTime;     // time has ran up
   break;                  // nothing else to do but stop
   }

  if( ( ceval < minceval ) || extrastep )  // too few components compute()-d
   continue;               // do not stop regardless of MPchgs

  if( MPchgs )             // the MP is guaranteed to change
   break;                  // happily stop

  }  // end( for( functions evaluation loop ) )

 return( ceval );

 }  // end( BundleSolver::InnerLoop )

/*--------------------------------------------------------------------------*/

bool BundleSolver::FiAndGi( Index wFi , bool getgi )
{
 // compute and set upper and lower cutoffs and the accuracy - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_log && ( LogVerb > 3 ) )
  *f_log << std::endl << "            Fi[ " << wFi;

 if( getgi )
  SetupFiLambda1( wFi );
 else
  SetupFiLambda( wFi );

 // compute the C05Function and retrieve upper and lower estimates - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto fwFi = v_c05f[ wFi ];

 auto start = std::chrono::system_clock::now();

 FiStatus[ wFi ] = fwFi->compute( ( FiStatus[ wFi ] == kUnEval ) );

 auto end = std::chrono::system_clock::now();
 std::chrono::duration< double > elapsed = end - start;

 if( ( FiStatus[ wFi ] <= kUnEval ) || ( FiStatus[ wFi ] >= kError ) ) {
  if( f_log && ( LogVerb > 3 ) )
   *f_log << " ] = Error #" <<  FiStatus[ wFi ] << ", stop";
  return( false );
  }

 auto ue = fwFi->get_upper_estimate();
 auto le = fwFi->get_lower_estimate();

 if( f_log && ( LogVerb > 3 ) ) {
  *f_log << " ]: UB = "<< def;
  pval( *f_log , ue );
  *f_log << ", LB = ";
  pval( *f_log , le );
  *f_log << " [" << fixd << elapsed.count() << "] " << def;
  }

 // very special case: the sub-problem is unbounded (value == -INF in convex
 // convention, or +INF in concave). Do *not* short-circuit here: even though
 // no finite subgradient exists, the oracle may have produced a feasibility
 // ray (vertical linearization) that the master needs in order to remain
 // bounded at the next iteration. We still propagate the unbounded marker
 // into UpFiLmb1[ wFi ] so GetGi() picks the prefer_vertical branch
 const bool unbounded = f_convex ? ( ue == -INFshift )
                                 : ( le ==  INFshift );

 // if getgi == false the method is actually being called on Lambda, hence
 // it is Lambda's estimates that need be updated, not Lambda1's
 if( ! getgi ) {
  update_UpFiLambd( wFi , f_convex ? ue : - le );
  update_LwFiLambd( wFi , f_convex ? le : - ue );

  // furthermore one immediately returns before getting any linearization
  return( false );
  }

 // update UpFiLambd1[ wFi ] (and possibly UpFiLambd1[ NrFi ])
 update_UpFiLambd1( wFi , f_convex ? ue : - le );

 if( unbounded )
  // skip the remaining accounting (Lipschitz upper-bound refresh and
  // LwFiLambd1 update both rely on a finite value) and go straight to
  // fetching a vertical certificate; GetGi() sees UpFiLmb1[ wFi ] = -INF
  // (resp. +INF in concave) and asks the oracle for a ray first
  return( GetGi( wFi ) );

  // if bit 4 of TrgtMng == 1, then compute the upper bound in Lambda
  // provided by the upper bound in Lambda1 and try to update UpFiLmb[ wFi ]
  // (and possibly UpFiLambd1[ NrFi ])
  // note that, even if this suceeds and therefore decreases UpFiLmb[ wFi ]
  // (and possibly UpFiLambd[ NrFi ], which would be a "rather big" decrease
  // from +INF to something finite), as the theory requires the upper target
  // is *not* changed
  if( ( TrgtMng & 16 ) && ( UpFiLmb1[ wFi ] < INFshift ) ) {
   c_VarValue LwFi = fwFi->get_Lipschitz_constant();
   if( LwFi < INFshift )
    update_UpFiLambd( wFi , UpFiLmb1[ wFi ] + LwFi * NrmD );
   }

 // update LwFiLambd1[ wFi ] (and possibly LwFiLambd1[ NrFi ])
 update_LwFiLambd1( wFi , f_convex ? le : - ue );

 // get new linearizations - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( GetGi( wFi ) );

 }  // end( BundleSolver::FiAndGi )

/*--------------------------------------------------------------------------*/

void BundleSolver::SetupFiStrPar( Index wFi )
{
 // set the component-specific string parameters

 if( v_C05_EI_SPAR_Names.empty() )  // ... if any
  return;

 if( v_C05_EI_SPAR_Names.size() > v_C05_EI_SPAR_Vals.size() )
  throw( std::logic_error(
             "BundleSolver::SetupFiStrPar: "
             "vstr_C05_SPAR_EI_Names.size() > vstr_C05_SPAR_EI_Vals.size()" ) );

 ComputeConfig CwFi;
 CwFi.set_diff( true );
 CwFi.set_relax( true );
 auto Vit = v_C05_EI_SPAR_Vals.begin();
 for( const auto & name : v_C05_EI_SPAR_Names ) {
  auto par = ps_insert( *(Vit++) ,
                        "_" + std::to_string( wFi ) +
                        "_" + std::to_string( get_elapsed_calls() ) +
                        "_" + std::to_string( get_elapsed_iterations() ) );
  if( ( name.size() > 4 ) && ( name.substr( 0 , 4 ) == "vstr" ) )
   CwFi.set_par( std::string( name ) ,
                 std::vector< std::string >( { par } ) );
  else
   CwFi.set_par( std::string( name ) , std::move( par ) );
  }

 v_c05f[ wFi ]->set_ComputeConfig( & CwFi );

 }  // end( BundleSolver::SetupFiStrPar )

/*--------------------------------------------------------------------------*/

void BundleSolver::SetupFiLambda1( Index wFi )
{
 auto fwFi = v_c05f[ wFi ];

 // start by setting a "time cutoff" with the remaining total time - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( MaxTime < INFshift )
  fwFi->set_par( dblMaxTime , MaxTime - get_elapsed_time() );

 // set the component-specific string parameters, if any - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 SetupFiStrPar( wFi );

 if( ! ( TrgtMng & 15 ) )  // if target management is not active
  return;                  // all done

 // compute upper and lower cutoffs and the accuracy - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 auto UpCutOff = INFshift;
 auto LwCutOff = -INFshift;

 if( ( TrgtMng & 7 ) &&
     ( UpFiLmb1.back() > UpTrgt ) && ( LwFiLmb1.back() < LwTrgt ) ) {
  // finite upper and lower cutoffs can (and make sense to) be set only if the
  // upper and lower targets are finite (they either both are or none of them
  // are), which depends on the fact that v* = vStar.back() (and its "cousin"
  // Delta* = D_t( z* ) + Sigma*) is finite, which means that *all* non-easy
  // components have diagonal linearizations in their bundle
  //
  // note that this implies that vStar[ wFi ] (for this particular component)
  // must also be well-defined (< INF)
  //
  // if some component does not have any linearization yet, any "finite"
  // information about them is better than what we currently have, and hence
  // the cutoffs should be "as weak as they can be"
  //
  // also, the upper cutoff only makes sense to be set if the SS condition
  //
  //  \bar{f}_{tot}( x ) <= \bar{tau}  \equiv  UpFiLmb1.back() <= UpTrgt
  //
  // it is not satisfied already, i.e., if \bar{f}_{tot}( x ) > \bar{tau}
  // (note that this can never hold if \bar{tau} = +INF, due to the ">",
  // which ensures that the upper target is finite when the if() is entered).
  // indeed, the second part of the UpCutOff computation corresponds to the
  // following argument: we want the SS condition
  //
  //  \bar{f}_{tot}( x ) = \sum_i \bar{f}_i( x ) <= \bar{tau}
  //
  // to hold, which means
  //
  //  \bar{f}_k( x ) <= \bar{tau} - \sum_{i \neq k} i \bar{f}_i( x )
  //
  // but
  //
  //  \sum_{i \neq k} i \bar{f}_i( x ) = \bar{f}_{tot}( x ) - \bar{f}_k( x )
  //
  // and hence we will be content if
  //
  //  \bar{f}_k( x ) <= \bar{tau} - ( \bar{f}_{tot}( x ) - \bar{f}_k( x ) )
  //
  // this is of course conditional to the fact that the function value is
  // available, i.e., \bar{f}_{tot}( x ) < INF ==> \bar{f}_k( x ) < INF, but
  // it is also conditional to the fact that the SS condition does not hold
  // already. in fact, in the above condition \bar{f}_k( x ) on the left means
  // "after having computed the function value at x", whereas \bar{f}_k( x )
  // on the right means "the current value it has". if the SS condition
  //
  //  \bar{f}_{tot}( x ) <= \bar{tau}
  //
  // holds already when the method is called, then
  //
  //  \bar{f}_k( x ) <= \bar{tau} - ( \bar{f}_{tot}( x ) - \bar{f}_k( x ) )
  //
  // also holds with \bar{f}_k( x ) meaning "the current value it has" in
  // both places, and therefore *a fortiori* after that the function is
  // computed. thus, if the SS condition holds already, then there is no
  // reason to put a finite upper cutoff, since any value would do.
  //
  // of course the symmetric argument holds for the lower cutoff and the
  // NS condition
  //
  //  \underline{f}_{tot}( x ) >= \underline{tau}
  //
  // and note that if either one of the two conditions holds already, then
  // no cutoff (be it upper or lower) need be set, and neither does the
  // accuracy, because the fate of the step is already sealed whatever the
  // information that the oracle provides

  // note: the reason why upper and lower cutoffs are not entirely separated
  // is the computation of BetaK(), which is currently very easy but it may
  // one day become more sophisticated and therefore costly
  const auto bk = BetaK( wFi );
  const auto LwFiK = UpRifFi[ wFi ] + vStar[ wFi ];

  UpCutOff = LwFiK - m2 * bk * vStar.back();
  LwCutOff = LwFiK + std::abs( m1 ) * bk * ( DST + Sigma );

  if( UpCutOff < LwCutOff )  // this should never happen, but account for
   LwCutOff = UpCutOff;      // numerical issues

  // note that here UpCutOff >= LwCutOff and the next step can only increase
  // UpCutOff and decrease LwCutOff, so the relationship will keep holding

  // note that UpFiLmb1.back() < INFshift ==> UpFiLmb1[ wFi ] < INFshift
  if( UpFiLmb1.back() < INFshift )
   UpCutOff = std::max( UpTrgt - ( UpFiLmb1.back() - UpFiLmb1[ wFi ] ) ,
                        UpCutOff );

  // note that LwFiLmb1.back() > -INFshift ==> LwFiLmb1[ wFi ] > -INFshift
  if( LwFiLmb1.back() > -INFshift )
   LwCutOff = std::min( LwTrgt - ( LwFiLmb1.back() - LwFiLmb1[ wFi ] ) ,
                        LwCutOff );
  }

 // assign the cutoff values to the C05Function - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ( TrgtMng & 1 ) && f_convex ) ||
     ( ( TrgtMng & 2 ) && ( ! f_convex ) ) ) {
  auto lt = f_convex ? LwCutOff : - UpCutOff;
  if( f_log && ( LogVerb > 3 ) )
   *f_log << " ~ lt = " << def << lt;
  fwFi->set_par( dblLwCutOff , lt );
  }

 if( ( ( TrgtMng & 2 ) && f_convex ) ||
     ( ( TrgtMng & 1 ) && ( ! f_convex ) ) ) {
  auto ut = f_convex ? UpCutOff : - LwCutOff;
  if( f_log && ( LogVerb > 3 ) )
   *f_log << " ~ ut = " << def << ut;
  fwFi->set_par( dblUpCutOff , ut );
  }

 if( TrgtMng & 12 ) {
  double EpsCurr = 100;  // 1e+2 relative error is "a finite INF"

 // set EpsCurr to the difference between upper and lower cutoffs if they
 // are finite, and leave it to "INF" == "anything goes" otherwise
  if( ( TrgtMng & 4 ) &&
      ( UpCutOff < INFshift ) && ( LwCutOff > -INFshift ) )
   EpsCurr = ( UpCutOff - LwCutOff ) /
                           std::max( 1.0 , std::abs( UpRifFi[ wFi ] ) );

  // set EpsCurr to (almost) the current achieved accuracy EpsU as
  // computed via tStar if this is larger than the one before (or no
  // value has been set before)
  if( ( TrgtMng & 8 ) &&
      ( ( EpsCurr < EpsU / Nearly ) || ( EpsCurr == 100 ) ) )
   EpsCurr = EpsU / Nearly;

  // EpsCurr never needs be (much) smaller than the required relative
  // accuracy for the overall computation,
  EpsCurr = std::max( EpsCurr , RelAcc / Nearly );

  if( f_log && ( LogVerb > 3 ) )
   *f_log << " ~ eps = " << shrt << EpsCurr;
  fwFi->set_par( dblRelAcc , EpsCurr );
  }

 }  // end( BundleSolver::SetupFiLambda1 )

/*--------------------------------------------------------------------------*/

void BundleSolver::SetupFiLambda( Index wFi )
{
 auto fwFi = v_c05f[ wFi ];

 // start by setting a "time cutoff" with the remaining total time - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( MaxTime < INFshift )
  fwFi->set_par( dblMaxTime , MaxTime - get_elapsed_time() );

 // set the component-specific string parameters, if any - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 SetupFiStrPar( wFi );

 if( ! ( TrgtMng & 15 ) )  // if target management is not active
  return;                  // all done

 // compute upper and lower cutoffs and the accuracy- - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this is just to recompute an already computed value, hence getting back
 // what one already had is OK

 double UpCutOff = UpFiLmb1[ wFi ];
 double LwCutOff = LwFiLmb1[ wFi ];

 double EpsCurr = ( UpCutOff - LwCutOff ) /
                                std::max( 1.0 , std::abs( UpRifFi[ wFi ] ) );

 // assign the cutoff values to the C05Function - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_convex ) {
  if( TrgtMng & 1 )
   fwFi->set_par( dblLwCutOff , LwCutOff );

  if( TrgtMng & 2 )
   fwFi->set_par( dblUpCutOff , UpCutOff );
  }
 else {
  if( TrgtMng & 1 )
   fwFi->set_par( dblUpCutOff , - LwCutOff );

  if( TrgtMng & 2 )
   fwFi->set_par( dblLwCutOff , - UpCutOff );
  }

 if( TrgtMng & 12 )
  fwFi->set_par( dblRelAcc , EpsCurr );

 }  // end( BundleSolver::SetupFiLambda )

/*--------------------------------------------------------------------------*/

bool BundleSolver::GetGi( Index wFi )
{
 bool insrtd = false;  // keep track if anything new at all was inserted
 auto fwFi = v_c05f[ wFi ];

 for( Index Ftchd = 0 ; Ftchd < aBP3 ; ++Ftchd ) {
  bool diagonal = true;
  bool HasLinearization;

  // when Fi at Lambda1 is *not* finite the sub-problem either is unbounded
  // (its solver has produced a recession direction rather than a primal
  // solution, signaled by UpFiLmb1 == -INFshift in convex-min convention
  // or +INFshift in the not-yet-computed flavour) or has not been touched
  // yet; in both cases we look for a vertical (feasibility) linearization
  // first and only fall back to a diagonal one if no ray is available
  const bool prefer_vertical = ( UpFiLmb1[ wFi ] == INFshift ) ||
                               ( UpFiLmb1[ wFi ] == -INFshift );
  if( ! Ftchd ) {  // first look for a constraint then for a sub-gradient
   if( prefer_vertical ) {
    HasLinearization = fwFi->has_linearization( diagonal = false );
    if( ! HasLinearization )
     HasLinearization = fwFi->has_linearization( diagonal = true );
    }
   else
    HasLinearization = fwFi->has_linearization( diagonal );
   }
  else {
   if( prefer_vertical ) {
    HasLinearization = fwFi->compute_new_linearization( diagonal = false );
    if( ! HasLinearization )
     HasLinearization = fwFi->compute_new_linearization( diagonal = true );
    }
   else
    HasLinearization = fwFi->compute_new_linearization( diagonal );
   }

  if( ! HasLinearization )  // no new linearization of either type available
   break;                   // nothing else to do

  if( ! diagonal )          // a vertical linearization changes the MP
   MPchgs = 2;              // no matter what else happens

  // check if aggregation has to be performed - - - - - - - - - - - - - - - -
  // doing this now could occasionally result in useless aggregations, but it
  // is necessary due to limitations in the MPSolver interface (there can be
  // only one "un-named item being inserted", so inserting Z[ wFi ] while
  // inserting the new item is complicated

  auto wh = BStrategy( wFi );

  // local buffer for the new linearization coefficients
  //
  // in sparse Lambda mode v_c05f[ wFi ] writes only its own loc_NV
  // active-var coefficients into G1k[ 0 .. loc_NV - 1 ]; the master will
  // read them through the v_local2global[ wFi ] map. In dense mode loc_NV
  // == NumVar and the map is empty.

  const Index loc_NV = f_sparse_lambda ? fwFi->get_num_active_var() : NumVar;
  std::vector< double > G1k_buf( loc_NV , 0.0 );
  double * G1k = G1k_buf.data();

  // fetch the item from the Oracle - - - - - - - - - - - - - - - - - - - - -

  fwFi->get_linearization_coefficients( G1k );
  if( ! f_convex )
   chgsign( G1k , loc_NV );

  auto Alfa1k = rs( fwFi->get_linearization_constant() );

  // The master receives the raw cut constant (no promotion to a Lambda-
  // dependent linearization error): each cut is shipped as ( g , alpha_raw )
  // and the master keeps that form. The translation to the linearization-
  // error form needed by the bundle stop test is done internally by
  // MasterProblemBlock::get_aggregated_alpha using the cached F_k( x_bar )
  // installed via set_reference. Alfa1k_for_master is therefore the
  // un-promoted alpha_raw (already in convex-min sign via rs()), while the
  // local Alfa1k is overwritten further below as the linearization error
  // in Lambda1, a separate form used by the heuristic accumulators Alfa1,
  // eps, ScPr1
  auto Alfa1k_for_master = Alfa1k;

  // eps is only meaningful in the diagonal branch below, where it gets
  // set to the linearization error in Lambda1 of the new subgradient.
  // The vertical-constraint branch does not touch it, and the only
  // downstream reader (the LogVerb > 2 log block) gates the access
  // behind the same diagonal flag, so eps would be safe even
  // uninitialised; default-initialise it to 0 to silence the static
  // analyser warning and keep the variable harmless if a future change
  // ever lifts the diagonal gate.
  double eps = 0;

  // compute ScPr1k and Alfa1k- - - - - - - - - - - - - - - - - - - - - - - -

  Index cp = InINF;
  double ScPr1k;
  bool is_rep = diagonal && ( ! Ftchd ) && ( ! CurrNrEvls[ wFi ] );
  // if it is the first subgradient of the first call to GetGi() for this
  // component, it may be the "representative subgradient"

  if( diagonal ) {  // it is a subgradient
   // compute the lower bound in Lambda provided by the subgradient
   double FikLmb;
   if( f_sparse_lambda ) {
    FikLmb = Alfa1k;
    const auto & m = v_local2global[ wFi ];
    for( Index li = 0 ; li < loc_NV ; ++li )
     FikLmb += Lambda[ m[ li ] ] * G1k[ li ];
    }
   else
    FikLmb = Alfa1k + std::inner_product( Lambda.begin() , Lambda.end() ,
                                          G1k , double( 0 ) );

   // try to update LwFiLmb[ wFi ] (and possibly LwFiLambd.back())
   // note that, even if this suceeds and therefore increases LwFiLmb[ wFi ]
   // (and possibly LwFiLambd.back(), which would be a "rather big" increase
   // from -INF to something finite), as the theory requires the lower target
   // is *not* changed
   update_LwFiLambd( wFi , FikLmb );

   // compute the linearization error in Lambda1
   if( f_sparse_lambda ) {
    double ip = 0;
    const auto & m = v_local2global[ wFi ];
    for( Index li = 0 ; li < loc_NV ; ++li )
     ip += Lambda1[ m[ li ] ] * G1k[ li ];
    Alfa1k = UpFiLmb1[ wFi ] - Alfa1k - ip;
    }
   else
    Alfa1k = UpFiLmb1[ wFi ] - Alfa1k -
     std::inner_product( Lambda1.begin() , Lambda1.end() , G1k ,
                         double( 0 ) );

   // this is the eps so that G1 is an eps-subgradient in Lambda1
   eps = Alfa1k;

   // (g, alpha) is built locally and pushed directly into
   // MasterPB->add_cut, so no in-place rewrite is needed; duplicate
   // detection is delegated to MasterPB further below
   if( f_sparse_lambda ) {
    ScPr1k = 0;
    const auto & m = v_local2global[ wFi ];
    for( Index li = 0 ; li < loc_NV ; ++li )
     ScPr1k += Lambda[ m[ li ] ] * G1k[ li ];
    }
   else
    ScPr1k = std::inner_product( Lambda.begin() , Lambda.end() , G1k ,
                                 double( 0 ) );

   #if CHECK_BAD_F & 2
    // if so required, check for "very negative Alfa1" and print a warning on
    // std::cerr if it is found; however, avoid the check if the true function
    // value of component wFi in Lambda has not been ever computed yet (or has
    // became invalid), as testified by the fact that UpRifFi[ wFi ] (which is
    // always finite) is different from UpFiLmb[ wFi ] (which is set to +INF
    // when unknown), since in this case a negative Alfa1 does not signify
    // that anything untowards has been done by the oracle
    if( UpRifFi[ wFi ] == UpFiLmb[ wFi ] ) {
     auto rel_error = ( Alfa1k / std::max( std::abs( UpFiLmb1[ wFi ] ) , 1.0 )
                        );
     if( rel_error < - CHECK_BAD_F_EPS )
      std::cerr << std::endl << "Warning[ " << SCalls << ", " << ParIter
                << ", " << wFi << " ]: Alfa1 = " << rel_error << std::endl;
     }
   #endif
   }
  else {             // it is a (vertical) constraint
   // the standard form of constraints in the master problem is
   //
   //       [ v / 0 ] >= g d - \alpha
   //
   // while the diagonal linearizations are
   //
   //       ( 1 , - g ) ( v , x ) >= \alpha
   //
   // The master receives the physical row pair (g, alpha) for both diagonal
   // and vertical linearizations; MasterProblemBlock converts it to the
   // active primal/dual storage convention. Keep the local sign used by the
   // old heuristic accumulators, but do not pre-convert the master copy.
   Alfa1k = - Alfa1k;
   Alfa1k_for_master = rs( fwFi->get_linearization_constant() );
   if( f_sparse_lambda ) {
    ScPr1k = 0;
    const auto & m = v_local2global[ wFi ];
    for( Index li = 0 ; li < loc_NV ; ++li )
     ScPr1k += Lambda[ m[ li ] ] * G1k[ li ];
    }
   else
    ScPr1k = std::inner_product( Lambda.begin() , Lambda.end() , G1k ,
                                 double( 0 ) );
   }

  // helper: materialise a dense (NumVar-sized) physical copy of G1k either by
  // scattering sparse coefficients into their global slots or, in dense mode,
  // by a straight copy. MasterProblemBlock owns the conversion from this
  // physical cut convention to the active primal/dual storage convention.
  auto make_dense_g1 = [ & ]() -> std::vector< double > {
   if( f_sparse_lambda ) {
    std::vector< double > g( NumVar , 0.0 );
    const auto & m = v_local2global[ wFi ];
    for( Index li = 0 ; li < loc_NV ; ++li )
     g[ m[ li ] ] = G1k[ li ];
    return( g );
    }
   std::vector< double > g( NumVar );
   for( Index j = 0 ; j < NumVar ; ++j )
    g[ j ] = G1k[ j ];
   return( g );
   };

  // helper: accumulate G1k into the representative-subgradient sum G1,
  // scattering through v_local2global[ wFi ] when in sparse mode.
  auto accumulate_G1 = [ & ]() {
   if( f_sparse_lambda ) {
    const auto & m = v_local2global[ wFi ];
    for( Index li = 0 ; li < loc_NV ; ++li )
     G1[ m[ li ] ] += G1k[ li ];
    }
   else
    vect_sum( G1 , G1k );
   };

  auto dense_g = make_dense_g1();
  if( MasterPB ) {
   const int identical =
    MasterPB->find_identical_cut( int( wFi ) , dense_g , ! diagonal );
   if( identical >= 0 )
    cp = Index( identical );
   }

  Index gpp = InINF;  // position in the global pool where to put it

  if( f_log && ( LogVerb > 2 ) ) {
   *f_log << std::endl << "            New " << shrt;
   if( diagonal ) {
    if( eps >= std::max( std::abs( UpRifFi[ wFi ] ) , double( 1 ) )
               * RelAcc / 10 )
     *f_log << "eps-subgradient with eps = " << eps;
    else
     *f_log << "subgradient";
    *f_log << " for Fi[ " << wFi << " ] ~ Alfa1 = " << Alfa1k_for_master
           << " ~ gd = " << rs( ScPr1k );
    }
   else
    // note: we don't print wh here because it has not been finalized yet
    // (BStrategy() may have returned InINF, the actual slot is decided
    // later by FindAPlace()); the "stored in <wh> (<gpp>)" log line below
    // shows the real slot
    *f_log << "constraint for Fi[ " << wFi << " ] ~ rhs = " << Alfa1k;
   }

  bool to_insert = true;  // if it has to be inserted

  if( cp < InINF ) {  // the item is a copy - - - - - - - - - - - - - - - - -
   BLOG( 2 , " is copy of " << cp << " (" << ItemVcblr[ cp ].second << ")" );

   wh = cp;  // we have it already

   const auto OldA1k = read_alpha_global( cp );
   const auto NewA1k =
    MasterPB->get_stored_constant( int( wFi ) , dense_g ,
                                   Alfa1k_for_master , ! diagonal );

   assert( ( ItemVcblr[ cp ].first == wFi ) &&
           ( ItemVcblr[ cp ].second < vBPar2[ wFi ] ) &&
           ( InvItemVcblr[ wFi ][ ItemVcblr[ cp ].second ] == cp ) );

   const auto AlfaTol = std::max( std::abs( NewA1k ) , double( 1 ) )
                        * RelAcc / 10;
   if( MasterPB->is_stronger_constant( int( wFi ) , OldA1k , NewA1k ,
                                       AlfaTol ) ) {
    // replace the original when the copy is substantially stronger; the
    // direction of the comparison depends on the PolyhedralFunction sense

    BLOG( 2 , " with stronger Alfa" );

    gpp = ItemVcblr[ cp ].second;
    if( ( BPar7 & 3 ) < 3 ) {
     // BundleSolver does not immediately replace the copy unless necessary,
     // but clearly if one linearization in the global pool has to be
     // sacrificed, it'll be the copy
     auto ngpp = find_place_in_global_pool( wFi );
     if( ngpp < InINF ) {       // a free place has been found
      // although the old linearization is kept, it is removed from the
      // bundle: the position cp is now associated with ngpp, which
      // means that position gpp is now free
      remove_from_global_pool( wFi , gpp , false );
      gpp = ngpp;                      // store the copy there
      BLOG( 2 , " (" << gpp << ")" );  // print the chosen place
      }
     }

    // if it is the "representative subgradient", add its contribution to
    // the required ones of Alfa1, ScPr1 and G1 (if any); do this before
    // the call to modify_cut() because the state of the G1k memory after
    // the call is unclear
    if( is_rep ) {
     whisG1[ wFi ] = cp;
     if( NeedsAlfa1() )
      Alfa1 += Alfa1k;
     if( NeedsScPr1() )
      ScPr1 += ScPr1k;
     if( NeedsG1() )
      accumulate_G1();
     }

    if( MasterPB ) {
     // replace the cut at global bundle slot cp of HardCmps[ wFi ]
     // with the new ( G1k , Alfa1k ) pair
     MasterPB->modify_cut( int( wFi ) , int( cp ) ,
                           std::move( dense_g ) , Alfa1k_for_master );
     }
    // note that the number of items of component wFi in the master problem
    // is unchanged
    }
   else                 // the item is a copy, not better than the original
    to_insert = false;  // do nothing
   }
  else {           // the item is not a copy- - - - - - - - - - - - - - - - -
   // insert the item, if there is space

   if( wh == InINF )  // the position has not been selected in BStrategy()
    wh = FindAPlace( wFi );  // find a free spot in the bundle

   if( wh == InINF ) {  // no space found ...
    if( ! Ftchd ) {     // ... and this was the first item
     BLOG( 1 , std::endl << " ERROR: No space in the bundle" << std::endl );
     Result = kError;   // signal an error to end the outer Fi-cycle
     }
    else
     BLOG( 1 , std::endl << " WARNING: No space in the bundle" << std::endl );
    break;              // the cycle ends
    }

   if( ItemVcblr[ wh ].second < vBPar2[ wFi ] )
    // the place is occupied already: this happens if the bundle was full
    // (and, possibly aggregation has been performed for safety)
    remove_cut_global( wh );  // the old item has to be removed first
   else {                   // the place is unoccupied
    ++NrItems[ wFi ];       // one more item in the bundle (otherwise the
    ++NrItems[ NrFi ];      // number remains the same as one is replaced)
    }

   // if it is the "representative subgradient", add its contribution to
   // the required ones of Alfa1, ScPr1 and G1 (if any); do this before
   // the call to add_cut() because the state of the G1k memory after
   // the call is unclear
   if( is_rep ) {
    whisG1[ wFi ] = wh;
    if( NeedsAlfa1() )
     Alfa1 += Alfa1k;
    if( NeedsScPr1() )
     ScPr1 += ScPr1k;
    if( NeedsG1() )
     accumulate_G1();
    }

   // ( G1k , Alfa1k_for_master ) is the cut at slot wh of HardCmps[ wFi ]
   if( MasterPB )
    MasterPB->add_cut( int( wFi ) , int( wh ) , std::move( dense_g ) ,
                       Alfa1k_for_master , ! diagonal );

   // now find a position in the global pool of component wFi where to store
   // the new linearization
   gpp = find_place_in_global_pool( wFi );

   if( gpp == InINF ) {     // there is none
    // this means that not only the global pool is full, but also the bundle
    // is full: one can therefore put it in the very place of the item it
    // replaces, which must be an item of the same component because
    // BStrategy() ensures this
    assert( ItemVcblr[ wh ].first == wFi );
    gpp = ItemVcblr[ wh ].second;
    assert( gpp < vBPar2.back() );
    }

   BLOG( 2 , " stored in " << wh << " (" << gpp << ")"  );
   }

  // if something was inserted, bookkeeping is needed - - - - - - - - - - - -

  if( to_insert ) {
   inhibit_Modification( true );
   v_c05f[ wFi ]->store_linearization( gpp );
   inhibit_Modification( false );

   insrtd = true;
   add_to_global_pool( wFi , gpp , wh );

   if( diagonal )       // it is a subgradient
    OOBase[ wh ] = -1;  // ensure it won't be touched again this round
   else                 // it is a constraint
    // mark it as permanently fixed: this may be a bad choice in practice,
    // although it is required by the theory (we'll see ...)
    OOBase[ wh ] = -Inf< SIndex >();
   }

  #if CHECK_DS & 1
   CheckBundle();
  #endif

  }  // end( items-collecting loop )- - - - - - - - - - - - - - - - - - - - -

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( insrtd );  // returns true if at least one item was inserted

 }  // end( BundleSolver::GetGi )

/*--------------------------------------------------------------------------*/

void BundleSolver::GotoLambda1( void )
{
 // compute DeltaFi - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector< VarValue > DF( NrFi + 1 );
 /* DeltaFi = UpFiLmb1 - UpRifFi; note that the code one may expect
  *
  * std::transform( UpFiLmb1.begin() , UpFiLmb1.end() , UpRifFi.begin() ,
  *                 DF.begin() , std::minus< double >() );
  *
  * is wrong since the format of DeltaFi expected by ChangeCurrPoint() is
  * different from the one used in BundleSolver; in particular, the total
  * value need be in DF.front() rather than in DF.back(), and the value for
  * component i need be in DF[ i + 1 ] rather than in DF[ i ]. */

 DF.front() = UpFiLmb1.back() - UpRifFi.back();
 std::transform( UpFiLmb1.begin() , --(UpFiLmb1.end()) , UpRifFi.begin() ,
                 ++(DF.begin()) , std::minus< double >() );

 // do the move - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // Lambda = Lambda1 and all associated data structures

 Lambda.swap( Lambda1 );
 UpFiLmb.swap( UpFiLmb1 );
 LwFiLmb.swap( LwFiLmb1 );
 UpRifFi = UpFiLmb;
 RifeqFi = true;
 UpFiLmbdef = UpFiLmb1def;
 LwFiLmbdef = LwFiLmb1def;
 Fi0Lmb = Fi0Lmb1;

 // change the current point in the Master Problem - - - - - - - - - - - - -

 if( MasterPB ) {
  // atomic refresh of the master reference (x_bar + per-hard F_k(x_bar)):
  // the cache fed here will eventually drive the on-demand translation of
  // raw cut alphas to linearization errors inside MasterPB, removing the
  // shift loop that follows. For now both sides run: the cache is poked
  // and the legacy promoted alphas are still refreshed below
  {
   std::vector< double > F_hard;
   F_hard.reserve( NrFi );
   for( Index k = 0 ; k < NrFi ; ++k ) {
    if( NrEasy && IsEasy[ k ] )
     continue;
    F_hard.push_back( UpRifFi[ k ] );
    }
   MasterPB->set_reference( Lambda , F_hard );
   }
  // set_x_bar is intentionally NOT called separately: set_reference
  // already forwards x_bar to set_x_bar internally

  // Every cut in the bundle stores its raw native constant alpha_raw_i
  // (immutable post-insertion). The translation to the bundle-method
  // linearization-error form is rebuilt on demand by
  // MasterProblemBlock::get_aggregated_alpha from the cached
  // F_k( x_bar ) installed by set_reference above and the cached
  // x_bar. The reference change Lambda_old -> Lambda_new is therefore
  // *fully* expressed by the set_reference call: no per-cut alpha
  // refresh is needed (vertical cuts were already invariant under
  // reference moves and stay that way)
  }

 #if CHECK_DS & 4
  CheckAlpha();
 #endif
 #if CHECK_DS & 8
  CheckLBs();
 #endif

 }  // end( BundleSolver::GotoLambda1 )

/*--------------------------------------------------------------------------*/

void BundleSolver::GotoLambda( void )
{
 UpRifFi = UpFiLmb;        // set UpFiLmb as the reference values
 RifeqFi = true;

 // refresh the master reference cache: only the per-component
 // F_k( x_bar ) values change (the stability centre x_bar itself is the
 // same as before, since GotoLambda is reached after a NS that does not
 // move the centre). Stored alpha_raw values are immutable, so no
 // per-cut shift is needed: the linearization-error form is rebuilt on
 // demand by MasterProblemBlock::get_aggregated_alpha
 if( MasterPB ) {
  std::vector< double > F_hard;
  F_hard.reserve( NrFi );
  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( NrEasy && IsEasy[ k ] )
    continue;
   F_hard.push_back( UpRifFi[ k ] );
   }
  MasterPB->set_reference( Lambda , F_hard );
  }

 #if CHECK_DS & 4
  CheckAlpha();
 #endif
 #if CHECK_DS & 8
  CheckLBs();
 #endif

 }  // end( BundleSolver::GotoLambda )

/*--------------------------------------------------------------------------*/

void BundleSolver::ResetAlfa( Index k )
{
 std::vector< VarValue > Gi( NumVar );
 // size the Alfa buffer to the global upper bound on names; the
 // slot-indexed MasterPB does not expose a per-cmp MaxName separately
 std::vector< VarValue > Alfa( Index( vBPar2.back() ) );
 if( k == NrFi ) {  // all components need be reset
  for( Index i = 0 ; i < get_max_name() ; ++i )
   if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) {
    auto kk = ItemVcblr[ i ].first;
    auto nm = ItemVcblr[ i ].second;
    auto Ai = rs( v_c05f[ kk ]->get_linearization_constant( nm ) );
    #ifndef NDEBUG
     if( std::isnan( Ai ) )  // linearization no longer valid
      throw( std::logic_error(
                 "BundleSolver::ResetAlfa: "
                 "inconsistent ItemVcblr" ) );
    #endif

    // recover the linearization coefficients
    v_c05f[ kk ]->get_linearization_coefficients( Gi.data() ,
                                                  Range( 0 , NumVar ) , nm );
    if( ! f_convex )
     chgsign( Gi.data() , NumVar );

    // MasterPB owns the raw -> linearization-error translation ( b = F_k -
    // alpha + g . x_bar, using the cut's stored g ), so feed it the RAW
    // constant, not the pre-translated lin-error: passing the lin-error here
    // would make modify_alpha translate a second time. Vertical cuts store
    // the raw constant directly ( sign-flipped ), diagonal cuts the raw alpha.
    if( v_c05f[ kk ]->is_linearization_vertical( nm ) )
     Alfa[ i ] = - Ai;
    else
     Alfa[ i ] = Ai;
    }
  }
 else {             // only that specific component need be reset
  for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
   if( InvItemVcblr[ k ][ i ] < vBPar2.back() ) {
    auto Ai = rs( v_c05f[ k ]->get_linearization_constant( i ) );

    #ifndef NDEBUG
     if( std::isnan( Ai ) )  // linearization no longer valid
      throw( std::logic_error(
                 "BundleSolver::ResetAlfa: "
                 "inconsistent ItemVcblr" ) );
    #endif

    // recover the linearization coefficients
    v_c05f[ k ]->get_linearization_coefficients( Gi.data() ,
                                                 Range( 0 , NumVar ) , i );
    if( ! f_convex )
     chgsign( Gi.data() , NumVar );

    // feed MasterPB the RAW constant ( it owns the raw -> lin-error
    // translation ); vertical cuts sign-flipped, diagonal cuts the raw alpha
    if( v_c05f[ k ]->is_linearization_vertical( i ) )
     Alfa[ InvItemVcblr[ k ][ i ] ] = - Ai;
    else
     Alfa[ InvItemVcblr[ k ][ i ] ] = Ai;
   }
  }

 if( MasterPB ) {
  if( k == NrFi ) {
   // global reset: every hard component
   for( Index kk = 0 ; kk < NrFi ; ++kk )
    if( ! ( NrEasy && IsEasy[ kk ] ) )
     MasterPB->set_alphas_bulk( int( kk ) , Alfa );
   }
  else
   MasterPB->set_alphas_bulk( int( k ) , Alfa );
  }

 }  // end( BundleSolver::ResetAlfa )

/*--------------------------------------------------------------------------*/

void BundleSolver::SimpleBStrat( void )
{
 if( ( BPar7 & 3 ) == 3 ) {  // "eager" deletion
  std::vector< Subset > tbdltd( NrFi );
  for( Index i = 0 ; i < get_max_name() ; ++i )
   if( ( OOBase[ i ] < Inf< SIndex >() ) &&
       ( OOBase[ i ] > SIndex( BPar1 ) ) ) {
    tbdltd[ ItemVcblr[ i ].first ].push_back( ItemVcblr[ i ].second );
    Delete( i );
    }

  inhibit_Modification( true );
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ! tbdltd[ k ].empty() )
    v_c05f[ k ]->delete_linearizations( std::move( tbdltd[ k ] ) , false );
  inhibit_Modification( false );
  }
 else                        // "lazy" deletion
  for( Index i = 0 ; i < get_max_name() ; ++i )
   if( ( OOBase[ i ] < Inf< SIndex >() ) && ( OOBase[ i ] > SIndex( BPar1 ) ) )
    Delete( i );

 #if CHECK_DS & 1
  CheckBundle();
 #endif

 }  // end( BundleSolver::SimpleBStrat )

/*--------------------------------------------------------------------------*/

double BundleSolver::BetaK( Index wFi ) {
 return( 1.0 / double( NrFi - NrEasy ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::Log1( void )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "{" << SCalls << "-" << ParIter << "-"
        << NrItems.back() << "-" << fixd << get_elapsed_time() << "} ";

 if( UsesPureLevelStabilization() )
  *f_log << "Lvl = " << shrt << f_level_value;
 else
  *f_log << "t = " << shrt << t;

 if( ( tStar > 0 ) || ( NrmZFctr == INFshift ) )
  *f_log << " ~ D*_1( z* ) = " << read_DStart( 1 );
 else
  *f_log << " ~ || z* || = " << NrmZ / NrmZFctr;

 *f_log << " ~ Sigma = " << Sigma << std::endl << "           ";

 if( UpFiLmb.back() == INFshift )
  *f_log << " Fi undefined";
 else
  *f_log << " Fi = " << def << rs( UpFiLmb.back() ) << " ~ eU = "
         << shrt << EpsU;

 if( BPar6 )
  *f_log << " ~ BP3 = " << aBP3;

 }  // end( BundleSolver::Log1 )

/*--------------------------------------------------------------------------*/

void BundleSolver::Log2( double ft )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "            [" << fixd << ft << "] " << def;

 if( LowerBound.back() > -INFshift ) {
  if( f_convex )
   *f_log << "LB = " << LowerBound.back() << " ~ ";
  else
   *f_log << "UB = " << - LowerBound.back() << " ~ ";
  }

 *f_log << "Fi1 = ";

 if( f_convex ) {
  if( UpFiLmb1.back() == -INFshift ) {
   *f_log << "-INF => STOP." << std::endl;
   return;
   }
  else
   if( UpFiLmb1.back() == INFshift ) {
    *f_log << "+INF" << std::endl;
    return;
    }
   else
    *f_log << UpFiLmb1.back() << shrt;
  }
 else
  if( UpFiLmb1.back() == -INFshift ) {
   *f_log << "+INF => STOP." << std::endl;
   return;
   }
  else
   if( UpFiLmb1.back() == INFshift ) {
    *f_log << "-INF" << std::endl;
    return;
    }
   else
    *f_log << - UpFiLmb1.back() << shrt;

 if( NeedsAlfa1() )
  *f_log << " ~ Alfa1 = " << Alfa1;

 if( NeedsScPr1() )
  *f_log << " ~ Gi1xd = " << ScPr1;

 *f_log << std::endl;

 }  // end( BundleSolver::Log2 )

/*--------------------------------------------------------------------------*/

void BundleSolver::compute_NrmZFctr( void )
{
 auto wf = ( WZNorm << 2 );
 // if we need to sum but some component has no linearization, return:
 // NrmZFctr remains undefined
 if( wf > 1 )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ( ( ! NrEasy ) || ( ! IsEasy[ k ] ) ) && ( NrItems[ k ] == 0 ) )
    return;

 // now we sum: note that if ! f_convex one should change the sign, but
 // the norm is invariant w.r.t. the sign
 std::vector< VarValue > tg( NumVar , 0 );
 if( f_lf ) {      // the linear 0-th component is there
  auto & cf = f_lf->get_v_var();
  for( Index i = 0 ; i < NumVar ; ++i )
   tg[ i ] = cf[ i ].second;
  }
 else              // there is no 0-th component
  if( wf <= 1 ) {  // and we just wanted is subgradient
   NrmZFctr = 1;   // ... which is all-0, so use 1
   return;
   }

 if( wf > 1 ) {  // also need to sum a subgradient for each component
  std::vector< VarValue > tg1( NumVar );

  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( NrEasy && IsEasy[ k ] )
    continue;

   Index i = 0;
   while( InvItemVcblr[ k ][ i ] == InINF )
    ++i;

   v_c05f[ k ]->get_linearization_coefficients( tg1.data() ,
                                                Range( 0 , NumVar ) , i );

   std::transform( tg.begin() , tg.end() , tg1.begin() , tg.begin() ,
                   std::plus< VarValue >() );
   }
  }

 NrmZFctr = ::norm( tg , WZNorm & 3 );

 }  // end( compute_NrmZFctr )

/*--------------------------------------------------------------------------*/

bool BundleSolver::FindNext( void )
{
 Index InitwFi = f_wFi;
 do {
  f_wFi = ( f_wFi + 1 ) % NrFi;    // next patient, please
  if( NrEasy && IsEasy[ f_wFi ] )  // skip easy components
   continue;
  // A fictitious lower bound only keeps the master well-defined while a
  // hard component has no cuts; it is not a substitute for its model.
  // Force that component to be evaluated again so a real linearization can
  // refill the bundle.
  if( ( NrItems[ f_wFi ] == 0 ) &&
      ( CurrNrEvls[ f_wFi ] < MaxNrEvls ) )
   return( true );
  if( ( FiStatus[ f_wFi ] == kUnEval ) ||
      ( ( FiStatus[ f_wFi ] < kError ) && ( FiStatus[ f_wFi ] > kOK ) &&
        ( CurrNrEvls[ f_wFi ] < MaxNrEvls ) ) )
   return( true );

  } while( f_wFi != InitwFi );

 return( false );

 }  // end( BundleSolver::FindNext )

/*--------------------------------------------------------------------------*/
/* Front-end for four different heuristics for short-term t management.
 *
 * Three of the heuristics (1, 2, and 3) are based on three slightly different
 * variants of the same idea: considering the objective f( x ) from the
 * current stability center \bar{x} along direction d* seen as a function of
 * t, assuming that d* = - t z* (which is, strictly speaking, only true in
 * the pure quadratic proximal case and therefore does not cleanly generalise
 * to the generalised one; yet these are heuristics).
 *
 * That is, we consider the translated function along z*
 *
 *    q( v ) = f( \bar{x} - v z* ) - f( \bar{x} )
 *
 * Assuming (only for notational simplicity) differentiability, we thus have
 *
 *    q'( v ) = < - z* , f'( \bar{x} - v z* ) >
 *
 * After that f( \bar{x} + d* ) = f( \bar{x} - t z* ) = q( t ) has been
 * computed, we know:
 *
 * - the aggregated subgradient z*, which is a Sigma*-subgradient in \bar{x}
 *
 * - the newly obtained subgradient g, which is an Alfa1-subgradient in 0 and
 *   eps-subgradient in t, where
 *
 *     eps = DeltaFi - ( Alfa1 + < g , d* > )
 *
 * We thus assume:
 *
 * - q( 0 ) = 0
 *
 * - q'( 0 ) = < - z* , z* > = - NrmZ^2
 *
 * Note, however, that z* is a Sigma*-subgradient in \bar{x}, and therefore
 * the value of the linearization there is rather -Sigma*; thus, we could
 * alternatively assume q( 0 ) = - Sigma*.
 *
 * - q( t ) = f( \bar{x} - t z* ) - f( \bar{x} ) = DeltaFi
 *
 * - q'( t ) = < - z* , g >; since we have ScPr1 = < d* , g > =
 *   < - t z* , g > (note again that this only holds in the quadratic case),
 *   we conclude q'( t ) = ScPr1 / t
 *
 * Note, however, that g* is a Alfa1-subgradient of \bar{x}, and therefore we
 * could alternatively take the value in t as that of the corresponding
 * linearization, i.e., q( t ) = ScPr1 - Alfa1.
 *
 * We can then consider the quadratic function
 *
 *    m( v ) = a v^2 + b v + c
 *
 * and construct different forms of it corresponding to different choices of
 * three of the four information we have, then use its minimum
 *
 *   v* = - b / ( 2 a )
 *
 * as the suggested new value for t. Since v* only makes sense if a > 0,
 * when a <= 0 we use the best possible convex approximation of a concave
 * function by setting a = 0, in which case the minimum is the extreme of
 * the interval [ tMinor , tMaior ] dictated by the sign of b.
 *
 * The fourth heuristic is based on an entirely different idea related to the
 * Moreau-Yoshida regularization, called "reversal form of the poorman's
 * quasi-Newton update". */

double BundleSolver::Heuristic( Index whch )
{
 switch( whch & 3 ) {
  case( 0 ): return( Heuristic1() );
  case( 1 ): return( Heuristic2() );
  case( 2 ): return( Heuristic3() );
  }

 return( Heuristic4() );
 }

/*--------------------------------------------------------------------------*/
/* With the notation above, in the first case we impose
 *
 *    m( 0 )  = c = - Sigma*
 *    m( t )  = a t^2 + b t + c = DeltaFi
 *    m'( 0 ) = [ 2 a 0 ] + b = - NrmZ^2
 *
 * which yields c = - Sigma*, b = - NrmZ^2,
 * a = ( DeltaFi + NrmZ^2 t + Sigma*  ) / t^2.
 * Note that z* is a Sigma*-subgradient in \bar{x}, and therefore
 *
 *    f( \bar{x} + d* ) >= f( \bar{x} ) + < d* , z* > - Sigma*
 *                       = f( \bar{x} ) - t < z* , z* > - Sigma*
 *    ==> DeltaFi = f( \bar{x} + d* ) - f( \bar{x} ) >= - NrmZ^2 t - Sigma*
 *    ==> a = DeltaFi + NrmZ^2 t + Sigma* >= 0
 *
 * which guarantees that m() is convex and therefore the minimum of m() is
 *
 *   v* = - ( - NrmZ^2 ) / ( 2 ( DeltaFi + NrmZ^2 t + Sigma* ) / t^2 )
 *      =   NrmZ^2 / ( 2 ( DeltaFi + NrmZ^2 t + Sigma* ) / t^2 )
 *      =   t^2 NrmZ^2 / ( 2 ( DeltaFi + NrmZ^2 t + Sigma* ) )
 *
 * Note that, conveniently, v* >= 0 always holds. This corresponds to the
 * fact that m'( 0 ) = b < 0, i.e., m() is surely decreasing in 0.
 *
 * However, this formula has a serious issue: we need to know DeltaFi,
 * which may well not be defined when a NS is performed and multiple
 * components are present since the incremental approach may stop the
 * inner loop before having computed them all. The obvious solution is
 * to replace DeltaFi with
 *
 *    \underline{f}( \bar{x} + d* ) - \bar{f}( \bar{x} ) =
 *    LwFiLmb1.back() - UpRifFi.back()
 *
 * which is always well-defined since a finite lower bound is always
 * available (unless some component evaluates to -INF, in which case
 * the algorithm stops and this method is not invoked). */

double BundleSolver::Heuristic1( void )
{
 auto DF = DeltaFi < INFshift ? DeltaFi : LwFiLmb1.back() - UpRifFi.back();
 auto NZ2 = NrmZ * NrmZ;
 if( DF + NZ2 * t + Sigma > 1e-16 )  // this should always be >=
  return( t * t * NZ2 / ( 2 * ( DF + NZ2 * t + Sigma ) ) );
 else                 // a == 0, all depends on the sign of b
  /* there is no "if" here, NrmZ >= 0 by definition, a fortiori NZ2
  if( - NZ2 <= 0 )    // b < 0  */
   return( tMaior );  // ==> tMaior
  /* there is no "else" here, - NZ2 > 0 cannot happen
  else                // b > 0
   return( tMinor );  // ==> tMinor */
 }

/*--------------------------------------------------------------------------*/
/* With the notation above, in the second case we rather impose
 *
 *    m( 0 )  = c = 0
 *    m( t )  = a t^2 + b t [ + 0 ] = ScPr1 - Alfa1
 *    m'( t ) = 2 a t + b = ScPr1 / t
 *
 * which yields c = 0, a = Alfa1 / t^2, b = ( ScPr1 - 2 Alfa1 ) / t
 *
 * Since Alfa1 >= 0, a >= 0 which implies that m() is surely convex and the
 * minimum is
 *
 *   v* = - [ ( ScPr1 - 2 Alfa1 ) / t ] / [ 2 Alfa1 / t^2 ]
 *      = t ( 2 Alfa1 - ScPr1 ) / ( 2 Alfa1 )
 *
 * Note that, unlike in the first case, there is no guarantee that v* >= 0,
 * because we fix the derivative in t and therefore m'( 0 ) = b may turn up
 * to be positive (m() in increasing in 0).
 *
 * Since this formula does not really use DeltaFi, it being undefined is not
 * an issue here. However, this formula has a somewhat similar issue with NS
 * (and SS alike) in that not all components may have been evaluated, and
 * therefore only a "partial" g may be available. ScPr1 and Alfa1 are
 * computed for all non-"easy" components using the "representative
 * subgradients" out of the previous iterations, *provided they have not by
 * chance been deleted* (which should not happen unless the bundle is very
 * very small). Yet, "easy" components are left out. There may be some way
 * put of this, e.g. by using z*_i in place of g_i for the "easy" components,
 * but this is nontrivial and therefore avoided for now. */

double BundleSolver::Heuristic2( void )
{
 if( Alfa1 > 1e-16 )            // it is always >= 0, but it may be ==
  return( t * ( 2 * Alfa1 - ScPr1 ) / ( 2 * Alfa1 ) );
 else                           // a == 0,  all depends on the sign of b
  if( ScPr1 - 2 * Alfa1 <= 0 )  // b < 0
   return( tMaior );            // ==> tMaior
  else                          // b > 0
   return( tMinor );            // ==> tMinor
 }

/*--------------------------------------------------------------------------*/
/* With the notation above, in the third case we rather impose
 *
 *    m'( 0 ) = [ 2 a 0 ] + b = - NrmZ^2
 *    m( t )  = a t^2 + b t + c = < something >
 *    m'( t ) = 2 a t + b = ScPr1 / t
 *
 * which yields b = - NrmZ^2, a  = ( ScPr1 / t + NrmZ^2 ) / ( 2 t ), and c
 * ... something depending on which value we choose for m( t ), but we
 * don't care about because c does not appear in the computation of v*.
 * If m() is convex, i.e.
 *
 *    ScPr1 / t + NrmZ^2 > 0
 *
 * yields
 *
 *   v* = - [ - NrmZ^2 ] / [ 2 ( ScPr1 / t + NrmZ^2 ) / ( 2 t ) ]
 *      = t NrmZ^2 / ( ScPr1 / t + NrmZ^2 )
 *
 * Note that, if m() is convex, then v* >= 0 holds because, as usual, we have
 * fixed m'( 0 ) = b < 0 and therefore m() is decreasing in 0.
 *
 * See above for the "issue" about ScPr1 having been computed with a
 * "partial" g; however, since this formula does not really use DeltaFi,
 * it being undefined is not an issue here.
 *
 * Note also that the possible fourth case
 *
 *    m( 0 )  = c = 0 [ or - Sigma* ]
 *    m'( 0 ) = [ 2 a 0 ] + b = - NrmZ^2
 *    m'( t ) = 2 a t + b = ScPr1 / t
 *
 * only changes c w.r.t. the current one, hence it does not change v*, and
 * therefore need not be separately considered. */

double BundleSolver::Heuristic3( void )
{
 auto NZ2 = NrmZ * NrmZ;
 if( ScPr1 / t + NZ2 > 1e-16 )
  return( t * NZ2 / ( ScPr1 / t + NZ2 ) );
 else                 // a == 0, all depends on the sign of b
  /* there is no "if" here, NrmZ >= 0 by definition, a fortiori NZ2
  if( - NZ2 <= 0 )    // b < 0  */
   return( tMaior );  // ==> tMaior
  /* there is no "else" here, - NZ2 > 0 cannot happen
  else                // b > 0
   return( tMinor );  // ==> tMinor */
 }

/*--------------------------------------------------------------------------*/
/* This heuristic is instead based on a completely different approach. It is
 * called "reversal form of the poorman's quasi-Newton update" and its
 * nontrivial rationale is described in details in
 *
 *  C. Lemarechal and C. Sagastizabal. Variable metric bundle methods: from
 *  conceptual to implementable forms. Mathematical Programming,
 *  76(3):393-410, 1997
 *
 * A more refined version of the same is proposed in
 *
 *  P.A. Rey and C. Sagastizabal. Dynamical adjustment of the prox-parameter
 *  in variable metric bundle methods. Optimization, 51(2):423-447, 2002
 *
 * It should be noted that this heuristic is explicitly developed for being
 * used at SS only.
 *
 * The proposed new value is
 *
 *   t = < v , u > / || v ||^2
 *
 * where
 *
 *   v = g - z*
 *
 *   u = ( \bar{x} - v z* ) - \bar{x} + t v = d* + t v
 *
 * although v would in general be g_{i+1} - g_i, hence the choice of z* as
 * g_i is somewhat arbitrary; but in general z* is considered to be "the best
 * (approximate) subgradient we have at \bar{x}".
 *
 * Hence
 *
 *   v = < v , d* + t v > / || v ||^2
 *     = [ < v , d* > + t < v , v > ] / || v ||^2
 *     = < g - z* , d* > / || g - z* ||^2 + t
 *     = t + [ < g , d* > + t || z* ||^2 ] /
 *           [ || g ||^2  - 2 < g , z* > + || z* ||^2 ]
 *     = t + [ < g , d* > + t || z* ||^2 ] /
 *           [ || g ||^2  + 2 < g , d* > / t + || z* ||^2 ]
 *
 * The issue with this formula is the || g || term. This is the same issue as
 * with ScPr1 = < g , d* >, i.e., what to do with the easy components which
 * do not explicitly compute a subgradient. For the scalar product we could
 * use < z*_i , d* > that should be available "for free" out of the Master
 * Problem but it currently isn't; in theory z*_i is also available, and we
 * could use it to compute < z*_i , d* >, but reading per-component z*
 * from MasterPB is not yet exposed as a cheap query. */

double BundleSolver::Heuristic4( void )
{
 auto NZ2 = NrmZ * NrmZ;
 if( G1Norm == INFshift )
  G1Norm = norm( G1 , 2 );
 return( t + ( ScPr1 + t * NZ2 ) / ( G1Norm + 2 * ScPr1 / t + NZ2 ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::CreateMPB( void )
{
 // allocate (or re-cycle) the MasterProblemBlock 
 if( MasterPB )
  MasterPB->clear();
 else
  MasterPB = new MasterProblemBlock();

 }  // end( BundleSolver::CreateMPB )

/*--------------------------------------------------------------------------*/

void BundleSolver::InitMPB( void )
{
 // configure() constructs the primal variables and Objective, so select the
 // storage frame before calling it. In primal raw mode the optimization
 // variable is the absolute point x rather than the displacement d.
 MasterPB->set_v2_form( MPV2Form );

 // pick primal vs dual: easy components force the dual form
 const bool want_primal = IsMPPrimal && ! ( DoEasy && ( NrEasy > 0 ) );
 if( IsMPPrimal && ! want_primal )
  BLOG( 1 , std::endl << "Warning: intMPPrimal is set to 1 but DoEasy is "
        "non-zero and NrEasy = " << NrEasy << "; the dual MP will be "
        "used anyway." );

 // one-shot MasterProblemBlock configuration: sizes, primal/dual form,
 // hard-component count and easy components are handed over in a single
 // call. configure() dispatches on the concrete type of each easy
 // C05Function and grafts the resulting sub-Block under MasterPB, so no
 // separate per-component wiring is needed.
 //
 // Any exclusion list installed on *this* via Solver::set_excluded_blocks()
 // is propagated to MasterPB so its inner Solver skips the same subtrees.
 //
 // build the easy_components vector by walking v_c05f and keeping only
 // the components flagged as easy (IsEasy is sized only when DoEasy != 0,
 // so the trivial "all hard" case yields an empty vector). Here we also
 // store the local2global map for easy components, as this may be needed
 // by MPB to build the coupling constraints.
 std::vector< C05Function * > easy_cmps;
 std::vector< std::vector< Index > > easy_local2global;
 if( IsEasy.size() == NrFi )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( IsEasy[ k ] ){
    easy_cmps.push_back( v_c05f[ k ] );

    if( f_sparse_lambda )
     easy_local2global.push_back( v_local2global[ k ] );
   }

 const int n_hard = int( NrFi ) - int( easy_cmps.size() );

 // The MasterPB MP sense is dictated by the easy sub-Block grafted by
 // configure(), whose Objective sense is the OPPOSITE of f_convex
 // (LagBFunction::IsConvex is set to inner_sense == eMax, i.e. the inner
 // primal of a convex LBF is eMax and of a concave LBF is eMin). To
 // match, MasterPB receives !f_convex.
 //
 // MaxBSize passed to MasterPB is the *global* pool size, not the
 // per-component cap BPar2. BundleSolver picks slot names through
 // FindAPlace() from a single global pool (with size vBPar2.back() =
 // sum_k vBPar2[ k ]), so add_cut() sees slot indices that can reach
 // up to that total even though each component never owns more than
 // vBPar2[ k ] = BPar2 of them simultaneously. Routing BPar2 (or even
 // max_k vBPar2[ k ]) to MasterPB would under-allocate slot_to_local[]
 // and trip "slot out of range" the first time a slot >= BPar2 is
 // legitimately allocated by FindAPlace
 const int max_bsize = vBPar2.empty()
                       ? int( BPar2 )
                       : int( vBPar2.back() );
 MasterPB->configure( want_primal ,
                      std::max( int( BPar2 ) , max_bsize ) ,
                      int( NumVar ) ,
                      n_hard ,
                      easy_cmps ,
                      get_excluded_blocks() ,
                      MPStbl ,
                      ! f_convex ,
                      IsEasy ,
                      MPHScaling ,
                      easy_local2global );

 // The storage frame has already been selected before configure(): 0 is the
 // translated/displacement form, 1 the raw/iterate form.
 // optional "lazy reference" for the displacement form: defer the per-cut
 // g . ( x_bar - x_ref ) shift to the lin-z and re-align x_ref only when the
 // centre drifts past the tolerance. 0 ( the default ) keeps the plain
 // displacement frame; > 0 is opt-in via the environment.
 if( const char * env = std::getenv( "BUNDLE_MPB_XREFTOL" ) ;
     env && env[ 0 ] )
  MasterPB->set_xref_tol( std::atof( env ) );

 tHasChgd = true;

 }  // end( BundleSolver::InitMPB )

/*--------------------------------------------------------------------------*/

Index BundleSolver::BStrategy( Index wFi )
{
 // this method implements the B-strategies of the code, i.e., which "old"
 // items are discarded if the bundle is full and a new item belonging to
 // component wFi has to be inserted
 //
 // this is called *before* that we know if the place will actually be
 // required, so it has a "loose attitude"; in particular, it will return
 // InINF under two opposite set of conditions:
 //
 // - there is plenty of space left in the bundle, so that no B-strategy (no
 //   removal or aggregation) is required;
 //
 // - there is no way in which space can be found, i.e., the bundle (for this
 //   component) is full, and removal/aggregation are not successful (a very
 //   strange occurrence due to an extremely small bundle or it being chock
 //   full of constraints)
 //
 // Picking a specific spot in the free space is the task of FindAPlace(),
 // which however is not called right away because the place may end up not
 // being needed. If BStrategy() return InINF because there is plenty of
 // space then FindAPlace() will suceed, if BStrategy() return InINF because
 // there is no way space can be found then FindAPlace() will fail and it
 // will be clear that disaster looms
 //
 // important note: if the returned wFi is not InINF, *and* it is the name
 // of an item still in the bundle, then:
 //
 // - the item belongs to the component wFi, so as to keep the number of
 //   items of that component constant
 //
 // - the item is *not* deleted, since it is not 100% sure this will need
 //   to be done, so deletion will be responsibility of the caller

 // there are "free" items in the global pool: there is "plenty of space"
 if( FrFItem[ wFi ] < vBPar2[ wFi ] )
  return( InINF );

 // there is not plenty of space, take 1- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // there are no "free" items in the gloal pool, but there are items in the
 // global pool that are not in the bundle; these will go first

 if( NrItems[ wFi ] < vBPar2[ wFi ] )
  return( InINF );

 // there is not plenty of space, take 2- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the bundle for component wFi is full, which implies that the global pool
 // is also full: among the items in the bundle for component wFi, find the
 // removable one with largest OOBase[], snd among these with the largest
 // OOBase[] select that with largest Alfa[]
 // note: the Z[ wFi ] in the bundle (if any) is not removable and
 //       therefore cannot be selected, which in particular happens if
 //       wFi has only *one* subgradient in base

 Index wh = InINF;
 SIndex OOwh = -Inf< SIndex >();
 double Awh = -Inf< double >();
 for( auto i : InvItemVcblr[ wFi ] ) {
  assert( i < vBPar2.back() );
  const auto ai = read_alpha_global( i );
  if( ( OOBase[ i ] > OOwh ) ||
      ( ( OOBase[ i ] == OOwh ) && ( ai > Awh ) ) ) {
   wh = i;
   OOwh = OOBase[ i ];
   Awh = ai;
   }
  }

 if( wh == InINF )         // InvItemVcblr[ wFi ] was empty: nothing to
  return( InINF );         // pick, signal the caller
 if( OOBase[ wh ] < 0 )    // all items are non-removable: nothing else to
  return( InINF );         // do (except maybe complaining very loudly)
 else
  if( OOBase[ wh ] > 0 )   // a place is found
   return( wh );           // there are no problems, all done

 // wh is a basic item- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this means that *all* items of component wFi are either in base or not
 // removable, for otherwise we would have selected an item with OOBase > 0;
 // we cannot discard anything before having performed aggregation, but in
 // order to do so we also need to free some space for the Z[ wFi ]
 //
 // note: there is an easy case where z is "naturally" in the base without
 //       any aggregation: there is only one subgradient in base (for the
 //       component wFi), and therefore, its Mlt[] is == 1. however this
 //       can never happen here because it would mean vBPar2[ wFi ] == 1
 //       which is not permitted. the other case in which this could happen
 //       is that all the items in the bundle for component wF are
 //       constraints, but then aggregation would not be useful (in fact
 //       this is checked and reported as a failure)
 //
 // note: this also means that wh is the item in base with largest Alpha
 //
 // note: that since *all* items of this component are either in base or
 //       not removable, we can scan MBse[] for the items to be removed,
 //       possibly ignoring items with Mlt[] == 0 -- but in fact not doing
 //       it because there is not any

 if( Zvalid[ wFi ] )  // a valid Z[ wFi ] is in the bundle: it is safe to
  return( wh );       // replace wh with anything the oracle provides us

 // a valid Z[ wFi ] is not already in: aggregation has to be performed.
 // Walk InvItemVcblr[ wFi ] (the per-component slot list) collecting
 // ( name , theta ) pairs from MasterPB; the aggregation coefficients
 // are the optimal multipliers theta^k_i of the bundle B^k.

 if( ! MasterPB )
  return( wh );  // no master to query: fall back to the slot picked above

 const bool pure_level_aggregate = MasterPB->uses_pure_level_aggregation();
 const double aggregate_mass = pure_level_aggregate
                               ? MasterPB->get_level_multiplier()
                               : MasterPB->get_lambda();
 const double aggregate_mass_eps =
  1e-12 * std::max( { std::abs( aggregate_mass ) , double( 1 ) } );
 if( aggregate_mass <= aggregate_mass_eps )
  return( InINF );

 LinearCombination coeff;
 coeff.reserve( InvItemVcblr[ wFi ].size() );
 for( Index slot = 0 ; slot < InvItemVcblr[ wFi ].size() ; ++slot ) {
  const auto name = InvItemVcblr[ wFi ][ slot ];
  if( name >= vBPar2.back() )
   continue;  // empty slot in the per-cmp pool
  const auto th = MasterPB->get_theta( int( wFi ) , int( slot ) );
  if( th == 0 )
   continue;
  coeff.emplace_back( slot , th / aggregate_mass );
  }

 Index whZ = InINF;  // the position where Z[ wFi ] has to go
 if( ( whisZ[ wFi ] < InINF ) && is_subgradient_global( whisZ[ wFi ] ) )
  whZ = whisZ[ wFi ];  // preferably re-use the last position
 else
  // pick the item in base with min theta different from wh
  for( const auto & p : coeff )
   if( p.first != ItemVcblr[ wh ].second &&
       OOBase[ InvItemVcblr[ wFi ][ p.first ] ] >= 0 ) {
    whZ = InvItemVcblr[ wFi ][ p.first ];
    break;
    }

 if( whZ == InINF )  // there is no removable item apart from wh
  return( InINF );   // nothing else to do except complaining very loudly

 // tell the C05Function what is going to happen
 // note that this only happens when the bundle (for component wFi) is
 // "very full", and therefore also the global pool (for component wFi)
 // is such. The natural choice is to put the new aggregate linearization
 // in the same position in the global pool where whZ was.
 inhibit_Modification( true );
 v_c05f[ wFi ]->store_combination_of_linearizations( coeff ,
                                                   ItemVcblr[ whZ ].second );
 inhibit_Modification( false );

 remove_cut_global( whZ );  // remove the old item in position whZ

 // materialise the V2 aggregate cut:
 //   g* = (1/mass) sum_i theta_i g_i
 //   a* = (1/mass) (sum_i theta_i a_i + gamma LB)
 // In pure level, get_aggregated_subgradient() already returns g*.
 std::vector< double > tZ_buf =
                              MasterPB->get_aggregated_subgradient( wFi );
 if( tZ_buf.empty() )
  tZ_buf.assign( NumVar , 0.0 );
 if( ! pure_level_aggregate )
  for( auto & v : tZ_buf )
   v /= aggregate_mass;

 const double Ai =
  MasterPB->get_raw_aggregated_alpha_with_LB( wFi ) / aggregate_mass;

 // First remove any pre-existing cut at slot whZ, then re-add the
 // aggregated one
 MasterPB->remove_cut( int( wFi ) , int( whZ ) );
 MasterPB->add_cut( int( wFi ) , int( whZ ) , std::move( tZ_buf ) , Ai );

 whisZ[ wFi ] = whZ;      // Z[ wFi ] is in the bundle in position whZ
 Zvalid[ wFi ] = true;    // ... and it is valid
 OOBase[ whZ ] = -1;      // ... and it won't be removed in this iteration

 BLOG( 2 , std::endl << "Aggregation performed into " << whZ );

 // at this point, Z[ wFi ] is in the bundle, hence it is safe to replace wh
 // with anything the oracle provides us

 return( wh );

 }  // end( BundleSolver::BStrategy )

/*--------------------------------------------------------------------------*/

Index BundleSolver::FindAPlace( Index wFi )
{
 // this method is used to return the index of an available position in the
 // bundle where to store a new item belonging to "component" wFi; if there
 // are no possible positions left, then InINF is returned
 //
 // note that positions in the bundle are not "reserved to components", so
 // the task of FindAPlace() is easy. it is the business of other parts of
 // the code to ensure that a component does not use more than its fair
 // share of positions in the bundle

 Index wh = InINF;

 if( ! FreList.empty() ) {       // there are deleted items
  wh = FreList.top();            // pick the one with smaller name
  if( wh >= get_max_name() )  // FreList has all items with "large" names
   FreList = {};                 // clear FreList
  else {                         // wh is a "small" name
   FreList.pop();                // take it away
   return( wh );                 // all done
   }
  }

 // if there are no deleted items (with suitably small names)
 if( get_max_name() < vBPar2.back() )  // ... but there is still space
  wh = get_max_name();                 // next name

 return( wh );

 }  // end( BundleSolver::FindAPlace )

/*--------------------------------------------------------------------------*/

bool BundleSolver::NeedsAlfa1( void )
{
 // Alfa1 is only used in Heuristic2
 return( ( ( ( tSPar1 & 1 ) && ( ( tSPar1 & tSPHMsk1 ) == 64 ) ) ) ||
         ( ( ( tSPar1 & 2 ) && ( ( tSPar1 & tSPHMsk2 ) == 256 ) ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

bool BundleSolver::NeedsScPr1( void )
{
 // ScPr1 is used by everyone save for Heuristic1
 return( ( ( ( tSPar1 & 1 ) && ( tSPar1 & tSPHMsk1 ) ) ) ||
         ( ( ( tSPar1 & 2 ) && ( tSPar1 & tSPHMsk2 ) ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

bool BundleSolver::NeedsG1( void )
{
 // G1 is only used in Heuristic4
 return( ( ( ( tSPar1 & 1 ) && ( ( tSPar1 & tSPHMsk1 ) == 172 ) ) ) ||
         ( ( ( tSPar1 & 2 ) && ( ( tSPar1 & tSPHMsk2 ) == 768 ) ) ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdateHeuristicInfo( void )
{
 // if required, update the "aggregated" Alfa1 and ScPr1, that are used in
 // the t heuristics, using the "representatives" of all components that
 // have not been compute()-d in the latest InnerLoop()

 if( NeedsAlfa1() && NeedsG1() ) {
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ( ! CurrNrEvls[ k ] ) && ( whisG1[ k ] < InINF ) ) {
    Alfa1 += read_alpha_global( whisG1[ k ] );
    ScPr1 += read_Gid_global( whisG1[ k ] );
    }

  return;
  }

 if( NeedsAlfa1() ) {
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ( ! CurrNrEvls[ k ] ) && ( whisG1[ k ] < InINF ) )
    Alfa1 += read_alpha_global( whisG1[ k ] );

  return;
  }

 if( NeedsG1() )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ( ! CurrNrEvls[ k ] ) && ( whisG1[ k ] < InINF ) )
    ScPr1 += read_Gid_global( whisG1[ k ] );

 }  // end( UpdateHeuristicInfo )

/*--------------------------------------------------------------------------*/

void BundleSolver::guts_of_destructor( void )
{
 whisG1.clear();
 vStar.clear();

 LowerBound.clear();
 LwFiLmb.clear();
 UpFiLmb.clear();
 LwFiLmb1.clear();
 UpFiLmb1.clear();
 UpRifFi.clear();

 FiStatus.clear();

 CurrNrEvls.clear();

 Zvalid.clear();
 whisZ.clear();
 FreList = {};

 MaxItem.clear();
 FrFItem.clear();
 NrItems.clear();

 ItemVcblr.clear();

 OOBase.clear();

 LmbdBst.clear();
 Lambda1.clear();
 Lambda.clear();

 InvItemVcblr.clear();
 vBPar2.clear();
 f_max_name = 0;

 if( NrEasy ) {
  // the easy-component sub-Blocks remain owned by their Function Blocks;
  // MasterPB only keeps non-owning registrations while it is configured
  IsEasy.clear();
  NrEasy = 0;
  }

 LamVcblr.clear();

 v_local2global.clear();
 v_local2global.shrink_to_fit();
 Lambda2Idx.clear();
 v_ref_count.clear();
 v_ref_count.shrink_to_fit();
 f_sparse_lambda = false;

 v_c05f.clear();

 }  // end( BundleSolver:guts_of_destructor )

/*--------------------------------------------------------------------------*/

void BundleSolver::ReSetAlg( unsigned char RstLvl )
{
 if( ! ( RstLvl & RstAlg ) ) {  // reset algorithmic parameters - - - - - - -
  ParIter = 0;             // reset iterations count
  CSSCntr = CNSCntr = 0;   // ... comprised consecutive NS/SS count

  if( t != tInit ) {       // reset t
   t = tInit;
   tHasChgd = true;
   }

  // reset the dynamic number of fetched items
  if( BPar6 && ( BPar5 > 0 ) )
   aBP3 = BPar4;
  else
   aBP3 = BPar3;
  }

 // reset function values in Lambda, since they are changing
 std::fill( UpFiLmb.begin() , UpFiLmb.end() ,  INFshift );  // upper
 std::fill( LwFiLmb.begin() , LwFiLmb.end() , -INFshift );  // lower
 UpFiLmbdef = LwFiLmbdef = 0;             // ... at the current point
 f_global_LB = -INFshift;         // algorithmic global LB

 if( RstLvl & RstCrr ) {  // get an initial point - - - - - - - - - - - - - -
  // note that the thusly constructed Master Problem assumes the stability
  // centre to be all-0, in particular when constructing the Lagrangian costs
  // of the easy components, so if that's not happening the right value has
  // to be explicitly passed

  for( Index i = 0 ; i < NumVar ; ++i )
   Lambda[ i ] = LamVcblr[ i ]->get_value();

  bool nonzero = false;
  for( Index i = 0 ; i < NumVar ; ++i )
   if( Lambda[ i ] ) { nonzero = true; break; }

  if( nonzero ) {
   // shift the stability centre to Lambda: oldLambda == 0 by
   // construction. UpRifFi is still all-0 at this point (no function
   // has been evaluated yet), so the per-component F_k( x_bar ) cache
   // installed inside MasterPB starts at 0; subsequent set_reference()
   // calls after the first SS will refresh both x_bar and F_k( x_bar )
   // atomically
   if( MasterPB ) {
    std::vector< double > F_hard;
    F_hard.reserve( NrFi );
    for( Index k = 0 ; k < NrFi ; ++k ) {
     if( NrEasy && IsEasy[ k ] )
      continue;
     F_hard.push_back( UpRifFi[ k ] );
     }
    MasterPB->set_reference( Lambda , F_hard );
    }
   Fi0Lmb = INFshift;  // the value of the linear part must be computed
   }
  else  // Lambda was all-0 anyway
   Fi0Lmb = 0;  // then the value of the linear part is quite obvious ...
  }
 else {                   // reset the current point to all-0 - - - - - - - -
  // note that the thusly constructed Master Problem precisely assumes
  // the stability centre to be all-0, in particular when constructing the
  // Lagrangian costs of the easy components, so that's fine as it is
  Lambda.assign( NumVar , 0 );
  // "tell" this to the ColVariable of the C05Function(s)
  for( Index i = 0 ; i < NumVar ; ++i )
   LamVcblr[ i++ ]->set_value( 0 );
  Fi0Lmb = 0;  // then the value of the linear part is quite obvious ...

  // Seed / refresh the master reference cache. Lambda has just been
  // reset to 0 so x_bar moves to 0 too; the per-component F_k( x_bar )
  // values, on the other hand, are *not* known at this point (Fi(.)
  // has not been re-evaluated since the previous compute()). To keep
  // the cuts still living in the bundle (from a previous compute())
  // consistent with the new x_bar, we re-use whatever F values
  // MasterPB has cached: this yields a pure x_bar shift (delta =
  // A . ( x_new - x_old )) on every stored b[ i ], with no spurious
  // dF jolt that would otherwise corrupt their lin-error meaning.
  // The "true" F_k( 0 ) will be installed at the next SS / NS via
  // GotoLambda1 / GotoLambda's set_reference call once Fi(.) has been
  // evaluated. At the very first ever compute() call the bundle is
  // empty and the cached values default to 0, so the same code path
  // also serves as initial seeding
  if( MasterPB ) {
   std::vector< double > F_hard;
   F_hard.reserve( NrFi );
   int hard_idx = 0;
   for( Index k = 0 ; k < NrFi ; ++k ) {
    if( NrEasy && IsEasy[ k ] )
     continue;
    F_hard.push_back( MasterPB->get_F_at_x_bar( hard_idx++ ) );
    }
   MasterPB->set_reference( Lambda , F_hard );
   }
  }
 }  // end( BundleSolver::ReSetAlg )

/*--------------------------------------------------------------------------*/

void BundleSolver::Delete( Index i , bool ModDelete )
{
 // deletes from the bundle the item in position i
 //
 // ModDelete == true means that this is called in response of a Modification
 // where the linearization has been removed from the global pool
 //
 // whatever BPar7 says, no linearization is physically deleted here inside,
 // this has to be done by the caller (if needed)

 Index k = ItemVcblr[ i ].first;

 // check if this item was the "representative" for its component - - - - - -

 if( whisG1[ k ] == i )  // it is the representative of k
  whisG1[ k ] = InINF;   // a new representative is needed

 // check if this item was the z* for its component - - - - - - - - - - - - -

 if( whisZ[ k ] == i ) {  // it is the aggregate subgradient of k
  whisZ[ k ] = InINF;     // no aggregate subgradient is in the bundle
  Zvalid[ k ] = false;    // a fortiori, no valid one
  }

 // delete the item from the MP - - - - - - - - - - - - - - - - - - - - - - -

 remove_cut_global( i );

 BLOG( 2 , std::endl << "Item " << i << " removed" );

 // bookkeeping of internal data structures - - - - - - - - - - - - - - - - -
 // note that any item whose name is >= get_max_name() is surely not in
 // the bundle (master problem), and therefore it need not be in FreList

 Index MxNm = get_max_name();
 if( i < MxNm )
  FreList.push( i );

 OOBase[ i ] = Inf< SIndex >();
 --NrItems[ k ];
 --NrItems[ NrFi ];

 // FormD() may delete an item as part of its emergency recovery and
 // immediately re-solve the master, without returning to the outer loop
 // where fictitious bounds are normally synchronized. If this was the last
 // cut of a hard component with no genuine lower bound, install its
 // fictitious bound now so the component normalization remains feasible.
 if( MasterPB && ( NrItems[ k ] == 0 ) &&
     ( ( ! NrEasy ) || ( ! IsEasy[ k ] ) ) &&
     ( LowerBound[ k ] <= - INFshift ) && ( ! FictLB[ k ] ) ) {
  MasterPB->set_fictitious_LB( int( k ) , true );
  FictLB[ k ] = true;
  }

 // remove from the global pool: the removal is "hard" if either BPar7 says
 // so, or the linearization had been deleted anyway

 remove_from_global_pool( k , ItemVcblr[ i ].second ,
                          ( ( BPar7 & 3 ) == 3 ) || ModDelete );
 ItemVcblr[ i ].second = InINF;

 // check if compacting FreList is appropriate- - - - - - - - - - - - - - - -
 // the issue with having indices of "free" position in the bundle stored in
 // a priority_queue is the following: if the bundle gets "full", but then is
 // "emptied", FreList may end up containing "many" elements, and in
 // particular elements that are >= get_max_name(), which therefore are
 // useless since they are obviously not in the bundle. the check above tries
 // to avoid that, but it may clearly fail (say, if small items are deleted
 // before large ones). checking if there are items with name >=
 // get_max_name() in FreList and deleting them is not cheap. the only
 // easy-to-check case is the one where FreList.size() > get_max_name():
 // if this happens, FreList is cleared and re-initialized

 if( FreList.size() > MxNm ) {
  FreList = {};
  for( Index h = 0 ; h < MxNm ; ++h )
   if( ItemVcblr[ h ].second == InINF )
    FreList.push( h );
  }
 }  // end( BundleSolver::Delete )

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdtaBP3( void )
{
 if( BPar5 == 0 )
  return;

 switch( BPar6 ) {
  case( 4 ):
   if( UpFiLmb[ NrFi ] > -INFshift )
    aBP3 = ( BPar5 > 0 ? BPar4 : BPar3 ) +
           Index( BPar5 / std::log10( EpsU / RelAcc ) );
   break;
  case( 3 ):
   if( UpFiLmb[ NrFi ] > -INFshift )
    aBP3 = ( BPar5 > 0 ? BPar4 : BPar3 ) +
           Index( BPar5 / std::sqrt( EpsU / RelAcc ) );
   break;
  case( 2 ):
   if( UpFiLmb[ NrFi ] > -INFshift )
    aBP3 = ( BPar5 > 0 ? BPar4 : BPar3 ) +
           Index( BPar5 * ( RelAcc / EpsU ) );
   break;
  case( 1 ):
   if( ! ( ParIter % Index( std::abs( BPar5 ) ) ) ) {
    if( BPar5 > 0 )
     aBP3++;
    else
     aBP3--;
    }
  }

 if( aBP3 > Index( BPar3 ) )
  aBP3 = BPar3;
 else
  if( aBP3 < BPar4 )
   aBP3 = BPar4;

 }  // end( BundleSolver::UpdtaBP3 )

/*--------------------------------------------------------------------------*/

bool BundleSolver::IsOptimal( double eps ) const
{
 if( ! RifeqFi )    // the linearization errors are not "properly computed"
  return( false );  // no way one can detect optimality

 if( eps <= 0 )
  eps = RelAcc;

 if( vStar.back() >= INFshift )  // some components have no subgradients
  return( false );               // no way one can detect optimality

 c_VarValue err = max_error( eps );
 if( err >= INFshift )
  return( false );

 // A significantly negative aggregate linearization error is inconsistent
 // with a valid cutting-plane model and must not trigger optimality. Let the
 // noise-reduction logic in compute() increase t and re-solve the master.
 if( Sigma < - err )
  return( false );

 if( ( tStar > 0 ) && ( DSTS + Sigma <= err ) )
  return( true );

 return( ( Sigma <= err ) &&
         ( NrmZFctr < INFshift ) && ( NrmZ <= NrmZFctr * NZEps ) );

 }  // end( BundleSolver::IsOptimal )

/*--------------------------------------------------------------------------*/

void BundleSolver::FModChg( VarValue shift , Index wFi )
{
 if( ( ! std::isnan( shift ) ) && ( ! f_convex ) )
  shift = - shift;

 if( shift == INFshift ) {      // function changed monotonically up
  if( UpFiLmb[ wFi ] < INFshift ) {
   UpFiLmb[ wFi ] = INFshift;   // reset upper function value for component
   --UpFiLmbdef;                // one less known
   }
  if( UpFiLmb.back() < INFshift ) {
   UpFiLmb.back() = INFshift;   // reset total upper function value
   --UpFiLmbdef;                // one less known
   }
  UpFiBest = INFshift;          // comprised best one
  return;
  }

 if( shift == -INFshift ) {     // function changed monotonically dn
  if( LwFiLmb[ wFi ] > -INFshift ) {
   LwFiLmb[ wFi ] = -INFshift;  // reset lower function value for component
   --LwFiLmbdef;                // one less known
   }
  if( LwFiLmb.back() > -INFshift ) {
   LwFiLmb.back() = -INFshift;  // reset total lower function value
   --LwFiLmbdef;                // one less known
   }
  f_global_LB = -INFshift;      // global LB no longer valid
  return;
  }

 if( std::isnan( shift ) ) {    // function changed unpredictably
  if( UpFiLmb[ wFi ] < INFshift ) {
   UpFiLmb[ wFi ] = INFshift;   // reset upper function value for component
   --UpFiLmbdef;                // one less known
   }
  if( UpFiLmb.back() < INFshift ) {
   UpFiLmb.back() = INFshift;   // reset total upper function value
   --UpFiLmbdef;                // one less known
   }
  UpFiBest = INFshift;          // and of course best one
  if( LwFiLmb[ wFi ] > -INFshift ) {
   LwFiLmb[ wFi ] = -INFshift;  // reset lower function value for component
   --LwFiLmbdef;                // one less known
   }
  if( LwFiLmb.back() > -INFshift ) {
   LwFiLmb.back() = -INFshift;  // reset total lower function value
   --LwFiLmbdef;                // one less known
   }
  f_global_LB = -INFshift;      // global LB no longer valid
  return;
  }

 // function changed by shift: just update everything

 if( UpFiLmb[ wFi ] < INFshift )
  UpFiLmb[ wFi ] += shift;

 if( UpFiLmb.back() < INFshift )
  UpFiLmb.back() += shift;

 if( UpRifFi[ wFi ] < INFshift )
  UpRifFi[ wFi ] += shift;

 if( UpRifFi.back() < INFshift )
  UpRifFi.back() += shift;

 if( UpFiBest < INFshift )
  UpFiBest += shift;

 if( LwFiLmb[ wFi ] > -INFshift )
  LwFiLmb[ wFi ] += shift;

 if( LwFiLmb.back() > -INFshift )
  LwFiLmb.back() += shift;

 f_global_LB += shift;

 }  // end( BundleSolver::FModChg )

/*--------------------------------------------------------------------------*/

void BundleSolver::remove_from_global_pool( Index k , Index i , bool hard )
{
 // update InvItemVcblr and all the associated fields to the fact that the
 // linearization currently in position i of the global pool of component k
 // is removed; actually, there may as well not be any linearization in
 // position i already
 //
 // the removal is "hard" (meaning the linearization is actually deleted) if
 // hard == true, and "soft" (the linearization is kept in the global pool
 // but deleted from the bundle) otherwise
 //
 // the actual removal from the global pool is not handled here

 InvItemVcblr[ k ][ i ] = hard ? InINF : vBPar2.back();
 while( MaxItem[ k ] && ( InvItemVcblr[ k ][ MaxItem[ k ] - 1 ] == InINF ) )
  --MaxItem[ k ];
 if( i < FrFItem[ k ] )   // creating a new "hole" before the FrFItem
  FrFItem[ k ] = i;       // this is the new FrFItem
 else                     // deleting something that may be FrFItem
  while( FrFItem[ k ] &&
         ( InvItemVcblr[ k ][ FrFItem[ k ] - 1 ] >=
                                ( ( BPar7 & 3 ) ? vBPar2.back() : InINF ) ) )
   --FrFItem[ k ];

 }  // end( BundleSolver::remove_from_global_pool )

/*--------------------------------------------------------------------------*/

Index BundleSolver::find_place_in_global_pool( Index k )
{
 // returns a suitable position in the global pool of component k, or InINF
 // if there is no free space

 if( FrFItem[ k ] < vBPar2[ k ] )  // there are free positions
  return( FrFItem[ k ] );          // return the first of them

 // if there are no free positions but there are less items in the bundle
 // (for component k) than the size of the global pool, there must be a
 // some linearization in the global pool that is not in the bundle: find
 // and return the first one (smallest position)
 Index gpp = InINF;
 if( NrItems[ k ] < vBPar2[ k ] )
  for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
   if( InvItemVcblr[ k ][ i ] >= vBPar2.back() ) {
    gpp = i;
    break;
    }

 return( gpp );

 }  // end( BundleSolver::find_place_in_global_pool )

/*--------------------------------------------------------------------------*/

void BundleSolver::add_to_global_pool( Index k , Index i , Index wh )
{
 // update ItemVcblr, InvItemVcblr and all the associated fields the to fact
 // that the linearization in position i in the global global_pool of
 // component k will be kept in the bundle at position wh; if wh == InINF,
 // this means the linearization is in the global pool but not in the bundle
 //
 // the actual addition to the global pool is not handled here

 if( wh < InINF ) {
  ItemVcblr[ wh ].first = k;
  ItemVcblr[ wh ].second = i;
  InvItemVcblr[ k ][ i ] = wh;
  if( wh + 1 > f_max_name )
   f_max_name = wh + 1;
  }
 else
  InvItemVcblr[ k ][ i ] = vBPar2.back();

 if( i >= MaxItem[ k ] )
  MaxItem[ k ] = i + 1;

 while( ( FrFItem[ k ] < MaxItem[ k ] ) &&
        ( InvItemVcblr[ k ][ FrFItem[ k ] ] <
                                ( ( BPar7 & 3 ) ? vBPar2.back() : InINF ) ) )
  ++FrFItem[ k ];

 }  // end( BundleSolver::add_to_global_pool( k , i , wh ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::reload_component_bundle( Index k )
{
 for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
  if( InvItemVcblr[ k ][ i ] < vBPar2.back() )
   add_to_bundle( k , i );
}

/*--------------------------------------------------------------------------*/

void BundleSolver::add_to_bundle( Index k , Index i )
{
 // add to the bundle (master problem) the item corresponding to the
 // linearization to be found at position i in the global pool of component
 // k; this assumes that the linearization is already there in the global
 // pool. if InvItemVcblr[ k ][ i ] < vBPar2.back(), i.e., the item is
 // already in the bundle, then it is replaced, otherwise it is added
 //
 // note that CheckSubG() or CheckCnst() need be called, but even if the
 // item is identical to some in the bundle already this information is
 // ignored and the item is inserted anyway; hence, if the check is active,
 // it is temporarily deactivated (and then re-activated)

 auto wh = InvItemVcblr[ k ][ i ];
 if( wh >= vBPar2.back() ) {  // the item is not there already
  wh = FindAPlace( k );       // find a "free" spot in the bundle
  if( wh == InINF )           // one must be there
   throw( std::logic_error(
              "BundleSolver::add_to_bundle: "
              "no space found in the bundle" ) );

  ++NrItems[ k ];              // keep count
  ++NrItems[ NrFi ];
  add_to_global_pool( k , i , wh );  // update dictionaries
  }
 else                          // the item is there already
  remove_cut_global( wh );       // remove it so that it can be replaced

 // local buffer for the new linearization coefficients
 std::vector< double > G1( NumVar , 0.0 );

 // recover the linearization from the C05Function
 v_c05f[ k ]->get_linearization_coefficients( G1.data() ,
                                              Range( 0 , NumVar ) , i );
 if( ! f_convex )
  chgsign( G1.data() , NumVar );

 // recover the physical raw constant. MasterPB owns the raw -> stored-b
 // translation/sign conversion for both diagonal and vertical cuts.
 auto Ai = rs( v_c05f[ k ]->get_linearization_constant( i ) );

 // append the physical ( G1 , Ai ) cut at slot wh of HardCmps[ k ];
 // MasterProblemBlock stores it in the requested primal/dual representation.
 if( MasterPB ) {
  MasterPB->add_cut( int( k ) , int( wh ) , std::move( G1 ) , Ai ,
                     v_c05f[ k ]->is_linearization_vertical( i ) );
  }

 }  // end( BundleSolver::add_to_bundle )

/*--------------------------------------------------------------------------*/

void BundleSolver::reset_bundle( void )
{
 // completely resets the bundle, because a (bunch of) Modification(s) saying
 // so has(ve) been received. this only affects the BundleSolver data
 // structures and the MPSolver, not the C05Function(s)

 OOBase.assign( vBPar2.back() , Inf< SIndex >() );

 ItemVcblr.assign( vBPar2.back() , std::make_pair( InINF , InINF ) );

 for( Index k = 0 ; k < NrFi ; ++k )
  InvItemVcblr[ k ].assign( vBPar2[ k ] , InINF );

 NrItems.assign( NrFi + 1 , 0 );
 FrFItem.assign( NrFi , 0 );
 MaxItem.assign( NrFi , 0 );

 f_max_name = 0;

 FreList = {};
 whisZ.assign( NrFi , InINF );
 Zvalid.assign( NrFi , false );

 whisG1.assign( NrFi , InINF );

 if( MasterPB )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ! ( NrEasy && IsEasy[ k ] ) )
    MasterPB->invalidate_subgradients( int( k ) );

 // reset the proximal parameter to its initial value: after a full reset the
 // bundle re-converges from an empty model, which takes a long run of null
 // steps; with t decreased on every null step ( e.g. intMnNSC == 1 ) t would
 // otherwise collapse ( observed ~1e-6 ) and, never reset, carry the collapsed
 // value across successive resets, so the trial point d* = -t z* stops moving
 // and every cut becomes a duplicate -> the bundle stalls in an endless
 // null-step loop. Starting each re-convergence from tInit avoids that
 t = tInit;
 Prevt = INFshift;
 if( MasterPB )
  MasterPB->set_t( t );

 }  // end( BundleSolver::reset_bundle )

/*--------------------------------------------------------------------------*/

Lst_sp_Mod::size_type BundleSolver::num_outstanding_Modification( void )
{
 return( v_mod.size() );
 }

/*--------------------------------------------------------------------------*/

bool BundleSolver::is_special_GroupMod( GroupModification & gmod )
{
 // recognise "special" GroupModification for changing the set of "active"
 // Variable of all the Objective at the same time; note that these
 // contain FunctionModVars* not necessarily C05FunctionModVars* because
 // the Modification may not be strongly quasi-additive

 if( gmod.sub_Modifications().size() != NrFi + ( f_lf ? 1 : 0 ) )
  return( false );

 auto smi = gmod.sub_Modifications().begin();
 auto sm0 = *(smi++);
 for( ; smi !=  gmod.sub_Modifications().end() ; ++smi )
  if( typeid( sm0 ) != typeid( *smi ) )
   return( false );

 smi = gmod.sub_Modifications().begin();
 ++smi;

 // check FunctionModVarsAddd
 if( const auto mod0 =
     std::dynamic_pointer_cast< FunctionModVarsAddd >( sm0 ) ) {
  for( ; smi != gmod.sub_Modifications().end() ; ++smi ) {
   auto modi = std::static_pointer_cast< FunctionModVarsAddd >( *smi );
   if( ( mod0->first() != modi->first() ) ||
       ( mod0->vars() != modi->vars() ) )
    throw( std::logic_error(
               "BundleSolver::is_special_GroupMod: "
               "different Variable change in components" ) );
   }

  return( true );
  }

 // check FunctionModVarsRngd
 if( const auto mod0 =
     std::dynamic_pointer_cast< FunctionModVarsRngd >( sm0 ) ) {
  for( ; smi != gmod.sub_Modifications().end() ; ++smi ) {
   auto modi = std::static_pointer_cast< FunctionModVarsRngd >( *smi );
   if( mod0->range() != modi->range() )
    throw( std::logic_error(
               "BundleSolver::is_special_GroupMod: "
               "different Variable change in components" ) );
   }

  return( true );
  }

 // check FunctionModVarsSbst
 if( const auto mod0 =
     std::dynamic_pointer_cast< FunctionModVarsSbst >( sm0 ) ) {
  for( ; smi != gmod.sub_Modifications().end() ; ++smi ) {
   auto modi = std::static_pointer_cast< FunctionModVarsSbst >( *smi );
   if( mod0->subset() != modi->subset() )
    throw( std::logic_error(
               "BundleSolver::is_special_GroupMod: "
               "different Variable change in components" ) );
   }

  return( true );
  }

 return( false );

 }  // end( BundleSolver::is_special_GroupMod )

/*--------------------------------------------------------------------------*/

void BundleSolver::flatten_Modification_list( Lst_sp_Mod & vmt , sp_Mod mod )
{
 const auto tmod = std::dynamic_pointer_cast< GroupModification >( mod );
 if( tmod && ( ! is_special_GroupMod( *tmod ) ) )
  for( auto submod : tmod->sub_Modifications() )
   flatten_Modification_list( vmt , submod );
 else
  vmt.push_back( mod );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::process_outstanding_Modification( void )
{
 // multiple-loop version, where several passes are done in order to gather
 // which kind of Modification have occurred and avoid doing costly work
 // more than once

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 0-th loop: "atomically flatten" v_mod into a temporary list to better
 // handle it, then clear it

 Lst_sp_Mod v_mod_tmp;

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;  // try to acquire lock, spin on failure

 for( auto mod : v_mod )
  flatten_Modification_list( v_mod_tmp , mod );

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

  #ifndef NDEBUG
  // high-verbosity diagnostic (LogVerb >= 7): account for every time the
  // Modification queue is drained, and how many Modification are in flight
  BLOG( 6 , std::endl << "process_outstanding_Modification: "
                      << v_mod_tmp.size() << " Modification(s)" );
 #endif


 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the 1st loop is made in reverse, from the latest Modification to the
 // earlies, and does the following:
 // - reset/change the upper/lower bounds that need to (check all shift() of
 //   all *FunctionMod*)
 // - check if some global pool has been "hard" reset, i.e., all the
 //   linearization in there have been deleted; this is brought about by
 //   a C05FunctionMod with type() == GlobalPoolRemoved and which().empty()
 //   or by any FunctionMod that does *not* imply strong quasi-additivity
 //   (i.e., that it is not a *C05*FunctionMod*)
 // - if a *FunctionMod* changing linearizations happens before a reset of the
 //   global pool (meaning it is found afterwards in the reverse order) it
 //   is deleted since it is useless (after having checked if it also impacts
 //   the upper/lower bounds)
 // - if there is more than one component, check that no "naked"
 //   *FunctionModVars* is there
 // - check that no ConstraintMod or VariableMod are there, since they are
 //   not handled (yet)
 //
 // Future optimization: during the 1st loop one could collect the set of
 // components that have been modified in any way and use it to skip the
 // allocation of NrFi-sized scratch vectors (such as reset[] below);
 // this would pay off when NrFi is very large and only a handful of
 // components actually change between successive compute()s

 std::vector< bool > reset( NrFi , false );

 // how many non-easy components are reset (easy never are)
 Index cntreset = 0;

 bool Fi0Chgd = false;  // true if the 0-th component changes

 bool to_delete;  // should have been defined inside, but there is not
                  // visible by the lambda

 for( auto rimod = v_mod_tmp.rbegin() ; rimod != v_mod_tmp.rend() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to rimod
      [ & to_delete , & v_mod_tmp ]( decltype( rimod ) & ri ) {
       if( to_delete )
        ri = std::reverse_iterator( v_mod_tmp.erase( std::next( ri ).base()
                                                     ) );
       else
        ++ri;
       }( rimod ) ) {
  to_delete = false;
  auto mod = *rimod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // first check if it is any kind of FunctionMod, since this gives immediate
  // access to the component, and any FunctionMod pertaining to an already
  // reset component can be almost immediately deleted
  if( const auto tmod = std::dynamic_pointer_cast< FunctionMod >( mod ) ) {
   if( tmod->function() == f_lf ) {
    const auto shift = tmod->shift();
    // special immediate treatment of the 0-th component, which is simple
    if( std::isnan( shift ) ) {  // is a C05FunctionModLin*
     Fi0Chgd = true;             // changing the coefficients
     Fi0Lmb = INFshift;          // the value is no longer known
     if( UpFiLmbdef == NrFi + 1 ) {  // the total upper bound was known
      UpFiLmb.back() = INFshift;     // it is no longer so
      UpFiLmbdef = NrFi;             // it will have to be recomputed
      }
     if( LwFiLmbdef == NrFi + 1 ) {  // the total lower bound was known
      LwFiLmb.back() = -INFshift;    // it is no longer so
      LwFiLmbdef = NrFi;             // it will have to be recomputed
      }
     UpFiBest = INFshift;            // and the best value as well
     }
    else {                       // a FunctionMod changing the constant
     #ifndef NDEBUG
      if( ( shift == INFshift ) || ( shift == -INFshift ) )
       throw( std::logic_error(
                  "BundleSolver::process_outstanding_Modification: "
                  "unexpected *FunctionMod* from LinearFunction" ) );
     #endif
     if( Fi0Lmb < INFshift )
      Fi0Lmb += shift;
     if( UpFiLmb.back() < INFshift )
      UpFiLmb.back() += shift;
     if( UpRifFi.back() < INFshift )
      UpRifFi.back() += shift;
     if( UpFiBest < INFshift )
      UpFiBest += shift;
     if( LwFiLmb.back() > -INFshift )
      LwFiLmb.back() += shift;
     f_global_LB += shift;
     }

    to_delete = true;  // in either case, all that had to be done
    continue;          // has been done
    }

   auto wFi = get_index_of_component( tmod->function() );

   // adjust or reset upper/lower values as needed
   // note that the list is scanned in reverse, hence these changes are
   // applied in reverse order. however, if the upper/lower values are reset
   // at any point in the list they stay reset forever. indeed, even if a
   // function has a finite shift after a reset, this says nothing because
   // there are no known values to shift. if, rather, the values are only
   // shifted by finite amounts, the total shift is the sum of the shift,
   // and the order of additions does not change the result
   if( wFi < NrFi )
    FModChg( tmod->shift() , wFi );
   else {
    // this is a FunctionMod coming from some unknown Function, not any
    // business of BundleSolver
    to_delete = true;
    continue;
    }

   if( NrEasy && IsEasy[ wFi ] ) {  // coming from an easy component
    // The easy component is represented exactly by its inner Block, registered
    // below MPBlock. The original physical Modification therefore reaches the
    // Solver attached to MPBlock independently and updates the master model.
    // The FunctionMod produced by the retained owner only tells BundleSolver
    // that its cached function values may be stale; FModChg() above has already
    // handled that. Easy components have no bundle/global-pool bookkeeping.
    to_delete = true;
    continue;
    }
   else                             // coming from a non-easy component
    if( reset[ wFi ] ) {
     // any FunctionMod after (before) one that completely resets the
     // component is useless, delete it and move forward (backward)
     to_delete = true;
     continue;
     }

   // if the component is not reset (yet), one must look in details what
   // exact type the *FunctionMod* is and react accordingly

   // a C05FunctionModRngd only changes existing linearizations, and
   // therefore is never a "hard" reset
   if( const auto ttmod =
       std::dynamic_pointer_cast< C05FunctionModRngd >( tmod ) ) {
    switch( ttmod->type() ) {
     case( C05FunctionMod::AllLinearizationChanged ):
     case( C05FunctionMod::AllEntriesChanged ):
      continue;
     default:
      throw( std::invalid_argument(
                 "BundleSolver::process_outstanding_Modification: "
                 "wrong type in C05FunctionModRngd" ) );
     }  // end( switch( ttmod->f_type ) )
    }  // end( if( ttmod == C05FunctionModRngd ) )

   // a C05FunctionModSbst only changes existing linearizations, and
   // therefore is never a "hard" reset; the only easy case is
   // NothingChanged, which by definition does nothing save for the
   // shift(), that has been dealt with already
   if( const auto ttmod =
       std::dynamic_pointer_cast< C05FunctionModSbst >( tmod ) ) {
    switch( ttmod->type() ) {
     case( C05FunctionMod::AllLinearizationChanged ):
     case( C05FunctionMod::AllEntriesChanged ):
      continue;
     default:
      throw( std::invalid_argument(
                 "BundleSolver::process_outstanding_Modification: "
                 "wrong type in C05FunctionModSbst" ) );
     }  // end( switch( ttmod->f_type ) )
    }  // end( if( ttmod == C05FunctionModSbst ) )

   // a C05FunctionMod of type GlobalPoolRemoved with which.empty() resets
   // all the component. NothingChanged by definition does nothing (save
   // for the shift(), that has been dealt with already) ... except
   // possibly doing a lot, i.e., signalling the switch from convex to
   // concave, or vice-versa, which is not allowed; all other cases
   // will have to be dealt with later
   if( const auto ttmod =
       std::dynamic_pointer_cast< C05FunctionMod >( tmod ) ) {
    switch( ttmod->type() ) {
     case( C05FunctionMod::GlobalPoolRemoved ):
      if( ttmod->which().empty() ) {
       reset[ wFi ] = true;
       ++cntreset;
       to_delete = true;
       }
      continue;
     case( C05FunctionMod::NothingChanged ):
      if( v_c05f[ wFi ]->is_convex() != f_convex )
       throw( std::logic_error(
                  "BundleSolver::process_outstanding_Modification: "
                  "convex/concave switch not allowed" ) );
      to_delete = true;
     case( C05FunctionMod::AllLinearizationChanged ):
     case( C05FunctionMod::AllEntriesChanged ):
     case( C05FunctionMod::AlphaChanged ):
     case( C05FunctionMod::GlobalPoolAdded ):
      continue;
     default:
      throw( std::invalid_argument(
                 "BundleSolver::process_outstanding_Modification: "
                 "wrong type in C05FunctionMod" ) );
     }  // end( switch( ttmod->f_type ) )
    }  // end( if( ttmod == C05FunctionMod ) )

   // a C05FunctionModLin* only changes existing linearizations, and
   // therefore is never a "hard" reset
   if( std::dynamic_pointer_cast< C05FunctionModLinRngd >( tmod ) )
    continue;

   // if control reaches this point, mod is (indistinguishable from) a base
   // FunctionMod, and in particular it is not a C05FunctionMod*. hence the
   // change in the Function is not quasi-additive, and therefore a fortiori
   // not strongly quasi-additive. as a result, this is a "hard" reset

   reset[ wFi ] = to_delete = true;
   ++cntreset;

   }  // end( if( tmod == FunctionMod ) )

  // a "naked" FunctionModVars is only allowed if there is only one
  // component (comprised the linear one). if it is allowed, it is
  // of no consequence here, except for the possible effect on the
  // function values, if it is a C05FunctionModVars*, meaning that it
  // represents a strongly quasi-additive variable change. if not, the
  // variable change also implies a reset
  // in no case, however, the Modification is removed from the list
  if( const auto tmod = std::dynamic_pointer_cast< FunctionModVars >( mod ) ) {
   if( ( NrFi > 1 ) || f_lf )
    throw( std::invalid_argument(
               "BundleSolver::process_outstanding_Modification: "
               "naked FunctionModVars not allowed" ) );

   auto wFi = get_index_of_component( tmod->function() );

   FModChg( tmod->shift() , wFi );  // change/reset upper/lower values

   if( std::dynamic_pointer_cast< C05FunctionModVarsAddd >( tmod ) )
    continue;

   if( std::dynamic_pointer_cast< C05FunctionModVarsRngd >( tmod ) )
    continue;

   if( std::dynamic_pointer_cast< C05FunctionModVarsSbst >( tmod ) )
    continue;

   // if control reaches here, this is a FunctionModVars* that is not a
   // C05FunctionModVars*, i.e., a non strongly quasi-additive variable
   // change, which implies a "hard" reset for the component
   reset[ wFi ] = true;
   ++cntreset;
   continue;

   }  // end( if( tmod == FunctionModVars ) )

  // a GroupModification here can only be a bunch of identical
  // *FunctionModVar*: pick the first one and act on it
  if( const auto tmod =
      std::dynamic_pointer_cast< GroupModification >( mod ) ) {
   auto fmod = tmod->sub_Modifications().front();

   if( std::dynamic_pointer_cast< C05FunctionModVarsAddd >( fmod ) )
    continue;

   if( std::dynamic_pointer_cast< C05FunctionModVarsRngd >( fmod ) )
    continue;

   if( std::dynamic_pointer_cast< C05FunctionModVarsSbst >( fmod ) )
    continue;

   // if control reaches here, this is a FunctionModVars* that is not a
   // C05FunctionModVars*, i.e., a non strongly quasi-additive variable
   // change, which implies a "hard" reset for *all* components
   if( NrEasy ) {
    for( Index k = 0 ; k < NrFi ; ++k )
     if( ! IsEasy[ k ] )
      reset[ k ] = true;
    }
   else
    reset.assign( NrFi , true );
   cntreset = NrFi - NrEasy;
   continue;

   }  // end( if( tmod == GroupModification ) )

  if( std::dynamic_pointer_cast< ConstraintMod >( mod ) )
   throw( std::invalid_argument(
              "BundleSolver::process_outstanding_Modification: "
              "ConstraintMod not handled (yet)" ) );

  if( std::dynamic_pointer_cast< VariableMod >( mod ) )
   throw( std::invalid_argument(
              "BundleSolver::process_outstanding_Modification: "
              "VariableMod not handled (yet)" ) );

  if( std::dynamic_pointer_cast< BlockMod >( mod ) )
   throw( std::invalid_argument(
              "BundleSolver::process_outstanding_Modification: "
              "BlockMod not handled (yet)" ) );

  if( std::dynamic_pointer_cast< BlockModAD >( mod ) )
   throw( std::invalid_argument(
              "BundleSolver::process_outstanding_Modification: "
              "BlockModAD not handled (yet)" ) );

  // if control reaches here, the Modification is "unknown", probably a
  // "physical" Modification that BundleSolver does not care about

  to_delete = true;

  }  // end( 1st loop, in reverse )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if any global pool has been reset/all of them have reset, then delete all
 // items in the corresponding bundle/reset all the bundle

 if( cntreset == NrFi - NrEasy )  // all (non-easy) components have been reset
  reset_bundle();
 else
  if( cntreset ) {           // at least a (non-easy) component has been reset
   for( Index k = 0 ; k < NrFi ; ++k )
    if( reset[ k ] ) {       // reset[ k ] ==> ! IsEasy[ k ]
     // if BundleSolver "plays nice" with other Solvers, it keeps track of
     // linearizations in the global pool even if they are not in the bundle
     // to avoid overriding them; but there is no longer anything in the
     // global pool, so reset any such information; do this first because
     // MaxItem[ k ] is zeroed during the next loop
     if( ( BPar7 & 3 ) < 3 )
      std::fill( std::next( InvItemVcblr[ k ].begin() , MaxItem[ k ] ) ,
                 InvItemVcblr[ k ].end() , InINF );

     // delete all linearizations in the bundle for this component (in
     // the sense of updating the BundleSolver data structures, since they
     // have been deleted already from the global pool)
     for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
      if( InvItemVcblr[ k ][ i ] < vBPar2.back() )
       Delete( InvItemVcblr[ k ][ i ] , true );
     }
   }

 // After this point, all the Modification adding, deleting or modifying
 // linearizations are significant: they either pertain to components that
 // have never been reset, or are the remaining ones after the (last) one
 // resetting the component

 if( v_mod_tmp.empty() ) {  // no more Modification to process
  // ordinarily, changes in the 0-th component would be dealt with at the
  // end, but if this is the only thing that happened the end is now, so
  // they have to be dealt with immediately
  if( Fi0Chgd ) {
   // the "0-th component" is the linear part, which has no
   // PolyhedralFunctionBlock counterpart -- it lives in the master
   // Objective itself (set_b()). No invalidation is needed; the
   // linear term changes are picked up by the next compute() through
   // the usual LinearFunctionMod issued by set_b().
   }

  return;                   // all done
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 2nd loop, again in reverse: check for "soft" reset of components, i.e.,
 // when all existing linearization changes. any Modification that changes
 // the linearizations happening before a "soft" reset of the global pool
 // (meaning it is found afterwards in the reverse order) is deleted since it
 // is useless.
 //
 // note that the linearization error of a linearization depends on both the
 // initial constant (\alpha), the linearization itself (g) and the current
 // stability centre (Lambda); thus, if any of those changes, the
 // linearization error need be recomputed. hence reset[ k ] == true
 // implies AlphaC[ k ] == true; the reverse implication does not hold, i.e.,
 // AlphaC[ k ] == true with reset[ k ] == false is possible and it means
 // that only the constants need be changed, but not all the rest. in this
 // loop we consider changes of \alpha and g, while in the 4th loop we will
 // consider changes of Lambda due to the removal of variables; additions
 // never create problems since new variables are always initialized to 0, and
 // therefore they never change the existing linearization error. in fact, not
 // all changes of g necessarily change the linearization error: if a
 // component g_i changes such that Lambda_i == 0, this has no impact. however,
 // in this reverse loop the map between the Lambda[] vector and the indices
 // in the Modification is nontrivial (if additions/removals happened), which
 // makes checking this too complicated. anyway, the issue will go away in the
 // version of BundleSolver that does not use MPSolver since there the
 // linearizations will (likely) be represented by means of their "naked"
 // constant \alpha rather than by their linearization error
 //
 // note that Modification changing the linearizations happening *after* a
 // "soft" reset of the global pool (meaning it is found *before* in the
 // reverse order) is also useless, since a reset forces the re-reading of all
 // linearizations, which by definition happens at their current (final)
 // state. yet, this is not done immediately
 //
 // note that we make no serious attempt at keeping track of the combined
 // effect of all changes, in order to detect if a large set of small
 // changes actually implies a reset. this is complicated for "horizontal"
 // changes (for all linearizations, a range/subset of entries) because the
 // names of the changed Variable may not be current (additions/deletions may
 // happen in the meantime), and keeping track is too burdensome. similarly
 // for "vertical" changes (a set of specific linearizations). some steps
 // in this direction will perhaps be done in later stages of development

 if( cntreset ) {                // if there was any hard reset
  reset.assign( NrFi , false );  // reset reset (couldn't resist)
  cntreset = 0;
  }

 std::vector< bool > AlphaC( NrFi , false );

 for( auto rimod = v_mod_tmp.rbegin() ; rimod != v_mod_tmp.rend() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to rimod
      [ & to_delete , & v_mod_tmp ]( decltype( rimod ) & ri ) {
       if( to_delete )
        ri = std::reverse_iterator( v_mod_tmp.erase( std::next( ri ).base()
                                                     ) );
       else
        ++ri;
       }( rimod ) ) {
  to_delete = false;
  auto mod = *rimod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // a C05FunctionModRngd only changes a range of the linearizations,
  // and therefore is not considered a "soft" reset even if which().empty()
  // in fact the range could be so large as to be (almost) all the
  // variables, which would count as a reset, but so far we don't attempt
  // at detecting this. the only easy case would be NothingChanged, but
  // any such C05FunctionModRngd has been deleted already. however, if the
  // component is "soft" reset already, it can be deleted
  //
  // AllEntriesChanged and AllLinearizationChanged are equivalent, since
  // even if only g changes, also the linearization error does (unless it
  // only changes in places where Lambda == 0, but this cannot be checked
  // efficiently)
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModRngd >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   switch( tmod->type() ) {
    case( C05FunctionMod::AllEntriesChanged ):
    case( C05FunctionMod::AllLinearizationChanged ):
     if( tmod->which().empty() ) {  // reset of the constants
      if( reset[ wFi ] )            // component reset already
       to_delete = true;            // nothing else to do
      else                          // component not reset
       AlphaC[ wFi ] = true;        // reset the constants
      }
     continue;
    default:  // this must not happen
     throw( std::invalid_argument(
                "BundleSolver::process_outstanding_Modification: "
                "wrong type() in C05FunctionModRngd" ) );
    }
   }  // end( if( tmod == C05FunctionModRngd ) )

  // a C05FunctionModSbst only changes a subset of the linearizations,
  // and therefore is not considered a "soft" reset even if which().empty()
  // in fact the subset could be so large as to be (almost) all the
  // variables, which would count as a reset, but so far we don't attempt
  // at detecting this. the only easy case would be NothingChanged, but
  // any such C05FunctionModSbst has been deleted already. however, if the
  // component is "soft" reset already, it can be deleted
  //
  // AllEntriesChanged and AllLinearizationChanged are equivalent, since
  // even if only g changes, also the linearization error does (unless it
  // only changes in places where Lambda == 0, but this cannot be checked
  // efficiently)
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModSbst >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   switch( tmod->type() ) {
    case( C05FunctionMod::AllEntriesChanged ):
    case( C05FunctionMod::AllLinearizationChanged ):
     if( tmod->which().empty() ) {  // reset of the constants
      if( reset[ wFi ] )            // component reset already
       to_delete = true;            // nothing else to do
      else                          // component not reset
       AlphaC[ wFi ] = true;        // reset the constants
      }
     continue;
    default:  // this must not happen
     throw( std::invalid_argument(
                "BundleSolver::process_outstanding_Modification: "
                "wrong type() in C05FunctionModSbst" ) );
    }
   }  // end( if( tmod == C05FunctionModSbst ) )

  // a C05FunctionMod of type AllLinearizationChanged or AllEntriesChanged
  // with which.empty() "soft" resets all the component and Alpha;
  // AlphaChanged only changes the constants (obviously)
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionMod >( mod ) ) {
   if( ! tmod->which().empty() )   // we only react to which().empty()
    continue;

   auto wFi = get_index_of_component( tmod->function() );

   switch( tmod->type() ) {
    case( C05FunctionMod::AlphaChanged ):
     AlphaC[ wFi ] = true;
     break;
    case( C05FunctionMod::AllEntriesChanged ):
    case( C05FunctionMod::AllLinearizationChanged ):
     if( ! reset[ wFi ] )
      ++cntreset;
     reset[ wFi ] = AlphaC[ wFi ] = true;
     break;
    default:  // this must not happen, as GlobalPoolRemoved with
              // which.empty() has been dealt with and deleted before
     throw( std::invalid_argument(
                "BundleSolver::process_outstanding_Modification: "
                "wrong type in C05FunctionMod with empty which()" ) );
    }

   to_delete = true;
   continue;

   }  // end( if( tmod == C05FunctionMod ) )

  // a C05FunctionModLinRngd implies that a specific range in all the
  // linearizations must be changed by adding; this is never considered a
  // "soft" reset even if in fact the range could be so large as to be
  // (almost) all the variables, which would count as a reset, but so far
  // we don't attempt at detecting this. however, if the component is "soft"
  // reset already, it can be deleted
  //
  // again, on the grounds that changing the linearization (presumably)
  // changes the linearization errors as well, this also resets all the
  // Alpha
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModLinRngd >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   if( reset[ wFi ] )            // component reset already
    to_delete = true;            // nothing else to do
   else                          // component not reset
    AlphaC[ wFi ] = true;        // reset the constants

   continue;

   }  // end( if( tmod == C05FunctionModLinRngd ) )

  // a C05FunctionModLinSbst implies that a specific subset in all the
  // linearizations must be changed by adding; this is never considered a
  // "soft" reset even if in fact the subset could be so large as to be
  // (almost) all the variables, which would count as a reset, but so far
  // we don't attempt at detecting this. however, if the component is "soft"
  // reset already, it can be deleted
  //
  // again, on the grounds that changing the linearization (presumably)
  // changes the linearization errors as well, this also resets all the
  // Alpha
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModLinSbst >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   if( reset[ wFi ] )            // component reset already
    to_delete = true;            // nothing else to do
   else                          // component not reset
    AlphaC[ wFi ] = true;        // reset the constants

   continue;

   }  // end( if( tmod == C05FunctionModLinSbst ) )

  // a C05FunctionModLin implies that *all* the linearizations must be
  // changed by adding them \delta; this may in principle be handled in
  // a specialised way by BundleSolver, but is currently not, and
  // therefore it is a "full" reset
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModLin >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   if( ! reset[ wFi ] )
    ++cntreset;
   AlphaC[ wFi ] = reset[ wFi ] = to_delete = true;

   }  // end( if( tmod == C05FunctionModLin ) )
  }  // end( 2nd loop, again in reverse )

 // note that even if there were no more Modification to process we could not
 // stop because this means that reset[ k ] == true and/or AlphaC[ k ] == true
 // for some k. in fact v_mod_tmp was not empty(), and elements can be removed
 // from it only if some component is "soft" reset.

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 3rd loop, forward: prepare for addition/removal/changes of individual
 // linearization for each component by computing the four sets
 // - linearizations that need be removed
 // - linearizations that need be added
 // - linearizations that need be changed
 // - constants that need be changed
 // Note that:
 // - due to limitations in the MPSolver interface, changing a linearization
 //   implies changing its constant; therefore, Cchg[] contains changes in
 //   the constants only and nothing else, which means that has empty
 //   intersection with all other three sets
 // - if a linearization that is added/changed is later removed, it is no
 //   longer added/changed
 // - if a linearization that is removed/changed is later added it is no
 //   longer removed/changed (adding over an existing linearization changes
 //   it anyway, no reason to remove it)
 // - thus, Addd[], Rmvd[], Chgd[] and Cchg[] all have empty intersection
 // - linearization changes to a reset component can be ignored
 // - constant changes when all constants change can be ignored
 // - due to limitations in the MPSolver interface, "horizontal" changes (of
 //   a given range/subset of entries) to a subset of linearizations are
 //   not supported. these must be either mapped in "horizontal" changes to
 //   *all* linearizations (of a given component), or to "vertical" changes
 //   of all components to a subset of linearization. somewhat arbitrarily,
 //   the second option is chosen here. as a consequence, C05FunctionModRngd
 //   and C05FunctionModSbst are considered C05FunctionMod. note that
 //   C05FunctionMod* with which().empty() still remain untreated, as well
 //   as C05FunctionModLinRngd and C05FunctionModLinSbst (that by definition
 //   always concern all the linearizations), while C05FunctionModLin have
 //   been dealt with already

 std::vector< Subset > Addd( NrFi );  // added linearizations
 std::vector< Subset > Rmvd( NrFi );  // removed linearizations
 std::vector< Subset > Chgd( NrFi );  // changed linearizations
 std::vector< Subset > Cchg( NrFi );  // changed constants

 for( auto imod = v_mod_tmp.begin() ; imod != v_mod_tmp.end() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to imod
      [ & to_delete , & v_mod_tmp ]( decltype( imod ) & it ) {
       if( to_delete )
        it =  v_mod_tmp.erase( it );
       else
        ++it;
       }( imod ) ) {
  to_delete = false;
  auto mod = *imod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly
  //
  // actually, only C05FunctionMod need be treated here. we do not need to
  // distinguish between the base and the derived classes because we take
  // the change of a range/subset of the entries for a change of the whole
  // linearization. however, we do in fact distinguish them because we
  // ignore the Modification if ttmod->which().empty() and type() ==
  // AllEntriesChanged or == AllLinearizationChanged. any such C05FunctionMod
  // has been deleted (at worst) in the 2nd loop, and therefore here it must
  // be a C05FunctionModRngd or C05FunctionModSbst
  //
  // also, note that NothingChanged and AlphaChanged cannot happen, since
  // again these cases have been dealt with already

  if( const auto tmod = std::dynamic_pointer_cast< C05FunctionMod >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );

   switch( tmod->type() ) {
    case( C05FunctionMod::AlphaChanged ):
     assert( ! tmod->which().empty() );
     // this cannot happen, the case has been dealt with and deleted
     if( AlphaC[ wFi ] )  // constants are reset already
      break;              // nothing else to do

     // note that reset[ wFi ] ==> AlphaC[ wFi ], hence surely reset[ wFi ]
     // is false here
     // add to Cchg[ wFi ] the names in tmod->which(), save those that are
     // in either Addd[ wFi ], or Rmvd[ wFi ], or Chgd[ wFi ]
     if( Addd[ wFi ].empty() && Rmvd[ wFi ].empty() && Chgd[ wFi ].empty() )
      set_union_in_place( Cchg[ wFi ] , tmod->which() );
     else {
      // only those constants corresponding to linearizations that are not
      // added or removed or changed need be changed
      Subset tmp( tmod->which() );

      if( ! Addd[ wFi ].empty() )
       set_difference_in_place( tmp , Addd[ wFi ] );

      if( ( ! Rmvd[ wFi ].empty() ) && ( ! tmp.empty() ) )
       set_difference_in_place( tmp , Rmvd[ wFi ] );

      if( ( ! Chgd[ wFi ].empty() ) && ( ! tmp.empty() ) )
       set_difference_in_place( tmp , Chgd[ wFi ] );

      set_union_in_place( Cchg[ wFi ] , std::move( tmp ) );
      }
     break;
    case( C05FunctionMod::AllLinearizationChanged ):
    case( C05FunctionMod::AllEntriesChanged ):
     // AllEntriesChanged is completely equivalent to AllLinearizationChanged
     // since the change in the linearization implies that of the constant
     // if tmod->which().empty(), this must actually be either a
     // C05FunctionModRngd or a C05FunctionModSbst: save it, for it
     // will be dealt with in the next loop
     if( tmod->which().empty() )
      continue;

     if( reset[ wFi ] )  // changes in reset components
      break;             // are ignored

     // add to Chgd[ wFi ] the names in tmod->which(), save those that are
     // in either Addd[ wFi ] or Rmvd[ wFi ]
     if( Addd[ wFi ].empty() && Rmvd[ wFi ].empty() )
      set_union_in_place( Chgd[ wFi ] , tmod->which() );
     else {
      // only those items that are not added and/or removed need be changed
      Subset tmp( tmod->which() );

      if( ! Addd[ wFi ].empty() )
       set_difference_in_place( tmp , Addd[ wFi ] );

      if( ( ! Rmvd[ wFi ].empty() ) && ( ! tmp.empty() ) )
       set_difference_in_place( tmp , Rmvd[ wFi ] );

      set_union_in_place( Chgd[ wFi ] , std::move( tmp ) );
      }

     // changed linearizations imply changed constants
     if( ! AlphaC[ wFi ] )
      set_difference_in_place( Cchg[ wFi ] , Chgd[ wFi ] );

     break;

    case( C05FunctionMod::GlobalPoolAdded ):
     // add to Addd[ wFi ] the names in tmod->which(), and remove them
     // from Rmvd[ wFi ], and from Chgd[ wFi ] if the component is not
     // reset, and from Cchg[ wFi ] is constants are not reset

     set_union_in_place( Addd[ wFi ] , tmod->which() );

     set_difference_in_place( Rmvd[ wFi ] , tmod->which() );

     if( ! reset[ wFi ] )
      set_difference_in_place( Chgd[ wFi ] , tmod->which() );

     if( ! AlphaC[ wFi ] )
      set_difference_in_place( Cchg[ wFi ] , tmod->which() );

     break;

    case( C05FunctionMod::GlobalPoolRemoved ):
     // add to Rmvd[ wFi ] the names in tmod->which(), and remove them
     // from Addd[ wFi ], and from Chgd[ wFi ] if the component is not
     // reset, and from Cchg[ wFi ] is constants are not reset

     set_union_in_place( Rmvd[ wFi ] , tmod->which() );

     set_difference_in_place( Addd[ wFi ] , tmod->which() );

     if( ! reset[ wFi ] )
      set_difference_in_place( Chgd[ wFi ] , tmod->which() );

     if( ! AlphaC[ wFi ] )
      set_difference_in_place( Cchg[ wFi ] , tmod->which() );

    }  // end( switch( tmod->type() ) )

   to_delete = true;   // delete it

   }  // end( if( tmod == C05FunctionMod ) )
  }  // end( 3rd loop, forward )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now delete all linearization that need to, if any

 if( std::find_if( Rmvd.begin() , Rmvd.end() ,
                   []( Subset & Rk ) { return( ! Rk.empty() ); }
                   ) != Rmvd.end() ) {
  // at least a component has had linearizations removed, but note that not
  // all linearizations need be in the bundle; if they are not they are
  // just removed from the global pool (if they are there)

  for( Index k = 0 ; k < NrFi ; ++k )
   for( auto i : Rmvd[ k ] ) {
    auto h = InvItemVcblr[ k ][ i ];
    if( h < vBPar2.back() )
     Delete( h , true );
    else
     remove_from_global_pool( k , i , true );
    }
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 4th loop: handle exclusively Variable addition and removal, in forward
 // order. all removals of existing Variable are immediately performed,
 // diminishing NumVar. additions of Variable are "cached", and only performed
 // once (if ever) after the end of the cycle. removals of Variable that have
 // not been added yet just decreases the number of new Variable to be added
 // at the end. all corresponding Modification are removed from the list
 //
 // this loop may also force some linearization errors to be reset since they
 // depend not only on the initial constant (\alpha), but also from the
 // linearization itself (g) and the current stability centre (Lambda); if
 // any of those changes, the linearization error need be recomputed. in this
 // loop we consider changes of Lambda due to the removal of variables, while
 // additions never create this problem since new variables are always
 // initialized to 0, and therefore they never change the existing
 // linearization error per-se (recall that changes in \alpha and g were
 // already handled by the 2nd loop). note, however, that not all changes of
 // Lambda necessarily change the linearization error: if a Lambda_i == 0 is
 // removed, this has no impact. because this loop (unlike the 2nd one) is
 // forward and does deletions immediately, Lambda is always "in synch" with
 // the indices of the Modification and therefore the check can be easily
 // done. this is all the more important since variables are removed in
 // parallel for all components, hence if a variable with Lambda_i != 0 is
 // removed then all the linearization errors must be reset
 Index to_add = 0;
 bool addd_vars = false;  // if any Variable has ever been added
 bool rmvd_vars = false;  // if any Variable has ever been removed

 // Sparse mode bookkeeping: list of LamVcblr global indices whose
 // v_ref_count reached 0 during this Modification batch (FunctionMod
 // VarsRngd / VarsSbst handlers append here). After the 4th loop, we
 // compact LamVcblr / Lambda* / each v_c05f's global-index map /
 // Lambda2Idx and call MasterPB->remove_vars on this set to reclaim
 // the master rows.
 std::vector< Index > globally_to_remove;

 for( auto imod = v_mod_tmp.begin() ; imod != v_mod_tmp.end() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to imod
      [ & to_delete , & v_mod_tmp ]( decltype( imod ) & it ) {
       if( to_delete )
        it =  v_mod_tmp.erase( it );
       else
        ++it;
       }( imod ) ) {
  to_delete = false;
  auto mod = *imod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // not that we do not distinguish C05FunctionModVars* from "plain"
  // FunctionModVars*, since the only difference is whether or not the
  // operation is strongly quasi-additive, i.e., it implies or not a "hard"
  // reset, but this has already been acted upon

  {
   // a "naked" FunctionModVars
   auto tmod = std::dynamic_pointer_cast< FunctionModVars >( mod );
   if( ! tmod ) {
    // if it is not a "naked" FunctionModVars, it can still be a group of
    // identical *FunctionModVars* "dressed" into a GroupModification
    if( const auto gmod =
        std::dynamic_pointer_cast< GroupModification >( mod ) )
     // if so, pick the first one and act on it
     tmod = std::static_pointer_cast< FunctionModVars >(
                                         gmod->sub_Modifications().front() );
    }

   if( tmod ) {
    // if we have a *FunctionModVars*, we have to distinguish its exact type
    // and add/delete Variable accordingly; in all cases, however, the
    // Modification is processed and can be deleted
    to_delete = true;

    // auto-promote dense → sparse on the first per-Function (i.e. NOT
    // arriving as a lockstep GroupModification) FunctionModVars*. The
    // dense path assumes a lockstep change across all components, so
    // applying a single-component Mod to it would corrupt the global
    // Lambda. Promotion materialises the identity local-to-global map
    // for every v_c05f[ h ] (and rebuilds Lambda2Idx + v_ref_count from
    // the current dense invariant) before letting the (sparse) handlers
    // below process the Mod. A naked FunctionModVars in the 4th loop
    // can come from either a single-Mod (per-Function by construction)
    // or a non-special GroupMod already flattened by
    // flatten_Modification_list — both signal divergence.
    if( ( ! f_sparse_lambda ) &&
        ( ! std::dynamic_pointer_cast< GroupModification >( mod ) ) ) {
     f_sparse_lambda = true;
     Lambda2Idx.clear();
     Lambda2Idx.reserve( LamVcblr.size() );
     for( Index i = 0 ; i < LamVcblr.size() ; ++i )
      Lambda2Idx.emplace( LamVcblr[ i ] , i );
     // in dense mode every v_c05f[ h ] (and f_lf, if any) sees the
     // full LamVcblr in identity order, so each global slot is
     // referenced by every component
     const Index refs = v_c05f.size() + ( f_lf ? 1 : 0 );
     v_ref_count.assign( LamVcblr.size() , refs );
     // identity map of size NumVar + 1 with trailing Inf< Index >()
     const Index lN = LamVcblr.size();
     v_local2global.assign( v_c05f.size() , {} );
     for( Index h = 0 ; h < v_c05f.size() ; ++h ) {
      auto & id_map = v_local2global[ h ];
      id_map.reserve( lN + 1 );
      for( Index i = 0 ; i < lN ; ++i )
       id_map.push_back( i );
      id_map.push_back( Inf< Index >() );
      }
     }

    if( const auto ttmod =
        std::dynamic_pointer_cast< FunctionModVarsAddd >( tmod ) ) {
     addd_vars = true;

     if( f_sparse_lambda ) {
      // Sparse Lambda Addd: ttmod->first() is local (= loc_NV[ h ] at
      // the time the Mod was issued) and ttmod->vars() is the subset
      // of new globals that v_c05f[ h ] actually couples to. We
      // translate each ColVariable * to its global LamVcblr index —
      // appending it to the global Lambda only the first time we
      // encounter it across all sparse Mods — and extend
      // v_local2global[ h ] accordingly.

      const auto h = get_index_of_component( ttmod->function() );

      // strip the Inf< Index >() terminator before appending; we will
      // re-append it once we're done with this h's add Mod
      auto & m = v_local2global[ h ];
      if( ( ! m.empty() ) && ( m.back() == Inf< Index >() ) )
       m.pop_back();

      for( auto v : ttmod->vars() ) {
       const auto p = static_cast< ColVariable * >( v );
       auto [ it , inserted ] =
                              Lambda2Idx.try_emplace( p , LamVcblr.size() );
       if( inserted ) {
        LamVcblr.push_back( p );
        v_ref_count.push_back( 1 );
        ++to_add;  // a genuinely new global variable was added
        }
       else
        // v_c05f[ h ] picks up an already-existing global; the global
        // slot is now referenced by one more component, which the
        // Rngd / Sbst handlers will decrement back on removal
        ++v_ref_count[ it->second ];
       m.push_back( it->second );
       }

      // re-append the terminator; the global slot for the new vars is
      // at the end of LamVcblr / v_local2global[ h ], so monotonicity
      // is preserved by construction
      m.push_back( Inf< Index >() );
      continue;
      }

     if( ! to_add ) {
      // Dense mode: every component sees the same active vars in the
      // same order, so first() must equal NumVar (the global position
      // of the next slot) on the very first Mod
      if( ttmod->first() != NumVar )
       throw( std::logic_error(
                  "BundleSolver::process_outstanding_Modification: "
                  "wrong Variable names in FunctionModVars" ) );
      }

     to_add += ttmod->vars().size();
     continue;

     } // end( if( tmod == FunctionModVarsAddd ) )

    if( const auto ttmod =
        std::dynamic_pointer_cast< FunctionModVarsRngd >( tmod ) ) {
     if( f_sparse_lambda ) {
      // Sparse Lambda Rngd: range is in v_c05f[ h ]'s LOCAL index
      // space. Drop the affected local slots from v_local2global[ h ]
      // and decrement the global refcount for each. Any slot whose
      // refcount reaches 0 is queued for global removal (LamVcblr /
      // Lambda* / MasterPB->remove_vars), to be applied in one shot
      // after the 4th loop. Also invalidate linearization errors on
      // any nonzero Lambda removed.
      const auto h = get_index_of_component( ttmod->function() );
      auto & m = v_local2global[ h ];
      // m has size loc_NV + 1 with trailing Inf< Index >(); the valid
      // range is [ 0 , loc_NV ).
      const Index loc_NV = m.size() - 1;
      const Index r0 = std::min( ttmod->range().first , loc_NV );
      const Index r1 = std::min( ttmod->range().second , loc_NV );
      for( Index l = r0 ; l < r1 ; ++l ) {
       const Index g = m[ l ];
       if( std::abs( Lambda[ g ] ) > 1e-12 )
        std::fill( AlphaC.begin() , AlphaC.end() , true );
       if( --v_ref_count[ g ] == 0 )
        globally_to_remove.push_back( g );
       }
      m.erase( m.begin() + r0 , m.begin() + r1 );
      rmvd_vars = true;
      continue;
      }

     rmvd_vars = true;
     Range rng = ttmod->range();

     // if any of the deleted Variable is nonzero, the linearization errors
     // will have to be recomputed for all components
     for( Index i = std::min( rng.first , NumVar ) ;
          i < std::min( rng.second , NumVar ) ; ++i )
      if( std::abs( Lambda[ i ] ) > 1e-12 ) {
       std::fill( AlphaC.begin() , AlphaC.end() , true );
       break;
       }

     if( ( rng.first == 0 ) && ( rng.second >= NumVar ) ) {
      NumVar = 0;                      // deleting *all* Variable
      Lambda.clear();
      Lambda1.clear();
      LmbdBst.clear();
      if( MasterPB )
       MasterPB->remove_vars( nullptr , 0 );  // remove from MP
      continue;                        // nothing else to do
      }

     if( rng.first >= NumVar ) {  // all the Variable are deleted already
      auto nr = rng.second - rng.first;
      if( nr > to_add )
       throw( std::logic_error(
                  "BundleSolver::process_outstanding_Modification: "
                  "removing non-existing Variable" ) );
      to_add -= nr;               // "virtually" remove them
      continue;                   // nothing else to do
      }
     if( rng.second >= NumVar ) {  // some of the Variable deleted already
      auto nr = rng.second - NumVar;
      if( nr > to_add )
       throw( std::logic_error(
                  "BundleSolver::process_outstanding_Modification: "
                  "removing non-existing Variable" ) );
      to_add -= nr;               // "virtually" remove them
      rng.second = NumVar;
      }
     if( rng.second < NumVar ) {
      // if deleting the last range of Variable nothing has to be done,
      // but deleting Variable "in the middle" rather requires moving
      // down the remaining range of values in Lambda
      std::copy( Lambda.begin() + rng.second , Lambda.end() ,
                 Lambda.begin() + rng.first );
      }
     Subset tdlt( rng.second - rng.first );
     NumVar -= tdlt.size();
     Lambda.resize( NumVar );  // adjust Lambda
     Lambda1.resize( NumVar );
     if( MaxSol > 1 )
      LmbdBst.resize( NumVar );
     std::iota( tdlt.begin() , tdlt.end() , rng.first );
     if( MasterPB )
      MasterPB->remove_vars( reinterpret_cast< const int * >( tdlt.data() ) ,
                             int( tdlt.size() ) );  // remove from MP
     continue;

     }  // end( if( ttmod == FunctionModVarsRngd ) )

    if( const auto ttmod =
        std::dynamic_pointer_cast< FunctionModVarsSbst >( tmod ) ) {
     if( f_sparse_lambda ) {
      // Sparse Lambda Sbst: subset() lists LOCAL indices in v_c05f[ h ];
      // mirror the Rngd path — drop them from v_local2global[ h ],
      // decrement refcounts, queue globally-dead slots for compaction.
      const auto h = get_index_of_component( ttmod->function() );
      auto & m = v_local2global[ h ];
      const Index loc_NV = m.size() - 1;  // exclude trailing Inf< Index >()
      const auto & sbst = ttmod->subset();

      if( sbst.empty() ) {
       // deleting *all* of v_c05f[ h ]'s local Lambda
       for( Index l = 0 ; l < loc_NV ; ++l ) {
        const Index g = m[ l ];
        if( std::abs( Lambda[ g ] ) > 1e-12 )
         std::fill( AlphaC.begin() , AlphaC.end() , true );
        if( --v_ref_count[ g ] == 0 )
         globally_to_remove.push_back( g );
        }
       m.assign( { Inf< Index >() } );  // keep only the terminator
       rmvd_vars = true;
       continue;
       }

      // bounded subset removal: drop subset() ∩ [ 0 , loc_NV ) from m
      Subset effective;
      effective.reserve( sbst.size() );
      for( auto l : sbst )
       if( l < loc_NV )
        effective.push_back( l );
      if( effective.empty() ) {
       rmvd_vars = true;
       continue;
       }
      for( auto l : effective ) {
       const Index g = m[ l ];
       if( std::abs( Lambda[ g ] ) > 1e-12 )
        std::fill( AlphaC.begin() , AlphaC.end() , true );
       if( --v_ref_count[ g ] == 0 )
        globally_to_remove.push_back( g );
       }
      // erase effective[] from m in descending order so earlier erases
      // don't invalidate later positions (effective is ordered ascending
      // per the FunctionModVarsSbst contract).
      for( auto it = effective.rbegin() ; it != effective.rend() ; ++it )
       m.erase( m.begin() + *it );
      rmvd_vars = true;
      continue;
      }

     rmvd_vars = true;

     if( ttmod->subset().empty() ) {  // deleting *all* Variable
      // if any Variable is nonzero, the linearization errors
      // will have to be recomputed for all components
      for( auto & el : Lambda )
       if( std::abs( el ) > 1e-12 ) {
        std::fill( AlphaC.begin() , AlphaC.end() , true );
        break;
        }

      NumVar = 0;
      Lambda.clear();
      Lambda1.clear();
      LmbdBst.clear();
      if( MasterPB )
       MasterPB->remove_vars( nullptr , 0 );  // remove from MP
      continue;                        // nothing else to do
      }

     if( ttmod->subset().front() >= NumVar ) {
      // all the Variable are deleted already
      if( ttmod->subset().size() > to_add )
       throw( std::logic_error(
                  "BundleSolver::process_outstanding_Modification: "
                  "removing non-existing Variable" ) );
      to_add -= ttmod->subset().size();  // "virtually" remove them
      continue;                          // nothing else to do
      }

     // if any of the deleted Variable is nonzero, the linearization errors
     // will have to be recomputed for all components
     for( auto el : ttmod->subset() ) {
      if( el >= NumVar )
       break;

      if( std::abs( Lambda[ el ] ) > 1e-12 ) {
       std::fill( AlphaC.begin() , AlphaC.end() , true );
       break;
       }
      }

     Subset tsbst;
     c_Subset * sbst = & tsbst;
     if( ttmod->subset().back() < NumVar )  // no Variable deleted already
      sbst = & ttmod->subset();             // delete them all
     else {                                 // construct the subset to delete
      // walk back from the end until the first index that already
      // refers to an existing Variable (< NumVar)
      auto sbstit = ttmod->subset().end();
      while( *(--sbstit) >= NumVar ) { }
      tsbst = Subset( ttmod->subset().begin() , ++sbstit );
      auto nr = ttmod->subset().size() - tsbst.size();
      if( nr > to_add )
       throw( std::logic_error(
                  "BundleSolver::process_outstanding_Modification: "
                  "removing non-existing Variable" ) );
      to_add -= nr;               // "virtually" remove them
      }

     Compact( Lambda , *sbst );  // adjust Lambda
     NumVar -= sbst->size();
     Lambda.resize( NumVar );
     Lambda1.resize( NumVar );
     if( MaxSol > 1 )
      LmbdBst.resize( NumVar );
     if( MasterPB )
      MasterPB->remove_vars( reinterpret_cast< const int * >( sbst->data() ) ,
                             int( sbst->size() ) );  // remove from MP
     continue;

     }  // end( if( ttmod == FunctionModVarsSbst ) )

    // if control reaches here, this is an unknown *FunctionModVars* (??)
    throw( std::logic_error(
               "BundleSolver::process_outstanding_Modification: "
               "unknown FunctionModVars" ) );

    }  // end( if( tmod == FunctionModVars ) )
   }  // end FunctionModVars
  }  // end( 4th loop, forward )

 // sparse-mode compaction: any LamVcblr slot whose refcount fell to 0
 // during the 4th loop is now truly dead and can be reclaimed.
 if( f_sparse_lambda && ( ! globally_to_remove.empty() ) ) {
  // sort + dedup (defensive — a slot can in principle be queued multiple
  // times if two h's removed their last reference)
  std::sort( globally_to_remove.begin() , globally_to_remove.end() );
  globally_to_remove.erase( std::unique( globally_to_remove.begin() ,
                                          globally_to_remove.end() ) ,
                            globally_to_remove.end() );

  // compact LamVcblr / v_ref_count in place; LamVcblr may have been
  // extended past NumVar by sparse FunctionModVarsAddd, those extra
  // entries (new vars, with positive refcount) survive and shift down
  Index w = 0;
  Index g_pos = 0;
  for( Index i = 0 ; i < LamVcblr.size() ; ++i ) {
   if( ( g_pos < globally_to_remove.size() ) &&
       ( globally_to_remove[ g_pos ] == i ) ) {
    ++g_pos;
    continue;
    }
   if( w != i ) {
    LamVcblr[ w ]    = LamVcblr[ i ];
    v_ref_count[ w ] = v_ref_count[ i ];
    }
   ++w;
   }
  LamVcblr.resize( w );
  v_ref_count.resize( w );

  // compact Lambda / Lambda1 / LmbdBst (sized to old NumVar pre-Adds;
  // only positions [ 0 , NumVar ) are affected — the new vars from
  // sparse Adds aren't in Lambda yet, they'll be appended by the
  // post-loop add_vars). globally_to_remove entries are all in
  // [ 0 , NumVar ) by construction.
  Index lw = 0;
  Index lg = 0;
  for( Index i = 0 ; i < NumVar ; ++i ) {
   if( ( lg < globally_to_remove.size() ) &&
       ( globally_to_remove[ lg ] == i ) ) {
    ++lg;
    continue;
    }
   if( lw != i ) {
    Lambda[ lw ]   = Lambda[ i ];
    Lambda1[ lw ]  = Lambda1[ i ];
    if( MaxSol > 1 )
     LmbdBst[ lw ] = LmbdBst[ i ];
    }
   ++lw;
   }
  NumVar = lw;
  Lambda.resize( NumVar );
  Lambda1.resize( NumVar );
  if( MaxSol > 1 )
   LmbdBst.resize( NumVar );

  // translate v_local2global[ h ] entries: every surviving global g
  // shifts down by the number of removed entries strictly below it
  for( auto & vmap : v_local2global )
   for( auto & e : vmap ) {
    if( e == Inf< Index >() )
     continue;
    const auto cnt = std::distance(
                       globally_to_remove.begin() ,
                       std::lower_bound( globally_to_remove.begin() ,
                                          globally_to_remove.end() , e ) );
    e -= cnt;
    }

  // rebuild Lambda2Idx with the new compacted indices
  Lambda2Idx.clear();
  for( Index i = 0 ; i < LamVcblr.size() ; ++i )
   Lambda2Idx.emplace( LamVcblr[ i ] , i );

  // tell the Master to drop the old global names (the Master still sees
  // the pre-compaction index space)
  if( MasterPB )
   MasterPB->remove_vars(
              reinterpret_cast< const int * >( globally_to_remove.data() ) ,
              int( globally_to_remove.size() ) );

  globally_to_remove.clear();
  }

 // at this point, the set of Variable in the BundleSolver/Master Problem
 // coincides with the set of Variable in the C05Function(s), save for the
 // Variable to be added: in other words, the positions from 0 no NumVar - 1
 // in the linearizations corresponds to what BundleSolver expects

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 5th loop: handle "horizontal" changes, i.e., changes of a given range
 // (subset) of entries in all the linearizations, i.e., C05FunctionMod* and
 // C05FunctionModLin* with which().empty(). note that due to limitations
 // of the MPSolver interface, subsets are anyway translated to a range,
 // thereby possibly requiring also entries that have not changed. if
 // Variable have been removed the original indices in the Modification need
 // be translated, and in fact if Variable in the range have been removed
 // the range can be shrank, up to disappearing altogether. the mapping is
 // more difficult if Variable have also been added. however we can exploit
 // the property that if a Variable ever has a name that is larger than its
 // index in the Modification, this can only mean that the Variable has been
 // removed and then re-added after that the Modification has been issued.
 // these Variable can therefore be ignored, since (unless they have been
 // re-removed) they will be added in the end
 //
 // this is the final loop, so the list must be empty at the end

 for( ; ! v_mod_tmp.empty() ; v_mod_tmp.pop_front() ) {
  auto mod = v_mod_tmp.front();  // pick the first Modification

  Index wFi;                     // the affected component
  Range range( NumVar , 0 );     // an empty range
  c_Subset * subset = nullptr;   // an empty subset
  c_Vec_p_Var * vars;            // the affected Variable

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // a C05FunctionModRngd, that at this point can only have which().empty()
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModRngd >( mod ) ) {
   if( ! tmod->which().empty() )
    throw( std::logic_error(
               "BundleSolver::process_outstanding_Modification: "
               "unexpected nonempty C05FunctionModRngd" ) );

   wFi = get_index_of_component( tmod->function() );
   vars = & tmod->vars();
   range = tmod->range();
   goto done;
   }

  // a C05FunctionModSbst, that at this point can only have which().empty()
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModSbst >( mod ) ) {
   if( ! tmod->which().empty() )
    throw( std::logic_error(
               "BundleSolver::process_outstanding_Modification: "
               "unexpected nonempty C05FunctionModSbst" ) );

   wFi = get_index_of_component( tmod->function() );
   vars = & tmod->vars();
   subset = & tmod->subset();
   goto done;
   }

  // a C05FunctionModLinRngd implies that a specific range in all the
  // linearizations must be changed (by adding something)
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModLinRngd >( mod ) ) {
   wFi = get_index_of_component( tmod->function() );
   vars = & tmod->vars();
   range = tmod->range();
   goto done;
   }

  // a C05FunctionModLinSbst implies that a specific subset in all the
  // linearizations must be changed (by adding something)
  if( const auto tmod =
           std::dynamic_pointer_cast< C05FunctionModLinSbst >( mod ) ) {
   wFi = get_index_of_component( tmod->function() );
   vars = & tmod->vars();
   subset = & tmod->subset();
   goto done;
   }

  // it is neither of the above: this should not happen
  throw( std::logic_error(
             "BundleSolver::process_outstanding_Modification: "
             "unexpected Modification slipped in" ) );

  // the range/subset (and component) have been identified: check if the
  // need to be translated due to addition/removals, and in case do it
  done:if( ! rmvd_vars ) {
   // Variable have never been removed, hence the names can be used directly
   if( subset ) {  // turn the subset into a range
    range.first = subset->front();
    range.second = subset->back() + 1;
    }
   }
  else {
   // Variable have been removed, and hence names need be actualised
   // this is done by directly checking vars() against the "active"
   // Variable of v_c05f[ 0 ], which is fairly taken as a representative
   // since all the C05Function have the same "active" Variable
   if( ! addd_vars ) {
    // ... but never added: names can have only decreased, but even more
    // importantly must have remained ordered, i.e., the first "active"
    // Variable in vars() is the first variable of the range, the last
    // "active" Variable vars() is the last variable of the range
    // note that we do not use subset and range here, as the range is
    // reconstructed from scratch using vars
    auto lit = vars->begin();
    for( ; lit != vars->end() ; ++lit ) {
     range.first = v_c05f[ 0 ]->is_active( *lit );
     if( range.first < v_c05f[ 0 ]->get_num_active_var() )
      break;
     }
    if( lit == vars->end() )  // no Variable in vars is still "active"
     continue;                // nothing else to do
    // since we know that here are some "active" Variable in vars(), this
    // second loop will necessarily end
     for( auto rit = vars->rbegin() ; ; ++rit ) {
      range.second = v_c05f[ 0 ]->is_active( *lit );
      if( range.second < v_c05f[ 0 ]->get_num_active_var() )
       break;
      }
     ++range.second;  // the range is [ first , second )
     }
   else {
    // the complicated case: Variable have both been removed and added
    // names can have changed in an almost arbitrary way, except that if
    // a name has increased then the Variable has been deleted and re-added
    // and therefore need not be included
    Subset newnames( vars->size() );
    auto lit = vars->begin();
    auto nni = newnames.begin();
    if( subset ) {
     auto sit = subset->begin();
     for( ; lit != vars->end() ; ++lit , ++sit ) {
      auto i = v_c05f[ 0 ]->is_active( *lit );
      if( ( i <= *sit ) && ( i < v_c05f[ 0 ]->get_num_active_var() ) )
       *(nni++) = i;
      }
     }
    else {
     for( ; lit != vars->end() ; ++lit , ++range.first ) {
      auto i = v_c05f[ 0 ]->is_active( *lit );
      if( ( i <= range.first ) && ( i < v_c05f[ 0 ]->get_num_active_var() ) )
       *(nni++) = i;
      }
     }
    if( nni == newnames.begin() )  // no Variable in vars is still "active"
     continue;                     // nothing else to do
    newnames.resize( std::distance( newnames.begin() , nni ) );
    range.first = newnames.front();
    range.second = newnames.back();
    ++range.second;  // the range is [ first , second )
    }
   }

  // now actually do it
  if( MasterPB )
   MasterPB->invalidate_subgradients( int( wFi ) );

  }  // end( 5th loop, forward )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are Variable to add, do it now in one blow.
 // There is a trade-off here: doing this now causes the master AddVars
 // pass to (indirectly) call get_linearization_coefficients() on a
 // smaller set of linearizations if additions are done, but on the
 // other hand increases NumVar and therefore the work done in later
 // stages. Hence, this is done here only if no additions are done.

 bool toadd = std::find_if( Addd.begin() , Addd.end() ,
                            []( Subset & Ak ) { return( ! Ak.empty() ); }
                            ) != Addd.end();
 if( to_add && ( ! toadd ) ) {
  NumVar += to_add;
  Lambda.resize( NumVar , 0 );
  Lambda1.resize( NumVar , 0 );
  if( MaxSol > 1 )
   LmbdBst.resize( NumVar , 0 );
  if( MasterPB )
   MasterPB->add_vars( int( to_add ) );
  to_add = 0;  // done already
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are linearization to add/change, do it now in one blow
 //
 // note that due to limitations in the MPSolver interface, changing a
 // linearization is identical to adding one, even if the change was limited
 // to a range/subset of the entries
 //
 // yet, handling of additions and changes differs depending on BPar7.
 // in particular, if ( BPar7 & 4 ), then additions to the global pool
 // imply additions to the bundle, otherwise just marking that the
 // linearization exists
 //
 // note that if ( BPar7 & 3 ) >= 2, BundleSolver will happily delete from
 // the global pool any linearization it deletes from the bundle; if also
 // ! ( BPar7 & 4 ), one could therefore think it appropriate to delete
 // from the global pool any linearization that is added. however, we do not
 // do that, just refraining to add them to the bundle. these will likely be
 // overwritten during the optimization, and if memory is a problem then the
 // global pools should just be sized accordingly

 if( toadd ||
     ( std::find_if( Chgd.begin() , Chgd.end() ,
                     []( Subset & Ck ) { return( ! Ck.empty() ); }
                     ) != Chgd.end() ) ) {
  // at least a component has had linearizations added or changed

  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( Addd[ k ].empty() && Chgd[ k ].empty() )  // nothing to see here
    continue;                                    // move on

   if( Addd[ k ].empty() && ( Chgd[ k ].size() >= NrItems[ k ] ) ) {
    // all existing items change and no new one is added
    if( ! reset[ k ] )
     ++cntreset;
    reset[ k ] = AlphaC[ k ] = true;             // this is a reset
    continue;
    }

   // first, cleanup Chgd[ k ] from linearizations not in the bundle
   if( ! Chgd[ k ].empty() ) {
    Subset tmp;
    for( auto i : Chgd[ k ] )
     if( InvItemVcblr[ k ][ i ] >= vBPar2.back() )
      tmp.push_back( i );

    if( ! tmp.empty() ) {
     if( tmp.size() >= Chgd[ k ].size() ) {
      Chgd[ k ].clear();
      if( Addd[ k ].empty() )  // nothing more to see here
       continue;               // move off
      }
     else
      set_difference_in_place( Chgd[ k ] , tmp );
     }
    }

   // now, if ! ( BPar7 & 4 ), also cleanup Addd[ k ] from linearizations
   // not in the bundle, but in doing so mark them into InvItemVcblr unless
   // ( BPar7 & 3 ) == 3, in which case BundleSolver does not care if there
   // are existing linearizations because it anyway freely overwrites them
   // note that for linearizations in the bundle, being in Addd[ k ] is the
   // same as being in Chgd[ k ]: the linearization has changed, and the
   // master problem must be changed to reflect this
   if( ( ! ( BPar7 & 4 ) ) && ( ! Addd[ k ].empty() ) ) {
    Subset tmp;
    for( auto i : Addd[ k ] )
     if( InvItemVcblr[ k ][ i ] >= vBPar2.back() ) {
      InvItemVcblr[ k ][ i ] = ( BPar7 & 3 ) == 3 ? InINF : vBPar2.back();
      tmp.push_back( i );
      }

    if( ! tmp.empty() ) {
     if( tmp.size() >= Addd[ k ].size() ) {
      Addd[ k ].clear();
      if( Chgd[ k ].empty() )  // nothing more to see here
       continue;               // move off
      }
     else
      set_difference_in_place( Addd[ k ] , tmp );
     }
    }

   // compute the union between Addd[ k ] and Chgd[ k ] into Addd[ k ]
   set_union_in_place( Addd[ k ] , std::move( Chgd[ k ] ) );

   // finally, add the resulting stuff to the bundle
   for( auto i : Addd[ k ] )
    add_to_bundle( k , i );
   }
  }  // end( if( additions or changes ) )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if some component need be reset, reset the linearizations: since
 // reset[ k ] ==> AlphaC[ k ], later on also the constants will be reset

 // MasterPB owns copies of the bundle rows, so merely invalidating its Solver
 // would reload the same stale rows. Always copy reset linearizations back
 // from the component global pools, also when every hard component is reset.
 if( cntreset )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( reset[ k ] )  // reset[ k ] ==> ! IsEasy[ k ]
    reload_component_bundle( k );

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are constants to change entirely, do it now in one blow
 // note: in case of a full reset, get_linearization_coefficients() is
 //       called twice, once in ChgSubG() (via GetGi()) and once in the
 //       loop below. this has the potential to be horribly inefficient,
 //       but the only clean way out is to do away with MPSolver entirely

 if( std::find( AlphaC.begin() , AlphaC.end() , false ) == AlphaC.end() ) {
  for( auto & cchk : Cchg )  // all components have been reset
   cchk.clear();             // no need to change them individually

  ResetAlfa( NrFi );
  }
 else
  for( Index k = 0 ; k < NrFi ; ++k )
   if( AlphaC[ k ] ) {  // all the constants of this component are reset
    Cchg[ k ].clear();  // no need to change them individually

    ResetAlfa( k );
    }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there subsets of Alphas to change, do it now

 if( std::find_if( Cchg.begin() , Cchg.end() ,
                   []( Subset & Ck ) { return( ! Ck.empty() ); }
                   ) != Cchg.end() ) {
  for( Index k = 0 ; k < NrFi ; ++k ) {
   const Index loc_NV = f_sparse_lambda ? v_c05f[ k ]->get_num_active_var()
                                        : NumVar;
   std::vector< VarValue > Gi( loc_NV );
   for( auto i : Cchg[ k ] )
    if( InvItemVcblr[ k ][ i ] < vBPar2.back() ) {
     auto Ai = rs( v_c05f[ k ]->get_linearization_constant( i ) );

     #ifndef NDEBUG
      if( std::isnan( Ai ) )  // linearization no longer valid
       throw( std::logic_error(
                  "BundleSolver::CheckBundle: "
                  "inexistent linearization" ) );
     #endif

     // recover the (possibly changed) subgradient of linearization i: an
     // AllEntriesChanged / AllLinearizationChanged oracle Modification is
     // routed here as a constant change, but it may also have moved g (the
     // PolyhedralFunction component case: a "modify rows" alters the active
     // piece). Re-installing only the constant via modify_alpha would leave
     // the master carrying a stale g and hence a cut that no longer
     // under-estimates the new function: the dual master would then pick a
     // direction along which the true Fi worsens. So push both g and the
     // constant through modify_cut. chgsign matches the convex-min space of
     // get_linearization_coefficients for concave problems
     v_c05f[ k ]->get_linearization_coefficients( Gi.data() ,
                                                  Range( 0 , loc_NV ) , i );
     if( ! f_convex )
      chgsign( Gi.data() , loc_NV );

     if( MasterPB ) {
      // scatter the physical subgradient to global Variable slots; MPB will
      // convert it to the active primal/dual storage convention.
      std::vector< double > g_master( NumVar , 0.0 );
      if( f_sparse_lambda ) {
       const auto & m = v_local2global[ k ];
       for( Index li = 0 ; li < loc_NV ; ++li )
        g_master[ m[ li ] ] = Gi[ li ];
       }
      else
       for( Index j = 0 ; j < NumVar ; ++j )
        g_master[ j ] = Gi[ j ];
      MasterPB->modify_cut( int( k ) , int( InvItemVcblr[ k ][ i ] ) ,
                            std::move( g_master ) , Ai );
      }
    }
   }
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are (still) Variable to add, do it now in one blow


 if( to_add ) {
  NumVar += to_add;
  Lambda.resize( NumVar , 0 );
  Lambda1.resize( NumVar , 0 );
  if( MaxSol > 1 )
   LmbdBst.resize( NumVar , 0 );
  if( MasterPB )
   MasterPB->add_vars( int( to_add ) );
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // and this, finally, is all! (save possibly some checks)

 //!! PrintBundle();
 #if CHECK_DS & 1
  CheckBundle();
 #endif
 #if CHECK_DS & 4
  CheckAlpha();
 #endif

 }  // end( BundleSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

void BundleSolver::CheckBundle( void )
{
 std::ostream * wlog = ( ( ! f_log ) || ( LogVerb <= 1 ) ) ? & std::cerr
                                                           : f_log;
 // check vBPar2
 for( Index k = 0 ; k < NrFi ; ++k )
  if( Index( v_c05f[ k ]->get_int_par( C05Function::intGPMaxSz ) )
      != vBPar2[ k ] )
   *wlog << "size of global pool " << k << " does not match" << std::endl;

 // check ItemVcblr against InvItemVcblr and Master
 Subset tmp( NrFi , 0 );
 for( Index i = 0 ; i < get_max_name() ; ++i )
  if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) {
   ++tmp[ ItemVcblr[ i ].first ];
   if( wcomponent_global( i ) != ItemVcblr[ i ].first + 1 ) {
    *wlog << "position " << i << " in the bundle should be of component "
          << ItemVcblr[ i ].first << " but Master says ";
    if( wcomponent_global( i ) == InINF )
     *wlog << "empty" << std::endl;
    else
     *wlog << wcomponent_global( i ) - 1 << std::endl;
    }

   if( InvItemVcblr[ ItemVcblr[ i ].first ][ ItemVcblr[ i ].second ] != i )
    *wlog << "position " << i << " in the bundle shoud be linearization "
          << ItemVcblr[ i ].second << " of component "
          << ItemVcblr[ i ].first << " but InvItemVcblr disagrees"
          << std::endl;
   }
  else
   if( wcomponent_global( i ) < InINF )
    *wlog << "position " << i
          << " in the bundle should be empty but Master says "
          << wcomponent_global( i ) << std::endl;

 // check NrItems
 if( int diff = ( NrItems.back() - std::accumulate( NrItems.begin() ,
                                                    NrItems.begin() + NrFi ,
                                                    Index( 0 ) ) ) != 0 )
  *wlog << " NrItems[ NrFi ] = " << NrItems.back() << " but the sum is "
        <<  NrItems.back() + diff  << std::endl;

 for( Index k = 0 ; k < NrFi ; ++k )
  if( tmp[ k ] != NrItems[ k ] )
   *wlog << "counted " << tmp[ k ] << " items in the bundle for component "
         << k << " but NrItems says " << NrItems[ k ] << std::endl;

 // check InvItemVcblr against ItemVcblr and C05Function
 for( Index k = 0 ; k < NrFi ; ++k )
  for( Index i = 0 ; i < vBPar2[ k ] ; ++i ) {
   if( ( InvItemVcblr[ k ][ i ] < InINF ) &&
       ( ! v_c05f[ k ]->is_linearization_there( i ) ) )
    *wlog << "linearization " << i << " in pool " << k
          << " does not exist" << std::endl;

   if( ( InvItemVcblr[ k ][ i ] == InINF ) &&
       v_c05f[ k ]->is_linearization_there( i ) )
    *wlog << "linearization " << i << " in pool " << k
          << " unaccounted for" << std::endl;

   if( ( InvItemVcblr[ k ][ i ] < vBPar2.back() ) &&
       ( ( ItemVcblr[ InvItemVcblr[ k ][ i ] ].first != k ) ||
         ( ItemVcblr[ InvItemVcblr[ k ][ i ] ].second != i ) ) )
    *wlog << "linearization " << i << " in pool " << k
          << " should be in bundle in position "
          << InvItemVcblr[ k ][ i ] << " but ItemVcblr disagrees"
          << std::endl;

   if( ( InvItemVcblr[ k ][ i ] < ( ( BPar7 & 3 ) ? vBPar2.back() : InINF ) )
       && ( i >= MaxItem[ k ] ) )
    *wlog << "item in position " << i << " of pool " << k
          << " but MaxItem says " << MaxItem[ k ] << std::endl;

   if( ( InvItemVcblr[ k ][ i ] >= ( ( BPar7 & 3 ) ? vBPar2.back() : InINF ) )
       && ( i < FrFItem[ k ] ) )
     *wlog << "free item in position " << i << " of pool " << k
           << " but FrFItem says " << FrFItem[ k ] << std::endl;

   }  // end( for( i ) )

 // check FreList (if there is anything to check)
 if( ! FreList.empty() ) {
  // copy FreList to a vector to check it (this only works since the
  // underlying container is a std::vector< Index >) and the topmost
  // element in a heap is the first element of the vector
  tmp.resize( FreList.size() );
  std::copy( &( FreList.top() ) , &( FreList.top() ) + FreList.size() ,
             tmp.begin() );
  std::sort( tmp.begin() , tmp.end() );

  for( auto i : tmp )
   if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] )
    *wlog << "item " << i << " in FreList is not free" << std::endl;

  for( Index i = 0 ; i < get_max_name() ; ++i )
   if( ItemVcblr[ i ].second >= vBPar2[ ItemVcblr[ i ].first ] ) {
    auto it = std::lower_bound( tmp.begin() , tmp.end() , i );
    if( ( it == tmp.end() ) || ( *it != i ) )
     *wlog << "free item " << i << " not in FreList" << std::endl;
    }
  }
 }  // end( BundleSolver::CheckBundle )

/*--------------------------------------------------------------------------*/

void BundleSolver::CheckAlpha( void )
{
 std::ostream * wlog = ( ( ! f_log ) || ( LogVerb <= 1 ) ) ? & std::cerr
                                                           : f_log;
 *wlog << def;
 std::vector< VarValue > G( NumVar );
 const double eps = 1e-8;

 for( Index i = 0 ; i < get_max_name() ; ++i )
  if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) {
   v_c05f[ ItemVcblr[ i ].first ]->get_linearization_coefficients( G.data() ,
                                                        Range( 0 , NumVar ) ,
                                                     ItemVcblr[ i ].second );
   if( ! f_convex )
    chgsign( G.data() , NumVar );
   double lin_cst = rs( v_c05f[ ItemVcblr[ i ].first
                           ]->get_linearization_constant(
                                                   ItemVcblr[ i ].second ) );
   double dotLG = std::inner_product( Lambda.begin() , Lambda.end() ,
                                     G.begin() , double( 0 ) );
   double ref = UpRifFi[ ItemVcblr[ i ].first ];
   double tAi = ref - lin_cst - dotLG;

   if( std::abs( tAi - read_alpha_global( i ) ) >= eps *
       std::max( std::max( std::abs( tAi ) ,
                           std::abs( UpRifFi[ ItemVcblr[ i ].first ] ) ) ,
                 double( 1 ) ) )
    *wlog << std::endl << "Alfa[ " << i << " ]: F = " << tAi << " ~ M = "
          << read_alpha_global( i ) << " (F-M = " << shrt4 << ( tAi - read_alpha_global( i ) ) << def << ")"
          << " | k=" << ItemVcblr[ i ].first
          << " UpRifFi=" << shrt4 << ref
          << " lin_const=" << lin_cst << " <L,G>=" << dotLG << def;
    }

 }  // end( BundleSolver::CheckAlpha )

/*--------------------------------------------------------------------------*/

void BundleSolver::CheckLBs( void )
{
 // diagnostic cross-check between the global lower bounds stored in
 // LowerBound / inside the C05Functions and those used by MasterPB.
 // Currently only the global bound is checked; per-component bounds
 // need a MasterPB-side getter that does not yet exist.
 std::ostream * wlog = ( ( ! f_log ) || ( LogVerb <= 1 ) ) ? & std::cerr
                                                           : f_log;
 *wlog << def;

 auto GLB = f_convex ?   f_Block->get_valid_lower_bound( false )
                     : - f_Block->get_valid_upper_bound( false );

 if( TrueLB ) {  // a finite global lower bound is set
  if( LowerBound.back() == -INFshift ) {
   *wlog << std::endl << "TrueLB but no stored global bound";
   if( GLB > -INFshift )
    *wlog << std::endl << "finite global bound " << rs( GLB )
          << " available but not set";
   }
  }
 else {
  if( GLB > -INFshift )
   *wlog << std::endl << "finite global bound " << rs( GLB )
         << " available but not set";

  GLB =  f_convex ?   f_Block->get_valid_lower_bound( true )
                  : - f_Block->get_valid_upper_bound( true );
  if( GLB != LowerBound.back() )
   *wlog << std::endl << "conditional global bound = " << rs( GLB )
         << " != from stored = " << rs( LowerBound.back() );
  }

 }  // end( BundleSolver::CheckLBs )

/*--------------------------------------------------------------------------*/

void BundleSolver::PrintBundle( void )
{
 std::ostream * wlog = ( ( ! f_log ) || ( LogVerb <= 1 ) ) ? & std::cerr
                                                           : f_log;
 *wlog << def << std::endl << "Lambda = [ ";
 for( Index h = 0 ; h < NumVar - 1 ; ++h )
  *wlog << Lambda[ h ] << ", ";
 *wlog << Lambda.back() << " ]";

  std::vector< VarValue > G( NumVar );

 *wlog << std::endl;
 for( Index i = 0 ; i < get_max_name() ; ++i ) {
  *wlog << i << "\t";
  if( ItemVcblr[ i ].second >= vBPar2[ ItemVcblr[ i ].first ]
          || ItemVcblr[ i ].second < 0 ) {
   *wlog << "[empty]" << std::endl;
   continue;
   }

  auto wFi = ItemVcblr[ i ].first;
  auto j = ItemVcblr[ i ].second;
  *wlog << wFi << "\t" << j << "\t[ ";

  v_c05f[ wFi ]->get_linearization_coefficients( G.data() ,
                                                 Range( 0 , NumVar ) , j );
  if( ! f_convex )
   chgsign( G.data() , NumVar );

  for( Index h = 0 ; h < NumVar - 1 ; ++h )
   *wlog << G[ h ] << ", ";

  *wlog << G.back() << " ]\t"
         << rs( v_c05f[ wFi ]->get_linearization_constant( j ) )
         << "\t" << read_alpha_global( i ) << std::endl;
  }
 }

#endif

/*--------------------------------------------------------------------------*/
/*------------------ CLASS BundleSolverState --------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolverState::deserialize( const netCDF::NcGroup & group )
{
 auto nv = group.getDim( "BundleSolver_NumVar" );
 if( nv.isNull() )
  throw( std::logic_error(
                      "BundleSolverState::deserialize: missing NumVar" ) );

 NumVar = nv.getSize();

 auto l = group.getVar( "BundleSolver_Lambda" );
 if( l.isNull() )
  throw( std::logic_error(
                      "BundleSolverState::deserialize: missing Lambda" ) );

 Lambda.resize( nv.getSize() );
 l.getVar( Lambda.data() );

 auto tt = group.getVar( "BundleSolver_t" );
 if( tt.isNull() )
  throw( std::logic_error( "BundleSolverState::deserialize: missing t" ) );

 tt.getVar( &t );

 auto nf = group.getDim( "BundleSolver_NrFi" );
 if( nf.isNull() )
  throw( std::logic_error(
                      "BundleSolverState::deserialize: missing NrFi" ) );

 NrFi = nf.getSize();
 --NrFi;

 auto nup = group.getDim( "BundleSolver_UpFiLmbdef" );
 if( nup.isNull() ) {
  UpFiLmbdef = 0;
  UpFiLmb.clear();
  }
 else {
  UpFiLmbdef = nup.getSize();
  auto vup = group.getVar( "BundleSolver_UpFiLmb" );
  if( vup.isNull() )
   throw( std::logic_error(
                      "BundleSolverState::deserialize: missing UpFiLmb" ) );

  UpFiLmb.resize( NrFi + 1 );
  vup.getVar( UpFiLmb.data() );
  }

 auto nlw = group.getDim( "BundleSolver_LwFiLmbdef" );
 if( nlw.isNull() ) {
  LwFiLmbdef = 0;
  LwFiLmb.clear();
  }
 else {
  LwFiLmbdef = nlw.getSize();
  auto vup = group.getVar( "BundleSolver_LwFiLmb" );
  if( vup.isNull() )
   throw( std::logic_error(
                      "BundleSolverState::deserialize: missing LwFiLmb" ) );

  LwFiLmb.resize( NrFi + 1 );
  vup.getVar( LwFiLmb.data() );
  }

 auto f0 = group.getVar( "BundleSolver_Fi0Lmb" );
 if( f0.isNull() )
  Fi0Lmb = 0;
 else
  f0.getVar( &Fi0Lmb );

 auto lb = group.getVar( "BundleSolver_global_LB" );
 if( lb.isNull() )
  global_LB = -BundleSolver::INFshift;
 else
  lb.getVar( &global_LB );

 v_comp_State.resize( nf.getSize() , nullptr );

 for( BundleSolver::Index i = 0 ; i < v_comp_State.size() ; ++i ) {
  auto gi = group.getGroup( "Component_State_" + std::to_string( i ) );
  if( ! gi.isNull() )
    v_comp_State[ i ] = State::new_State( gi );
  }
 }  // end( BundleSolverState::deserialize )

/*--------------------------------------------------------------------------*/

void BundleSolverState::serialize( netCDF::NcGroup & group ) const
{
 // always call the method of the base class first
 State::serialize( group );

 auto nv = group.addDim( "BundleSolver_NumVar" , NumVar );

 ( group.addVar( "BundleSolver_Lambda" , netCDF::NcDouble() , nv ) ).putVar(
                                       { 0 } , { NumVar } , Lambda.data() );

 ( group.addVar( "BundleSolver_t" , netCDF::NcDouble() ) ).putVar( &t );

 auto nfi = group.addDim( "BundleSolver_NrFi" , NrFi + 1 );

 if( UpFiLmbdef ) {
  group.addDim( "BundleSolver_UpFiLmbdef" , UpFiLmbdef );
  ( group.addVar( "BundleSolver_UpFiLmb" , netCDF::NcDouble() , nfi )
    ).putVar( { 0 } , { NrFi + 1 } , UpFiLmb.data() );
  }

 if( LwFiLmbdef ) {
  group.addDim( "BundleSolver_LwFiLmbdef" , LwFiLmbdef );
  ( group.addVar( "BundleSolver_LwFiLmb" , netCDF::NcDouble() , nfi )
    ).putVar( { 0 } , { NrFi + 1 } , LwFiLmb.data() );
  }

 if( Fi0Lmb != 0 )
  ( group.addVar( "BundleSolver_Fi0Lmb" , netCDF::NcDouble() )
    ).putVar( & Fi0Lmb );

 if( global_LB > - BundleSolver::INFshift )
  ( group.addVar( "BundleSolver_global_LB" , netCDF::NcDouble() )
    ).putVar( & global_LB );

 for( Index i = 0 ; i < NrFi ; ++i ) {
  if( ! v_comp_State[ i ] )
   continue;

  auto gi = group.addGroup( "Component_State_" + std::to_string( i ) );
  v_comp_State[ i ]->serialize( gi );
  }
 }  // end( BundleSolverState::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- End File BundleSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
