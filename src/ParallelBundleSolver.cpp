/*--------------------------------------------------------------------------*/
/*-------------------- File ParallelBundleSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ParallelBundleSolver class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ParallelBundleSolver.h"

#include <deque>

#include <iomanip>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define VERBOSE_LOG 1

#if VERBOSE_LOG
 #define BLOG( l , x ) if( f_log && ( LogVerb > l ) ) *f_log << x

 #define BLOG2( l , c , x ) if( f_log && ( LogVerb > l ) && c ) *f_log << x
#else
 #define BLOG( l , x )

 #define BLOG2( l , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
using std::ios;
using std::setprecision;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto InINF = SMSpp_di_unipi_it::Inf< Block::Index >();

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/
// set precision for long floats (10 digits)

static inline std::ostream & def( std::ostream & os ) {
 os.setf( ios::scientific, ios::floatfield );
 os << setprecision( 10 );
 return( os );
 }

/*--------------------------------------------------------------------------*/
// set precision for short floats (2 digits)

static inline std::ostream & shrt( std::ostream & os ) {
 os.setf( ios::scientific, ios::floatfield );
 os << setprecision( 2 );
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
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ParallelBundleSolver to the Solver factory
SMSpp_insert_in_factory_cpp_0( ParallelBundleSolver );

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF ParallelBundleSolver -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- OTHER PROTECTED METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

BundleSolver::Index ParallelBundleSolver::InnerLoop( bool extrastep )
{
 /* The inner loop is divided in three phases:
  *
  * - in the ramp-up phase, min( MaxThread , NrFi - NrEasy ) tasks are
  *   started, each one compute()-ing a different component;
  *
  * - in the cruise phase, existing tasks are checked every PoolingInt
  *   seconds; if one (or more) is found ready it is processed by extracting
  *   all the relevant information, and it is then substituted by another
  *   for a different component;
  *
  * - in the ramp-down phases, all existing tasks not "consumed" already are
  *   checked every PoolingInt seconds; if one (or more) is found ready it is
  *   processed by extracting all the relevant information, but no other task
  *   takes its place, so that eventually the process ends.
  *
  * Note that gathering function values and linearizations from the evaluated
  * components is done in the main thread, and therefore it is a part of the
  * sequential bottleneck; however, this is required since the master problem
  * and all the other BundleSolver data structures are not protected from
  * concurrent access. */
 
 // if there is nothing to parallelize, call the base class version - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ( NrFi == 1 ) || ( MaxThread == 0 ) )
  return( BundleSolver::InnerLoop( extrastep ) );

 // if a deterministic formulation is selected, dispatch to it - - - - - - - -
 // bit 0 of ParFrm picks fixed-order consumption; bits 1/2 select the
 // discard (faithful) and batch sub-behaviours
 if( ParFrm & 1 )
  return( InnerLoopOrdered( extrastep , ParFrm & 2 , ParFrm & 4 ) );

 // compute the minimum number of components to evaluate
 Index minceval = ( MinNrEvls >= 0 ? Index( MinNrEvls )
		                   : ( NrFi - NrEasy ) * ( - MinNrEvls ) );
 Index ceval = 0;  // how many components have been evaluated so far

 // define the vector of std::future
 using EvalEl = std::pair< Index , std::future< int > >;

 std::vector< EvalEl > EvalV( std::min( MaxThread , NrFi - NrEasy ) );

 // ramp-up phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // fill-in the vector of std::future
 // note: each time we start the evaluation of a component we provisionally
 // set FiStatus == kOK so that FindNext() will not produce it multiple
 // times before its computation is actually over. when it is, its FiStatus
 // may become, say, kStopTime or kStopIter and the component may be again
 // eligible to be re-evaluated

 for( auto & el : EvalV ) {
  if( ! FindNext() )
   throw( std::logic_error( "no component to evaluate in ramp-up phase" ) );

  BLOG( 6 , std::endl << "ramp-up: component " << f_wFi << " in position "
	    << ( & el - & EvalV.front() ) );

  el.first = f_wFi;
  if( extrastep )
   SetupFiLambda( f_wFi );
  else
   SetupFiLambda1( f_wFi );
  el.second = v_c05f[ f_wFi ]->compute_async(
					  ( FiStatus[ f_wFi ] == kUnEval ) );
  FiStatus[ f_wFi ] = kOK;
  }

 // cruise phase- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // this is defined outside so that we can see what the last one was
 std::vector< EvalEl >::iterator it;

 bool insrtd = false;
 for( ; ; ) {
  // check if any future is ready - - - - - - - - - - - - - - - - - - - - - -
  // note that we assume without checking that all future are valid()
  it = std::find_if( EvalV.begin() , EvalV.end() ,
		     []( auto & el ) { return( el.second.wait_for(
				       std::chrono::duration< int >::zero() )
				       == std::future_status::ready ); } );

  // if not, sleep over and retry later on- - - - - - - - - - - - - - - - - -
  if( it == EvalV.end() ) {
   std::this_thread::sleep_for( std::chrono::duration< double >(
							      PoolingInt ) );
   continue;
   }
  
  // if one future is ready, read it- - - - - - - - - - - - - - - - - - - - -
  Index wFi = it->first;
  if( ! CurrNrEvls[ wFi ] )  // not evaluated before
   ++ceval;                  // one more evaluated
  ++CurrNrEvls[ wFi ];       // evaluated once more

  FiStatus[ wFi ] = it->second.get();  // get() the status of compute()

  BLOG( 6 , std::endl << "cruise: component " << wFi << " in position "
	    << it - EvalV.begin() << " has status " << FiStatus[ wFi ] );

  // if an unrecoverable error happens, immediately start the ramp-down - - -
  if( ( FiStatus[ wFi ] <= kUnEval ) || ( FiStatus[ wFi ] >= kError ) ) {
   BLOG( 3 , std::endl << "            Component " << wFi
	     << " evaluated: Error" );
   Result = kError;
   break;
   }

  // collect function values (upper and lower bound)- - - - - - - - - - - - -
  auto fwFi = v_c05f[ wFi ];
  auto ue = fwFi->get_upper_estimate();
  auto le = fwFi->get_lower_estimate();

  #if VERBOSE_LOG
   if( f_log && ( LogVerb > 3 ) ) {
    *f_log << std::endl << "            Component " << wFi
	   << " evaluated: UB = "<< def;
    pval( *f_log , ue );
    *f_log << ", LB = ";
    pval( *f_log , le );
    }
  #endif

  // if any component evaluates to -INF, then the whole problem evaluates to
  // -INF and it is therefore unbounded below; due to convexity, it is "very
  // seriously unbounded" in that a convex function being -INF anywhere is
  // -INF everywhere, hence if this ever happens it will do it "very soon"
  // (the very first time the offending component is evaluated); immediately
  // start the ramp-down- - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( f_convex && ( ue == -INFshift ) ) ||
      ( ( ! f_convex ) && ( le == INFshift ) ) ) {
   UpFiLmb1[ wFi ] = UpFiLmb1.back() = -INFshift;
   break;
   }

  if( extrastep ) {
   // if extrastep == true the method is actually being called on Lambda,
   // hence it is Lambda's estimates that need be updated, not Lambda1's
   update_UpFiLambd( wFi , f_convex ? ue : - le );
   update_LwFiLambd( wFi , f_convex ? le : - ue );

   // furthermore one immediately goes to put in the new task
   goto PutInNewTask;
   }
  
  // update UpFiLambd1[ wFi ] (and possibly UpFiLambd1[ NrFi ])
  update_UpFiLambd1( wFi , f_convex ? ue : - le );

  // if bit 4 of TrgtMng == 1, then compute the upper bound in Lambda
  // provided by the upper bound in Lambda1 and try to update UpFiLmb[ wFi ]
  // (and possibly UpFiLambd1[ NrFi ])
  // note that, even if this succeeds and therefore decreases UpFiLmb[ wFi ]
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

  // get new linearizations - - - - - - - - - - - - - - - - - - - - - - - - -
  if( GetGi( wFi ) )
   insrtd = true;

  // check if the accrued information changes the MP- - - - - - - - - - - - -
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

  // check if we can/must wind down - - - - - - - - - - - - - - - - - - - - -
  if( ( MaxTime < INFshift ) && ( get_elapsed_time() > MaxTime ) ) {
   Result = kStopTime;     // time has ran up
   break;                  // start ramp-down
   }

  // the MP is guaranteed to change and enough components evaluated
  if( MPchgs && ( ceval >= minceval ) )
   break;                  // start ramp-down

  // run a new task in the same position- - - - - - - - - - - - - - - - - - -
  PutInNewTask:

  if( ! FindNext() )       // find next component
   break;                  // if none, start ramp-down

  it->first = f_wFi;
  if( extrastep )
   SetupFiLambda( f_wFi );
  else
   SetupFiLambda1( f_wFi );
  it->second = v_c05f[ f_wFi ]->compute_async(
				          ( FiStatus[ f_wFi ] == kUnEval ) );
  FiStatus[ f_wFi ] = kOK;

  BLOG( 6 , std::endl << "cruise: in component " << f_wFi );

  }  // end( for( cruise phase loop ) )

 BLOG( 6 , std::endl << "end of cruise" );

 it->first = InINF;  // mark the entry as invalid
 
 // ramp-down phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index cnt = EvalV.size() - 1 ; cnt ; ) {
  // check if any future is both valid() and ready- - - - - - - - - - - - - -
  it = std::find_if( EvalV.begin() , EvalV.end() ,
		     []( auto & el ) { return( ( el.first != InINF ) &&
					       el.second.wait_for(
				       std::chrono::duration< int >::zero() )
				       == std::future_status::ready ); } );

  // if not, sleep over and retry later on- - - - - - - - - - - - - - - - - -
  if( it == EvalV.end() ) {
   std::this_thread::sleep_for( std::chrono::duration< double >(
							      PoolingInt ) );
   continue;
   }

  // if one future is ready, read it- - - - - - - - - - - - - - - - - - - - -
  Index wFi = it->first;
  it->first = InINF;
  if( ! CurrNrEvls[ wFi ] )  // not evaluated before
   ++ceval;                  // one more evaluated
  ++CurrNrEvls[ wFi ];       // evaluated once more
  --cnt;                     // one std::future less to wait for

  FiStatus[ wFi ] = it->second.get();  // get() the status of compute()

  BLOG( 6 , std::endl << "ramp-down: component " << wFi << " in position "
 	    << it - EvalV.begin() << " has status " << FiStatus[ wFi ] );

  // if an unrecoverable error happens, do nothing else - - - - - - - - - - -
  if( ( FiStatus[ wFi ] <= kUnEval ) || ( FiStatus[ wFi ] >= kError ) ) {
   if( f_log && ( LogVerb > 3 ) )
    *f_log << std::endl << "            Component " << wFi
	   << " evaluated: Error";
   Result = kError;
   continue;
   }

  // collect function values (upper and lower bound)- - - - - - - - - - - - -
  auto fwFi = v_c05f[ wFi ];
  auto ue = fwFi->get_upper_estimate();
  auto le = fwFi->get_lower_estimate();

  #if VERBOSE_LOG
   if( f_log && ( LogVerb > 3 ) ) {
    *f_log << std::endl << "            Component " << wFi
	   << " evaluated: UB = "<< def;
    pval( *f_log , ue );
    *f_log << ", LB = ";
    pval( *f_log , le );
    }
  #endif

  // if any component evaluates to -INF, then the whole problem evaluates to
  // -INF and it is therefore unbounded below; due to convexity, it is "very
  // seriously unbounded" in that a convex function being -INF anywhere is
  // -INF everywhere, hence if this ever happens it will do it "very soon"
  // (the very first time the offending component is evaluated); immediately
  // start the ramp-down- - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( f_convex && ( ue == -INFshift ) ) ||
      ( ( ! f_convex ) && ( le == INFshift ) ) ) {
   UpFiLmb1[ wFi ] = UpFiLmb1.back() = -INFshift;
   continue;
   }

  if( extrastep ) {
   // if extrastep == true the method is actually being called on Lambda,
   // hence it is Lambda's estimates that need be updated, not Lambda1's
   update_UpFiLambd( wFi , f_convex ? ue : - le );
   update_LwFiLambd( wFi , f_convex ? le : - ue );

   continue;  // amd there is nothing left to do
   }
 
  // update UpFiLambd1[ wFi ] (and possibly UpFiLambd1[ NrFi ])
  update_UpFiLambd1( wFi , f_convex ? ue : - le );

  // if bit 4 of TrgtMng == 1, then compute the upper bound in Lambda
  // provided by the upper bound in Lambda1 and try to update UpFiLmb[ wFi ]
  // (and possibly UpFiLambd1[ NrFi ])
  // note that, even if this succeeds and therefore decreases UpFiLmb[ wFi ]
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

  // if an unrecoverable error had happened previously, or the problem had
  // already been found unbounded below, do nothing else
  if( ( Result == kError ) || ( UpFiLmb1.back() == -INFshift ) )
   continue;
  
  // get new linearizations - - - - - - - - - - - - - - - - - - - - - - - - -
  if( GetGi( wFi ) )
   insrtd = true;

  // check if the accrued information changes the MP- - - - - - - - - - - - -
  // see the cruise phase for detailed comments

  if( ( ! MPchgs ) && ( UpFiLmb1.back() < UpTrgt ) )
   MPchgs = 1;

  if( ( ! MPchgs ) && insrtd && RifeqFi && ( LwFiLmb1.back() > LwTrgt ) )
   MPchgs = 1;

  }  // end( for( ramp-down phase loop ) )

 return( ceval );

 }  // end( ParallelBundleSolver::InnerLoop )

/*--------------------------------------------------------------------------*/

BundleSolver::Index ParallelBundleSolver::InnerLoopOrdered( bool extrastep ,
							   bool discard ,
							   bool batch )
{
 /* Deterministic parallel inner loop. Up to MaxThread components are kept
  * in flight (each compute()-d by a std::async), but, unlike the legacy
  * formulation, their results are *consumed in the fixed round-robin order
  * in which they were launched*: the main thread always blocks on the head
  * of the in-flight queue, even if a later task finished first. This removes
  * any dependence on thread timing, so both the set and the order of the
  * processed components are reproducible. Blocking on the head also makes
  * the dampened active wait (PoolingInt) unnecessary.
  *
  * The two flags select the remaining behaviour:
  *
  * - discard: when the master problem is guaranteed to change, the still
  *   in-flight tasks are dropped without being used (their provisional
  *   FiStatus is restored), so that exactly the components the sequential
  *   BundleSolver would process, and no more, are processed; if false the
  *   in-flight tasks are consumed instead (work-conserving);
  *
  * - batch: every ( NrFi - NrEasy ) component is evaluated, with no early
  *   stop on MPchgs (non-incremental). */

 // with at most one "hard" component there is nothing to parallelize: defer
 // to the base version, which also has the right edge-case semantics
 if( NrFi - NrEasy <= 1 )
  return( BundleSolver::InnerLoop( extrastep ) );

 // minimum number of components to evaluate; batch forces them all
 Index minceval = ( MinNrEvls >= 0 ? Index( MinNrEvls )
		                   : ( NrFi - NrEasy ) * ( - MinNrEvls ) );
 if( batch )
  minceval = NrFi - NrEasy;

 Index ceval = 0;  // how many distinct components have been evaluated

 // an in-flight task: the component, its FiStatus before the (provisional)
 // overwrite (needed to restore it on discard), and the std::future
 struct Task { Index wFi; int prev; std::future< int > fut; };
 std::deque< Task > inflight;

 // launch the evaluation of component w and push it at the back of the queue
 // note: FiStatus is provisionally set to kOK so that FindNext() does not
 // produce w again while it is in flight (see the legacy InnerLoop())
 auto launch = [ & ]( Index w ) {
  if( extrastep )
   SetupFiLambda( w );
  else
   SetupFiLambda1( w );
  Task t { w , int( FiStatus[ w ] ) , v_c05f[ w ]->compute_async(
					    ( FiStatus[ w ] == kUnEval ) ) };
  FiStatus[ w ] = kOK;
  inflight.push_back( std::move( t ) );
  };

 // drain the still in-flight tasks without using their results, restoring
 // the FiStatus they had before being launched (so they are re-evaluated as
 // if never touched this round)
 auto drain = [ & ]( void ) {
  for( auto & t : inflight ) {
   t.fut.get();
   FiStatus[ t.wFi ] = t.prev;
   }
  inflight.clear();
  };

 // ramp-up: launch up to MaxThread components in round-robin order - - - - -
 for( Index i = std::min( MaxThread , NrFi - NrEasy ) ; i ; --i ) {
  if( ! FindNext() )
   break;
  launch( f_wFi );
  }

 bool insrtd = false;
 bool stoplaunch = false;  // once true no new task is started
 Index lastproc = f_wFi;   // last processed component, to restore round-robin

 // consume the queue in launch order - - - - - - - - - - - - - - - - - - - -
 while( ! inflight.empty() ) {
  // block on the head (fixed order), even if a later task is ready first
  Task t = std::move( inflight.front() );
  inflight.pop_front();
  Index wFi = t.wFi;

  FiStatus[ wFi ] = t.fut.get();  // get() the status of compute()

  if( ! CurrNrEvls[ wFi ] )  // not evaluated before
   ++ceval;                  // one more evaluated
  ++CurrNrEvls[ wFi ];       // evaluated once more
  lastproc = wFi;

  BLOG( 6 , std::endl << "ordered: component " << wFi << " has status "
	    << FiStatus[ wFi ] );

  // unrecoverable error: stop, the whole compute() aborts- - - - - - - - - -
  if( ( FiStatus[ wFi ] <= kUnEval ) || ( FiStatus[ wFi ] >= kError ) ) {
   BLOG( 3 , std::endl << "            Component " << wFi
	     << " evaluated: Error" );
   Result = kError;
   drain();
   break;
   }

  // collect function values (upper and lower bound)- - - - - - - - - - - - -
  auto fwFi = v_c05f[ wFi ];
  auto ue = fwFi->get_upper_estimate();
  auto le = fwFi->get_lower_estimate();

  #if VERBOSE_LOG
   if( f_log && ( LogVerb > 3 ) ) {
    *f_log << std::endl << "            Component " << wFi
	   << " evaluated: UB = " << def;
    pval( *f_log , ue );
    *f_log << ", LB = ";
    pval( *f_log , le );
    }
  #endif

  // a component evaluating to -INF makes the whole problem unbounded below;
  // by convexity this shows up immediately, so just stop- - - - - - - - - - -
  if( ( f_convex && ( ue == -INFshift ) ) ||
      ( ( ! f_convex ) && ( le == INFshift ) ) ) {
   UpFiLmb1[ wFi ] = UpFiLmb1.back() = -INFshift;
   drain();
   break;
   }

  if( extrastep ) {
   // the method is actually being called on Lambda: update Lambda's
   // estimates and move on to the next component (no early stop)
   update_UpFiLambd( wFi , f_convex ? ue : - le );
   update_LwFiLambd( wFi , f_convex ? le : - ue );
   }
  else {
   // update UpFiLambd1[ wFi ] (and possibly UpFiLambd1[ NrFi ])
   update_UpFiLambd1( wFi , f_convex ? ue : - le );

   // if bit 4 of TrgtMng == 1, try to tighten UpFiLmb[ wFi ] using the
   // Lipschitz constant; the upper target is *not* changed (see legacy)
   if( ( TrgtMng & 16 ) && ( UpFiLmb1[ wFi ] < INFshift ) ) {
    c_VarValue LwFi = fwFi->get_Lipschitz_constant();
    if( LwFi < INFshift )
     update_UpFiLambd( wFi , UpFiLmb1[ wFi ] + LwFi * NrmD );
    }

   // update LwFiLambd1[ wFi ] (and possibly LwFiLambd1[ NrFi ])
   update_LwFiLambd1( wFi , f_convex ? le : - ue );

   // get new linearizations
   if( GetGi( wFi ) )
    insrtd = true;

   // check if the accrued information changes the MP (see legacy for the
   // detailed rationale of the strict inequalities)
   if( ( ! MPchgs ) && ( UpFiLmb1.back() < UpTrgt ) )
    MPchgs = 1;

   if( ( ! MPchgs ) && insrtd && RifeqFi && ( LwFiLmb1.back() > LwTrgt ) )
    MPchgs = 1;

   // wind-down tests (skipped in the extrastep path, as in the legacy one).
   // "earlystop" means the MP is already guaranteed to change (or time is
   // up): we must stop launching new tasks, and in faithful (discard) mode
   // we also drop the still in-flight tasks so that exactly the components
   // the sequential run would process, and no more, are processed. This is
   // *different* from simply running out of components to launch (FindNext()
   // below returning false): in that case the in-flight tasks are the last
   // components and must still be consumed, never dropped.
   bool earlystop = false;
   if( ( MaxTime < INFshift ) && ( get_elapsed_time() > MaxTime ) ) {
    Result = kStopTime;
    earlystop = true;
    }

   if( ( ! batch ) && MPchgs && ( ceval >= minceval ) )
    earlystop = true;

   if( earlystop ) {
    if( discard ) {     // faithful: drop the in-flight tasks and stop
     drain();
     break;
     }
    stoplaunch = true;  // work-conserving: keep consuming the in-flight ones
    }
   }

  // refill: launch the next component unless we have stopped launching - - - -
  if( ! stoplaunch ) {
   if( ! FindNext() )
    stoplaunch = true;  // no component left; just drain what is in flight
   else
    launch( f_wFi );
   }
  }  // end( while( consume queue ) )

 // restore the round-robin pointer to the last processed component, so that
 // the next InnerLoop() resumes exactly where a sequential run would
 f_wFi = lastproc;

 return( ceval );

 }  // end( ParallelBundleSolver::InnerLoopOrdered )

/*--------------------------------------------------------------------------*/
/*------------------- End File ParallelBundleSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
