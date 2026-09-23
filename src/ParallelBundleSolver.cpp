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
  * - in the cruise phase, the main thread waits for any task to end; the
  *   one that has ended is processed by extracting all the relevant
  *   information, and it is then substituted by another for a different
  *   component;
  *
  * - in the ramp-down phases, the main thread waits for the tasks not
  *   "consumed" already to end, processing each one as it does, but no
  *   other task takes its place, so that eventually the process ends.
  *
  * Note that gathering function values and linearizations from the evaluated
  * components is done in the main thread, and therefore it is a part of the
  * sequential bottleneck; however, this is required since the master problem
  * and all the other BundleSolver data structures are not protected from
  * concurrent access. */
 
 /* The threads this solver does not spend on evaluating the components at
  * once are given to the groups, which spend them on their members: with
  * fewer components than threads the outer loop leaves some of them idle,
  * and a group of many members is exactly what can use them. A group is
  * told once per inner loop, the number depending on how many components
  * there are to evaluate together; the single component is the extreme
  * case, where the outer loop has nothing to spend the threads on and they
  * are all the group's, which is why this comes before the check below. */

  if( MaxThread && ( ! v_groups.empty() ) ) {
  const Index at_once = std::max( std::min( MaxThread , NrFi - NrEasy ) ,
                                  Index( 1 ) );
  const int each = int( MaxThread / at_once );

  /* The two levels run on the same threads. The pool is sized for both,
   * i.e. for the components evaluated at once times the members each of
   * them may evaluate at once, which is MaxThread by construction; a group
   * runs one of its members in the thread that is waiting for them anyway,
   * so that no thread of the pool is ever held doing nothing and the
   * members cannot be starved by the components. */

  if( each > 1 )
   get_pool( MaxThread );

  for( auto & group : v_groups ) {
   group->set_members_at_once( each );
   if( each > 1 )
    group->set_submitter( [ this ]( ThinComputeInterface * f , bool cv ) {
                           return( f_pool->submit( f , cv , InINF ) );
                           } );
   else
    group->set_submitter();
   }
  }

 // if there is nothing to parallelize here, call the base class version - - -
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

 // the threads, and the end of each evaluation is announced by the index of
 // its position in EvalV
 auto & pool = get_pool( EvalV.size() );

 // FindNext() does not know which components are being evaluated right now:
 // one of them must not be handed out again, or two compute() would run on
 // the same C05Function at once
 std::vector< char > inflight( NrFi , 0 );
 auto find_next = [ & ]( void ) {
  for( Index i = 0 ; i < NrFi ; ++i ) {
   if( ! FindNext() )
    return( false );
   if( ! inflight[ f_wFi ] )
    return( true );
   }
  return( false );
  };

 // ramp-up phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // fill-in the vector of std::future
 // note: each time we start the evaluation of a component it is marked in
 // flight, and FiStatus is provisionally set to kOK; when the computation is
 // over, its FiStatus may become, say, kStopTime or kStopIter and the
 // component may be again eligible to be re-evaluated

 for( auto & el : EvalV ) {
  if( ! find_next() )
   throw( std::logic_error( "no component to evaluate in ramp-up phase" ) );

  BLOG( 6 , std::endl << "ramp-up: component " << f_wFi << " in position "
	    << ( & el - & EvalV.front() ) );

  el.first = f_wFi;
  if( extrastep )
   SetupFiLambda( f_wFi );
  else
   SetupFiLambda1( f_wFi );
  el.second = pool.submit( v_c05f[ f_wFi ] ,
                           ( FiStatus[ f_wFi ] == kUnEval ) ,
                           Index( & el - & EvalV.front() ) );
  FiStatus[ f_wFi ] = kOK;
  inflight[ f_wFi ] = 1;
  }

 // cruise phase- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // this is defined outside so that we can see what the last one was
 std::vector< EvalEl >::iterator it;

 bool insrtd = false;
 for( ; ; ) {
  // wait for any future to be ready, and read it- - - - - - - - - - - - - -
  it = EvalV.begin() + pool.wait_any();

  Index wFi = it->first;
  if( ! CurrNrEvls[ wFi ] )  // not evaluated before
   ++ceval;                  // one more evaluated
  ++CurrNrEvls[ wFi ];       // evaluated once more

  FiStatus[ wFi ] = it->second.get();  // get() the status of compute()
  inflight[ wFi ] = 0;

  BLOG( 6 , std::endl << "cruise: component " << wFi << " in position "
	    << it - EvalV.begin() << " has status " << FiStatus[ wFi ] );

  // if an unrecoverable error happens, immediately start the ramp-down - - -
  // kLowPrecision is not one of them [see BundleSolver::InnerLoop()]
  if( ( FiStatus[ wFi ] <= kUnEval ) ||
      ( ( FiStatus[ wFi ] >= kError ) &&
	( FiStatus[ wFi ] != kLowPrecision ) ) ) {
   if( f_log && ( LogVerb > 0 ) )
    *f_log << std::endl << "            Component " << wFi
	   << " evaluated: Error, status " << FiStatus[ wFi ];
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
   update_Fi_estimates( wFi , false , ue , le );

   // furthermore one immediately goes to put in the new task
   goto PutInNewTask;
   }
  
  update_Fi_estimates( wFi , true , ue , le );

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

  if( ! find_next() )      // find next component
   break;                  // if none, start ramp-down

  it->first = f_wFi;
  if( extrastep )
   SetupFiLambda( f_wFi );
  else
   SetupFiLambda1( f_wFi );
  it->second = pool.submit( v_c05f[ f_wFi ] ,
                           ( FiStatus[ f_wFi ] == kUnEval ) ,
                           Index( it - EvalV.begin() ) );
  FiStatus[ f_wFi ] = kOK;
  inflight[ f_wFi ] = 1;

  BLOG( 6 , std::endl << "cruise: in component " << f_wFi );

  }  // end( for( cruise phase loop ) )

 BLOG( 6 , std::endl << "end of cruise" );

 it->first = InINF;  // mark the entry as invalid
 
 // ramp-down phase - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index cnt = EvalV.size() - 1 ; cnt ; ) {
  // wait for any of the still running futures to be ready, and read it - -
  // note that the one consumed last in the cruise phase is not running, and
  // therefore it cannot be announced
  it = EvalV.begin() + pool.wait_any();

  Index wFi = it->first;
  it->first = InINF;
  if( ! CurrNrEvls[ wFi ] )  // not evaluated before
   ++ceval;                  // one more evaluated
  ++CurrNrEvls[ wFi ];       // evaluated once more
  --cnt;                     // one std::future less to wait for

  FiStatus[ wFi ] = it->second.get();  // get() the status of compute()
  inflight[ wFi ] = 0;

  BLOG( 6 , std::endl << "ramp-down: component " << wFi << " in position "
 	    << it - EvalV.begin() << " has status " << FiStatus[ wFi ] );

  // if an unrecoverable error happens, do nothing else - - - - - - - - - - -
  // kLowPrecision is not one of them [see BundleSolver::InnerLoop()]
  if( ( FiStatus[ wFi ] <= kUnEval ) ||
      ( ( FiStatus[ wFi ] >= kError ) &&
	( FiStatus[ wFi ] != kLowPrecision ) ) ) {
   if( f_log && ( LogVerb > 0 ) )
    *f_log << std::endl << "            Component " << wFi
	   << " evaluated: Error, status " << FiStatus[ wFi ];
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
   update_Fi_estimates( wFi , false , ue , le );

   continue;  // amd there is nothing left to do
   }
 
  update_Fi_estimates( wFi , true , ue , le );

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
  * in flight (each compute()-d by a thread of the pool), but, unlike the legacy
  * formulation, their results are *consumed in the fixed round-robin order
  * in which they were launched*: the main thread always blocks on the head
  * of the in-flight queue, even if a later task finished first. This removes
  * any dependence on thread timing, so both the set and the order of the
  * processed components are reproducible.
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
 std::deque< Task > queue;

 // launch the evaluation of component w and push it at the back of the queue
 // note: w is marked in flight so that it is not produced again until it is
 // consumed, and FiStatus is provisionally set to kOK (see InnerLoop())
 auto & pool = get_pool( std::min( MaxThread , NrFi - NrEasy ) );

 // FindNext() does not know which components are in flight: one of them
 // must not be launched again [see InnerLoop()]
 std::vector< char > inflight( NrFi , 0 );
 auto find_next = [ & ]( void ) {
  for( Index i = 0 ; i < NrFi ; ++i ) {
   if( ! FindNext() )
    return( false );
   if( ! inflight[ f_wFi ] )
    return( true );
   }
  return( false );
  };

 auto launch = [ & ]( Index w ) {
  if( extrastep )
   SetupFiLambda( w );
  else
   SetupFiLambda1( w );
  Task t { w , int( FiStatus[ w ] ) ,
           pool.submit( v_c05f[ w ] , ( FiStatus[ w ] == kUnEval ) ,
                        InINF ) };
  FiStatus[ w ] = kOK;
  inflight[ w ] = 1;
  queue.push_back( std::move( t ) );
  };

 // drain the still in-flight tasks without using their results, restoring
 // the FiStatus they had before being launched (so they are re-evaluated as
 // if never touched this round)
 auto drain = [ & ]( void ) {
  for( auto & t : queue ) {
   t.fut.get();
   FiStatus[ t.wFi ] = t.prev;
   inflight[ t.wFi ] = 0;
   }
  queue.clear();
  };

 // ramp-up: launch up to MaxThread components in round-robin order - - - - -
 for( Index i = std::min( MaxThread , NrFi - NrEasy ) ; i ; --i ) {
  if( ! find_next() )
   break;
  launch( f_wFi );
  }

 bool insrtd = false;
 bool stoplaunch = false;  // once true no new task is started
 Index lastproc = f_wFi;   // last processed component, to restore round-robin

 // consume the queue in launch order - - - - - - - - - - - - - - - - - - - -
 while( ! queue.empty() ) {
  // block on the head (fixed order), even if a later task is ready first
  Task t = std::move( queue.front() );
  queue.pop_front();
  Index wFi = t.wFi;

  FiStatus[ wFi ] = t.fut.get();  // get() the status of compute()
  inflight[ wFi ] = 0;

  if( ! CurrNrEvls[ wFi ] )  // not evaluated before
   ++ceval;                  // one more evaluated
  ++CurrNrEvls[ wFi ];       // evaluated once more
  lastproc = wFi;

  BLOG( 6 , std::endl << "ordered: component " << wFi << " has status "
	    << FiStatus[ wFi ] );

  // unrecoverable error: stop, the whole compute() aborts- - - - - - - - - -
  // kLowPrecision is not one of them [see BundleSolver::InnerLoop()]
  if( ( FiStatus[ wFi ] <= kUnEval ) ||
      ( ( FiStatus[ wFi ] >= kError ) &&
	( FiStatus[ wFi ] != kLowPrecision ) ) ) {
   if( f_log && ( LogVerb > 0 ) )
    *f_log << std::endl << "            Component " << wFi
	   << " evaluated: Error, status " << FiStatus[ wFi ];
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
   update_Fi_estimates( wFi , false , ue , le );
   }
  else {
   update_Fi_estimates( wFi , true , ue , le );

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
   if( ! find_next() )
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
/*----------------------------- THE EvalPool -------------------------------*/
/*--------------------------------------------------------------------------*/

ParallelBundleSolver::EvalPool::EvalPool( Index n )
{
 v_thread.reserve( n );
 for( Index i = 0 ; i < n ; ++i )
  v_thread.emplace_back( [ this ]() {
   for( ; ; ) {
    Job job;
    {
     std::unique_lock< std::mutex > lock( f_mutex );
     f_job_cv.wait( lock , [ this ]() {
                              return( f_stop || ( ! q_job.empty() ) );
                              } );
     if( q_job.empty() )  // f_stop and nothing left to run
      return;
     job = std::move( q_job.front() );
     q_job.pop_front();
     }

    job.task();  // the outcome, exception included, goes in the future

    if( job.tag != InINF ) {
     {
      std::lock_guard< std::mutex > lock( f_mutex );
      q_done.push_back( job.tag );
      }
     f_done_cv.notify_one();
     }
    }
   } );
 }

/*--------------------------------------------------------------------------*/

ParallelBundleSolver::EvalPool::~EvalPool()
{
 {
  std::lock_guard< std::mutex > lock( f_mutex );
  f_stop = true;
  }
 f_job_cv.notify_all();
 for( auto & t : v_thread )
  t.join();
 }

/*--------------------------------------------------------------------------*/

std::future< int > ParallelBundleSolver::EvalPool::submit(
                  ThinComputeInterface * f , bool changedvars , Index tag )
{
 std::packaged_task< int() > task( [ f , changedvars ]() {
                                    return( f->compute( changedvars ) );
                                    } );
 auto fut = task.get_future();
 {
  std::lock_guard< std::mutex > lock( f_mutex );
  q_job.push_back( Job{ std::move( task ) , tag } );
  }
 f_job_cv.notify_one();
 return( fut );
 }

/*--------------------------------------------------------------------------*/

BundleSolver::Index ParallelBundleSolver::EvalPool::wait_any( void )
{
 std::unique_lock< std::mutex > lock( f_mutex );
 f_done_cv.wait( lock , [ this ]() { return( ! q_done.empty() ); } );
 const auto tag = q_done.front();
 q_done.pop_front();
 return( tag );
 }

/*--------------------------------------------------------------------------*/

void ParallelBundleSolver::EvalPool::clear_finished( void )
{
 std::lock_guard< std::mutex > lock( f_mutex );
 q_done.clear();
 }

/*--------------------------------------------------------------------------*/

ParallelBundleSolver::EvalPool & ParallelBundleSolver::get_pool( Index n )
{
 if( ( ! f_pool ) || ( f_pool->size() < n ) ) {
  f_pool.reset();  // the old threads end before the new ones start
  f_pool = std::make_unique< EvalPool >( n );
  }

 f_pool->clear_finished();
 return( *f_pool );
 }

/*--------------------------------------------------------------------------*/
/*------------------- End File ParallelBundleSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
