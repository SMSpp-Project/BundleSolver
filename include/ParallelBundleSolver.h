/*--------------------------------------------------------------------------*/
/*-------------------- File ParallelBundleSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ParallelBundleSolver class, which extends BundleSolver
 * to implement the computation of the C0Function in parallel, with a simple
 * master/slave scheme.
 *
 * Apart from that ParallelBundleSolver is completely equivalent to
 * BundleSolver, in particular regarding the class of :Block it can solve.
 * However, having ParallelBundleSolver solve a :Block with a single
 * C0Function hardly makes sense (although it also makes little difference,
 * in that for a single C0Function ParallelBundleSolver behaves exactly as
 * BundleSolver with an extremely limited overhead).
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ParallelBundleSolver
 #define __ParallelBundleSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BundleSolver.h"

#include <condition_variable>
#include <deque>
#include <future>
#include <memory>
#include <mutex>
#include <thread>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ParallelBundleSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a BundleSolver where multiple C05Function computations are parallelized
/** The ParallelBundleSolver class extends BundleSolver to implement the
 * computation of the C0Function in parallel, with a basic master/slave
 * scheme.
 *
 * Apart from that ParallelBundleSolver is completely equivalent to
 * BundleSolver, in particular regarding the class of :Block it can solve.
 * However, having ParallelBundleSolver solve a :Block with a single
 * C0Function hardly makes sense (although it also makes little difference,
 * in that for a single C0Function ParallelBundleSolver behaves exactly as
 * BundleSolver with an extremely limited overhead). */

class ParallelBundleSolver : public BundleSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 *  @{ */

/*----------------------------- CONSTANTS ----------------------------------*/

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that ParallelBundleSolver has in addition to these of BundleSolver (and
  * CDASolver, Solver, ThinComputeInterface). The value intLastPBndSlvPar is
  * provided so that the list can be easily further extended by derived
  * classes.
  *
  * In fact ParallelBundleSolver currently has no extra parameters w.r.t.
  * BundleSolver, except the fact that it does really react to the
  * intMaxThread one of ThinComputeInterface which BundleSolver instead
  * ignores (since it's completely sequential implementation).
 */

 enum int_par_type_PBndSlv {

 intParFrm = intLastBndSlvPar ,
 ///< bit-wise selector of the parallel InnerLoop() formulation
 /**< Bit-wise encoded choice of how the parallel inner loop evaluates and
  * consumes the function components:
  *
  * - bit 0 (1): if 0, the "legacy" formulation is used, whereby ready
  *   std::future are consumed in completion order; this is fast but the set
  *   and order of the evaluated components, hence the whole trajectory of
  *   the algorithm, depend on thread timing and is therefore not
  *   reproducible. If 1, futures are instead consumed in a fixed round-robin
  *   index order, which makes the run deterministic;
  *
  * - bit 1 (2): only relevant with bit 0 == 1 and bit 2 == 0. If 1, as soon
  *   as the master problem is guaranteed to change the still in-flight tasks
  *   are discarded rather than used, so that exactly the same components are
  *   processed, in the same order, as the sequential BundleSolver would
  *   (faithful to the sequential run, but not work-conserving);
  *
  * - bit 2 (4): only relevant with bit 0 == 1. If 1, all the ( NrFi - NrEasy
  *   ) components are always evaluated at each iteration (non-incremental
  *   "batch" mode), which maximises the available parallelism.
  *
  * The meaningful combinations are: 0 = legacy, 1 = work-conserving
  * deterministic, 3 = faithful-to-sequential, 5 = deterministic batch. */

 intLastPBndSlvPar  ///< first allowed new int parameter for derived classes
 /**< Convenience value for easily allow derived classes
  * to extend the set of int algorithmic parameters. */

 };  // end( int_par_type_PBndSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that ParallelBundleSolver has in addition to these of BundleSolver
  * (and CDASolver, Solver, ThinComputeInterface). The value dblLastPBndSlvPar
  * is provided so that the list can be easily further extended by derived
  * classes. */

 enum dbl_par_type_PBndSlv {
  dblPoolingInt = dblLastBndSlvPar ,
  ///< parameter for declaring the frequency of pooling of results

  dblLastPBndSlvPar ///< first new double parameter for derived classes
                    /**< Convenience value for easily allow derived classes
		     * to extend the set of double algorithmic parameters. */

  };  // end( dbl_par_type_BndSlv )

/*@} -----------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING ParallelBundleSolver ----------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ParallelBundleSolver
 *  @{ */

 /// constructor: ensure every field is initialized

 ParallelBundleSolver( void ) : BundleSolver() {
  // ensure all parameters are properly given their default value
  MaxThread = ThinComputeInterface::get_dflt_int_par( intMaxThread );
  PoolingInt = 1e-4;
  ParFrm = 0;  // legacy (completion-order) formulation
  }

/*--------------------------------------------------------------------------*/
 /// destructor: does nothing special (explicitly)

 virtual ~ParallelBundleSolver() {}

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */

 /// set the int parameters of ParallelBundleSolver
 /** Set the int parameters specific of ParallelBundleSolver, which actually
  * is a parameter of ThinComputeInterface that ParallelBundleSolver
  * "listens to" while BundleSolver does not:
  *
  * - intMaxThread [0]: maximum number of threads that compute() uses to
  *                     evaluate the components. The threads are started
  *   the first time InnerLoop() needs them and live as long as the
  *   ParallelBundleSolver does, waiting between one evaluation and the next,
  *   so that no thread is created per evaluation; their number is
  *   min( intMaxThread , number of hard components ). */

 void set_par( idx_type par , int value ) override {
  if( par == intMaxThread )
   MaxThread = value;
  else
   if( par == intParFrm )
    ParFrm = value;
   else
    BundleSolver::set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameters of ParallelBundleSolver
 /** Set the double parameters specific of ParallelBundleSolver, or calls the
  * BundleSolver version to deal with the rest:
  *
  * - dblPoolingInt [1e-4]: accepted for compatibility, and not used: the main
  *                         thread is woken up by the evaluation that ends,
  *   rather than checking for one at fixed intervals. */

 void set_par( idx_type par , double value ) override {
  if( par == dblPoolingInt ) {
   if( value < 0 )
    throw( std::invalid_argument( "PoolingInt must be >= 0" ) );
   PoolingInt = value;
   }
  else
   BundleSolver::set_par( par , value );
  }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the BundleSolver
 *
 *  @{ */

 idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastPBndSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastPBndSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 int get_dflt_int_par( idx_type par ) const override {
  if( par == intParFrm )
   return( 0 );
  else
   return( BundleSolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 double get_dflt_dbl_par( idx_type par ) const override {
  if( par == dblPoolingInt )
   return( 1e-4 );
  else
   return( BundleSolver::get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/
 
 int get_int_par( idx_type par ) const override {
  if( par == intMaxThread )
   return( MaxThread );
  if( par == intParFrm )
   return( ParFrm );
  return( BundleSolver::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dbl_par( idx_type par ) const override {
  if( par == dblPoolingInt )
   return( PoolingInt );
  else
   return( BundleSolver::get_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 idx_type int_par_str2idx( const std::string & name ) const override {
  if( name == "intParFrm" )
   return( intParFrm );

  return( BundleSolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type dbl_par_str2idx( const std::string & name ) const override {
  if( name == "dblPoolingInt" )
   return( dblPoolingInt );
 
  return( BundleSolver::dbl_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 const std::string & int_par_idx2str( idx_type idx ) const override {
  static const std::string __pfname = "intParFrm";
  if( idx == intParFrm )
   return( __pfname );

  return( BundleSolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & dbl_par_idx2str( idx_type idx ) const override {
  static const std::string __psname = "dblPoolingInt";
  if( idx == dblPoolingInt )
   return( __psname );

  return( BundleSolver::dbl_par_idx2str( idx ) );
  }

/*@} -----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
 /* Performs the parallel inner loop: keeps (at most) MaxThread evaluations
  * running, each one compute()-ing a different component, up until the
  * conditions to stop are satisfied or there no longer are available
  * components to evaluate. */

 Index InnerLoop( bool extrastep = false ) override;

/*--------------------------------------------------------------------------*/
 /* Deterministic variant of the parallel inner loop: ready std::future are
  * consumed in a fixed round-robin index order rather than in completion
  * order, which makes the whole run reproducible. Selected by bit 0 of
  * intParFrm. The two boolean arguments encode the remaining bits:
  *
  * - discard (bit 1): when the master problem is guaranteed to change, the
  *   still in-flight tasks are discarded instead of being used, so that the
  *   processed components match exactly those of the sequential BundleSolver;
  *
  * - batch (bit 2): all the ( NrFi - NrEasy ) components are evaluated at
  *   each iteration, with no early stop (non-incremental). */

 Index InnerLoopOrdered( bool extrastep , bool discard , bool batch );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index MaxThread;    ///< maximum number of different threads (tasks)

 double PoolingInt;  ///< waiting time between each pooling round

 int ParFrm;         ///< bit-wise selector of the InnerLoop() formulation

 // generic fields- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// a fixed set of threads that compute() the components
 /** The threads wait for an evaluation to run and, once it has run, for the
  * next one. An evaluation submitted with a tag other than Inf< Index >()
  * puts it in the queue of the finished ones as soon as it ends, and
  * wait_any() blocks until that queue is not empty; the outcome itself
  * (the value returned by compute(), or the exception it has thrown) is in
  * the std::future that submit() returns. */

 class EvalPool {
 public:

  explicit EvalPool( Index n );

  ~EvalPool();

  Index size( void ) const { return( Index( v_thread.size() ) ); }

  /// runs f->compute( changedvars ) on one of the threads
  std::future< int > submit( ThinComputeInterface * f , bool changedvars ,
                             Index tag );

  /// blocks until a tagged evaluation has ended, and returns its tag
  Index wait_any( void );

  /// forgets the tags of the evaluations ended and not yet waited for
  void clear_finished( void );

 private:

  struct Job {
   std::packaged_task< int() > task;
   Index tag;
   };

  std::vector< std::thread > v_thread;
  std::deque< Job > q_job;
  std::deque< Index > q_done;
  std::mutex f_mutex;
  std::condition_variable f_job_cv;
  std::condition_variable f_done_cv;
  bool f_stop = false;
  };

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 std::unique_ptr< EvalPool > f_pool;  ///< the threads, once they are needed

 /// returns f_pool, (re)started if it has fewer than n threads
 EvalPool & get_pool( Index n );

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class ParallelBundleSolver )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ParallelBundleSolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File ParallelBundleSolver.h --------------------*/
/*--------------------------------------------------------------------------*/
