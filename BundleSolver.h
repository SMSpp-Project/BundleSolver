/*--------------------------------------------------------------------------*/
/*---------------------- File BundleFSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the BunldeSolver class, which implements the Solver
 * interface, in particular in its CDASolver version. The class contains
 * a NonDifferentiable Optimization Solver using a "Generalized Bundle"
 * algorithm.
 *
 * The user is assumed to be familiar with the algorithm: refer to
 *
 *  A. Frangioni "Generalized Bundle Methods"
 *  SIAM Journal on Optimization 13(1), p. 117 - 156, 2002
 *
 * available at
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#SIOPT02
 * \endlink
 *
 * The class requires that the function to be minimized be available under
 * the form of a FiOracle object, as described in FiOracle.h.
 *
 * The class is parametric over the type of Master Problem used: it just
 * relies over an object of class MPSolver to solve it, see MPSolver.h.
 *
 * \version 0.01
 *
 * \date 19 - 05 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 2019 by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BundleSolver
 #define __BundleSolver  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"

#include "C05Function.h"
#include "LinearFunction.h"

#include "Block.h"
#include "ColVariable.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"

#include "MILPSolver.h"
#include "MPSolver.h"

#include "NDOSlver.h"

/*------------------------------- LOG_BND ----------------------------------*/

#define LOG_BND 1

/* If LOG_BND > 0, the Bundle class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetBLog() [see below]. */

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 using namespace NDO_di_unipi_it;
#endif

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup BundleSolver_CLASSES Classes in BundleSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BundleSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// CDASolver
/** The BunldeSolver implements the Solver interface for the (Generalized)
    Bundle algorithm. */

class BundleSolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

  typedef unsigned int Index;  /// type for the indices
  typedef std::vector<Index> Vec_Index;  ///< a std::vector of Index

  typedef int SIndex;  /// type for the signed indices
  typedef std::vector<SIndex> Vec_SIndex;  ///< a std::vector of SIndex

  typedef double VarValue;     ///< type of the value of the ColVariable
  typedef std::vector<VarValue> Vec_VarValue; ///< a std::vector of VarValue

  typedef std::vector<bool> Vec_Bool;  ///< a std::vector of bool
  typedef std::vector<OFValue> Vec_OFValue; ///< a std::vector of OFValue

  typedef unsigned int LinearizationName;
  ///< type used to define names of linearizations

  typedef double FunctionValue;
  ///< type of the returned value by Function

  typedef std::vector< std::pair < LinearizationName , FunctionValue > >
  LinearCombination;
  ///< type used to define linear combinations of linearizations

/*--------------------------------------------------------------------------*/

 /// public enum for the int algorithmic parameters
 /** Public enum describing the different types of algorithmic parameters
  * of "int" type that BundleSolver should have. The value
  * intLastAlgPar is provided so that the list can be easily further extended
  * by derived classes. */

 enum int_par_type_BndSlv {

 intBPar1 = CDASolver::intLastParCDAS ,    ///<
       /**< If an item has had a zero multiplier for the last BPar1 steps,
        * it is eliminated; if BPar1 is "too small" precious information
 *  may be lost, but keeping the "bundle" small obviously makes the Master
 * Problem cheaper. */

 intBPar6 ,   ///<
         /**< These parameters control how the actual number of
         * subgradients/constraints (items) that are requested
 * to the C05Function, between BPar4 and BPar3, as the
 * algorithm proceeds; note that what varies in practice is
 * the maximum number, as it is always legal for the FiOracle to refuse
 * giving other items, although the BundleSolver code will complain and
 * stop if less than BPar4 are given. In the Bundle code, the number
 *
 *      EpsU = Sigma + D_{tStar}*( z* ) / max( | ReadFiVal() | , 1 ) ,
 *
 * where Sigma = Sum_i Fi[ i ]_{B,Lambda}*( z[ i ]* ) + \sigma_L( w ) and
 * z* = - Sum_i z[ i ]* is the optimal solution of the stabilized Dual
 * Master Problem [see MPSolver.h], is used as an estimate of the relative
 * gap between the current and the optimal solution; that is, IsOptimal()
 * returns true if EpsU <= RAccSol. Thus, the number RAccSol / EpsU is always
 * smaller than one, and typically increases as the algorithm proceeds.
 * Depending on the value of BPar6, the following formulae for the actual
 * value of BPar3, aBP3, are used:
 *  0: aBP3 is set to (the "finalized" value of) BPar3 and never changed;
 *  1: if BPar5 > 0 then aBP3 is initialized to (...) BPar4 and increased
 *        every BPar5 iterations, while if BPar5 <= 0 then aBP3 is
 *        initialized to BPar3 and decreased every - BPar5 iterations;
 *  2: aBP3 is set to
 *        ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 * ( RAccSol / EpsU )
 *  3: aBP3 is set to
 *         ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 / sqrt( EpsU / RAccSol )
 *  4: aBP3 is set to
 *        ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 / log10( EpsU / RAccSol ) */

 intMnSSC ,   ///< minimum number of consecutive serious steps
       /**< minimum number of consecutive SS with the same t that
       * have to be performed before t is allowed to grow  */

 intMnNSC ,   ///< minimum number of consecutive null steps
       /**< minimum number of consecutive NS with the same t that
       * have to be performed before t is allowed to diminish  */

 inttSPar1 ,    ///< first t-strategy parameter
       /**<  Select the t-strategy used. This field is coded
       * bit-wise in the following way.
 * The first two bits control which heuristics are used to
 * compute a new value of t when increasing/decreasing it.
 * There are two heuristics avaliable, H1 and H2, both based on a quadratic
 * interpolation of the function but differing in which derivative is used:
 * H1 uses the derivative in the new tentative point, and it is guaranteed
 * to produce a value greater than the current one if and only if the scalar
 * product between the direction and the newly obtained subgradient is < 0
 * (indicating that a longer step along the same direction could have been
 * advantageous), while H2 uses the derivative in the current point and it
 * does not possess this property. The value of the first two bits of
 * tSPar1 has the following meaning:
 *   bit 0:  which heuristic is used to increase t: 0 = H1, 1 = H2
 *   bit 1:  which heuristic is used to decrease t: 0 = H2, 1 = H1
 *           The following bits of tSPar1 tell which long-term t-strategy is
 *           used, with the following values:
 *    0 (+ 0):  none, only the heuristics are used
 *    1 (+ 4):  the "soft" long-term t-strategy is used: an optimality
 *              estimate EpsU is mantained which estimates how far from the
 *       optimal value one currently is, and decreases of t are
 *       inhibited whenever v < tSPar2 * EpsU * | Fi |
 *    2 (+ 8):  the "hard" long-term t-strategy is used: an optimality
 *              estimate EpsU is mantained as above and t is increased
 *       whenever v < tSPar2 * EpsU * | Fi |
 *    3 (+12):  the "balancing" long-term t-strategy is used, where the two
 *              terms D*_t( -z* ) and Sigma* are kept of "roughly the same
 *       size": if D*_1( -z* ) <= tSPar2 * Sigma* then t increases
 *       are inhibited (increasing t causes a decrease of D*_1( -z* )
 *	     that is already small), if tSPar2 * D*_1( -z* ) >= Sigma*
 *       then t decreases are inhibited (decreasing t causes an
 *       increase of D*_1( -z* )	that is already big).
 *       Still later bits of tSPar1 activate "special cases" t-strategies:
 *    4 (+16):  the "endgame" t-strategy is used, where if D*_1( -z* ) is
 *               "small" (~ 1/10 of the current absolute epsilon) t is
 *       decreased no matter what the other strategies dictated.
 *       The rationale is that we are "towards the end" of the
 *       optimization and here t needs decrease. However, note that
 *       having D*_1( -z* ) "small" is no guarantee that we actually
 *       are at the end, especially if the FiOracle dynamically
 *       generates its variables, so use with caution.   */

 intMaxNrEvls , ///< max number of function evaluation for each iteration
 intKpBstL , ///< true if LmbdBst has to be kept, false otherwise

 intMPName, ///< true if MP solver is a QPPenalty, otherwise MP is a OSiMPSolver
 intMPlvl , ///< log verbosity of Master Problem

 intQPmp1, ///< MxAdd parameter for QPPenaltyMP solver only
 intQPmp2, ///< MxRmv parameter for QPPenaltyMP solver only

 intOSImp1 , ///< algorithm type for OsiMP solver only
 intOSImp2 , ///< reduction parameter for OsiMP solver only
 intOSImp3 , ///< threads parameter for OsiMP solver only
 intOSImp4 , ///< stabilization type for OsiMP solver only
 intOSImp5 , ///< which OsiXXXSolverInterface is used in OsiMP solver

 intLastBndSlvPar ///< first allowed new int parameter for derived classes
       /**< Convenience value for easily allow derived classes
 * to extend the set of int algorithmic parameters. */

 };  // end( int_par_type_S )

 /*--------------------------------------------------------------------------*/
  /// public enum for the double algorithmic parameters
  /** Public enum describing the different types of algorithmic parameters
   * of "double" type that any Solver should reasonably have. The value
   * dblLastAlgPar is provided so that the list can be easily further extended
   * by derived classes. */

  enum dbl_par_type_S {
  dbltStar = dblLastParCDAS ,    ///<
       /**< Optimality related parameters. Proving that some point Lambda is
       * optimal for a NonDifferentiable Optimization problem involves finding
 * an all-0 subgradient of the function at Lambda. If an all-0 vector is found in
 * the epsilon-subdifferential of Lambda, then the point is epsilon-optimal.
 * Note that if the minimization problem is subject to constraints, i.e.,
 * Fi() has to be minimized only on the points Lambda \in L, the latter
 * being a convex set, then the above is referred to a subgradient of the
 * "actual function" ( Fi + I_L )( Lambda ), where I_L is the indicator
 * function of L (evaluating to 0 inside L and to +INF otherwise). In other
 * words, one has to show that there exists a( enspilon-)subgradient of
 * Fi() at Lambda that, *after projection on the frontier of L*, is all-0.
 * A general stopping condition requires that, if RAccSol is the *relative*
 * precision required, a solver can stop if it finds an epsilon-subgradient
 * g at Lambda such that
              tStar * || g || + epsilon <= RAccSol * | MaxFi |
 * where MaxFi is an estimate of the optimal solution value of the NDO
 * problem, tStar is an estimate of the longest step that can be performed
 * and || || is a norm-like function. tStar is related to the "scaling" of
 * Fi(), and it can be seen as an estimate of the actual decrease that can be
 * obtained by moving of an unitary step in the direction of any subgradient.
 * Alternatively, the above condition can be seen as a weaker form of
        epsilon <= RAccSol * | MaxFi | / 2
        || g || <= RAccSol * | MaxFi | / ( 2 * tStar )
 * which says that the solver stops when both epsilon and || g || are
 * "small", with tStar dictating what "small" means for || g ||. Note that
 * each derived class can use different norm-like functions to evaluate one
 * or the other of the above conditions. Also, different estimates can be
 * used for MaxFi, although using the bast Fi-value found so far is pretty
 * common. */

 dblEInit ,    ///< *maximum* precision required to the C05Function
       /**< ABS( EInit ) is the initial, and *maximum*, precision
                   required to the C05Function, but the sign tells how the
		   "emergency mechanism" alluded to above interacts with
 * the "regular mechanism" controlled by these parameters. In particular,
 * if EInit > 0 then the accuracy is monotonically non-increasing: if the
 * "emergency mechanism" reduces it, then it will remain "at least as small"
 * in all the following iterations, even if the value computed by the
 * "regular mechanism" would be larger. If EInit < 0 instead, then each time
 * a "regular step" is computed the precision is reset to that dictated by
 * the "regular mechanism", even if it is larger than the current one (for
 * instance, if EStps == 0, see below, then the precision is set to a "fixed"
 * value). */

 dblBPar3 ,    ///< maximum number of new subgradients/constraints
 dblBPar4 ,    ///< minimum number of new subgradients/constraints
       /**< Maximum and minimum number of new subgradients/constraints (items)
        * to be fetched from the FiOracle for
 * each function evaluation. Two different ways are given for specifying
 * these numbers: positive values are (rounded up and) taken as
 * absolute values, while negative numbers are first multiplied by
 * the number of components of the C05Function (and then
 * rounded up); thus, the default vale "-1" stands for "one for each of the
 * components of Fi()". Clearly, the "finalized" value of BPar3 has to be
 * <= BPar2, and the "finalized" value of BPar4 has to be <= than that.  */

 dblBPar5 ,    ///<
       /**< see intbPar6 above */


 dblm1 ,    ///< m1 factor
       /**< SS condition: if DeltaFi >= | m1 | * Deltav, then a
       * SS is done. What is taken as Deltav depends on the sign of
 * m1: if m1 > 0 then Deltav = - v* (the decrease predicted
 * by the model), while if m1 < 0 then Deltav = - ( v* + D_t( d* ) ), i.e.
 * the optimal objective function value the dual Master Problem. Since
 * - v* >= - ( v* + D_t( d* ) ), the second condition is weaker and may
 * lead to a larger number of SS (and therefore possibly a fater convergence)
 * while still ensuring global convergence. The value m1 = 0, i.e., perform a
 * SS for whatever small improvement in the objective function, can only be
 * used, at least in theory, for some classes of functions (the polyhedral
 * ones) and some under assumptions on the MP. */

 dblm3 ,    ///< m3 factor
       /**< A nevly obtained subgradient is deemed "useless" if
       * Alfa >= m3 * Sigma; in this case, if a NS has to be done,
 * t is decreased. This parameter is mostly critical: if no
 * "long-term" t-strategy [see tSPar1 below] is used, values < 2/3 usually
 * make t to decrease rather fast to tMinor [see below], possibly making the
 * algorithm to perform very short steps and therefore to converge very
 * slowly. When a "long-term" t-strategy is used, .9 may be a good value.*/

 dblmxIncr ,    ///< maximum increasing t-factor
 dblmnIncr ,    ///< minimum increasing t-factor
       /**< each time t grows, the new value of t is chosen in
       * the interval [t * mnIncr, t * mxIncr] (t is the
 * previous value) */

 dblmxDecr ,    ///< maximum decreasing t-factor
 dblmnDecr ,    ///< minimum decreasing t-factor
       /**<  each time t diminishes, the new value of t is chosen
       * in the interval [t * mxDecr, t * mnDecr] (t is the
 * previous value) */

 dbltMaior ,    ///< maximum value of t
 dbltMinor ,    ///< minimum value of t
 dbltInit ,     ///< initial value of t
       /**< Maximum, minimum, and initial value of t. These parameters may be
       * critical, but they are not very difficult to set. Usually,
 * there is a "right" order of magnitude for t, that is the
 * one that is guessed by the t-heuristics during most of the run, even
 * though the starting value is very different. Hence, a good setting for
 * tInit is in that order of magnitude, while tMinor should be set small
 * enough to never enter into play. Note that t is always kept <= tStar
 * [see NDOSolver.h], and that a "good" value for tStar (i.e., one that
 * actually ensures that the stopping point is RAccSol-optimal) is usually
 * one or two orders of magnitude larger than a "good" tInit. */

 dbltSPar2 ,    ///< Numerical parameter for the long-term t-strategies
       /**< see inttSPar1 above. */

 dblQPmp1 , ///< cut-off value for QPPenaltyMP solver only

  dblLastBndSlvPar ///< first allowed new double parameter for derived classes
        /**< Convenience value for easily allow derived classes
  * to extend the set of double algorithmic parameters. */
  };

/*@} -----------------------------------------------------------------------*/
/*----------------- CONSTRUCTING AND DESTRUCTING BundleSolver --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing BundleSolver
 *  @{ */

 /// constructor: does nothing special
 /** Void constructor: does nothing special, except verifying that the
  * template argument derives from MCFClass. */

 BundleSolver( ) : CDASolver() , FakeFi( this )  {
    
  }

/*--------------------------------------------------------------------------*/
 /// destructor: it has to release all the Modifications

 virtual ~BundleSolver() {

  // delete[] Fi;
  }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve

 virtual void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 // set the ostream for the Solver log
 // not really, MCFClass objects are remarkably silent
 //
 // virtual void set_log( std::ostream *log_stream = nullptr ) override;

/*--------------------------------------------------------------------------*/

 virtual void set_par( const idx_type par , const int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual void set_par( const idx_type par , const double value ) override;

/*--------------------------------------------------------------------------*/

 virtual void set_log( std::ostream *log_stream = nullptr ) override;

/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the MCF encoded by the current MCFBlock
 *  @{ */

 /// (try to) solve the MCF encoded in the MCFBlock

 virtual int compute( bool changedvars = true ) override;

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 virtual OFValue get_lb( void ) override { return(UpRifFi[ NrFi ]); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual OFValue get_ub( void )  override { return(UpRifFi[ NrFi ]); }

/*--------------------------------------------------------------------------*/

 virtual bool has_var_solution( void ) override { return( true ); }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual bool has_dual_solution( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*
 virtual bool is_var_feasible( void ) override { return( true ); }

 virtual bool is_dual_feasible( void ) override { return( true ); }
*/
/*--------------------------------------------------------------------------*/
 /// write the "current" solution
 virtual void get_var_solution( Configuration *solc = nullptr ) override
 {
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the "current" dual solution
 virtual void get_dual_solution( Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 virtual bool new_var_solution( void ) override
 {
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual bool new_dual_solution( void )  override
 {
  return( false );
  }

/*--------------------------------------------------------------------------*/
/*
 virtual void set_unbounded_threshold( const OFValue thr ) override { }
*/

/*--------------------------------------------------------------------------*/

 virtual bool has_var_direction( void ) override { return( true ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual bool has_dual_direction( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
 /// write the current direction
 virtual void get_var_direction( Configuration *dirc = nullptr ) override
 {
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the current dual direction
 virtual void get_dual_direction( Configuration *dirc = nullptr ) override
 {
  }

/*--------------------------------------------------------------------------*/
/*
 virtual bool new_var_direction( void ) override { return( false ); }

 virtual bool new_dual_direction( void ) override{ return( false ); }
*/
/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/

/*
 virtual bool is_dual_exact( void ) const override { return( true ); }
*/
 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the BundleSolver
 *
 * Each BundleSolver< MCFC > may have its own extra int / double parameters. If
 * this is the case, it will have to specialize the following methods to
 * handle them. The general definition just handles the case of the
 *
 * intLastParCDAS ==> kReopt             whether or not to reoptimize
 *
 * extra (int) parameter and otherwise issues the method of the base
 * CDASolver class, which is OK for each MCFC that does *not* have any extra
 * parameter of the corresponding type (apart from that). The get_*_par()
 * methods exploit the same two const static arrays Solver_2_MCFClass_int and
 * Solver_2_MCFClass_dbl as the set_*_par(), with a negative entry meaning
 * "there is no such parameter in BundleSolver".
 *  @{ */

 virtual idx_type get_num_int_par( void ) const override
 {
  return( idx_type( intLastBndSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual idx_type get_num_dbl_par( void ) const override
 {
  return( idx_type( dblLastBndSlvPar ) );
  }

/*--------------------------------------------------------------------------*/
 
 virtual int get_dflt_int_par( const idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 virtual double get_dflt_dbl_par( const idx_type par ) const override;

/*--------------------------------------------------------------------------*/
 
 virtual int get_int_par( const idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 virtual double get_dbl_par( const idx_type par ) const override;

/*--------------------------------------------------------------------------*/

 virtual idx_type int_par_str2idx( const std::string & name ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual idx_type dbl_par_str2idx( const std::string & name ) const override;

/*--------------------------------------------------------------------------*/

 virtual const std::string & int_par_idx2str( const idx_type idx )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual const std::string & dbl_par_idx2str( const idx_type idx )
  const override;

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

 /*
 virtual void add_Modification( sp_Mod &mod ) {
  v_mod.push_back( mod );
  }
 */

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

 void FormD( void );

/* When no variables generation is done (PPar2 == 0), FormD() just calls
 SolveMP() once and calculates the direction d: however, it also implements
 some strategies to survive to "fatal" failures in the subproblem solver,
 typically eliminating some of the items in the bundle.

 Set the protected field Result to kOK if (evenctually after some "fatal"
 failure) a tentative descent direction could be found, to kUnfsbl if the
 MP is dual unfeasible and to kError if this was returned by SolveMP(): in
 the latter cases, the whole algorithm must abort.

 If variables generation is done (PPar2 > 0), this is where the
 corresponding strategies are implemented: in this case, SolveMP() can be
 called more than once within the same call to FormD(), since the resulting
 direction has to be optimal w.r.t. all the current "active set" of
 variables. */

/*--------------------------------------------------------------------------*/
 ///< Updates the out-of-base counters for all items in the Bundle.
 void UpdtCntrs( void );

/*--------------------------------------------------------------------------*/

 void FormLambda1( HpNum Tau );
 /* After a (succesfull) call to FormD(), sets the new tentative point Lambda1
    (a protected field of type LMRow) as Lambda1 = Lambda + ( Tau / t ) * d. */

/*--------------------------------------------------------------------------*/

 bool FiAndGi( Index wFi );

 /* Computes Fi( Lambda1 ), inserting the obtained items (subgradients or
    constraints) in the bundle. Returns true <=> the newly obtained information
    changes the solution of the MP. */

/*--------------------------------------------------------------------------*/

 void GotoLambda1( void );

 /* Move the current point to Lambda1. */

/*--------------------------------------------------------------------------*/
 ///< Eliminate outdated items, i.e., these with "large" out-of-base counter.

 void SimpleBStrat( void );

/*--------------------------------------------------------------------------*/

 void UpdtLowerBound( void );

/*--------------------------------------------------------------------------*/

 double BetaK( Index wFi );

/*--------------------------------------------------------------------------*/

    void Log1( void );

    void Log2( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

int MaxSol;         ///< maximum number of different solutions to report
double RelAcc;      ///< relative accuracy for declaring a solution optimal
double AbsAcc;      ///< absolute accuracy for declaring a solution optimal
double RAccSol;     ///< maximum relative error in any reported solution
double AAccSol;     ///< maximum absolute error in any reported solution
double FAccSol;     ///< maximum constraint violation in any reported solution

Index MaxIter;      ///< maximum number of iterations
double MaxTime;     ///< maximum time (in seconds) for each call to Solve()

double tStar;       ///< optimality related parameter: "scaling" of Fi

double EInit;      ///< precision-related parameter: initial precision
int Result;        ///< result of the latest call to Solve()

Index NumVar;      ///< (current) number of variables
Index NrFi;        ///< number of components of Fi()

Index SCalls;      ///< nuber of calls to Solve() (the current included)
Index ParIter;     ///< nuber of iterations in this run
Index FiEvaltns;   ///< total number of Fi() calls
Index GiEvaltns;   ///< total number of Gi() calls

int LogVerb;       ///< "verbosity" of the log

int BPar1;           // parameter for removal of items (B-strategy)
int BPar2;           // max Bundle size
double BPar3;        // max number of items fetched from Fi() at each call
double BPar4;        // min number of items fetched from Fi() at each call
double BPar5;        // control how the actual BPar3 changes over time
int BPar6;           // same as above

double mxIncr;       // max increase/decrease t parameters:
double mnIncr;       // see the description in the constructor
int MnSSC;
double mxDecr;
double mnDecr;
int MnNSC;

double m1;            // parameters for deciding if a SS/NS is done:
double m2;
double m3;            // see the description in the constructor

double tMaior;        // max value for t
double tMinor;        // min value for t
double tInit;         // initial value for t

int tSPar1;          // long-term t-strategy parameters
double tSPar2;        // see the description in the constructor

Vec_Bool IsEasy;     // tells whether any component of Fi is "easy"
Index NrEasy;

std::vector<VarValue> Lambda;  // the current point
std::vector<VarValue> Lambda1; // the tentative point
std::vector<VarValue> LmbdBst; // the best point found so far


bool KpBstL;         // if LmbdBst has to be kept
bool LHasChgd;       // true if Lambda has changed since the latest call
                     // to FiAndGi(): allows repeated calls in the same
                     // Lambda, e.g. with increasing precision
bool tHasChgd;       // true if t has changed since the last MP

Vec_Index whisZ;     // the position in the bundle where the "aggregate
                     // subgradient" Z[ k ] of "component" k is kept in
                     // whisZ[ k ]; Inf<Index>() means it is not in the
                     // bundle
Vec_Index whisG1;    // name of the "representative subgradient" for each
                     // component of Fi()
Vec_OFValue ScPr1;   // ScalarProduct( dir , G[ WhIsG1[ k ] ] )
Vec_OFValue Alfa1;   // linearization error of G[ WhIsG1[ k ] ] w.r.t. the
                     // current point Lambda
Vec_OFValue DeltaAlfa; // correction of Fi-values due to inexactness

OFValue LowerBound;// Lower Bound over (the various components of) Fi

double t;             // the (tremendous) t parameter
double Prevt;         // what t were before being changed for funny reasons

double Sigma;         // Sigma*: convex combination of the Alfa's
double DSTS;          // D*_{t*}( -z* ), the other part of the dual objective
Vec_OFValue vStar;         // v*, the predicted improvement
double Deltav;        // the "desired improvement" in the Fi-value

double DeltaFi;       // FiLambda - FiLambda1
double EpsU;          // precison required by the long-term t-strategy

double EpsFi;         // the last precision passed to the FiOracle (can be
                     // different from EpsCurr)

int ParSS;           // number of SS within the present call to Solve()

int CSSCntr;         // counter of consecutive SS
int CNSCntr;         // counter of consecutive NS

Index FreDim;        // number of free positions in the Bundle
Vec_Index FreList;   // list (heap) of free positions in the Bundle

Vec_SIndex OOBase;   // Out-Of-Base counters:
                     // = Inf<SIndex>() means no item is there
                     // = k > 0 means out of base since k iterations
                     // = 0 means in the current base but potentially
                     //   removable
                     // = a *finite* negative value - k means not
                     //   removable for the next k iterations: note that
                     //   some items in base may be such
                     // = - Inf<SIndex>() means unremovable

bool TrueLB;         // true if LowerBound is a "true" lower bound rather
                     // than just the "minus infinity"
bool LBHasChgd;      // true some LowerBound has changed
bool SSDone;         // true if the laste step was a SS

Vec_Index FiStatus;

/*--------------------------------------------------------------------------*/

 std::vector< std::pair < LinearizationName , LinearCombination > > zA;
 /* the vector of the pairs  important linearization name and the
    linear combination used to form it */

 std::vector< C05Function * > v_c05f; /* the vector of the components of the
                                         sum function */
 LinearFunction * linear_function; ///< the 0-th component of the sum function
 MPSolver *Master;    ///< (pointer to) the Master Problem Solver
 std::vector<MILPSolver*> MILP_s; /* MILP solver used to read the
                  easy part of the sub-problems */

 std::vector<ColVariable *> LamVcblr;    // the set of indices of Lambda
 bool MPName; // true if MP solver is a QPPenalty, otherwise MP is a OSiMPSolver
 double UpTrgt; // upper target
 double LwTrgt; // lower target

 Vec_OFValue UpFiBest;   // Fi best value vector
 Vec_OFValue UpRifFi;    /* The value of Fi[ k ]() where the zero of the Cutting
                            Plane models are fixed: it is == FiLambda[ k ]() but
                            when FiLambda[ k ]() == HpINF */

 Vec_OFValue UpFiLmb1;   ///< upper function value vector at the tentative point
 Vec_OFValue LwFiLmb1;   ///< lower function value vector at the tentative point

 Vec_OFValue UpFiLmb;    ///< upper function value vector at the current point
 Vec_OFValue LwFiLmb;    ///< lower function value vector at the current point

 Index MaxNrEvls;
 std::vector<Index> CurrNrEvls;

 double DeltaStar;
 double NrmD;

 Index MxAdd;
 Index MxRmv; /* How many variables can be "moved" at each iteration of the
        "upper-level method" for solving bound-constrained MPs
        implemented into the base class BMinQuad (if any, see
		HV_NNVAR above, otherwise the parameters are just ignored);
		see SetMaxVar[Add/Rmv]() in BMinQuad.h. If 0 [the default]
		is used, no bound is given. */

 double CtOff; /* The "break" value for the pricing in the base class
                 MinQuad [see SetPricing() in MinQuad.h]: positive values are
                 passed untouched, while any negative value is turned into
                 Inf<HpNum>() */

 Index algo;      ///< algorithm type ( for OSIMPSolver only )
 Index reduction; ///< pre-processing (reduction) ( for OSIMPSolver only )
 Index threads;   ///< number of threads ( for OSIMPSolver only )
 Index stblztn;   ///< stabilization type ( for OSIMPSolver only )
 bool osi_type;   /* which OsiXXXSolverInterface (for OsiMPSolver):
                      0 = Clp, 1 = Cplex */

 Index MPlvl; // log verbosity of master problem

 bool log_chgd ;

/*--------------------------------------------------------------------------*/

 const static std::vector<int> dflt_int_par;
 ///< the (static const) vector of int parameters default values

 const static std::vector<double> dflt_dbl_par;
 ///< the (static const) vector of double parameters default values

 const static std::vector< std::string > int_pars_str;
 ///< the (static const) vector of int parameters names

 const static std::vector< std::string > dbl_pars_str;
 ///< the (static const) vector of double parameters names

 const static std::map< std::string , idx_type > int_pars_map;
  ///< the (static const) map for int parameters names

 const static std::map< std::string , idx_type > dbl_pars_map;
 ///< the (static const) map for double parameters names

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS FakeFiOracle  ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

class FakeFiOracle : public FiOracle
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
    @{ */


/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

/** Constructor of the class: takes no arguments, since everything that
    concerns the real evaluation of the function must be done in derived
    classes, which will have their parameters. */

   FakeFiOracle( BundleSolver *solver );

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

   virtual void SetNDOSolver( NDOSolver *NwSlvr = 0 ) override;

/*--------------------------------------------------------------------------*/

   virtual void SetFiLog( ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

   virtual void SetFiTime( const bool TimeIt = true ) override;

/*--------------------------------------------------------------------------*/

   virtual void SetMaxName( cIndex MxNme = 0 ) override;

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

/// get the number of Variable
/** Variable cannot be changed. This means that is used the default
 *  implementation of GetMaxNumVar(). The maximum number of variables is
 *  equal to the current number of variable*/

  virtual Index GetNumVar( void ) const override;

/*--------------------------------------------------------------------------*/

  virtual Index GetNrFi( void ) const override;

/*--------------------------------------------------------------------------*/

  virtual Index GetMaxName( void ) const override;

/*--------------------------------------------------------------------------*/

  virtual HpNum GetMinusInfinity( void ) override;

/*--------------------------------------------------------------------------*/

  virtual Index GetMaxNZ( cIndex wFi = Inf<Index>() ) const override;

/*--------------------------------------------------------------------------*/

  virtual Index GetMaxCNZ( cIndex wFi = Inf<Index>() ) const override;

/*--------------------------------------------------------------------------*/

  virtual bool GetUC( cIndex i ) override;

/*--------------------------------------------------------------------------*/

  virtual LMNum GetUB( cIndex i ) override;

/*--------------------------------------------------------------------------*/

  virtual LMNum GetBndEps( void ) override;

/*--------------------------------------------------------------------------*/

  virtual HpNum GetGlobalLipschitz( cIndex wFi = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

  virtual Index GetBNC( cIndex wFi ) override;

/*--------------------------------------------------------------------------*/

  virtual Index GetBNR( cIndex wFi ) override;

/*--------------------------------------------------------------------------*/

  virtual Index GetBNZ( cIndex wFi ) override;

/*--------------------------------------------------------------------------*/

  virtual void GetBDesc( cIndex wFi , int *Bbeg , int *Bind , double *Bval ,
 			  double *lhs , double *rhs , double *cst ,
 			  double *lbd , double *ubd ) override;

/*--------------------------------------------------------------------------*/

  virtual Index GetANZ( cIndex wFi , cIndex strt = 0 ,
  			 Index stp = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

  virtual void GetADesc( cIndex wFi , int *Abeg , int *Aind , double *Aval ,
 			  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

  virtual NDOSolver *GetNDOSolver( void ) override;

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SETTING LAMBDA ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Setting Lambda
   @{ */

  virtual void SetLambda( cLMRow Lmbd = 0 ) override;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  virtual void SetLamBase( cIndex_Set LmbdB = 0 , cIndex LmbdBD = 0 ) override;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  virtual bool SetPrecision( HpNum Eps ) override;

/*@} -----------------------------------------------------------------------*/
/*------------------------ METHODS FOR COMPUTING Fi() ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Computing Fi()
   @{ */

   virtual HpNum Fi( cIndex wFi = Inf<Index>() ) override;

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/


   virtual bool NewGi( cIndex wFi = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   virtual Index GetGi( SgRow SubG , cIndex_Set &SGBse ,
			cIndex Name = Inf<Index>() ,
			cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   virtual HpNum GetVal( cIndex Name = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   virtual void SetGiName( cIndex Name ) override;

/*@} -----------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading other results
   @{ */


   virtual HpNum GetLowerBound( cIndex wFi = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   virtual FiStatus GetFiStatus( Index wFi = Inf<Index>() ) override;

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
   @{ */

   virtual void Deleted( cIndex i = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   virtual void Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm ,
			   cIndex NwNm ) override;

/*@} -----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~FakeFiOracle()
   {
	GiNameVcblr.clear();
    }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

  BundleSolver *bslv;
  std::vector< std::tuple< Index , Index , bool > >  GiNameVcblr; /* vocabulary
     for handling the items name; this is done to map the item name
     from the FiOracle to that of C05Function.  */

  Index last_c05;

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class FakeFiOracle )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  void InitMP( void );

  bool FindNext( Index &wFi );

/*--------------------------------------------------------------------------*/

  Index BStrategy( cIndex wFi );

/*--------------------------------------------------------------------------*/

  Index FindAPlace( cIndex wFi );

/*--------------------------------------------------------------------------*/

  void AggregateZ( cHpRow Mlt , cIndex_Set MBse , Index MBDm ,
  			  cIndex wFi , cIndex whr );

/*--------------------------------------------------------------------------*/

  HpNum Heuristic1( void );

  HpNum Heuristic2( void );

/*--------------------------------------------------------------------------*/
 /**< Remove all the items from the bundle, except the (sub)gradient of the
     linear 0-th component of Fi). */

  void RemoveItems( void );

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  void guts_of_destructor( );

/*--------------------------------------------------------------------------*/

  void ReSetAlg( unsigned char RstLvl = 0 );

  /**< Resets the internal state of the Bundle algorithm. Since several
     different things can be reset independently, RstLvl is coded bit-wise:

     - bit 0: if 0, all the algorithmic parameters are reset to the default
       values read by the stream/set by SetPar(), while if 1 they are left
       untouched;

     - bit 1: if 0 the current point is reset to the all-0 vector, while if
       1 it is left untouched;

     - bit 2: if 0, all the subgradients are removed from the bundle, except
       the constant (sub)gradient of the linear 0-th component, while if 1
       the subgradients are left there;

     - bit 3: if 0, all the constraints are removed from the bundle, while
       if 1 the constraints are left there.

     - bit 4: if 0 the value of Fi() in the current point is reset to HpINF
       (i.e., unknown), while if 1 it is left untouched; note that resetting
       the current point [see bit 1] has this as a side-effect, regardless to
       the value of bit 4. */

/*--------------------------------------------------------------------------*/

  void Delete( cIndex i );

/*--------------------------------------------------------------------------*/

  void UpdtaBP3( void );

  void CmptaBPX( void );

/*--------------------------------------------------------------------------*/

  bool IsOptimal( double eps = 0 ) const;

/*--------------------------------------------------------------------------*/

  bool CheckAlfa( const bool All = false );

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

  Index MBDim;      // number of items in the optimal multiplier base
  Index aBP3;       // current max number of items to be fetched
  Index aBP4;       // min number of items to be fetched

  FakeFiOracle FakeFi;  ///< a pointer to a FakeFiOracle object

SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class BundleSolver )

/*@}  end( group( Solver_CLASSES ) ) ---------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

/* template< class MCFC >
void BundleSolver< MCFC >::process_outstanding_Modification( void )
{
 }  // end( BundleSolver::process_outstanding_Modification ) */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BundleSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File BundleSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
