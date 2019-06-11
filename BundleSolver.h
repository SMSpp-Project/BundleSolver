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

// #include <Eigen/Dense>
// #include <Eigen/Sparse>

#include "MPSolver.h"
#include "FakeFiOracle.h"

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

 class FakeFiOracle;     // forward declaration of class FakeFiOracle

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
    Bundle algorithm.


    - aggiungere Incrementale (approssimato)

    - aggiungere doppia stabilizzazione, verificare formula per norme
      diverse dalla norma 2

    - parametro tolleranza assoluta subgradiente aggregato

       max |g_i| \leq \epsilon

    - test stopping con or tStar = 0, epsilon = 0

    - specificare proprieta' blocco a cui si attacca il solver

     */

class BundleSolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

  friend class FakeFiOracle;

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

  typedef Eigen::SparseVector<FunctionValue> SparseVector;
  ///< type used to store a sparse vector

/*--------------------------------------------------------------------------*/
 /// Public enum "extending" sol_type to a specific case of CDASolvers
 enum sol_type {
  kEILoopNow = kBothInfeasible + 1 ,

  kEIAbort,
  };

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different types of algorithmic parameters
  * of "int" type that BundleSolver should have. The value
  * intLastAlgPar is provided so that the list can be easily further extended
  * by derived classes. */

 enum int_par_type_BndSlv {
 intBPar1 = CDASolver::intLastParCDAS ,    ///<
       /**< If an item has had a zero multiplier [see ReadMult()
       * in NDOSolver.h] for the last BPar1 steps, it is eliminated;
 * if BPar1 is "too small" precious information may be lost,
 * but keeping the "bundle" small obviously makes the Master
 * Problem cheaper. */

 // intBPar2 ,    ///<
        /**< Maximum dimension of the bundle: has more or less the
        * Same "problems" as BPar1, but if the latter is well chosen
 * then BPar2 can be kept big while the "B-strategy" keeps the
 * actual number of items low. A small BPar2 can affect the convergence of
 * the algorithm, in theory as well as in practice, if aggregation is not
 * allowed. However, an unnecessarily large BPar2 may force the MP Solver to
 * allocate a large amount of memory without a real need.. */

 intBPar6 ,   ///<
         /**< These parameters control how the actual number of
         * subgradients/constraints (items) that are requested
 * to the FiOracle varies, between BPar4 and BPar3, as the
 * algorithm proceeds; note that what varies in practice is
 * the maximum number, as it is always legal for the FiOracle to refuse
 * giving other items, although the Bundle code will complain and stop if
 * less than BPar4 are given. In the Bundle code, the number
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

 intEStps ,    ///<
       /**< The evaluation of function to be minimized may be a
       * costly task: in many cases, it requires the solution
 * of a - possibly hard - optimization problem. Often,
 * time can be saved if the function is only approximately computed at the
 * beginning of the optimization process; of course, the computation should
 * become more and more "exact" as the optimization proceeds. The relative
 * precision required to the FiOracle [see SetPrecision() in FiOracle.h] is
 * initially set to EInit, and decreased down to EFnal by multiplying it by
 * EDcrs every EStps iterations. If the NDO algorithm allows it, EStps can be set
 * to 0 meaning that the precision is decreased only if necessary, i.e., when it is
 * impossible to proceed otherwise (the algorithm must have some way to
 * detect this). The precision is kept fixed if EInit == EFnal. A NDO solver
 * which is *not* capable of handling approximate computation of the function
 * can ignore the values of these parameters. */

 intMnSSC ,   ///<
       /**< minimum number of consecutive SS with the same t that
       * have to be performed before t is allowed to grow  */

 intMnNSC ,   ///<
       /**< minimum number of consecutive NS with the same t that
       * have to be performed before t is allowed to diminish  */

 inttSPar1 ,    ///<
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

 // intPPar1 ,    ///<
 // intPPar2 ,    ///<
 // intPPar3 ,    ///<
       /**< Parameters controlling the variables generator: "price
       * in" (discover if new variables have to be added) is done
 * all the first PPar1 iterations and then every PPar2 iterations; note
 * that the price in is done anyway each time convergence is detected.
 * If PPar2 == 0, all the variables are present from the
 * beginning (PPar1 is ignored if PPar2 == 0). A variable that has
 * been inactive for the last PPar3 pricings (this one included) is
 * eliminated: note that the "price out" operation is done every PPar2
 * iterations, so that a variable that is eliminated is likely to have been
 * inactive for (about) PPar2 * PPar3 iterations. For PPar3 == 1, a variable
 * is eliminated in the very pricing in which it is discovered to be zero
 * (and the direction saying that it would stay zero). If PPar3 == 0, variables
 *  are *never* removed. PPar3 is ignored if PPar2 == 0.  */

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
       * optimal for a NonDifferentiable
 * Optimization problem involves finding an all-0 subgradient
 * of the function at Lambda. If an all-0 vector is found in
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

 dblEInit ,    ///<
       /**< ABS( EInit ) is the initial, and *maximum*, precision
                   required to the oracle, but the sign tells how the
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

 dblEFnal ,    ///<
 dblEDcrs ,    ///<
       /**<  The other three parameters define a general formula that
       * sets the "regular mechanism" for changing the precision
 * along iterations. Note that EFnal has a completely
 * different meaning as the one postulated by the base class (smallest
 * precision) because that makes no sense: the "final" precision clearly
 * has to be RAccSol. In fact, a value larger than RAccSol would make it
 * impossible (in theory) to reach RAccSol-accuracy for the overall
 * optimization, and a value smaller than RAccSol is wasteful as a higher
 * precision than RAccSol is not required. The idea is that the precision
 * should improve along the iterations, and the "speed" at which this
 * happens is dictated by EStps and EFnal; however, one can also keep the
 * precision "fixed" by setting EStps == 0. In this case, having defined
 *
 *  EpsU = ( ReadDStart( tStar ) + ReadSigma() ) / ReadFiVal()
 *
 * the current estimate of the optimality measure, the formula is
 *
 *   precision = / ABS( EInit )          if EDcrs >= 0
 *               \ ABS( EDcrs ) * EpsU   if EDcrs < 0
 *
 * and this is kept fixed along all the iterations (except, EpsU is not
 * really fixed); only the "emergency mechanism" will increase it if this
 * is absolutely needed. If instead EStps != 0, having defined
 *
 *   opt = / ABS( EInit )   if EDcrs >= 0
 *         \ EpsU           if EDcrs < 0
 *
 *   h = / NrIter()   if EStps > 0    ,     k = ceil( h / ABS( EStps ) )
 *       \ NrSSs()    if EStps < 0
 *
 * the formula is:
 *
 *   precision = opt * / ABS( EDcrs )^{ EFnal * k }   if EFnal >= 0
 *                     \ ABS( EDcrs ) * k^{ EFnal }   if EFnal < 0
 *
 * (while ensuring precision <= ABS( EInit )). */

 dblBPar3 ,    ///< maximum number of new subgradients/constraints
 dblBPar4 ,    ///< minimum number of new subgradients/constraints
       /**< Maximum and minimum number of new subgradients/constraints (items)
        * to be fetched from the FiOracle for
 * each function evaluation. Two different ways are given for specifying
 * these numbers: positive values are (rounded up and) taken as
 * absolute values, while negative numbers are first multiplied by
 * FiOracle::GetNrFi()---the number of components of Fi()---(and then
 * rounded up); thus, the default vale "-1" stands for "one for each of the
 * components of Fi()". Clearly, the "finalized" value of BPar3 has to be
 * <= BPar2, and the "finalized" value of BPar4 has to be <= than that.  */

 dblBPar5 ,    ///<
       /**< see intbPar6 above */


 dblm1 ,    ///<
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

 dblm3 ,    ///<
       /**< A nevly obtained subgradient is deemed "useless" if
       * Alfa >= m3 * Sigma; in this case, if a NS has to be done,
 * t is decreased. This parameter is mostly critical: if no
 * "long-term" t-strategy [see tSPar1 below] is used, values < 2/3 usually
 * make t to decrease rather fast to tMinor [see below], possibly making the
 * algorithm to perform very short steps and therefore to converge very
 * slowly. When a "long-term" t-strategy is used, .9 may be a good value.*/

 dblmxIncr ,    ///<
 dblmnIncr ,    ///<
       /**< each time t grows, the new value of t is chosen in
       * the interval [t * mnIncr, t * mxIncr] (t is the
 * previous value) */

 dblmxDecr ,    ///<
 dblmnDecr ,    ///<
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

 dblMPEFsb ,    ///<
         /**< (relative) precision required to the MP Solver as
         * far as constraints satisfaction is concerned. */

 dblMPEOpt ,    ///<
         /**< (relative) precision required to the MP Solver as
         * far as optimality of the solution is concerned*/

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

 BundleSolver( ) : CDASolver()  {
    
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

 void SetMPSolver( MPSolver *MPS = 0 );

/**< Gives to the BundleSolver object a pointer to an object of class MPSolver
 that will be used as Master Problem Solver during the Bundle algorithm. */

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

 virtual OFValue get_lb( void ) override {  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual OFValue get_ub( void )  override { }

/*--------------------------------------------------------------------------*/

 virtual bool has_var_solution( void ) override
 {
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual bool has_dual_solution( void ) override
 {
  }

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
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual bool new_dual_solution( void )  override
 {
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

 virtual int EveryIteration( void );

 /**< This method is an "hook" for derived classes: it is called at Every
    Iteration, between the computation of the tentative direction and the
    computation of Fi(). It can serve to various purposes, primarly checking
    extra stoping conditions or interfering with the usual stopping conditions
    of the Bundle code: however, any kind of operation can be performed here
    inside, e.g. adding or removing variables from the problem [see
    [Add/Remove]Variable() above].
    More in general, this method can be used to merge the main cycle of the
    Bundle method within any other however complex code: the Bundle gives out
    the control at this time, and resumes its operations when EveryIteration()
    returns. The returned value influences the behaviour of the Bundle for the
    current iteration:

    kEINorm        the current iteration is continued normally;

    kEIAbort       the whole algorithm is aborted, and Solve() is immediately
                   terminated returning kAbort: this is useful for instance
                   to enforce new termination criteria;

    kEILoopNow     the current iteration is aborted, i.e. the stopping
                   condition is *not* checked, and Fi() is *not* called: the
 		  next iteration is immediately started, but the iterations
 		  count is *not* increased. This is useful e.g. if something
 		  has been changed in the data of the problem that suggests
 		  to try a new direction, like a new "active" or variable
 		  [see AddVariable() above] to be inserted;

    kEIContAnyway  the current iteration is continued normally but for the
                   fact that the stopping condition is *not* checked. */

/*--------------------------------------------------------------------------*/

 void FormLambda1( HpNum Tau );
 /* After a (succesfull) call to FormD(), sets the new tentative point Lambda1
    (a protected field of type LMRow) as Lambda1 = Lambda + ( Tau / t ) * d. */

/*--------------------------------------------------------------------------*/

 bool FiAndGi( void );

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

    void Log1( void );

    void Log2( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

int MaxSol;         ///< maximum number of different solutions to report
double RelAcc;      ///< relative accuracy for declaring a solution optimal
double AbsAcc;      ///< absolute accuracy for declaring a solution optimal
double UpCutOff;    ///< upper cutoff for stopping the algorithm
double LwCutOff;    ///< lower cutoff for stopping the algorithm
double RAccSol;     ///< maximum relative error in any reported solution
double AAccSol;     ///< maximum absolute error in any reported solution
double FAccSol;     ///< maximum constraint violation in any reported solution

Index MaxIter;      ///< maximum number of iterations
double MaxTime;     ///< maximum time (in seconds) for each call to Solve()

double tStar;       ///< optimality related parameter: "scaling" of Fi

double EInit;      ///< precision-related parameter: initial precision
double EFnal;      ///< precision-related parameter: final precision
double EDcrs;      ///< precision-related parameter: rate of decrease
int EStps;         ///< precision-related parameter: number of steps

int Result;        ///< result of the latest call to Solve()

Index NumVar;      ///< (current) number of variables
Index NrFi;        ///< number of components of Fi()

Index SCalls;      ///< nuber of calls to Solve() (the current included)
Index ParIter;     ///< nuber of iterations in this run
Index FiEvaltns;   ///< total number of Fi() calls
Index GiEvaltns;   ///< total number of Gi() calls

int LogVerb;       ///< "verbosity" of the log
// ostream *NDOLog;   ///< the output stream object for log purposes
// OPTtimers *NDOt;   ///< OPTtimer for timing purposes


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
double m3;            // see the description in the constructor

double tMaior;        // max value for t
double tMinor;        // min value for t
double tInit;         // initial value for t

int tSPar1;          // long-term t-strategy parameters
double tSPar2;        // see the description in the constructor

int PPar1;           // pricing related parameters
int PPar2;           // see the description in the constructor
int PPar3;

double MPEFsb;        // precision required to the MP Solver (feasibility)
double MPEOpt;        // precision required to the MP Solver (optimality)

Index MaxNumVar;     // maximum number of variables
Vec_Bool IsEasy;     // tells whether any component of Fi is "easy"

SparseVector Lambda;  // the current point
SparseVector Lambda1; // the tentative point
SparseVector LmbdBst; // the best point found so far

/*
Vec_Index LamBase;   // the set of indices of Lambda
Vec_Index Lam1Bse;   // the set of indices of Lambda1
Index LamDim;        // dimension of LamBase */

bool KpBstL;         // if LmbdBst has to be kept
bool BHasChgd;       // true if LamBase has changed during the latest
                     // pricing (never set to true if PPar2 == 0, unless
                     // at the very first call to the oracle)
bool LHasChgd;       // true if Lambda has changed since the latest call
                     // to FiAndGi(): allows repeated calls in the same
                     // Lambda, e.g. with increasing precision
bool tHasChgd;       // true if t has changed since the last MP

Vec_OFValue FiLambda;      // Fi[ k ]( Lambda )
Vec_OFValue FiBest;        // best value(s) of Fi found so far
Vec_OFValue FiLambda1;     // Fi[ k ]( Lambda1 )
Vec_OFValue RfrncFi;       // the value of Fi[ k ]() where the zero of the Cutting
                     // Plane models are fixed: it is == FiLambda[ k ]() but
                     // when FiLambda[ k ]() == HpINF
// double b0;         // the constant in the affine 0-th component of Fi
// we could get to know it, if it was useful (which it is not)

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

Vec_OFValue LowerBound;// Lower Bound over (the various components of) Fi

double t;             // the (tremendous) t parameter
double Prevt;         // what t were before being changed for funny reasons

double Sigma;         // Sigma*: convex combination of the Alfa's
double DSTS;          // D*_{t*}( -z* ), the other part of the dual objective
double vStar;         // v*, the predicted improvement
double Deltav;        // the "desired improvement" in the Fi-value

double DeltaFi;       // FiLambda - FiLambda1
double EpsU;          // precison required by the long-term t-strategy

double EpsCurr;       // the precision currently required to the FiOracle
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
// Vec_Index InctvCtr;  // "out of base" counter for variables
// Vec_Index nBase;     // temporary

bool TrueLB;         // true if LowerBound is a "true" lower bound rather
                     // than just the "minus infinity"
bool LBHasChgd;      // true some LowerBound has changed
bool SSDone;         // true if the laste step was a SS

bool ItemsChgd;      // true if no itmes have been added to MP

Vec_Index FiStatus;

#if( NONMONOTONE )
 HpRow FiVals;       // Fi-values for the last NONMONOTONE SSs
#endif

/*--------------------------------------------------------------------------*/

 std::vector< std::pair < LinearizationName , LinearCombination > > zA;
 /* the vector of the pairs  important linearization name and the
    linear combination used to form it */

 std::vector< C05Function * > v_c05f; /* the vector of the components of the
                                         sum function */
 LinearFunction * linf; ///< the 0-th component of the sum function
 MPSolver *Master;    // (pointer to) the Master Problem Solver

 FakeFiOracle * FakeFi;  ///< a pointer to a FakeFiOracle object

 std::vector<ColVariable *> LamVcblr;    // the set of indices of Lambda

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
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  void InitMP( void );

/*--------------------------------------------------------------------------*/

  bool FindNextSG( Index &wFi );

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

  bool DoSS( void );

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

  void StrongCheckAlfa( void );

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

  Index MBDim;      // number of items in the optimal multiplier base
  Index aBP3;       // current max number of items to be fetched
  Index aBP4;       // min number of items to be fetched

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





