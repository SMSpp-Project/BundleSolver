/*--------------------------------------------------------------------------*/
/*---------------------- File BundleFSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the BunldeSolver class, which implements the Solver
 * interface, in particular in its CDASolver version, using a "Generalized
 * Bundle" algorithm.
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
 * or
 *
 *  A. Frangioni "Standard Bundle Methods: Untrusted Models and Duality"
 *  Technical Report, Dipartimento di Informatica, Università di Pisa, 2018
 *
 * available at
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#NDOB18
 * \endlink
 *
 * BunldeSolver is capable of solving any Block such that:
 *
 * - only has "continuous" ColVariable (is_integer() == false);
 *
 * - has no Constraint, except possibly BoxConstraint (upper and/or lower
 *   bounds on the ColVariable);
 *
 * - either its Objective is a FRealObjective containing a C05Function, and
 *   it has no sub-Block;
 *
 * - or its Objective is a is a FRealObjective containing a LinearFunction,
 *   and each of its sub-Block has no Constraint and Variable, and its
 *   Objective is a FRealObjective containing a C05Function.
 *
 * A special treatment is given to the case where some of the C05Function
 * actually are LagBFunction whose inner Block only contains ColVariable,
 * whose Objective is a is a FRealObjective containing a LinearFunction,
 * ans whose Constraint are linear (either FRowConstraint containing a
 * LinearFunction, or BoxConstraint). These are passes to the Master Problem
 * of the Bundle as "easy components" see
 *
 *   A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
 *   with ``Easy'' Components: Applications to Multicommodity Network Design"
 *   Mathematical Programming 145(1), 133–161, 2014
 *
 * available at
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#MP11c
 * \endlink
 *
 * In that case, the LagBFunction is never evaluated, which means that there
 * is no need for a Solver to be attached to the inner Block.
 *
 * \version 0.10
 *
 * \date 23 - 11 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Matematica ed Informatica \n
 *         Universita' di Cagliari \n
 *
 * Copyright &copy 2019 by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BundleSolver
 #define __BundleSolver
                      /* self-identification: #endif at the end of the file */

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
#include <queue>

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
/// A CDASolver using a (Generalized) Bundle algorithm
/** BunldeSolver implements the Solver interface, in particular in its
 * CDASolver version, using a "Generalized Bundle" algorithm.
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
 * or
 *
 *  A. Frangioni "Standard Bundle Methods: Untrusted Models and Duality"
 *  Technical Report, Dipartimento di Informatica, Università di Pisa, 2018
 *
 * available at
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#NDOB18
 * \endlink
 *
 * BunldeSolver is capable of solving any Block such that:
 *
 * - only has "continuous" ColVariable (is_integer() == false);
 *
 * - has no Constraint, except possibly BoxConstraint (upper and/or lower
 *   bounds on the ColVariable);
 *
 * - either its Objective is a FRealObjective containing a C05Function, and
 *   it has no sub-Block;
 *
 * - or its Objective is a is a FRealObjective containing a LinearFunction,
 *   and each of its sub-Block has no Constraint and Variable, and its
 *   Objective is a FRealObjective containing a C05Function.
 *
 * A special treatment is given to the case where some of the C05Function
 * actually are LagBFunction whose inner Block only contains ColVariable,
 * whose Objective is a is a FRealObjective containing a LinearFunction,
 * ans whose Constraint are linear (either FRowConstraint containing a
 * LinearFunction, or BoxConstraint). These are passes to the Master Problem
 * of the Bundle as "easy components" see
 *
 *   A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
 *   with ``Easy'' Components: Applications to Multicommodity Network Design"
 *   Mathematical Programming 145(1), 133–161, 2014
 *
 * available at
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#MP11c
 * \endlink
 *
 * In that case, the LagBFunction is never evaluated, which means that there
 * is no need for a Solver to be attached to the inner Block. */

class BundleSolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from Function and C05Function.
 *
 *  @{ */

 using Index = Function::Index;
 using c_Index = Function::c_Index;

 using Range = Function::Range;
 using c_Range = Function::c_Range;

 using Subset = Function::Subset;
 using c_Subset = Function::c_Subset;

 using VarValue = Function::FunctionValue;
 using c_VarValue = Function::c_FunctionValue;

 using Vec_VarValue = Function::Vec_FunctionValue;
 using c_Vec_VarValue = Function::c_Vec_FunctionValue;

 using LinearCombination = C05Function::LinearCombination;
 using c_LinearCombination = C05Function::c_LinearCombination;

/*----------------------------- CONSTANTS ----------------------------------*/

 static constexpr Function::FunctionValue NaNshift
                             = std::numeric_limits<Function::FunctionValue>::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

 static constexpr Function::FunctionValue INFshift
                              = std::numeric_limits<Function::FunctionValue>::infinity();
 ///< convenience constexpr for "Infty"

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different types of algorithmic parameters
  * of "int" type that BundleSolver has in addition to these of CDASolver.
  * The value intLastBndSlvPar is provided so that the list can be easily
  * further extended by derived classes. */

 enum int_par_type_BndSlv {

 intBPar1 = CDASolver::intLastParCDAS ,
 ///< remove linearizations unused for more than this consecutive iterations

 intBPar2 ,  ///< max number linearizations per component

 intBPar3 ,  ///< max number of new linearizations per iteration per component

 intBPar4 ,  ///< min number of new linearizations per iteration per component

 intBPar6 ,  ///< control how the min/max number of new linearizations changes

 intMnSSC ,  ///< minimum number of consecutive Serious Steps

 intMnNSC ,  ///< minimum number of consecutive Null Steps

 inttSPar1 ,  ///< first t-strategy parameter

 intMaxNrEvls ,  ///< max number of function evaluation for each iteration

 intMPName,  ///< true == MP solver is QPPenalty, false == MP is OSiMPSolver

 intMPlvl ,  ///< log verbosity of Master Problem

 intQPmp1 ,  ///< MxAdd parameter for QPPenaltyMP solver only

 intQPmp2 ,  ///< MxRmv parameter for QPPenaltyMP solver only

 intOSImp1 , ///< algorithm type for OsiMP solver only

 intOSImp2 ,  ///< reduction parameter for OsiMP solver only

 intOSImp3 ,  ///< threads parameter for OsiMP solver only

 intRstAlg ,  ///< reset parameter

 intLastBndSlvPar  ///< first allowed new int parameter for derived classes
                   /**< Convenience value for easily allow derived classes
		    * to extend the set of int algorithmic parameters. */

 };  // end( int_par_type_BndSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different types of algorithmic parameters
  * of "double" type that BundleSolver has in addition to these of CDASolver.
  * The value intLastBndSlvPar is provided so that the list can be easily
  * further extended by derived classes. */

 enum dbl_par_type_BndSlv {
  dbltStar = dblLastParCDAS ,
  ///< optimality parameter: "scaling" of the linearizations

  dblRelMPAcc , ///< relative optimality accuracy for the Master Problem

  dblRMPAccSol , ///< relative feasibility accuracy for the Master Problem

  dblBPar5 ,   ///< see intBPar6 above

  dblm1 ,      ///< m1 factor in NS/SS decision

  dblm2 ,      ///< m2 factor in NS/SS decision

  dblm3 ,      ///< m3 factor in NS/SS decision

  dblmxIncr ,  ///< maximum increasing t-factor

  dblmnIncr ,  ///< minimum increasing t-factor

  dblmxDecr ,  ///< maximum decreasing t-factor

  dblmnDecr ,  ///< minimum decreasing t-factor

  dbltMaior ,  ///< maximum value of t

  dbltMinor ,  ///< minimum value of t

  dbltInit ,   ///< initial value of t

  dbltSPar2 ,  ///< numerical parameter for the long-term t-strategies

  dblCtOff ,   ///< cut-off value for QPPenaltyMP solver only

  dblLastBndSlvPar ///< first allowed new double parameter for derived classes
                   /**< Convenience value for easily allow derived classes
		    * to extend the set of double algorithmic parameters. */

  };  // end( dbl_par_type_BndSlv )

/*@} -----------------------------------------------------------------------*/
/*----------------- CONSTRUCTING AND DESTRUCTING BundleSolver --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing BundleSolver
 *  @{ */

 /// constructor: ensure every field is initialized
 /** Void constructor: does nothing special, except verifying that the
  * template argument derives from MCFClass. */

 BundleSolver( void ) : CDASolver() , FakeFi( this ) , Result( kUnEval ) ,
  NumVar( 0 ) , NrFi( 0 ) , SCalls( 0 ) , ParIter( 0 ) , FiEvaltns( 0 ) ,
  GiEvaltns ( 0 ) , NrEasy( 0 ) , LHasChgd( true ) ,
  tHasChgd( true ) , LowerBound( - Inf< VarValue >() ) , t( 0 ) ,
  Prevt( 0 ) , Sigma( 0 ) , DSTS( 0 ) , vStar( 0 ) , DeltaFi( 0 ) ,
  EpsU( 0 ) , CSSCntr( 0 ) , CNSCntr( 0 ) , TrueLB( false ) ,
  LBHasChgd( false ) , SSDone( true ) , MBDim( 0 ) , aBP3( 0 ) , 
  f_lf( nullptr ) , Master( nullptr ) , UpTrgt( 0 ) , LwTrgt( 0 ) ,
  UpFiBest( Inf< VarValue >() ) , MaxNrEvls( 0 ) , DeltaStar( 0 ) , NrmD( 0 )
 {
  // ensure all parameters are properly given their default value
  MaxIter = CDASolver::get_dflt_int_par( intMaxIter );
  MaxSol = CDASolver::get_dflt_int_par( intMaxSol );
  LogVerb = CDASolver::get_dflt_int_par( intLogVerb );
  BPar1 = dflt_int_par[ intBPar1 - intLastParCDAS ];
  BPar2 = dflt_int_par[ intBPar2 - intLastParCDAS ];
  BPar3 = dflt_int_par[ intBPar3 - intLastParCDAS ];
  BPar4 = dflt_int_par[ intBPar4 - intLastParCDAS ];
  BPar6 = dflt_int_par[ intBPar6 - intLastParCDAS ];
  MnSSC = dflt_int_par[ intMnSSC - intLastParCDAS ];
  MnNSC = dflt_int_par[ intMnNSC - intLastParCDAS ];
  tSPar1 = dflt_int_par[ inttSPar1 - intLastParCDAS ];
  MaxNrEvls = dflt_int_par[ intMaxNrEvls - intLastParCDAS ];
  MPName = dflt_int_par[ intMPName - intLastParCDAS ];
  MPlvl = dflt_int_par[ intMPlvl - intLastParCDAS ];
  MxAdd = dflt_int_par[ intQPmp1 - intLastParCDAS ];
  MxRmv = dflt_int_par[ intQPmp2 - intLastParCDAS ];
  algo = dflt_int_par[ intOSImp1 - intLastParCDAS ];
  reduction = dflt_int_par[ intOSImp2 - intLastParCDAS ];
  threads = dflt_int_par[ intOSImp3 - intLastParCDAS ];

  MaxTime = dflt_dbl_par[ dblMaxTime - dblLastParCDAS ];
  RelAcc = dflt_dbl_par[ dblRelAcc - dblLastParCDAS ];
  AbsAcc = dflt_dbl_par[ dblAbsAcc - dblLastParCDAS ];
  RAccSol = dflt_dbl_par[ dblRAccSol - dblLastParCDAS ];
  AAccSol = dflt_dbl_par[ dblAAccSol - dblLastParCDAS ];
  tStar = dflt_dbl_par[ dbltStar - dblLastParCDAS ];
  RelMPAcc = dflt_dbl_par[ dblRelMPAcc - dblLastParCDAS ];
  RMPAccSol = dflt_dbl_par[ dblRMPAccSol - dblLastParCDAS ];
  m1 = dflt_dbl_par[ dblm1 - dblLastParCDAS ];
  m2 = dflt_dbl_par[ dblm2 - dblLastParCDAS ];
  m3 = dflt_dbl_par[ dblm3 - dblLastParCDAS ];
  mxIncr = dflt_dbl_par[ dblmxIncr - dblLastParCDAS ];
  mnIncr = dflt_dbl_par[ dblmnIncr - dblLastParCDAS ];
  mxDecr = dflt_dbl_par[ dblmxDecr - dblLastParCDAS ];
  mnDecr = dflt_dbl_par[ dblmnDecr - dblLastParCDAS ];
  tMaior = dflt_dbl_par[ dbltMaior - dblLastParCDAS ];
  tMinor = dflt_dbl_par[ dbltMinor - dblLastParCDAS ];
  tInit = dflt_dbl_par[ dbltInit - dblLastParCDAS ];
  tSPar2 = dflt_dbl_par[ dbltSPar2 - dblLastParCDAS ];
  CtOff = dflt_dbl_par[ dblCtOff - dblLastParCDAS ];
  }

/*--------------------------------------------------------------------------*/
 /// destructor: cleanly detaches the BundleSolver from the Block

 virtual ~BundleSolver() { set_Block( nullptr ); }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// set the "int" paramaters of BundleSolver
 /** Set the "int" paramaters specific of BundleSolver, together with the
  * paramaters of CDASolver that BundleSolver actually "listens to":
  *
  * - intMaxIter [Inf<int>]: maximum iterations for the next call to solve()
  *
  * - intMaxSol [1]: maximum number of different solutions to report. Since
  *                  Bundle methods are "almost" monotone ones, they typically
  *   produce only one solution, which is the stability center at termination.
  *   However, this is "almost" true: in fact, that solution may not be the
  *   best one ever found. If intMaxSol > 1, also the best solution wll to
  *   be kept and separately reported (assuming the stability center at
  *   termination is not it).
  *
  * - intLogVerb [0]: "verbosity" of the BundleSolver log
  *
  * - intBPar1 [10]: if an item has had a zero multiplier for the last
  *                  intBPar1 steps, it is eliminated; if intBPar1 is "too
  *   small" precious information may be lost, but keeping the "bundle" small
  *   obviously makes the Master Problem cheaper
  *
  * - intBPar2 [100]: maximum number of linearizations that are kept *for
  *                   each C05Function*; hence, the maximum total number
  *   is intBPar2 * < number of C05Function >. Note that intBPar2 must be
  *   >= 2, with the "poorman's" case intBPar2 == 2 forcing the
  *   BundleSolver to perform aggregation at every iteration (for every
  *   C05Function). Of course, keeping the "bundle" small makes the Master
  *   Problem cheaper, but on the other hand acquiring enough first-order
  *   information is typically the name of the game, hence keeping this
  *   value too low can have a dramatc effect on convergence speed that
  *   can easily counterbalance any improvement in Master Problem cost.
  *
  * - intBPar3 [1]: maximum number of new linearizations to be fetched from
  *                 each (non-easy) C05Function at each function evaluation
  *
  * - intBPar4 [1]: minmum number of new linearizations to be fetched from
  *                 each (non-easy) C05Function at each function evaluation
  * - intBPar6 [0]: together with the double parameters dblBPar3, dblBPar4
  *                 and dblBPar5, controls how the actual number of
  *   linearization that are requested to the C05Function evolves as the
  *   algorithm proceeds; note that what varies in practice is the maximum
  *   number, as it is always legal for the C05Function to refuse giving
  *   other items, although the BundleSolver will complain and stop if less
  *   than BPar4 are given. In the Bundle code, the number
  *
  *      EpsU = Sigma + D_{tStar}*( z* ) / max( | FiVal | , 1 ) ,
  *
  *   where Sigma = \sum_i Fi[ i ]_{B,Lambda}*( z[ i ]* ) + \sigma_L( w ) and
  *   z* = - Sum_i z[ i ]* is the optimal solution of the stabilized Dual
  *   Master Problem, and FiVal =  \sum_i Fi[ i ]( Lambda ) is the function
  *   value in the current point, is used as an estimate of the relative
  *   gap between the current and the optimal solution; that is, IsOptimal()
  *   returns true if EpsU <= RAccSol. Thus, the number RAccSol / EpsU is
  *   always smaller than one, and typically increases as the algorithm
  *   proceeds. Depending on the value of intBPar6, the following formulae
  *   for the actual value of BPar3, aBP3, are used:
  *
  *    0: aBP3 is set to BPar3 and never changed;
  *
  *    1: if BPar5 > 0 then aBP3 is initialized to BPar4 and increased
  *        every BPar5 iterations, while if BPar5 <= 0 then aBP3 is
  *        initialized to BPar3 and decreased every - BPar5 iterations;
  *
  *    2: aBP3 is set to
  *       ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 * ( RAccSol / EpsU )
  *
  *    3: aBP3 is set to
  *       ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 / sqrt( EpsU / RAccSol )
  *
  *    4: aBP3 is set to
  *       ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 / log10( EpsU / RAccSol )
  *
  * - intMnSSC [0]: minimum number of consecutive SS with the same t that
  *                 have to be performed before t is allowed to grow
  *
  * - intMnNSC [3]: minimum number of consecutive NS with the same t that
  *                 have to be performed before t is allowed to diminish
  *
  * - inttSPar1 [12]: select the t-strategy used. This field is coded
  *                   bit-wise in the following way.
  *   The first two bits control which heuristics are used to compute a new
  *   value of t when increasing/decreasing it. There are two heuristics
  *   avaliable, H1 and H2, both based on a quadratic interpolation of the
  *   function but differing in which derivative is used: H1 uses the
  *   derivative in the new tentative point, and it is guaranteed to produce
  *   a value greater than the current one if and only if the scalar product
  *   between the direction and the newly obtained subgradient is < 0
  *   (indicating that a longer step along the same direction could have been
  *   advantageous), while H2 uses the derivative in the current point and it
  *   does not possess this property. The value of the first two bits of
  *   tSPar1 has the following meaning:
  *
  *   bit 0:  which heuristic is used to increase t: 0 = H1, 1 = H2
  *
  *   bit 1:  which heuristic is used to decrease t: 0 = H2, 1 = H1
  *
  *   The following bits of inttSPar1 tell which long-term t-strategy is used,
  *   with the following values:
  *
  *    0 (+ 0):  none, only the heuristics are used
  *
  *    1 (+ 4):  the "soft" long-term t-strategy is used: the value EpsU is
  *              mantained (cf. intBPar6) which estimates the current
  *              relative error, and decreases of t are inhibited whenever
  *              v < tSPar2 * EpsU * | FiVal |
  *
  *    2 (+ 8):  the "hard" long-term t-strategy is used: the value EpsU is
  *              mantained as above, and t is increased whenever
  *              v < tSPar2 * EpsU * | FiVal |
  *
  *    3 (+12):  the "balancing" long-term t-strategy is used, where the two
  *              terms D*_t( -z* ) and Sigma* are kept of "roughly the same
  *              size": if D*_1( -z* ) <= tSPar2 * Sigma* then t increases
  *              are inhibited (increasing t causes a decrease of D*_1( -z* )
  *	         that is already small), if tSPar2 * D*_1( -z* ) >= Sigma*
  *              then t decreases are inhibited (decreasing t causes an
  *              increase of D*_1( -z* ) that is already big)
  *
  *    4 (+16):  the "endgame" t-strategy is used, where if D*_1( -z* ) is
  *              "small" (~ 1/10 of the current absolute epsilon) t is
  *              decreased no matter what the other strategies dictated.
  *              The rationale is that we are "towards the end" of the
  *              optimization and here t needs decrease. However, note that
  *              having D*_1( -z* ) "small" is no guarantee that we actually
  *              are at the end, especially if the oracle dynamically
  *              generates its variables, so use with caution
  *
  * - intMaxNrEvls [2]: max number of function evaluation for each iteration
  *
  * - intMPName [1]: bit-wise encoding of which MPSolver is used:
  *                  bit 0: 0 = QPPenalty, 1 = OSiMPSolver
  *                  bit 1 = 1 OsiCpxInterface, bit 1 = 0 OsiCLPInterface
  *                  bit 2 = 1 Quadratic, bit 2 = 0 BoxStep
  *
  * - intMPlvl [0]: log verbosity of Master Problem solver
  *
  * - intQPmp1 [0]: MxAdd parameter ( for QPPenaltyMP solver only )
  *
  * - intQPmp2 [0]: MxRmv parameter ( for QPPenaltyMP solver only )
  *
  * - intOSImp1 [4]: algorithm type ( for OsiMP solver only )
  *
  * - intOSImp2 [0]: reduction parameter ( for OsiMP solver only )
  *
  * - intOSImp3 [1]: number of threads ( for OsiMP solver only )
  *
  * - intRstAlg [ ]: parameter to handle the reset of the algorithm
  *                  bit-wise coded:
  *
  *                  0 bit -> true if don't reset algorithmic parameters
  *                  1 bit -> true if don't reset current point
  *                  2 bit -> true if don't reset subgradients
  *                  3 bit -> true if don't reset constraints
  *                  4 bit -> true if don't reset FiVals
  *                  5 bit -> true if don't get an initial point  */

 void set_par( const idx_type par , const int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the "double" paramaters of BundleSolver
 /** Set the "double" paramaters specific of BundleSolver, together with the
  * paramaters of CDASolver that BundleSolver actually "listens to":
  *
  * - dblMaxTime [Inf<double>()]: maximum CPU time for the next call to
  *                               compute(), in seconds
  *
  * - dblRelAcc [1e-6]: relative accuracy for declaring a solution optimal
  *                     (the "easy part", see dbltStar below for the
  *                     "complicated part")
  *
  * - dblAbsAcc [Inf<double>()]: absolute accuracy for declaring a solution
  *                              optimal; if INF, it is disabled; see dbltStar
  *   below for a detailed account on how this is used
  *
  * - dblRAccSol [Inf<double>()]: maximum relative error in any reported
  *                               solution; it is used in the stopping
  *   condition if dbltStar < 0 (see below for details);
  *
  * - dblAAccSol [Inf<double>()]: maximum absolute error in any reported
  *                               solution; it is used in the stopping
  *   condition if dbltStar < 0 (see below for details);
  *
  * - dbltStar [1e+2]: optimality parameter related to subgradient scaling.
  *                    Proving that some point Lambda is optimal for a
  *   NonDifferentiable Optimization problem involves finding an all-0
  *   subgradient of the function at Lambda. If an all-0 vector is found in
  *   the epsilon-subdifferential of Lambda, then the point is
  *   epsilon-optimal. Note that if the minimization problem is subject to
  *   constraints, i.e., Fi() has to be minimized only on the points Lambda
  *   \in L, the latter being a convex set, then the above is referred to a
  *   subgradient of the "actual function" ( Fi + I_L )( Lambda ), where I_L
  *   is the indicator function of L (evaluating to 0 inside L and to +INF
  *   otherwise). In other words, one has to show that there exists a(n
  *   epsilon-)subgradient of Fi() at Lambda that, *after projection on the
  *   frontier of L*, is all-0. A general stopping condition requires that
  *   it finds an epsilon-subgradient g at Lambda such that
  *
  *      tStar * || g || + epsilon <= min( dblAbsAcc , dblRelAcc * | Fi | )
  *
  *   where Fi is the current estimate of the optimal solution value (the
  *   value of the objective at the current stability center), tStar is an
  *   estimate of the longest step that can be performed, and || || is a
  *   norm-like function. tStar is related to the "scaling" of  Fi(), and it
  *   can be seen as the longest possible step that one can perform along
  *   (the opposite of) any (epsilon-)subgradient and still achieve a
  *   decrease; this means that < - tStar * g , g > is an estimate of the
  *   maximum decrease in objective value that one can expect to achieve
  *   because of the mere existence of the non-0 (epsion-)subgradient g.
  *   Estimating tStar is clearly nontrivial, although in proximal and
  *   trust-region Bundle methods tStar should be "just consistently larger,
  *   but not too much, of the value of t that actually results in SS". In
  *   some cases, however, a different way to estimate tStar exists. In
  *   particular, in the Lagrangian case || g || is tied to the norm of the
  *   unfeasibility of the current primal (convexified) solution, and
  *   therefore an alternative stopping criterion can be
  *   
  *        epsilon <= min( dblAbsAcc , dblRelAcc * | Fi | )
  *
  *        || g || <= min( dblAbsAcc , dblRAccSol * | tStar | )
  *
  *   The above is implemented if tStar < 0 (whence the absolute value).
  *   Thus, one can directly provide the required absolute accuracy, or
  *   still use tStar to provide the order-of-magnitude of || g ||, and
  *   then use the relative accuracy dblRAccSol.
  *
  * - dblRelMPAcc [1e-8]: relative optimality accuracy for the Master Problem
  *
  * - dblRMPAccSol [1e-8]: relative feasibility accuracy for the Master
  *                        Problem
  *
  * - dblBPar5 [30]: parameter controlling the dynamic number of 
  *                  linearizations to be fetched from each oracle at each
  *   iteration, see intbPar6 for details  *
  *
  * - dblm1 [0.01]: m1 factor in the all-important NS/SS decision. This
  *                 factor sets the lower target for the objective function
  *   value at the new iterate Lambda1, with the following formula: if
  *   mi > 0, then
  *
  *        lower_target = UpFi + v* + m1 * Delta*
  *
  *   otherwise
  *
  *        lower_target = UpFi + ( 1 + m1 ) * v*
  *
  *   where UpFi is the current upper estimate of the value of the objective
  *   function at the current stability center Lambda, v* < 0 is the
  *   predicted decrease at the new iterate Lambda1, and Delta* > -v* > 0 is
  *   the optimal dual value of the Master Problem. Typically, lower_target
  *   should be "just a bit above UpFi + v*", which is the lowest possible
  *   value of Fi( Lambda1 ). Whenever a lower bound on Fi( Lambda1 ) is
  *   found that is >= lower_target, a Null Step can be safely declared.
  *
  * - dblm2 [0.99]: m2 factor in the all-important NS/SS decision. This
  *                 factor sets the upper target for the objective function
  *   value at the new iterate Lambda1, with the formula
  *
  *        upper_target = UpFi + ( 1 - m2 ) * v*
  *
  *   where UpFi is the current upper estimate of the value of the objective
  *   function at the current stability center Lambda and v* < 0 is the
  *   predicted decrease at the new iterate Lambda1. Typically, upper_target
  *   should be "just a bit below UpFi": whenever an upper bound on
  *   Fi( Lambda1 ) is found that is <= upper_target (which menas that a
  *   "sizable decrease" has been achieved), a Serious Step can be safely
  *   declared. Note that it must be 0 < m1 <= m2 < 1; m1 = m2 = 0 is
  *   theoretically possible, but not practically advisable, for
  *   polyhedral functions provided that both function values and v* are
  *   computed without numerical errors (which is typically impossible).
  *   Also, note that whenever m1 < m2 both a SS and a NS may be possible
  *   at the same time, in which case the BundleSolver will typically favor
  *   the SS. 
  *
  * - dblm3 [0.99]: factor governing the Noise Reduction for "unfaithful"
  *                 oracles that pretend to provide information with the
  *   required accuracy but in fact they do not. This results in negative
  *   linearization errors, and therefore possibly in delecting directions
  *   that are non-decreasing even for the model (hence even less so for
  *   the real functon). To avoid this, if the aggregate linearization
  *   error \sigma* is "too negative", i.e.,
  *
  *      \sigma* < - m3 * t * || z* ||^2
  *
  *   then a NR step is performed by increasing t (if this is still possible,
  *   otherwise error is given). Traditionally m3 < 0.5 was required, but
  *   the latest developments have shown what m3 ~= 1 is possible, and it
  *   would appear that keeping m3 close to one could be preferable in
  *   practice.
  *
  * - dblmxIncr [10]: maximum increasing t-factor: each time t is increased
  *                   (for whatever reason), the new value of t must anyway
  *   be <= t * dblmxIncr (t is the previous value). Clearly, dblmxIncr must
  *   be > 1.
  *
  * - dblmnIncr [1.5]: minimum increasing t-factor: each time t is increased
  *                    (for whatever reason), the new value of t must anyway
  *   be >= t * dblmnIncr (t is the previous value). Clearly, dblmnIncr must
  *   be > 1.
  *
  * - dblmxDecr [0.1]: maximum decreasing t-factor: each time t is decreased
  *                    (for whatever reason), the new value of t must anyway
  *   be >= t * dblmxDecr (t is the previous value). Clearly, dblmnIncr must
  *   be < 1.
  *
  * - dblmnDecr [0.66]: maximum decreasing t-factor: each time t is decreased
  *                     (for whatever reason), the new value of t must anyway
  *   be <= t * dblmnDecr (t is the previous value). Clearly, dblmnDecr must
  *   be < 1.
  *
  * - dbltMaior [1e+6]: maximum value of t
  *
  * - dbltMinor [1e-6]: minimum value of t
  *
  * - dbltInit [1]: initial value of t. Choosing the "right" initial value
  *                 of t can clearly help the BundleSolver to perform
  *   better in the initial iterations, although the t-strategies should see
  *   to the fact that blatantly wrong t values are rapidly corrected.
  *   Giving a reasonanble value to this parameter (and, consequently, to
  *   dbltMaior and dbltMinor) is in general nontrivial, but in practice a
  *   minor amount of tuning suffices to find reasonable values. Usually,
  *   there is a "right" order of magnitude for t, that is the one that is
  *   guessed by the t-heuristics during most of the run, even though the
  *   starting value is very different. Hence, a good setting for dbltInit
  *   is in that order of magnitude, while dbltMinor and dbltMaior should be
  *   set small/large enough to never enter into play. Also, a reasonable
  *   value for dbltInit typically provides good cues to a reasonable value
  *   for dbltStar (and vice-versa): a "good" value for dbltStar (i.e., one
  *   that actually ensures that the stopping point is RAccSol-optimal) is
  *   usually one or two orders of magnitude larger than a "good" tInit.
  *
  * - dbltSPar2 [1e-3]: numerical parameter for the long-term t-strategies;
  *                     see inttSPar1 for details
  *
  * - dblCtOff [1e-1]: cut-off value for pricing in QPPenaltyMP solver only
  */

 void set_par( const idx_type par , const double value ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the BundleSolver log

 void set_log( std::ostream *log_stream = nullptr ) override;

/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the MCF encoded by the current MCFBlock
 *  @{ */

 /// (try to) solve the MCF encoded in the MCFBlock

 int compute( bool changedvars = true ) override;

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 VarValue get_lb( void ) override { return( LwFiLmb[ NrFi ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 VarValue get_ub( void ) override
 {
  if( ( MaxSol > 1 ) && ( UpFiBest < UpRifFi[ NrFi ] ) )
   return( UpFiBest );
  else
   return( UpRifFi[ NrFi ] );
  }

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
  // TODO: do it!!
  if( ( MaxSol > 1 ) && ( UpFiBest < UpRifFi[ NrFi ] ) )
   throw( std::logic_error( "writing LmbdBst non implemented yet" ) );
  else
   throw( std::logic_error( "writing Lambda non implemented yet" ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the "current" dual solution
 virtual void get_dual_solution( Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 virtual bool new_var_solution( void ) override
 {
  if( ( MaxSol > 1 ) && ( UpFiBest < UpRifFi[ NrFi ] ) ) {
   // dirty trick: pretend that LmbdBst is not better than Lambda by
   // re-defining UpFiBest, so that next time Lambda is given
   // anyway, the value of UpFiBest is no longer used after that the
   // corresponding solution has been given and a new one if "generated"
   UpFiBest = UpRifFi[ NrFi ];
   return( true );
   }
  else
   return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual bool new_dual_solution( void )  override
 {
  return( false );
  }

/*--------------------------------------------------------------------------*/
/*
 virtual void set_unbounded_threshold( const VarValue thr ) override { }
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
 
 int get_dflt_int_par( const idx_type par ) const override
 {
  if( ( par >= intLastParCDAS ) && ( par < intLastBndSlvPar ) )
   return( dflt_int_par[ par - intLastParCDAS ] );
  else
   return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dflt_dbl_par( const idx_type par ) const override
 {
  if( ( par >= dblLastParCDAS ) && ( par < dblLastBndSlvPar ) )
   return( dflt_dbl_par[ par - dblLastParCDAS ] );
  else
   return( CDASolver::get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/
 
 int get_int_par( const idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dbl_par( const idx_type par ) const override;

/*--------------------------------------------------------------------------*/

 idx_type int_par_str2idx( const std::string & name ) const override
 {
  const auto it = int_pars_map.find( name );
  if( it != int_pars_map.end() )
   return( it->second );
  else
   return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type dbl_par_str2idx( const std::string & name ) const override
 {
  const auto it = dbl_pars_map.find( name );
  if( it != dbl_pars_map.end() )
   return( it->second );
  else
   return( CDASolver::dbl_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 const std::string & int_par_idx2str( const idx_type idx ) const override
 {
  if( ( idx >= intLastParCDAS ) && ( idx < intLastBndSlvPar ) )
   return( int_pars_str[ idx - intBPar1 ] );
  else
   return( CDASolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & dbl_par_idx2str( const idx_type idx ) const override
 {
  if( ( idx >= dblLastParCDAS ) && ( idx < dblLastBndSlvPar ) )
   return( dbl_pars_str[ idx - dblLastParCDAS ] );
  else
   return( CDASolver::dbl_par_idx2str( idx ) );
  }

/*@} -----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 using SIndex = int;                        ///< type for "signed" indices

 using Vec_SIndex = std::vector< SIndex >;  ///< a std::vector of SIndex

 using Vec_Bool = std::vector< bool >;      ///< a std::vector of bool

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
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

 void FormD( void );

/*--------------------------------------------------------------------------*/
 // Updates the out-of-base counters for all items in the Bundle.

 void UpdtCntrs( void );

/*--------------------------------------------------------------------------*/
 /* After a (succesfull) call to FormD(), sets the new tentative point Lambda1
  * as Lambda1 = Lambda + ( Tau / t ) * d. */

 void FormLambda1( HpNum Tau );

/*--------------------------------------------------------------------------*/
 /* Computes Fi( Lambda1 ), inserting the obtained items (subgradients or
  * constraints) in the bundle. Returns true <=> the newly obtained
  * information changes the solution of the MP. */

 bool FiAndGi( Index wFi );

/*--------------------------------------------------------------------------*/
 /* Move the current point to Lambda1. */

 void GotoLambda1( void );

/*--------------------------------------------------------------------------*/
 /// Eliminate outdated items, i.e., these with "large" out-of-base counter.

 void SimpleBStrat( void );

/*--------------------------------------------------------------------------*/

 void UpdtLowerBound( void );

/*--------------------------------------------------------------------------*/

 double BetaK( Index wFi );

/*--------------------------------------------------------------------------*/

 void Log1( void );

 void Log2( void );

/*--------------------------------------------------------------------------*/

 VarValue max_error( VarValue releps ) const {
  VarValue FiL = UpFiLmb[ NrFi ];
  if( ( FiL >= Inf< VarValue >() ) || ( FiL <= - Inf< VarValue >() ) )
   return( Inf< VarValue >() );
  if( FiL < 0 ) FiL = - FiL;
  if( FiL < 1 ) FiL = 1;
  return( std::min( releps * FiL , RelAcc ) );
  }

 VarValue max_error( void ) const { return( max_error( RelAcc ) ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorthmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -

 int MaxSol;       ///< maximum number of different solutions to report

 double RelAcc;    ///< relative accuracy for declaring a solution optimal
 double AbsAcc;    ///< absolute accuracy for declaring a solution optimal
 double RAccSol;   ///< maximum relative error in any reported solution
 double AAccSol;   ///< maximum absolute error in any reported solution
 double RelMPAcc;  ///< relative optimality accuracy for the Master Problem
 double RMPAccSol; ///< relative feasibility accuracy for the Master Problem

 Index MaxIter;    ///< maximum number of iterations
 double MaxTime;   ///< maximum time (in seconds) for each call to Solve()

 double tStar;     ///< optimality related parameter: "scaling" of Fi

 int LogVerb;       ///< "verbosity" of the log

 int BPar1;         ///< parameter for removal of items (B-strategy)
 int BPar2;         ///< max Bundle size
 int BPar3;         ///< max number of items fetched from Fi() at each call
 int BPar4;         ///< min number of items fetched from Fi() at each call
 double BPar5;      ///< control how the actual BPar3 changes over time
 int BPar6;         ///< control how the actual BPar3 changes over time

 double mxIncr;     ///< max increase t parameter
 double mnIncr;     ///< min increase t parameter
 int MnSSC;         ///< min good iterations to do a SS 
 double mxDecr;     ///< max decrease t parameter
 double mnDecr;     ///< min decrease t parameter
 int MnNSC;         ///< max bad iterations to do a NS 

 double m1;         ///< m1 parameter for deciding if a SS/NS
 double m2;         ///< m2 parameter for deciding if a SS/NS
 double m3;         ///< m3 parameters for deciding if a SS/NS

 double tMaior;     ///< max value for t
 double tMinor;     ///< min value for t
 double tInit;      ///< initial value for t

 int tSPar1;        ///< int parameter for long-term t-strategy
 double tSPar2;     ///< double parameter for long-term t-strategy

 Index MxAdd;       ///< max variables added per iteration in QPPenaltyMP
 Index MxRmv;       ///< max variables added per iteration in QPPenaltyMP

 double CtOff;      ///< "break" value for the pricing in MinQuad
 
 Index algo;        ///< algorithm type ( for OSIMPSolver only )
 Index reduction;   ///< pre-processing (reduction) ( for OSIMPSolver only )
 Index threads;     ///< number of threads ( for OSIMPSolver only )

 Index MPlvl;       ///< log verbosity of master problem

 // generic fields- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int Result;        ///< result of the latest call to Solve()

 Index NumVar;      ///< (current) number of variables
 Index NrFi;        ///< number of components of Fi()

 Index SCalls;      ///< nuber of calls to Solve() (the current included)
 Index ParIter;     ///< nuber of iterations in this run
 Index FiEvaltns;   ///< total number of Fi() calls
 Index GiEvaltns;   ///< total number of Gi() calls

 Vec_Bool IsEasy;   ///< tells which component of Fi is "easy"
 Index NrEasy;      ///< number of "easy" component of Fi

 std::vector<VarValue> Lambda;   ///< the current point

 std::vector<VarValue> Lambda1;  ///< the tentative point

 std::vector<VarValue> LmbdBst;  ///< the best point found so far

 bool LHasChgd;       /**< true if Lambda has changed since the latest call
		       * to FiAndGi(): allows repeated calls in the same
		       * Lambda, e.g. with increasing precision */
 bool tHasChgd;       ///< true if t has changed since the last MP

 Subset whisZ;     /**< the position in the bundle where the "aggregate
		       * subgradient" Z[ k ] of "component" k is kept in
		       * whisZ[ k ]; Inf<Index>() means it is not in the
		       * bundle */
 Subset whisG1;    ///< "representative subgradient" for each component
 Vec_VarValue ScPr1;  ///< ScalarProduct( dir , G[ WhIsG1[ k ] ] )
 Vec_VarValue Alfa1;  /**< linearization error of G[ WhIsG1[ k ] ] w.r.t. the
		       * current point Lambda. */
 Vec_VarValue DeltaAlfa;  ///< correction of Fi-values due to inexactness

 VarValue LowerBound;  ///< Lower Bound over (the various components of) Fi

 double t;             ///< the (tremendous) t parameter
 double Prevt;         ///< what t were before being changed for funny reasons

 double Sigma;         ///< Sigma*: convex combination of the Alfa's
 double DSTS;          /**< D*_{t*}( -z* ), the other part of the dual
			* objective */
 Vec_VarValue vStar;   ///< v*, the predicted improvement

 double DeltaFi;       ///< FiLambda - FiLambda1
 double EpsU;          ///< precison required by the long-term t-strategy

 int CSSCntr;          ///< counter of consecutive SS

 int CNSCntr;          ///< counter of consecutive NS

 std::priority_queue<Index> FreList;  ///< list of free positions

 Vec_SIndex OOBase;   /**< Out-Of-Base counters:
		       * = Inf<SIndex>() means no item is there
		       * = k > 0 means out of base since k iterations
		       * = 0 means in the current base but potentially
		       *   removable
		       * = a *finite* negative value - k means not
		       *    removable for the next k iterations: note that
		       *  some items in base may be such
		       * = - Inf<SIndex>() means unremovable */

 bool TrueLB;         /**< true if LowerBound is a "true" lower bound rather
		       * than a "conditional" one */
 bool LBHasChgd;      ///< true some LowerBound has changed
 bool SSDone;         ///< true if the laste step was a SS

 Subset FiStatus;

 int RstAlgPrm; // reset parameter bt-wise coded

/*--------------------------------------------------------------------------*/

 std::vector< std::pair < Index , LinearCombination > > zA;

 /* the vector of the pairs  important linearization name and the
    linear combination used to form it */

 std::vector< C05Function * > v_c05f; /* the vector of the components of the
                                         sum function */
 LinearFunction * f_lf;  ///< the 0-th component of the sum function

 MPSolver * Master;      ///< (pointer to) the Master Problem Solver

 std::vector<MILPSolver*> MILP_s; /* MILP solver used to read the
                                     easy part of the sub-problems */

 std::vector<ColVariable *> LamVcblr;    // the set of indices of Lambda

 int MPName;       // 0 MP solver == QPPenalty
                   // otherwise MP == OSiMPSolver
                   // bit 1 = 1 Cplex, bit 1 = 0 CLP
                   // bit 2 = 1 Quadratic, bit 2 = 0 BoxStep

 VarValue UpTrgt; // upper target
 VarValue LwTrgt; // lower target

 VarValue UpFiBest;   // Fi best value vector

 Vec_VarValue UpRifFi;    /* The value of Fi[ k ]() where the zero of the
			   * Cutting Plane models are fixed: it is ==
			   * FiLambda[ k ]() except when FiLambda[ k ]()
			   * == INF */

 Vec_VarValue UpFiLmb1;   ///< upper function value vector at the Lambda1
 Vec_VarValue LwFiLmb1;   ///< lower function value vector at the Lambda1

 Vec_VarValue UpFiLmb;    ///< upper function value vector at Lambda
 Vec_VarValue LwFiLmb;    ///< lower function value vector at Lambda

 Index MaxNrEvls;
 std::vector<Index> CurrNrEvls;

 double DeltaStar;
 double NrmD;

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

   virtual ~FakeFiOracle() { GiNameVcblr.clear();  }

/*--------------------------------------------------------------------------*/

   void initialize( void );

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

  BundleSolver *bslv;

  std::vector< std::tuple< Index , Index , bool > >  GiNameVcblr;

  /* vocabulary
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
 /** Remove all the items from the bundle, except the (sub)gradient of the
     linear 0-th component of Fi). */

 void RemoveItems( void );

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

 void guts_of_destructor( void );

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

/*--------------------------------------------------------------------------*/

 bool IsOptimal( double eps = 0 ) const;

/*--------------------------------------------------------------------------*/

 bool CheckAlfa( const bool All = false );

/*--------------------------------------------------------------------------*/

 Index get_index_of_component( Function * f )
 {
  if( f == f_lf )
   return( Inf< Index >() );

  const auto fit = std::find( v_c05f.begin() , v_c05f.end() , f );
  if( fit != v_c05f.end() )
   return( std::distance( v_c05f.begin() , fit ) );

  throw( std::logic_error( "Modifiction from unknonw Function" ) );
  }

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 Index MBDim;      // number of items in the optimal multiplier base

 Index aBP3;       // current max number of items to be fetched

 FakeFiOracle FakeFi;  ///< the FakeFiOracle object

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class BundleSolver )

/*@}  end( group( Solver_CLASSES ) ) ---------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BundleSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File BundleSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
