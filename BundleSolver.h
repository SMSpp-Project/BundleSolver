/*--------------------------------------------------------------------------*/
/*---------------------- File BundleFSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the BunldeSolver class, which implements the Solver
 * interface, in particular in its CDASolver version, using a "Generalized
 * Bundle" algorithm for the solution of convex nondifferentiable problems.
 *
 * The user is assumed to be familiar with the algorithm: refer to
 *
 *  A. Frangioni "Generalized Bundle Methods"
 *  SIAM Journal on Optimization 13(1), p. 117 - 156, 2002
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#SIOPT02
 * \endlink
 *
 * or
 *
 *  A. Frangioni "Standard Bundle Methods: Untrusted Models and Duality"
 *  in Numerical Nonsmooth Optimization: State of the Art Algorithms,
 *  A.M. Bagirov, M. Gaudioso, N. Karmitsa, M. Mäkelä, S. Taheri (Eds.),
 *  61--116, Springer, 2020
 *
 * available at
 *
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
 * whose Objective is linear (a FRealObjective containing a LinearFunction)
 * and whose Constraint are linear (either FRowConstraint containing a
 * LinearFunction, or BoxConstraint). These can be passed to the Master
 * Problem of the bundle algorithm as "easy components", see
 *
 *   A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
 *   with ``Easy'' Components: Applications to Multicommodity Network Design"
 *   Mathematical Programming 145(1), 133–161, 2014
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#MP11c
 * \endlink
 *
 * In that case, the LagBFunction is never evaluated, which means that there
 * is no need for a Solver to be attached to the inner Block.
 *
 * If the Block has multiple Objective (that is, it has sub-Block whose
 * Objective is a FRealObjective containing a C05Function), a very strong
 * assumption is required on them:
 *
 *     ALL THE Function IN THE Objective HAVE EXACTLY THE SAME SET OF
 *     "ACTIVE" Variable, ORDERED IN THE SAME WAY, AT ALL TIMES; THIS MEANS
 *     THAT IF THE SET OF "ACTIVE" Variable IS MODIFIED FOR ONE OF THE
 *     Function, IT MUST BE MODIFIED FOR ALL OF THEM AT THE SAME TIME
 *
 * The only exception is that
 *
 *     THE Objective OF THE Block CAN BE EMPTY, I.E., EITHER THERE IS NO
 *     Objective, OR THE FRealObjective HAS NO Function, OR THE
 *     LinearFunction IN THE FRealObjective HAS EXACTLY ZERO "ACTIVE"
 *     Variable; IN THE LATTER CASE, THE SET OF "ACTIVE" Variable IN THE
 *     LinearFunction MUST NEVER CHANGE
 *
 * To ensure that the rule about the list of "ACTIVE" Variable in the
 * (multiple) Objective is respected, an analogous very strong assumption is
 * made on the Modification that change that:
 *
 *     ALL THE Modification THAT CHANGE THE "ACTIVE" Variable MUST BE
 *     BUNCHED TOGETHER IN A SINGLE GroupModification. THIS MUST CONTAIN
 *     EXACTLY AS MANY Modification AS THERE ARE sub-Block (AND, THEREFORE,
 *     DIFFERENT OBJECTIVE), PLUS ONE IF THE (LinearFunction IN THE)
 *     Objective OF THE Block IS NOT EMPTY. ALL Modification MUST BE OF
 *     THE VERY SAME TYPE, I.E., EITHER ALL C05FunctionModVarsAddd, OR ALL
 *     C05FunctionModVarsRngd, OR ALL C05FunctionModVarsSbst, AND THEY MUST
 *     CHANGE THE "ACTIVE" Variable IN PRECISELY THE SAME WAY.
 *
 * \version 0.50
 *
 * \date 14 - 10 - 2020
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
 * Copyright &copy by Antonio Frangioni, Enrico Gorgone
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

#include <ctime>
#include <queue>

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

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 using namespace NDO_di_unipi_it;

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
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#SIOPT02
 * \endlink
 *
 * or
 *
 *  A. Frangioni "Standard Bundle Methods: Untrusted Models and Duality"
 *  in Numerical Nonsmooth Optimization: State of the Art Algorithms,
 *  A.M. Bagirov, M. Gaudioso, N. Karmitsa, M. Mäkelä, S. Taheri (Eds.),
 *  61--116, Springer, 2020
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#NDOB18
 * \endlink
 *
 * BunldeSolver is capable of solving any Block such that:
 *
 * - only has "continuous" ColVariable (is_integer() == false) that are
 *   either unconstrained below (lower bound == - INF) or non-negative
 *   (lower bound == 0);
 *
 * - has no Constraint, except possibly NNConstraint (bounds >= 0) or
 *   BoxConstraint (upper and/or lower bounds) on some of the ColVariable,
 *   but for BoxConstraint the lower bound can only be either 0 or -INF
 *   (general linear constraints would in principle be handled by the
 *   current master problem solver but the interface is not implemented);
 *
 * - the Objective of the Block and of all the sub-Block must be all
 *   FRealObjective containing a C05Function, and:
 *
 *   = either the Objective of the Block is a "generic" C05Function and the
 *     Block has no sub-Block;
 *
 *   = or the Objective of the Block is a LinearFunction (or is empty, which
 *     is taken to be the constantly 0 LinearFunction), each of its sub-Block
 *     has no Constraint and Variable and its Objective is a "generic"
 *     C05Function;
 *
 * - all the C05Function in all the Objective are either all convex or all
 *   concave (note that a linear function is both convex and concave and
 *   therefore is fine in both cases).
 *
 * A special treatment is given to the case where some of the C05Function
 * actually are LagBFunction whose inner Block only contains ColVariable,
 * whose Objective is linear (a FRealObjective containing a LinearFunction)
 * and whose Constraint are linear (either FRowConstraint containing a
 * LinearFunction, or OneVarConstraint). These can be passed to the Master
 * Problem of the bundle algorithm as "easy components", see
 *
 *   A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
 *   with ``Easy'' Components: Applications to Multicommodity Network Design"
 *   Mathematical Programming 145(1), 133–161, 2014
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#MP11c
 * \endlink
 *
 * In that case, the LagBFunction is never evaluated, which means that there
 * is no need for a Solver to be attached to the inner Block.
 *
 * If the Block has multiple Objective (that is, it has sub-Block whose
 * Objective FRealObjective containing a C05Function), a very strong
 * assumption is required on them:
 *
 *     ALL THE Function IN THE Objective HAVE EXACTLY THE SAME SET OF
 *     "ACTIVE" Variable, ORDERED IN THE SAME WAY, AT ALL TIMES; THIS MEANS
 *     THAT IF THE SET OF "ACTIVE" Variable IS MODIFIED FOR ONE OF THE
 *     Function, IT MUST BE MODIFIED FOR ALL OF THEM AT THE SAME TIME
 *
 * The only exception is that
 *
 *     THE Objective OF THE Block CAN BE EMPTY, I.E., EITHER THERE IS NO
 *     Objective, OR THE FRealObjective HAS NO Function, OR THE
 *     LinearFunction IN THE FRealObjective HAS EXACTLY ZERO "ACTIVE"
 *     Variable; IN THE LATTER CASE, THE SET OF "ACTIVE" Variable IN THE
 *     LinearFunction MUST NEVER CHANGE
 *
 * To ensure that the rule about the list of "ACTIVE" Variable in the
 * (multiple) Objective is respected, an analogous very strong assumption is
 * made on the Modification that change that:
 *
 *     IF THE Block HAS MORE THAN ONE C05Function, THAT IS, IT HAS A
 *     NON-EMPTY SET OF sub-Block AND THE (LinearFunction IN THE)
 *     Objective OF THE Block IS NOT EMPTY, THEN THE FunctionModVar THAT
 *     CHANGE THE "ACTIVE" Variable MUST BE BUNCHED TOGETHER IN A SINGLE
 *     GroupModification. THIS MUST CONTAIN EXACTLY AS MANY Modification AS
 *     THERE ARE C05Function, I.E., THE NUMBER OF sub-Block PLUS ONE IF 
 *     THE (LinearFunction IN THE) Objective OF THE Block IS NOT EMPTY. ALL
 *     Modification MUST BE OF THE VERY SAME TYPE, I.E., EITHER ALL
 *     C05FunctionModVarsAddd, OR ALL C05FunctionModVarsRngd, OR ALL
 *     C05FunctionModVarsSbst, AND THEY MUST CHANGE THE "ACTIVE" Variable IN
 *     PRECISELY THE SAME WAY.
 *
 * Failure to comply with the above rules will result in an exception being
 * thrown, either at set_Block() time (if the rules are violated from the
 * start), or when the offending Modification is processed. */

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

 static constexpr VarValue NaNshift
                              = std::numeric_limits< VarValue >::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

 static constexpr VarValue INFshift
                               = std::numeric_limits< VarValue >::infinity();
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

 intBPar7 ,  ///< how well-behaved BundleSolver is w.r.t. other Solver

 intMnSSC ,  ///< minimum number of consecutive Serious Steps

 intMnNSC ,  ///< minimum number of consecutive Null Steps

 inttSPar1 ,  ///< first t-strategy parameter

 intMaxNrEvls ,  ///< max number of function evaluation for each iteration

 intNoEasy,  ///< whether "easy" components are ignored

 intWZNorm,  ///< how to compute the norm of z*

 intMPName,  ///< whether the MP solver is QPPenalty or OSIMPSolver

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
  dblNZEps = dblLastParCDAS ,
               ///< parameter for declaring z* "almost 0"

  dbltStar ,   ///< optimality parameter: "scaling" of the linearizations

  dblMinNrEvls ,  ///< min fraction of components to evaluate at each iter.

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

 BundleSolver( void ) : CDASolver() , Result( kUnEval ) , NumVar( 0 ) ,
  NrFi( 0 ) , SCalls( 0 ) , ParIter( 0 ) , NrEasy( 0 ) , LHasChgd( true ) ,
  tHasChgd( true ) , MPchgs( 0 ) , t( 0 ) , f_global_LB( -INFshift ) ,
  Prevt( 0 ) , Sigma( 0 ) , DSTS( 0 ) , DeltaFi( 0 ) , EpsU( 0 ) ,
  CSSCntr( 0 ) , CNSCntr( 0 ) , TrueLB( false ) , SSDone( true ) ,
  f_wFi( 0 ) , f_lf( nullptr ) , f_convex( true ) , Master( nullptr ) ,
  UpTrgt( 0 ) , LwTrgt( 0 ) , RifeqFi( false ) , UpFiBest( INFshift ) ,
  UpFiLmb1def( 0 ) , LwFiLmb1def( 0 ) , UpFiLmbdef( 0 ) , LwFiLmbdef( 0 ) ,
  Fi0Lmb( 0 ) , Fi0Lmb1( 0 ) , DST( 0 ) , NrmD( 0 ) , NrmZ( 0 ) ,
  NrmZFctr( 1 ) , c_start( 0 ) , aBP3( 0 ) , FakeFi( this ) 
 {
  // ensure all parameters are properly given their default value
  MaxIter = CDASolver::get_dflt_int_par( intMaxIter );
  MaxSol = CDASolver::get_dflt_int_par( intMaxSol );
  EverykIt = CDASolver::get_dflt_int_par( intEverykIt );
  LogVerb = CDASolver::get_dflt_int_par( intLogVerb );
  BPar1 = dflt_int_par[ intBPar1 - intLastParCDAS ];
  BPar2 = dflt_int_par[ intBPar2 - intLastParCDAS ];
  BPar3 = dflt_int_par[ intBPar3 - intLastParCDAS ];
  BPar4 = dflt_int_par[ intBPar4 - intLastParCDAS ];
  BPar6 = dflt_int_par[ intBPar6 - intLastParCDAS ];
  BPar7 = dflt_int_par[ intBPar7 - intLastParCDAS ];
  MnSSC = dflt_int_par[ intMnSSC - intLastParCDAS ];
  MnNSC = dflt_int_par[ intMnNSC - intLastParCDAS ];
  tSPar1 = dflt_int_par[ inttSPar1 - intLastParCDAS ];
  MaxNrEvls = dflt_int_par[ intMaxNrEvls - intLastParCDAS ];
  NoEasy = bool( dflt_int_par[ intNoEasy - intLastParCDAS ] );
  WZNorm = bool( dflt_int_par[ intWZNorm - intLastParCDAS ] );
  MPName = dflt_int_par[ intMPName - intLastParCDAS ];
  MPlvl = dflt_int_par[ intMPlvl - intLastParCDAS ];
  MxAdd = dflt_int_par[ intQPmp1 - intLastParCDAS ];
  MxRmv = dflt_int_par[ intQPmp2 - intLastParCDAS ];
  algo = dflt_int_par[ intOSImp1 - intLastParCDAS ];
  reduction = dflt_int_par[ intOSImp2 - intLastParCDAS ];
  threads = dflt_int_par[ intOSImp3 - intLastParCDAS ];

  MaxTime = CDASolver::get_dflt_dbl_par( dblMaxTime );
  RelAcc = CDASolver::get_dflt_dbl_par( dblRelAcc );
  AbsAcc = CDASolver::get_dflt_dbl_par( dblAbsAcc );
  EveryTTm = CDASolver::get_dflt_dbl_par( dblEveryTTm );
  NZEps = dflt_dbl_par[ dblNZEps - dblLastParCDAS ];
  tStar = dflt_dbl_par[ dbltStar - dblLastParCDAS ];
  MinNrEvls = dflt_dbl_par[ dblMinNrEvls - dblLastParCDAS ];
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

  v_events.resize( max_event_number() );
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
 /** Gives the BundleSolver access to the Block it has to solve; note that
  * this does not register the BundleSolver among the Solver of the Block,
  * because the converse happens. Extensive checks are performed during
  * set_Block() to ensure that the Block does satify the requirements of
  * BundleSolver, and all the nontrivial internal data structures of
  * BundleSolver are set up.
  *
  * If \p block == nullptr, the BundleSolver is completely cleaned up and
  * prepared for either destruction or receiving an entirely unrelated Block
  * to solve. */

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
  *   best one ever found. If intMaxSol > 1, also the best solution will be
  *   kept and separately reported (assuming it is not the stability center
  *   at termination).
  *
  * - intEverykIt [0]: after how many iteration call the eEverykIteration
  *                    events
  *
  * - intLogVerb [0]: "verbosity" of the BundleSolver log
  *                   0 = no log
  *                   1 = only final state of the call and errors
  *                   2 = detailed step-by-step log
  *                   3 = as 2 + print every linearization added/removed
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
  *
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
  *   returns true if EpsU <= RelAcc. Thus, the number RelAcc / EpsU is
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
  *       ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 * ( RelAcc / EpsU )
  *
  *    3: aBP3 is set to
  *       ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 / sqrt( EpsU / RelAcc )
  *
  *    4: aBP3 is set to
  *       ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 / log10( EpsU / RelAcc )
  *
  * - intBPar7 [2]: This parameter, coded bit-wise, controls if BundleSolver
  *                 "tries to play nice" with any other Solver that may
  *   concurrently be using the same C05Function. The point is that each of
  *   these Solver is producing new linearizations, and possibly storing
  *   them in, or removing them from, the "finite resource" of the global
  *   pool(s) of the C05Function(s). Hence, what the BundleSolver does to the
  *   global pool may have an impact on the other Solver, if any. This
  *   parameter controls whether BundleSolver tries as hard as possible to
  *   avoid impacting the other Solver operations, or if it rather assumes to
  *   be "the only one" working with the C05Function, and therefore "treats
  *   the global pool as its exclusive property". To do so, BundleSolver
  *   handles the slot of the global pool in different ways according to the
  *   value in the first two bits of intBPar7 ( intBPar7 & 3 ):
  *
  *   = 0 means that BundleSolver will never override any position in the
  *     global pool unless it strictly needs to. This means that even if a
  *     linearization is removed from the bundle (the master problem), it is
  *     kept in the global pool of the corresponding component until the
  *     latter is completely full. Only then linearizations are removed,
  *     when necessary to make space for newly generated ones. Note that
  *     BundleSolver always "proceeds from left to right", i.e., selects the
  *     linearization in the global pool with smallest "name". This creates
  *     a sort of FIFO order whereby the oldest linearizations are removed
  *     first, which makes general sense.
  *
  *   = 1 means that BundleSolver will not immediately delete from the global
  *     pool a linearization that it removes from the bundle (the master
  *     problem). While the lineariztion is kept there, BundleSolver
  *     considers it "free", and can immediately after re-use that position
  *     to store a newly computed linearization. Again, the order is that if
  *     smaller names first, so if a linearization with "large name" is
  *     removed from the global pool it may take some time before it is
  *     actually overwritten by BundleSolver, thereby leaving it available
  *     to other Solver.
  *
  *   = 2 means that BundleSolver will immediately delete from the global
  *     pool any linearization that it removes from the bundle (the master
  *     problem). This makes sense if BundleSolver is the only Solver
  *     producing and consuming linearizations in these C05Function(s),
  *     since it allwas them to immediately delete all the memory (which may
  *     be significant) associated with that linearization in the global pool.
  *     However, if a linearization is found to be a "better copy" of a known
  *     one (the new linearization has the same linear part but a larger
  *     constant, and therefore provides a tighter constraint on the epigraph
  *     of the convex function, so that no Solver should complain if the
  *     weaker constraint is removed provided that the better one is added),
  *     still the old linearization is kept in the global pool (but not in the
  *     bundle) unless it is structly necessary to do so.
  *
  *   = 3 means that BundleSolver will immediately delete from the global
  *     pool any linearization that it removes from the bundle (the master
  *     problem); furthermore, if it finds a "better copy" of an existing
  *     linearization the new one immediately replaces the old one, in the
  *     global pool as well as in the bundle.
  *
  *   The bit 2 ( intBPar7 & 4 ) rather decides how BundleSolver reacts to
  *   Modification telling that some other Solver have generated a new
  *   linearization. If the bit is 0, then BundleSolver plainly ignores it,
  *   which is likely the best strategy if producing linearizations is
  *   "cheap". However, BundleSolver does record that a linearization is
  *   there: if ( intBPar7 & 3 ) == 0, it will avoid to touch it unless
  *   strictly necessary. If the bit is 1 instead, then BundleSolver will
  *   right away add the linearization to its bundle (the master problem),
  *   which is likely the best strategy if producing linearizations is
  *   "costly" and therefore it makes sense to profit from the effort that
  *   the C05Function(s) have done on behalf of the other Solver(s).
  *
  *   The bit 3 ( intBPar7 & 8 ) has a similar role for the initialization
  *   phase: if it is == 1, then BundleSolver will also scan the global pool
  *   of each component when it is attached to the Block, and immediately
  *   add to the bundle every linearization it finds there.
  *
  *   Of course, setting these bits to 1 has no impact if no other Solver is
  *   attached to the same Block, which is why the default value is 0 (to
  *   pair with the default value of 3 for the first two bits, indicating
  *   exclusive ownership).
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
  * - intMaxNrEvls [2]: max number of function evaluations for each
  *                (non-easy) C05Function for each iteration; multiple
  *   iterations may be needed if the C05Function has something like time
  *   or resource limits which causes it to stop its computation before
  *   having reached the required thresholds/accuracy. computation may be
  *   resumed multiple time to try to reach the required results, and this
  *   is the limit on how many times this will be attempted (for each
  *   non-easy C05Function) before giving up for good
  *
  * - intNoEasy [0]: if nonzero, instructs BundleSolver to disregard potential
  *                  "easy" components and to treat each as a "difficult" one
  *
  * - intWZNorm [2]: Proving that some point Lambda is epsilon-optimal for a
  *                  NonDifferentiable Optimization problem involves finding
  *   an all-0 epsilon-subgradient of the function at Lambda; see dbltStar for
  *   more comments. This is (tentatively) done at each iteration by, roughly
  *   speaking computing a convex combination (called z*) of the currently
  *   available epsilon-subgradients, and the corresponding convex combination
  *   of the linearization errors (called Sigma*). If z* is "almost 0" and
  *   Sigma* is "small", the algorithm can stop. While defining what "small"
  *   means for Sigma* is easy, since it is directly tied to the magnitude of
  *   the optimal function value (see comments to dblRelAcc, dblAbsAcc and
  *   again dbltStar), defining what "almost 0" means it is much less so. The
  *   typical format is that some norm of z* is smaller than some given
  *   threshold, but both the norm and the threshold may be nontrivial to set.
  *   This field is coded bit-wise and control these aspects. The first two
  *   bits control the choice of the norm, with the meaning
  *
  *    0 = INF-norm, 1 = 1-norm, 2 = 2-norm (Euclidean norm)
  *
  *   While the numerical value of the threshold is specified by the different
  *   parameter dblZNEps, the following two bits control how this is used,
  *   with the following meaning:
  *   
  *    0 = the parameter is taken as an absolute value (norm <= dblZNEps)
  *
  *    1 = the parameter is taken as a scaling factor of the corresponding
  *        norm of the fixed gradient of the linear 0-th component of the
  *        objective; if that is not present or all-0 the setting is
  *        equivalent to 0 (absolute value)
  *
  *    2 = the parameter is taken as a scaling factor of the corresponding
  *        norm of the "the first available full subgradient of the whole
  *        objective": the first time that one diagonal linearization is
  *        available for each component (which should happen quite early on)
  *        these are summed (if more then one is available for a component,
  *        the choice is arbitrary among them) and the norm of the obtained
  *        vector is used as the scaling factor for dblZNEps.
  *
  *   Besides for the stopping condition, these choices are crucial for the
  *   capability of BundleSolver to produce global valid lower bounds (for a
  *   minimization problem, upper bounds for a maximization one). Indeed,
  *   these can only be produced when z* is "0"; this is taken to mean
  *   "almost 0" in the specific sense dictated by this parameter together
  *   with dblZNEps.
  *
  * - intMPName [1]: bit-wise encoding of which MPSolver is used:
  *                  bit 0: 0 = QPPenalty, 1 = OSiMPSolver
  *                  bit 1: 1 = OsiCpxInterface, 0 = OsiCLPInterface
  *                  bit 2: 1 = Quadratic, 0 = BoxStep
  *                  bit 3: 1 = CheckIdentical( true ) is called, 0 = not
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
  * - intRstAlg [2]: parameter to handle the reset of the algorithm when
  *                  a new Block is set, bit-wise coded:
  *                  0 bit == 1 -> don't reset algorithmic parameters
  *                  1 bit == 1 -> set current point to using current values
  *                                of the Variable (otherwise reset to all-0)
  */

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
  *   below for a detailed recount on how this is used
  *
  * - dblEveryTTm [0]: periodicity of eEveryTTime events
  *
  * - dblNZEps [0]:    parameter controlling when the norm of the aggregated
  *                    subgradient z* is declared to be "almost 0". See
  *   intWZNorm for the details of how this is done in terms of which norm is
  *   used and how this constant is treated, as well as on the impact it has
  *   on the ability of BundleSolver to declare globally valid lower bounds
  *   (for a minimization problem, upper bounds for a maximization one).
  *   Choosing a very small value for dblZNEps may result in BundleSolver not
  *   being able to declare any global valid lower bound (especially if the
  *   alternative stopping criterion is used, see dbltStar), but on the other
  *   hand using a loose tolerance may result in declaring invalid global
  *   upper bound.
  *
  *   A relevant case where choosing a fair value for dblNZEps can be easier
  *   is that when the C05Function(s) is (are) Lagrangian function(s), since
  *   then z* is the violation of the relaxed constraints. Hence, the right
  *   value for dblNZEps (and its actual form, see intWZNorm) is the one
  *   corresponding to the tolerance required for declaring a solution of the
  *   original problem feasibile.
  *
  * - dbltStar [1e+2]: optimality parameter related to subgradient scaling.
  *                    Proving that some point Lambda is optimal for a convex
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
  *   frontier of L*, is all-0. Similar relationships hold for maximising
  *   concave NonDifferentiable functions with obvious changes.
  *
  *   A significant issue with this is that one typically cannot expect a true
  *   all-0 subgradient to be produced in practice; this eventually should
  *   happen but at the very least numerical errors have to be taken into
  *   account. Even worse, often the algorithm may find a "good" solution
  *   Lambda rapidly enough, but then struggle a long time to produce the
  *   corresponding all-0 subgradient that certifies its (almost) optimality.
  *   Of course, accepting a (projected) subgradient as being "almost 0" too
  *   early incurs the risk of stopping way before the optimal value is
  *   approached.
  *
  *   A stopping condition that may offer a good compromise between
  *   reliability (not stopping too far from the true optimum) and efficiency
  *   is
  *
  *     tStar * || z* ||^2_2 + Sigma* <= min( dblAbsAcc , dblRelAcc * | Fi | )
  *
  *   where Fi is the current estimate of the optimal solution value (the
  *   value of the objective at the current stability center), tStar is an
  *   estimate of the longest step that can be performed, z* is the current
  *   aggregated Sigma*-subgradient (both z* and Sigma* being produced by the
  *   master problem) and || ||_2 is the Euclidean norm. tStar is related to
  *   the "scaling" of  Fi(), and it can be seen as the longest possible step
  *   that one can perform along (the opposite of) any (epsilon-)subgradient
  *   and still achieve a decrease; this means that
  *
  *      < - tStar * z* , z* > - Sigma*
  *
  *   an estimate of the maximum decrease in objective value that one can
  *   expect to achieve because of the mere existence of the non-0
  *   Sigma*-subgradient z*. Estimating tStar is nontrivial, although in
  *   proximal and trust-region Bundle methods tStar should be "just
  *   consistently larger, but not too much, of the value of t that actually
  *   results in SS". This may make it possible to find good values for the
  *   parameter by experiments, which are often stable enough between
  *   different instances of the same kind. This is why this stopping
  *   criterion is offered. Furthermore, a proper choice of tStar may help in
  *   the on-line adjustment of the crucial proximal parameter t, see
  *   inttSPar1, dbltSPar2 and all the numerous corresponding parameters.
  *
  *   In some cases, estimating tStar is not easy, while it may be easier to
  *   come up with a direct estimate of "when the norm is small enough"; see
  *   dblNZEps and intWZNorm. Thus, the alternative stopping criterion
  *   
  *        Sigma* <= min( dblAbsAcc , dblRelAcc * | Fi | )
  *
  *        || z* || <= dblNZEps * < scaling factor >
  *
  *   is *always* used, where the choice of the norm and the scaling factor
  *   are controlled by dblNZEps. If tStar < 0, this is actually the *only*
  *   stopping criterion employed, with (the absolute value of) tStar then
  *   only playing a role in the t-strategies (if any). This may be
  *   appropriate e.g. if the value of tStar providing the best performances
  *   turns out not to give a reliable stopping test. If, instead, tStar > 0,
  *   then both stopping tests are employed in parallel; note that it is
  *   always easy to ensure that the second stopping condition never "wrongly
  *   fires" by just setting dblNZEps == 0 (although this will make it very
  *   difficult to ever generate a valid global upper bound).
  *
  * - dblMinNrEvls [0]: min number/fraction of non-easy C05Function evaluated
  *                     at each iteration. The solver can stop computing
  *   function values (and linearizations) as soon as the conditions required
  *   to declare either a SS or a NS are satisfied. If there are many non-easy
  *   C05Function, this may lead to many master problems being solved, which
  *   may not be convenient depending on the relative cost of the master
  *   problem and of the oracle. This parameter specifies the number/fraction
  *   of the total number of non-easy components that need be evaluated before
  *   the conditions for NS/SS conditions are even checked and the function
  *   values. If the parameter is >= 0, then the minimum number of evaluated
  *   components is just int( dblMinNrEvls ). If, instead, dblMinNrEvls < 0,
  *   then the minimum number of evaluated components is
  *   < number of non-easy C05Function > * ( - dblMinNrEvls ), i.e.,
  *   (-) dblMinNrEvls indicates the fraction of components that necessarily
  *   have to be evaluated.
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
  *   that actually ensures that the stopping point is RelAcc-optimal) is
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

/**@} ----------------------------------------------------------------------*/
/*----------------- METHODS FOR ACCESSING THE DATA OF THE Block ------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the data of the Block
 *
 * These methods provide convenient shortcuts for directly asking to the
 * BundleSolver some relevant data about the Block it is solving.
 *  @{ */

 /// returns the number of "components", i.e., C05Function in the objective
 /** Returns the total number of "components", i.e., the C05Function whose
  * sum (possibly together with one LinearFunction) makes up the objective of
  * the Block that the BundleSolver is solving. This method should not be
  * called if set_Block() has not been called, or has last been called with
  * nullptr argument. */

 Index n_components( void ) const { return( NrFi ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns a pointer to the i-th C05Function in the objective
 /** Returns a pointer to the i-th C05Function, i.e., the i-th term of the
  * sum of C05Function that (possibly together with one LinearFunction) 
  * makes up the objective of the Block that the BundleSolver is solving. 
  * It must ve 0 <= \p i <= n_components(). All returned pointers are not
  * nullptr provided that set_Block() has last been called with not nullptr
  * argument (otherwise this method should not be called). */

 C05Function * component( Index i ) const {
  #ifndef NDEBUG
   if( i >= NrFi )
    throw( std::invalid_argument( "wrong component number" ) );
  #endif
  return( v_c05f[ i ] );
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns a pointer to the LinearFunction in the objective
 /** Returns a pointer to the single LinearFunction that is summed with the
  * C05Function in the objective, if any. If the method returns nullptr,
  * there is no such LinearFunction (it is constantly 0, i.e., all its
  * coefficients are 0). This method should not be called if set_Block() has
  * not been called, or has last been called with nullptr argument. */

 LinearFunction * l_component( void ) const { return( f_lf ); }
  
/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR EVENTS HANDLING -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Set event handlers
 *
 *  BundleSolver heeds to all three "basic" types of events:
 *
 * - eBeforeTermination, called just before optimality stop (but not all
 *   other kinds of stop), with return action eForceContinue forcing one
 *   new iteration (master problem solution) to be performed;
 *
 * - eEverykIteration, called every intEverykIt iterations (if intEverykIt
 *   != 0), with possible return actions eStopOK and eStopError;
 *
 * - eEveryTTime, called every dblEveryTTm seconds (if dblEveryTTm != 0),
 *   with possible return actions eStopOK and eStopError;
 *
 * Of course, events have been set with set_event_handler() for them to be
 * called.
 *  @{ */

 EventID max_event_number( void ) const override { return( 3 ); }

/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the MCF encoded by the current MCFBlock
 *  @{ */

 /// (try to) solve the MCF encoded in the MCFBlock

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
 /// returns the number of calls to compute() (the current included)

 Index n_calls( void ) const { return( SCalls ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the number of iterations in the current call to compute()
 /** Returns the number of iterations in the current call to compute(), the
  * last one included. Note that this is the number of master problem
  * solutions, as clearly the number of function evaluations may be rather
  * different. */

 Index n_iter( void ) const { return( ParIter ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the number of evaluations of a component in the current compute()
 /** Returns the number of evaluations of the component \p i in the current
  * call to compute(). */

 Index n_f_eval( Index i ) const {
  #ifndef NDEBUG
   if( i >= NrFi )
    throw( std::invalid_argument( "wrong component number" ) );
  #endif
  return( CurrNrEvls[ i ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the elapsed CPU time since the start of the last compute()

 double elapsed( void ) const {
  return( ( std::clock() - c_start ) / CLOCKS_PER_SEC );
  }
 
/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 VarValue get_lb( void ) override {
  if( f_convex ) {
   if( f_global_LB > - INFshift )
    return( f_global_LB );

   return( TrueLB ? LowerBound.back() : - INFshift );
   }
  else
   if( ( MaxSol > 1 ) && ( UpFiBest < UpFiLmb.back() ) )
    return( - UpFiBest );
   else
    return( - UpFiLmb.back() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 VarValue get_ub( void ) override {
  if( f_convex )
   if( ( MaxSol > 1 ) && ( UpFiBest < UpFiLmb.back() ) )
    return( UpFiBest );
   else
    return( UpFiLmb.back() );
  else {
   if( f_global_LB > - INFshift )
    return( - f_global_LB );

   return( TrueLB ? - LowerBound.back() : INFshift );
   }
  }

/*--------------------------------------------------------------------------*/

 bool has_var_solution( void ) override { return( true ); }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_dual_solution( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*
 virtual bool is_var_feasible( void ) override { return( true ); }

 virtual bool is_dual_feasible( void ) override { return( true ); }
*/
/*--------------------------------------------------------------------------*/
 /// write the "current" solution

 void get_var_solution( Configuration *solc = nullptr ) override
 {
  if( ( MaxSol > 1 ) && ( UpFiBest < UpRifFi.back() ) ) {
   for( Index i = 0 ; i < NumVar ; i++ )
    LamVcblr[ i ]->set_value( LmbdBst[ i ] );
   }
  else
   for( Index i = 0 ; i < NumVar ; i++ )
    LamVcblr[ i ]->set_value( Lambda[ i ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the "current" dual solution

 void get_dual_solution( Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 bool new_var_solution( void ) override
 {
  if( ( MaxSol > 1 ) && ( UpFiBest < UpRifFi.back() ) ) {
   // dirty trick: pretend that LmbdBst is not better than Lambda by
   // re-defining UpFiBest, so that next time Lambda is given
   // anyway, the value of UpFiBest is no longer used after that the
   // corresponding solution has been given and a new one if "generated"
   UpFiBest = UpRifFi.back();
   return( true );
   }
  else
   return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool new_dual_solution( void )  override { return( false ); }

/*--------------------------------------------------------------------------*/
/*
  void set_unbounded_threshold( const VarValue thr ) override { }

  bool has_var_direction( void ) override { return( true ); }

  bool has_dual_direction( void ) override { return( true ); }

  void get_var_direction( Configuration *dirc = nullptr ) override {}

  void get_dual_direction( Configuration *dirc = nullptr ) override {}

  bool new_var_direction( void ) override { return( false ); }
  
  bool new_dual_direction( void ) override{ return( false ); }
*/

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/
/*
 virtual bool is_dual_exact( void ) const override { return( true ); }
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 c_Vec_VarValue & get_current_point( void ) const { return( Lambda ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 c_Vec_VarValue & get_tentative_point( void ) const { return( Lambda1 ); }
 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the BundleSolver
 *
 *  @{ */

 idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastBndSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastBndSlvPar ) );
  }

/*--------------------------------------------------------------------------*/
 
 int get_dflt_int_par( const idx_type par ) const override {
  if( ( par >= intLastParCDAS ) && ( par < intLastBndSlvPar ) )
   return( dflt_int_par[ par - intLastParCDAS ] );
  else
   return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dflt_dbl_par( const idx_type par ) const override {
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

 idx_type int_par_str2idx( const std::string & name ) const override {
  const auto it = int_pars_map.find( name );
  if( it != int_pars_map.end() )
   return( it->second );
  else
   return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type dbl_par_str2idx( const std::string & name ) const override {
  const auto it = dbl_pars_map.find( name );
  if( it != dbl_pars_map.end() )
   return( it->second );
  else
   return( CDASolver::dbl_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 const std::string & int_par_idx2str( const idx_type idx ) const override {
  if( ( idx >= intLastParCDAS ) && ( idx < intLastBndSlvPar ) )
   return( int_pars_str[ idx - intBPar1 ] );
  else
   return( CDASolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & dbl_par_idx2str( const idx_type idx ) const override {
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
 /* Performs the inner loop: repeatedly compute components up until the
  * conditions for either a SS or a NS (or both) are satisfied, or something
  * very bad happens (errors, out of time, ...). Sets MPchgs > 0 if the SS
  * and/or NS conditions are satisfied.
  *
  * It is virtual because this is precisely the point where a sequential and
  * a "basic" asynchronous implementation of the approach differ, and
  * therefore this is the obvious hook for a derived AsynchBundleSolver. */

 virtual void InnerLoop( void );

/*--------------------------------------------------------------------------*/
 /* Computes Fi( Lambda1 ), inserting the obtained items (subgradients or
  * constraints) in the bundle. Returns true <=> at least one item was
  * inserted. It also "sneakily" sets MPchgs if appropriate. */

 bool FiAndGi( Index wFi );

/*--------------------------------------------------------------------------*/

 void update_UpFiLambd1( Index wFi , VarValue nval )
 {
  if( UpFiLmb1[ wFi ] <= nval )
   return;

  if( UpFiLmb1[ wFi ] == INFshift )
   ++UpFiLmb1def;

  if( UpFiLmb1def == NrFi ) {
   ++UpFiLmb1def;  // all components + the sum computed
   UpFiLmb1[ wFi ] = nval;
   UpFiLmb1.back() = std::accumulate( UpFiLmb1.begin() , --(UpFiLmb1.end()) ,
				      Fi0Lmb1 );
   }
  else {
   if( UpFiLmb1def > NrFi )
    UpFiLmb1[ NrFi ] += nval - UpFiLmb1[ wFi ];
   UpFiLmb1[ wFi ] = nval;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void update_LwFiLambd1( Index wFi , VarValue nval )
 {
  if( LwFiLmb1[ wFi ] >= nval )
   return;

  if( LwFiLmb1[ wFi ] == -INFshift )
   ++LwFiLmb1def;

  if( LwFiLmb1def == NrFi ) {
   ++LwFiLmb1def;  // all components + the sum computed
   LwFiLmb1[ wFi ] = nval;
   LwFiLmb1.back() = std::max( LwFiLmb1.back() ,
			       std::accumulate( LwFiLmb1.begin() ,
						--(LwFiLmb1.end()) ,
						Fi0Lmb1 ) );
   }
  else {
   if( LwFiLmb1def > NrFi )
    LwFiLmb1[ NrFi ] += nval - LwFiLmb1[ wFi ];
   LwFiLmb1[ wFi ] = nval;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void update_UpFiLambd( Index wFi , VarValue nval )
 {
  if( UpFiLmb[ wFi ] <= nval )
   return;

  if( UpFiLmb[ wFi ] == INFshift )
   ++UpFiLmbdef;

  if( UpFiLmbdef == NrFi ) {
   ++UpFiLmbdef;  // all components + the sum computed
   UpFiLmb[ wFi ] = nval;
   UpFiLmb.back() = std::accumulate( UpFiLmb.begin() , --(UpFiLmb.end()) ,
				     Fi0Lmb );
   }
  else {
   if( UpFiLmbdef > NrFi )
    UpFiLmb[ NrFi ] += nval - UpFiLmb[ wFi ];
   UpFiLmb[ wFi ] = nval;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void update_LwFiLambd( Index wFi , VarValue nval )
 {
  if( LwFiLmb[ wFi ] >= nval )
   return;

  if( LwFiLmb[ wFi ] == -INFshift )
   ++LwFiLmbdef;

  if( LwFiLmbdef == NrFi ) {
   ++LwFiLmbdef;  // all components + the sum computed
   LwFiLmb[ wFi ] = nval;
   LwFiLmb.back() = std::max( LwFiLmb.back() ,
			      std::accumulate( LwFiLmb.begin() ,
					       --(LwFiLmb.end()) ,
					       Fi0Lmb ) );
   }
  else {
   if( LwFiLmbdef > NrFi )
    LwFiLmb[ NrFi ] += nval - LwFiLmb[ wFi ];
   LwFiLmb[ wFi ] = nval;
   }
  }

/*--------------------------------------------------------------------------*/
 /* Move the current point to Lambda1. */

 void GotoLambda1( void );

/*--------------------------------------------------------------------------*/
 /* Ensure that the linearization errors agree with the current point.  */

 void GotoLambda( void );

/*--------------------------------------------------------------------------*/
 /* Ensure that the linearization errors of the component k; if k >= NrFi,
  * do it for all. */

 void ResetAlfa( Index k );

/*--------------------------------------------------------------------------*/
 /* Eliminate outdated items, i.e., these with "large" out-of-base counter. */

 void SimpleBStrat( void );

/*--------------------------------------------------------------------------*/

 double BetaK( Index wFi );

/*--------------------------------------------------------------------------*/

 void Log1( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void Log2( void );

/*--------------------------------------------------------------------------*/

 VarValue eps_fi( VarValue fi , VarValue releps ) const {
  if( fi < 0 ) fi = - fi;
  if( fi < 1 ) fi = 1;
  return( releps * fi );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 VarValue max_error( VarValue fi , VarValue releps ) const {
  return( std::min( eps_fi( fi , releps ) , AbsAcc ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 VarValue max_error( VarValue releps ) const {
  c_VarValue FiL = UpFiLmb[ NrFi ];
  if( ( FiL >= Inf< VarValue >() ) || ( FiL <= - Inf< VarValue >() ) )
   return( Inf< VarValue >() );
  else
   return( max_error( FiL , releps ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 VarValue max_error( void ) const { return( max_error( RelAcc ) ); }

/*--------------------------------------------------------------------------*/

 void compute_NrmZFctr( void );

/*--------------------------------------------------------------------------*/
 /* Concave functions to be maximised are sneakily turned into convex
  * functions to be minimized inside by changing the sign of funcion values
  * and linearizations, but they have to be output with the right sign. */
 
 VarValue rs( const VarValue fv ) {
  return( f_convex ? fv : - fv );
  }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorthmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -

 int MaxSol;        ///< maximum number of different solutions to report
 Index MaxNrEvls;   ///< maximum total number of function evaluations

 double RelAcc;     ///< relative accuracy for declaring a solution optimal
 double AbsAcc;     ///< absolute accuracy for declaring a solution optimal
 double EveryTTm;   ///< periodicity of eEveryTTime events

 Index MaxIter;     ///< maximum number of iterations
 double MaxTime;    ///< maximum time (in seconds) for each call to Solve()

 double NZEps;      ///< parameter for declaring that the norm of z* is "0"

 double tStar;      ///< optimality related parameter: "scaling" of Fi

 double MinNrEvls;  ///< min fraction of components to evaluate at each iter.
 int LogVerb;       ///< "verbosity" of the log
 int EverykIt;      ///< periodicity of eEverykIteration events

 int BPar1;         ///< parameter for removal of items (B-strategy)
 int BPar2;         ///< max Bundle size
 int BPar3;         ///< max number of items fetched from Fi() at each call
 int BPar4;         ///< min number of items fetched from Fi() at each call
 double BPar5;      ///< control how the actual BPar3 changes over time
 int BPar6;         ///< control how the actual BPar3 changes over time
 int BPar7;         ///< if BundleSolver "plays nice" with other Solver

 double mxIncr;     ///< max increase t parameter
 double mnIncr;     ///< min increase t parameter
 Index MnSSC;       ///< min good iterations to do a SS 
 double mxDecr;     ///< max decrease t parameter
 double mnDecr;     ///< min decrease t parameter
 Index MnNSC;       ///< max bad iterations to do a NS 

 double m1;         ///< m1 parameter for deciding if a SS/NS
 double m2;         ///< m2 parameter for deciding if a SS/NS
 double m3;         ///< m3 parameters for deciding if a SS/NS

 double tMaior;     ///< max value for t
 double tMinor;     ///< min value for t
 double tInit;      ///< initial value for t

 int tSPar1;        ///< int parameter for long-term t-strategy
 double tSPar2;     ///< double parameter for long-term t-strategy

 bool NoEasy;       ///< true if easy components are ignored

 char WZNorm;       ///< how to compute the norm of z*
 
 int MPName;        /**< bit 0 = 0: MP solver == QPPenalty
		     * bit 0 = 1: MP == OSiMPSolver
		     * bit 1 = 1: Cplex, bit 1 = 0 CLP
		     * bit 2 = 1: Quadratic, bit 2 = 0 BoxStep
		     * + bit 3 = 1 (+8) = check for duplicates. */

 Index MxAdd;       ///< max variables added per iteration in QPPenaltyMP
 Index MxRmv;       ///< max variables added per iteration in QPPenaltyMP

 double CtOff;      ///< "break" value for the pricing in MinQuad
 
 Index algo;        ///< algorithm type ( for OSIMPSolver only )
 Index reduction;   ///< pre-processing (reduction) ( for OSIMPSolver only )
 Index threads;     ///< number of threads ( for OSIMPSolver only )

 Index MPlvl;       ///< log verbosity of master problem

 int RstAlgPrm;     ///< reset parameter, bit-wise coded

 // generic fields- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int Result;        ///< result of the latest call to Solve()

 Index NumVar;      ///< (current) number of variables
 Index NrFi;        ///< number of components of Fi()

 Index SCalls;      ///< number of calls to compute() (the current included)
 Index ParIter;     ///< number of iterations in this call to compute() 

 Vec_Bool IsEasy;   ///< tells which component of Fi is "easy"
 Index NrEasy;      ///< number of "easy" component of Fi

 Vec_VarValue Lambda;   ///< the current point

 Vec_VarValue Lambda1;  ///< the tentative point

 Vec_VarValue LmbdBst;  ///< the best point found so far

 bool LHasChgd;       /**< true if Lambda has changed since the latest call
		       * to FiAndGi(): allows repeated calls in the same
		       * Lambda, e.g. with increasing precision */
 bool tHasChgd;       ///< true if t has changed since the last MP

 char MPchgs;         ///< nonzero if we can prove no cycling will occur
                      /**< MPchgs == 1 means that the conditions for
		       * ensuring that no cycle will occur have been found
  * due to the function value (a SS can be done) or a diagonal linearization
  * (a NS can be done); MPchgs == 2 means that a vertical linearization
  * (cutting off Lambda1) has been found, and this by itself ensures no
  * cycling. */
 
 Subset whisZ;     /**< the position in the bundle where the "aggregate
		    * subgradient" Z[ k ] of component k is kept in
		    * whisZ[ k ]; Inf<Index>() == it is not in the bundle */
 std::vector< bool > Zvalid;  /**< Zvalid[ k ] == true if the item in position
			       * whisZ[ k ] is exactly Z[ k ] as computed by
 * the last master problee. Zvalid[ k ] == true ==> whisZ[ k ] < INF.
 * if Zvalid[ k ] == false and whisZ[ k ] < INF, then Z[ k ] had been
 * previously stored in position whisZ[ k ], but the master problem has
 * been re-solved since and therefore Z[ k ] is no longer current. */
 
 Subset whisG1;    ///< "representative subgradient" for each component
 Vec_VarValue ScPr1;  ///< ScalarProduct( dir , G[ WhIsG1[ k ] ] )
 Vec_VarValue Alfa1;  /**< linearization error of G[ WhIsG1[ k ] ] w.r.t. the
		       * current point Lambda. */

 Vec_VarValue LowerBound;  ///< Lower Bound over (each component of) Fi
 VarValue f_global_LB;     ///< an algorithmically discovered global LB
 
 VarValue t;           ///< the (tremendous) t parameter
 VarValue Prevt;       ///< what t were before being changed for funny reasons

 VarValue Sigma;       ///< Sigma*: convex combination of the Alfa's
 VarValue DSTS;        /**< D*_{t*}( -z* ), the other part of the dual
			* objective */
 Vec_VarValue vStar;   ///< v*, the predicted improvement

 VarValue DeltaFi;     ///< FiLambda - FiLambda1
 VarValue EpsU;        ///< precison required by the long-term t-strategy

 Index CSSCntr;        ///< counter of consecutive SS

 Index CNSCntr;        ///< counter of consecutive NS

 Subset vBPar2;        ///< size of the global pools of each component

 std::priority_queue< Index , std::vector< Index > ,
                      std::greater< Index > > FreList;
 ///< list of free positions in the bundle
 /**< FreList is the priority queue of free positions in the bundle, where
  * the position with higher priority is that with smaller name. This is why
  * the comparison parameter has to be explicitly set to
  * std::greater< Index >, since a priority queue with the default
  * std::less< Index > spits out the element with larger value. */

 /** NrItems[ k ] contains the number of items in the bundle (master problem)
  * for component k. If ( BPar7 & 3 ) < 3, this number may be strictly less
  * than the number of linearizations in the global pool of component k,
  * since removals from the bundle do not imply removals from the global pool
  */

 Subset NrItems;  ///< number of items in the bundle for each component

 /** FrFItem[ k ] contains the index of the first position in the global
  * pool of component k where BundleSolver can put a new linearization when
  * the corresponding item is added to the bundle (master problem). Note
  * that whether a position is suitable to this depends on BPar2: if
  * ( BPar7 & 3 == 0 ), then the position must be completely empty
  * (InvItemVcblr[ k ][ i ] == INF), while if ( BPar7 & 3 != 0 ) then
  * another linearization can be in that position already provided that
  * there is no corresponding item in the bundle
  * (vBPar2[ k ] <= InvItemVcblr[ k ][ i ] < INF). */

 Subset FrFItem;  ///< the first free item in each global pool

 /** MaxItem[ k ] contains 1 + the maximum index of a position in the
  * global pool of component k where a linearization is stored; if
  * MaxItem[ k ] == 0, then the global pool is empty. Note that this
  * ignores the fact that a position in the global pool corresponds or
  * not to an item in the bundle: all linearizations count. As a
  * consequence it must always be MaxItem[ k ] >= FrFItem[ k ]. */

 Subset MaxItem;  ///< the first unused item in each global pool

 /** Vocabulary of items: ItemVcblr[ i ].first is the component name
  * and ItemVcblr[ i ].second is the name in the global pool of that
  * component for item in position i of the bundle (master problem).
  * ItemVcblr[ i ].second == INF means that position i in the bundle is
  * not used. */

 std::vector< std::pair< Index , Index > > ItemVcblr;

 /** Inverse vocabulary of items. InvItemVcblr[ k ] is a Subset of size
  * vBPar2[ k ] and describes the global pool of component k. With
  * p = InvItemVcblr[ k ][ i ], if p < vBPar2[ NrFi ], then the linearization
  * with name i in the global pool of h is in the bundle at position p.
  * If p == INF, then there is no linearization with name i in the global
  * pool of k. If vBPar2[ NrFi ] <= p < INF, then there is a inearization with
  * name i in the global pool of h, but it is not in the bundle.
  *
  * NOTE: THE GLOBAL POOL OF SOME C05Function CAN BE LARGER THAN vBPar2[ k ],
  * BUT ALL ELEMENTS WITH NAME LARGER THAN vBPar2[ k ] ARE NEVER USED OR
  * CHANGED BY BundleSolver. */

 std::vector< Subset > InvItemVcblr;

  /** Out-Of-Base counters: if OOBase[ i ]
   * = Inf<SIndex>() then there is no item in position i of the bundle
   * = k > 0 means that the item in position i is out of base since k
   *   iterations
   * = 0 means in the current base but potentially removable
   * = a *finite* negative value - k means not removable for the next k
   *   iterations: note that some items in base may be such
   * = - Inf<SIndex>() means unremovable */

 Vec_SIndex OOBase;

 bool TrueLB;         /**< true if LowerBound is a "true" lower bound rather
		       * than a "conditional" one */
 bool SSDone;         ///< true if the last step was a SS

 Index f_wFi;         ///< which component was evaluated last

 Subset FiStatus;     ///< status of last computation of each component

 std::vector< C05Function * > v_c05f;
 ///< the vector of (pointers to) the components of the sum function

 LinearFunction * f_lf;  ///< the 0-th component of the sum function

 bool f_convex;          ///< true if all objectives are convex
 
 MPSolver * Master;      ///< (pointer to) the Master Problem Solver

 std::vector< MILPSolver * > MILP_s; /**< MILP solvera used to read the
				      * easy components */

 std::vector< ColVariable * > LamVcblr;  ///< map Lambda -> ColVariable

 VarValue UpTrgt;        ///< upper target
 VarValue LwTrgt;        ///< lower target

 VarValue UpFiBest;      ///< Fi best value vector

 Vec_VarValue UpRifFi;   /** The value of Fi[ k ]() where the zero of the
			  * translated Cutting Plane models are fixed */
 bool RifeqFi;           ///< true if UpRifFi == UpFiLmb

 Vec_VarValue UpFiLmb1;  ///< upper function values at Lambda1
 Vec_VarValue LwFiLmb1;  ///< lower function values at Lambda1
 Index UpFiLmb1def;      ///< how many entries of UpFiLmb1 are < INF
 Index LwFiLmb1def;      ///< how many entries of LwFiLmb1 are > -INF

 Vec_VarValue UpFiLmb;   ///< upper function value vector at Lambda
 Vec_VarValue LwFiLmb;   ///< lower function value vector at Lambda
 Index UpFiLmbdef;       ///< how many entries of UpFiLmb are < INF
 Index LwFiLmbdef;       ///< how many entries of LwFiLmb are > -INF

 VarValue Fi0Lmb;        ///< value of the linear 0-th component in Lambda
 VarValue Fi0Lmb1;       ///< value of the linear 0-th component in Lambda1
 
 Subset CurrNrEvls;      /**< how many times compute() has been called for
			  * each component in the current iteration */

 double DST;             ///< D_t( z* ), used to comupute the crucial Delta*
 double NrmD;            ///< Euclidean norm of the current direction d*
 double NrmZ;            ///< some norm of the aggregated sungradient z*
 double NrmZFctr;        ///< scaling factor to declare NrmZ "small"

 // fields for events - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::clock_t c_start;   ///< starting instant of last call to compute()

 // static fields - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
/*--------------------------- PRIVATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- CLASS FakeFiOracle  ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** FakeFiOracle implements the part of the FiOracle interface that is
 * strictly necessary to use a MPSolver inside BundleSolver. This hack will
 * one day be replaced with a native implementation of the master problem
 * solver, but until then, there you go. */

class FakeFiOracle : public FiOracle
{

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/** Constructor of the class: takes the pointer to the BundleSolver it has
 * to "serve". */

 FakeFiOracle( BundleSolver *solver ) : FiOracle() {
  bslv = solver;
  }

/*-------------------------- OTHER INITIALIZATIONS -------------------------*/

 void SetNDOSolver( NDOSolver *NwSlvr = 0 ) override {
  throw( std::logic_error( "this method cannot be called" ) );
  }

/*--------------------------------------------------------------------------*/

 void SetFiLog( ostream *outs = 0 , const char lvl = 0 ) override {
  throw( std::logic_error( "this method cannot be called" ) );
  }

/*--------------------------------------------------------------------------*/

 void SetFiTime( const bool TimeIt = true ) override {
  throw( std::logic_error( "this method cannot be called" ) );
  }

/*--------------------------------------------------------------------------*/

 void SetMaxName( cIndex MxNme = 0 ) override {
  throw( std::logic_error( "this method cannot be called" ) );
  }

/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/// get the number of Variable
/** Variable cannot be changed. This means that is used the default
 *  implementation of GetMaxNumVar(). The maximum number of variables is
 *  equal to the current number of variable*/

 Index GetNumVar( void ) const override;

/*--------------------------------------------------------------------------*/

 Index GetNrFi( void ) const override;

/*--------------------------------------------------------------------------*/

 Index GetMaxName( void ) const override;

/*--------------------------------------------------------------------------*/

 bool GetUC( cIndex i ) override;

/*--------------------------------------------------------------------------*/

 LMNum GetUB( cIndex i ) override;

/*--------------------------------------------------------------------------*/

 Index GetBNC( cIndex wFi ) override;

/*--------------------------------------------------------------------------*/

 Index GetBNR( cIndex wFi ) override;

/*--------------------------------------------------------------------------*/

 Index GetBNZ( cIndex wFi ) override;

/*--------------------------------------------------------------------------*/

 void GetBDesc( cIndex wFi , int *Bbeg , int *Bind , double *Bval ,
		double *lhs , double *rhs , double *cst ,
		double *lbd , double *ubd ) override;

/*--------------------------------------------------------------------------*/

 Index GetANZ( cIndex wFi , cIndex strt = 0 , Index stp = Inf<Index>() )
  override;

/*--------------------------------------------------------------------------*/

 void GetADesc( cIndex wFi , int *Abeg , int *Aind , double *Aval ,
		cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

 HpNum Fi( cIndex wFi = Inf<Index>() ) override {
  throw( std::logic_error( "this method cannot be called" ) );
  }

/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/

 bool NewGi( cIndex wFi = Inf<Index>() ) override { return( true ); }

/*--------------------------------------------------------------------------*/

 Index GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name = Inf<Index>() ,
	      cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*------------------------------ DESTRUCTOR --------------------------------*/

 virtual ~FakeFiOracle() { }

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

 BundleSolver * bslv;  ///< the BundleSolver that I "serve"

 };  // end( class FakeFiOracle )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void InitMP( void );

 bool FindNext( void );

/*--------------------------------------------------------------------------*/

 Index BStrategy( cIndex wFi );

/*--------------------------------------------------------------------------*/

 Index FindAPlace( Index wFi );

/*--------------------------------------------------------------------------*/

 HpNum Heuristic1( void );

 HpNum Heuristic2( void );

/*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );

/*--------------------------------------------------------------------------*/

 void ReSetAlg( unsigned char RstLvl = 0 );

 /* Resets the internal state of the Bundle algorithm. Since several
    different things can be reset independently, RstLvl is coded bit-wise:

    - bit 0: if 0, all the algorithmic parameters are reset to the default
      values read by the stream/set by SetPar(), while if 1 they are left
      untouched;

    - bit 1: if 0 the current point is reset to the all-0 vector, while if
      1 it is left untouched;

    - bit 2: if 0, the current point is reset to the value currently in the
      active Variable of the C05Function, while if 1 it is left untouched. */

/*--------------------------------------------------------------------------*/

 void Delete( cIndex i , bool ModDelete = false );

/*--------------------------------------------------------------------------*/

 void UpdtaBP3( void );

/*--------------------------------------------------------------------------*/

 bool IsOptimal( double eps = 0 ) const;

/*--------------------------------------------------------------------------*/

 Index get_index_of_component( Function * f )
 {
  const auto fit = std::find( v_c05f.begin() , v_c05f.end() , f );
  if( fit != v_c05f.end() )
   return( std::distance( v_c05f.begin() , fit ) );

  return( Inf< Index >() );
  }

/*--------------------------------------------------------------------------*/

 void remove_from_global_pool( Index k , Index i , bool hard );

 Index find_place_in_global_pool( Index k );

 void add_to_global_pool( Index k , Index i , Index wh = Inf<Index>() );

 void add_to_bundle( Index k , Index i );

 void reset_bundle( void );

/*--------------------------------------------------------------------------*/

 bool is_special_GroupMod( GroupModification & gmod );

 void flatten_Modification_list( Lst_sp_Mod & vmt , sp_Mod mod );

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/

 void FModChg( VarValue shift , Index wFi );

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

 void CheckBundle( void );

 void CheckAlpha( void );

 void CheckLBs( void );

 void PrintBundle( void );

#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 Index aBP3;       // current max number of items to be fetched

 FakeFiOracle FakeFi;  ///< the FakeFiOracle object

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class BundleSolver )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BundleSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File BundleSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
