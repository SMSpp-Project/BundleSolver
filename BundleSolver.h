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
 *     ALL THE Modification THAT CHANGE THE "ACTIVE" Variable MUST BE
 *     BUNCHED TOGETHER IN A SINGLE GroupModification. THIS MUST CONTAIN
 *     EXACTLY AS MANY Modification AS THERE ARE sub-Block (AND, THEREFORE,
 *     DIFFERENT OBJECTIVE), PLUS ONE IF THE (LinearFunction IN THE)
 *     Objective OF THE Block IS NOT EMPTY. ALL Modification MUST BE OF
 *     THE VERY SAME TYPE, I.E., EITHER ALL C05FunctionModVarsAddd, OR ALL
 *     C05FunctionModVarsRngd, OR ALL C05FunctionModVarsSbst, AND THEY MUST
 *     CHANGE THE "ACTIVE" Variable IN PRECISELY THE SAME WAY.
 *
 * \version 0.30
 *
 * \date 07 - 07 - 2020
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

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 using namespace NDO_di_unipi_it;

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
 * is no need for a Solver to be attached to the inner Block.
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
  NumVar( 0 ) , NrFi( 0 ) , SCalls( 0 ) , ParIter( 0 ) , NrEasy( 0 ) ,
  LHasChgd( true ) , tHasChgd( true ) , t( 0 ) , Prevt( 0 ) , Sigma( 0 ) ,
  DSTS( 0 ) , vStar( 0 ) , DeltaFi( 0 ) , EpsU( 0 ) , CSSCntr( 0 ) ,
  CNSCntr( 0 ) , TrueLB( false ) , SSDone( true ) , aBP3( 0 ) ,
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
  BPar7 = dflt_int_par[ intBPar7 - intLastParCDAS ];
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

  MaxTime = CDASolver::get_dflt_dbl_par( dblMaxTime );
  RelAcc = CDASolver::get_dflt_dbl_par( dblRelAcc );
  AbsAcc = CDASolver::get_dflt_dbl_par( dblAbsAcc );
  RAccSol = CDASolver::get_dflt_dbl_par( dblRAccSol );
  AAccSol = CDASolver::get_dflt_dbl_par( dblAAccSol );
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
  * - intMaxNrEvls [2]: max number of function evaluation for each iteration
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
  if( ( MaxSol > 1 ) && ( UpFiBest < UpFiLmb[ NrFi ] ) )
   return( UpFiBest );
  else
   return( UpFiLmb[ NrFi ] );
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
  if( ( MaxSol > 1 ) && ( UpFiBest < UpRifFi[ NrFi ] ) ) {
   for( Index i = 0 ; i < NumVar ; i++ )
    LamVcblr[ i ]->set_value( LmbdBst[ i ] );
   }
  else {
   for( Index i = 0 ; i < NumVar ; i++ )
    LamVcblr[ i ]->set_value( Lambda[ i ] );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the "current" dual solution

 void get_dual_solution( Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 bool new_var_solution( void ) override
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

 bool new_dual_solution( void )  override { return( false ); }

/*--------------------------------------------------------------------------*/
/*
  void set_unbounded_threshold( const VarValue thr ) override { }

  bool has_var_direction( void ) override { return( true ); }

  bool has_dual_direction( void ) override { return( true ); }

  void get_var_direction( Configuration *dirc = nullptr ) override {}

  void get_dual_direction( Configuration *dirc = nullptr ) override {}

  virtual bool new_var_direction( void ) override { return( false ); }
  
  virtual bool new_dual_direction( void ) override{ return( false ); }
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
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
 /* Computes Fi( Lambda1 ), inserting the obtained items (subgradients or
  * constraints) in the bundle. Returns true <=> the newly obtained
  * information changes the solution of the MP. */

 bool FiAndGi( Index wFi );

/*--------------------------------------------------------------------------*/
 /* Move the current point to Lambda1. */

 void GotoLambda1( void );

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
 int BPar7;         ///< if BundleSolver "plays nice" with other Solver

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
		    * subgradient" Z[ k ] of component k is kept in
		    * whisZ[ k ]; Inf<Index>() == it is not in the bundle */
 std::vector< bool > Zvalid;  /**< Zvalid[ k ] == true if the item in position
			       * whisZ[ k ] is exactly Z[ k ] as computed by
		    * the last master problee. Zvalid[ k ] == true ==>
		    * whisZ[ k ] < INF. if Zvalid[ k ] == false and
		    * whisZ[ k ] < INF, then Z[ k ] had been previously
		    * stored in position whisZ[ k ], but the master problem has
		    * been re-solved since and therefore Z[ k ] is no longer
		    * current. */
 
 Subset whisG1;    ///< "representative subgradient" for each component
 Vec_VarValue ScPr1;  ///< ScalarProduct( dir , G[ WhIsG1[ k ] ] )
 Vec_VarValue Alfa1;  /**< linearization error of G[ WhIsG1[ k ] ] w.r.t. the
		       * current point Lambda. */

 Vec_VarValue LowerBound;  ///< Lower Bound over (all the components of) Fi

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

 Subset vBPar2;  ///< dimension of the global pools of each component

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
  * and ItemVcblr[ i ].second field is the name in the global pool of
  * that component for item in position i of the bundle (master problem).
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
 Subset CurrNrEvls;

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

   void SetNDOSolver( NDOSolver *NwSlvr = 0 ) override;

/*--------------------------------------------------------------------------*/

   void SetFiLog( ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

   void SetFiTime( const bool TimeIt = true ) override;

/*--------------------------------------------------------------------------*/

   void SetMaxName( cIndex MxNme = 0 ) override;

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

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

   HpNum GetMinusInfinity( void ) override;

/*--------------------------------------------------------------------------*/

   Index GetMaxNZ( cIndex wFi = Inf<Index>() ) const override;

/*--------------------------------------------------------------------------*/

   Index GetMaxCNZ( cIndex wFi = Inf<Index>() ) const override;

/*--------------------------------------------------------------------------*/

   bool GetUC( cIndex i ) override;

/*--------------------------------------------------------------------------*/

   LMNum GetUB( cIndex i ) override;

/*--------------------------------------------------------------------------*/

   LMNum GetBndEps( void ) override;

/*--------------------------------------------------------------------------*/

   HpNum GetGlobalLipschitz( cIndex wFi = Inf<Index>() ) override;

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

   NDOSolver *GetNDOSolver( void ) override;

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SETTING LAMBDA ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Setting Lambda
   @{ */

   void SetLambda( cLMRow Lmbd = 0 ) override;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

   void SetLamBase( cIndex_Set LmbdB = 0 , cIndex LmbdBD = 0 ) override;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

   bool SetPrecision( HpNum Eps ) override;

/*@} -----------------------------------------------------------------------*/
/*------------------------ METHODS FOR COMPUTING Fi() ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Computing Fi()
   @{ */

   HpNum Fi( cIndex wFi = Inf<Index>() ) override;

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/

   bool NewGi( cIndex wFi = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   Index GetGi( SgRow SubG , cIndex_Set &SGBse ,
		cIndex Name = Inf<Index>() ,
		cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   HpNum GetVal( cIndex Name = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   void SetGiName( cIndex Name ) override;

/*@} -----------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading other results
   @{ */

   HpNum GetLowerBound( cIndex wFi = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   FiStatus GetFiStatus( Index wFi = Inf<Index>() ) override;

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
   @{ */

   void Deleted( cIndex i = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   void Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm , cIndex NwNm )
    override;

/*@} -----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~FakeFiOracle() { }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

  BundleSolver *bslv;

  /* vocabulary
     for handling the items name; this is done to map the item name
     from the FiOracle to that of C05Function.  */

  Index last_c05;

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class FakeFiOracle )

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS  -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void InitMP( void );

 bool FindNext( Index &wFi );

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
