/*--------------------------------------------------------------------------*/
/*----------------------- File BundleSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the BundleSolver class, which implements the
 * Solver interface, in particular in its CDASolver version, using a
 * "Generalized Bundle" algorithm for the solution of convex nondifferentiable
 * problems.
 *
 * The user is assumed to be familiar with the algorithm: refer to
 *
 *  A. Frangioni "Generalized Bundle Methods"
 *  SIAM Journal on Optimization 13(1), 117--156, 2002
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
 *  61 - 116, Springer, 2020
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#NDOB18
 * \endlink
 *
 * In particular, BundleSolver implements the Incremental version of
 * the (Generalized) Proximal Bundle approach using upper models (for all the
 * components that provide a Lipschitz constant) described in
 *
 *  W. van Ackooij, A. Frangioni "Incremental Bundle Methods Using Upper
 *  Models" SIAM Journal on Optimization 28(1), 379 – 410, 2018
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#SIOPT16
 * \endlink
 *
 * BundleSolver is capable of solving any Block such that:
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
 * actually are LagBFunction whose inner Block only contains ColVariable, whose
 * Objective is linear (a FRealObjective containing a LinearFunction) and whose
 * Constraint are linear (either FRowConstraint containing a LinearFunction, or
 * BoxConstraint). These can be passed to the Master Problem of the bundle
 * algorithm as "easy components", see
 *
 *   A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
 *   with ``Easy'' Components: Applications to Multicommodity Network Design"
 *   Mathematical Programming 145(1), 133 – 161, 2014
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#MP11c
 * \endlink
 *
 * In that case, the LagBFunction is never evaluated, which means that there is
 * no need for a Solver to be attached to the inner Block.
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
 * To ensure that the rule about the list of "ACTIVE" Variable in the (multiple)
 * Objective is respected, an analogous very strong assumption is made on the
 * Modification that change that:
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BundleSolver
 #define __BundleSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <chrono>
#include <queue>
#include <unordered_map>

#include "CDASolver.h"

#include "C05Function.h"
#include "LinearFunction.h"

#include "Block.h"
#include "ColVariable.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"

#include "FakeSolver.h"
#include "MILPSolver.h"

#include "MasterProblemBlock.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class BundleSolverState;     // forward declaration of the state

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup LagBFunction_CLASSES Classes in BundleSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS BundleSolver -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A CDASolver using a (Generalized) Bundle algorithm
/** BundleSolver implements the Solver interface, in particular
 * in its CDASolver version, using a "Generalized Bundle" algorithm.
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
 *  61 - 116, Springer, 2020
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#NDOB18
 * \endlink
 *
 * In particular, BundleSolver implements the Incremental version of
 * the (Generalized) Proximal Bundle approach using upper models (for all the
 * components that provide a Lipschitz constant) described in
 *
 *  W. van Ackooij, A. Frangioni "Incremental Bundle Methods Using Upper
 *  Models" SIAM Journal on Optimization 28(1), 379 – 410, 2018
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#SIOPT16
 * \endlink
 *
 * BundleSolver is capable of solving any Block such that:
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
 * actually are LagBFunction whose inner Block only contains ColVariable, whose
 * Objective is linear (a FRealObjective containing a LinearFunction) and whose
 * Constraint are linear (either FRowConstraint containing a LinearFunction, or
 * OneVarConstraint). These can be passed to the Master Problem of the bundle
 * algorithm as "easy components", see
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
 * In that case, the LagBFunction is never evaluated, which means that there is
 * no need for a Solver to be attached to the inner Block.
 *
 * If the Block has multiple Objective (that is, it has sub-Block whose
 * Objective FRealObjective containing a C05Function), a very strong assumption
 * is required on them:
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
 * To ensure that the rule about the list of "ACTIVE" Variable in the (multiple)
 * Objective is respected, an analogous very strong assumption is made on the
 * Modification that change that:
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
 * "Import" basic types from Function, C05Function.
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

 static constexpr auto NaNshift
                              = std::numeric_limits< VarValue >::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

 static constexpr auto INFshift = Inf< VarValue >();
 ///< convenience constexpr for "Infty"

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that BundleSolver has in addition to these of CDASolver. The
  * value intLastBndSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

 enum int_par_type_BndSlv {

 intBPar1 = intLastParCDAS ,
 ///< remove linearizations unused for more than this consecutive iterations

 intBPar2 ,  ///< max number linearizations per component

 intBPar3 ,  ///< max number of new linearizations per iteration per component

 intBPar4 ,  ///< min number of new linearizations per iteration per component

 intBPar6 ,  ///< control how the min/max number of new linearizations changes

 intBPar7 ,  ///< how well-behaved GBS is w.r.t. other Solver

 intMnSSC ,  ///< minimum number of consecutive Serious Steps

 intMnNSC ,  ///< minimum number of consecutive Null Steps

 inttSPar1 ,  ///< first t-strategy parameter

 intMaxNrEvls ,  ///< max number of function evaluation for each iteration

 intDoEasy ,  ///< whether "easy" components are considered

 intWZNorm ,  ///< how to compute the norm of z*

 intFrcLstSS ,  ///< whether to force the last step to be a SS

 intTrgtMng ,   ///< how to manage targets and accuracy in the functions

 intMPStbl ,  ///< type of stabilization for the master problem

 intMPPrimal ,  ///< whether the MP is in its primal or dual form

 intRstAlg ,  ///< reset parameter

 intMPV2Form ,  ///< dual MP storage frame: displacement or iterate form

 intMPHScaling ,  ///< bit-wise scaling of hard-component PFBs

 intMaxLevelNR ,  ///< maximum number of NR steps allowed

 intLastBndSlvPar  ///< first allowed new int parameter for derived classes
                   /**< Convenience value for easily allow derived classes
                    * to extend the set of int algorithmic parameters. */

 };  // end( int_par_type_BndSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that BundleSolver has in addition to these of CDASolver. The
  * value dblLastBndSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

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

  dbltSPar3 ,  ///< numerical parameter for "small" heuristic-based t changes

  dblLStabM ,  ///< m_l factor for level stabilization updates

  dblLStabDlt ,   ///< exogenous Delta fraction when no reliable LB is known

  dblLStabIncr ,  ///< exogenous Delta increase factor after consecutive SS

  dblLStabSmall , ///< small model-error threshold for level Delta increases

  dblLastBndSlvPar ///< first allowed new double parameter for derived classes
                   /**< Convenience value for easily allow derived classes
                    * to extend the set of double algorithmic parameters. */

  };  // end( dbl_par_type_BndSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the string algorithmic parameters
 /** Public enum describing the different algorithmic parameters of string
  * type that BundleSolver has in addition to these of CDASolver. The
  * value strLastBndSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

 enum str_par_type_BndSlv {
  strEasyCfg = strLastParCDAS ,  ///< string name for "easy" Configurations

  strHardCfg ,                   ///< string name for not-easy Configurations

  strMPBSolverCfg ,                 ///< string name for Configuration of the
                                 ///  solver associated with the MPBlock

  strLastBndSlvPar ///< first allowed new string parameter for derived classes
                   /**< Convenience value for easily allow derived classes
                    * to extend the set of string algorithmic parameters. */
  };  // end( str_par_type_BndSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the vector-of-int parameters
 /** Public enum describing the different algorithmic parameters of
  * vector-of-int type that BundleSolver has in addition to these of
  * CDASolver. The value vintLastBndSlvPar is provided so that the list can be
  * easily further extended by derived classes. */

 enum vint_par_type_BndSlv {
  vintNoEasy = vintLastParCDAS ,
  ///< parameter for excluding certain components from being "easy"
  /**< The vector vintNoEasy is assumed to contain the indices (numbers in
   * 0, ..., total number of components - 1, ordered in increasing sense and
   * therefore not repeated) of the components of the problem that must not be
   * treated as "easy" even if they could. This clearly only applies if
   * intDoEasy == 1, for otherwise no component is ever treated as "easy".
   */

  vintLastBndSlvPar ///< first allowed new vector-of-int parameter
                    /**< Convenience value for easily allow derived classes
                     * to extend the set of vector-of-int parameters. */

  };  // end( vint_par_type_BndSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the vector-of-string parameters
 /** Public enum describing the different parameters of vector-of-string type
  * that BundleSolver has in addition to these of CDASolver. The
  * value vstrLastBndSlvPar is provided so that the list can be easily further
  * extended by derived classes. */

 enum vstr_par_type_BndSlv {
  vstrCmpCfg = vstrLastParCDAS ,
  ///< parameter for configuring (possibly) each component individually

  vstr_C05_SPAR_Names ,
  ///< string parameters names that are set differently to each C05Function

  vstr_C05_SPAR_Vals ,
  ///< baseline values for different string parameters to each C05Function

  vstr_C05_EI_SPAR_Names ,
  /**< string parameters names that are changed at every iteration and call
   * differently to each C05Function */

  vstr_C05_EI_SPAR_Vals ,
  /**< baseline values for string parameters changed at every iteration and
   * call differently to each C05Function */

  vstrLastBndSlvPar ///< first allowed new vector-of-string parameter
                    /**< Convenience value for easily allow derived classes
                     * to extend the set of vector-of-string parameters. */

  };  // end( vstr_par_type_BndSlv )

/** @} ---------------------------------------------------------------------*/
/*--------------- CONSTRUCTING AND DESTRUCTING BundleSolver ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing BundleSolver
 *  @{ */

 /// constructor: ensure every field is initialized

 BundleSolver( void ) : CDASolver() , Result( kUnEval ) ,
  NumVar( 0 ) , NrFi( 0 ) , SCalls( 0 ) , ParIter( 0 ) , NrEasy( 0 ) ,
  LHasChgd( true ) , tHasChgd( true ) , MPchgs( 0 ) , G1Norm( 0 ) ,
  ScPr1( 0 ) , Alfa1( 0 ) , f_global_LB( -INFshift ) , t( 0 ) , Prevt( 0 ) ,
  Sigma( 0 ) , DSTS( 0 ) , DeltaFi( 0 ) , EpsU( 0 ) , CSSCntr( 0 ) ,
  CNSCntr( 0 ) , TrueLB( false ) , SSDone( true ) , f_wFi( 0 ) ,
  f_lf( nullptr ) , f_convex( true ) , MasterPB( nullptr ) ,
  UpTrgt( 0 ) , LwTrgt( 0 ) , RifeqFi( false ) ,
  CmptdinL( false ) , UpFiBest( INFshift ) , UpFiLmb1def( 0 ) ,
  LwFiLmb1def( 0 ) , UpFiLmbdef( 0 ) , LwFiLmbdef( 0 ) , Fi0Lmb( 0 ) ,
  Fi0Lmb1( 0 ) , DST( 0 ) , NrmD( 0 ) , NrmZ( 0 ) , NrmZFctr( 1 ) ,
  c_start() , aBP3( 0 ) , LevelNRCntr( 0 )
 {
  // ensure all parameters are properly given their default value
  MaxIter = CDASolver::get_dflt_int_par( intMaxIter );
  MaxSol = CDASolver::get_dflt_int_par( intMaxSol );
  EverykIt = CDASolver::get_dflt_int_par( intEverykIt );
  LogVerb = CDASolver::get_dflt_int_par( intLogVerb );
  BPar1 = Index( get_dflt_int_par( intBPar1 ) );
  BPar2 = Index( get_dflt_int_par( intBPar2 ) );
  BPar3 = Index( get_dflt_int_par( intBPar3 ) );
  BPar4 = Index( get_dflt_int_par( intBPar4 ) );
  BPar6 = get_dflt_int_par( intBPar6 );
  BPar7 = get_dflt_int_par( intBPar7 );
  MnSSC = get_dflt_int_par( intMnSSC );
  MnNSC = get_dflt_int_par( intMnNSC );
  tSPar1 = get_dflt_int_par( inttSPar1 );
  MaxNrEvls = get_dflt_int_par( intMaxNrEvls );
  DoEasy = get_dflt_int_par( intDoEasy );
  WZNorm = char( get_dflt_int_par( intWZNorm ) );
  FrcLstSS = get_dflt_int_par( intFrcLstSS );
  TrgtMng = Index( get_dflt_int_par( intTrgtMng ) );
  MPStbl = static_cast< MasterProblemBlock::stabilization_type >(
                                              get_dflt_int_par( intMPStbl ) );
  IsMPPrimal = bool( get_dflt_int_par( intMPPrimal ) );
  MPV2Form = get_dflt_int_par( intMPV2Form );
  MPHScaling = get_dflt_int_par( intMPHScaling );

  MaxTime = CDASolver::get_dflt_dbl_par( dblMaxTime );
  RelAcc = CDASolver::get_dflt_dbl_par( dblRelAcc );
  AbsAcc = CDASolver::get_dflt_dbl_par( dblAbsAcc );
  EveryTTm = CDASolver::get_dflt_dbl_par( dblEveryTTm );
  NZEps = get_dflt_dbl_par( dblNZEps );
  tStar = get_dflt_dbl_par( dbltStar );
  MinNrEvls = get_dflt_dbl_par( dblMinNrEvls );
  BPar5 = get_dflt_dbl_par( dblBPar5 );
  m1 = get_dflt_dbl_par( dblm1 );
  m2 = get_dflt_dbl_par( dblm2 );
  m3 = get_dflt_dbl_par( dblm3 );
  mxIncr = get_dflt_dbl_par( dblmxIncr );
  mnIncr = get_dflt_dbl_par( dblmnIncr );
  mxDecr = get_dflt_dbl_par( dblmxDecr );
  mnDecr = get_dflt_dbl_par( dblmnDecr );
  tMaior = get_dflt_dbl_par( dbltMaior );
  tMinor = get_dflt_dbl_par( dbltMinor );
  tInit = get_dflt_dbl_par( dbltInit );
  tSPar2 = get_dflt_dbl_par( dbltSPar2 );
  tSPar3 = get_dflt_dbl_par( dbltSPar3 );
  LStabM = get_dflt_dbl_par( dblLStabM );
  LStabDlt = get_dflt_dbl_par( dblLStabDlt );
  LStabIncr = get_dflt_dbl_par( dblLStabIncr );
  LStabSmall = get_dflt_dbl_par( dblLStabSmall );
  MaxLevelNR = get_dflt_int_par( intMaxLevelNR );

  v_events.resize( max_event_number() );
  }

/*--------------------------------------------------------------------------*/
 /// destructor: cleanly detaches the BundleSolver from the Block

 virtual ~BundleSolver() { set_Block( nullptr ); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve
 /** Gives the BundleSolver access to the Block it has to solve;
  * note that this does not register the BundleSolver among the
  * Solver of the Block, because the converse happens. Extensive checks are
  * performed during set_Block() to ensure that the Block does satisfy the
  * requirements of BundleSolver, and all the nontrivial internal
  * data structures of BundleSolver are set up.
  *
  * If \p block == nullptr, the BundleSolver is completely cleaned up
  * and prepared for either destruction or receiving an entirely unrelated Block
  * to solve. */

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// set the int parameters of BundleSolver
 /** Set the int parameters specific of BundleSolver, together
  * with the parameters of CDASolver that BundleSolver actually
  * "listens to":
  *
  * - intMaxIter [Inf< int >]: maximum iterations for the next call to solve()
  *
  * - intMaxSol [1]: maximum number of different solutions to report. Since
  *                  Bundle methods are "almost" monotone ones, they typically
  *   produce only one solution, which is the stability center at termination.
  *   However, this is "almost" true: in fact, that solution may not be the
  *   best one ever found. If intMaxSol > 1, also the best solution will be
  *   kept and separately reported (assuming it is not the stability center
  *   at termination).
  *
  * - intEverykIt [0]: after how many iterations call the
  *                    eEverykIteration events
  *
  * - intLogVerb [0]: "verbosity" of the BundleSolver log
  *                   0 = no log
  *                   1 = only final state of the call and errors
  *                   2 = detailed step-by-step log
  *                   3 = as 2 + print every linearization added/removed
  *                   4 = as 3 + print every function value computed
  *                   5 = as 4 + print the tentative point at every iteration
  *                   6 = as 5 + print all the UB/LB at the tentative point
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
  * BundleSolver to perform aggregation at every iteration (for every
  *   C05Function). Of course, keeping the "bundle" small makes the Master
  *   Problem cheaper, but on the other hand acquiring enough first-order
  *   information is typically the name of the game, hence keeping this
  *   value too low can have a dramatic effect on convergence speed that
  *   can easily counterbalance any improvement in Master Problem cost.
  *
  * - intBPar3 [1]: maximum number of new linearizations to be fetched from
  *                 each (non-easy) C05Function at each function evaluation
  *
  * - intBPar4 [1]: minimum number of new linearizations to be fetched from
  *                 each (non-easy) C05Function at each function evaluation
  *
  * - intBPar6 [0]: together with the double parameters dblBPar3, dblBPar4
  *                 and dblBPar5, controls how the actual number of
  *   linearization that are requested to the C05Function evolves as the
  *   algorithm proceeds; note that what varies in practice is the maximum
  *   number, as it is always legal for the C05Function to refuse giving
  *   other items, although the BundleSolver will complain and
  *   stop if less than BPar4 are given. In BundleSolver,
  *   the number
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
  * - intBPar7 [2]: This parameter, coded bit-wise, controls if
  *                 BundleSolver "tries to play nice" with any
  *   other Solver that may concurrently be using the same C05Function. The
  *   point is that each of these Solver is producing new linearizations,
  *   and possibly storing them in, or removing them from, the "finite
  *   resource" of the global pool(s) of the C05Function(s). Hence, what
  *   the BundleSolver does to the global pool may have an
  *   impact on the other Solver, if any. This parameter controls whether
  *   BundleSolver tries as hard as possible to avoid impacting
  *   the other Solver operations, or if it rather assumes to be "the only
  *   one" working with the C05Function, and therefore "treats the global
  *   pool as its exclusive property". To do so, BundleSolver
  *   handles the slot of the global pool in different ways according to
  *   the value in the first two bits of intBPar7 ( intBPar7 & 3 ):
  *
  *   = 0 means that BundleSolver will never override any
  *       position in the global pool unless it strictly needs to. This
  *       means that even if a linearization is removed from the bundle
  *       (the master problem), it is kept in the global pool of the
  *       corresponding component until the latter is completely full.
  *       Only then linearizations are removed, when necessary to make
  *       space for newly generated ones. Note that BundleSolver
  *       always "proceeds from left to right", i.e., selects the
  *       linearization in the global pool with smallest "name". This
  *       creates a sort of FIFO order whereby the oldest linearizations
  *       are removed first, which makes general sense.
  *
  *   = 1 means that BundleSolver will not immediately delete
  *       from the global pool a linearization that it removes from the
  *       bundle (the master problem). While the linearization is kept
  *       there, BundleSolver considers it "free", and can
  *       immediately after re-use that position to store a newly computed
  *       linearization. Again, the order is that if smaller names first,
  *       so if a linearization with "large name" is removed from the
  *       global pool it may take some time before it is actually
  *       overwritten by BundleSolver, thereby leaving it
  *       available to other Solver.
  *
  *   = 2 means that BundleSolver will immediately delete from
  *       the global pool any linearization that it removes from the
  *       bundle (the master problem). This makes sense if
  *       BundleSolver is the only Solver producing and
  *       consuming linearizations in these C05Function(s), since it
  *       allows them to immediately delete all the memory (which may be
  *       significant) associated with that linearization in the global
  *       pool. However, if a linearization is found to be a "better copy"
  *       of a known one (the new linearization has the same linear part
  *       but a larger constant, and therefore provides a tighter
  *       constraint on the epigraph of the convex function, so that no
  *       Solver should complain if the weaker constraint is removed
  *       provided that the better one is added), still the old
  *       linearization is kept in the global pool (but not in the bundle)
  *       unless it is strictly necessary to do so.
  *
  *   = 3 means that BundleSolver will immediately delete from
  *       the global pool any linearization that it removes from the
  *       bundle (the master problem); furthermore, if it finds a "better
  *       copy" of an existing linearization the new one immediately
  *       replaces the old one, in the global pool as well as in the
  *       bundle.
  *
  *   The bit 2 ( intBPar7 & 4 ) rather decides how BundleSolver
  *   reacts to Modification telling that some other Solver have generated
  *   a new linearization. If the bit is 0, then BundleSolver
  *   plainly ignores it, which is likely the best strategy if producing
  *   linearizations is "cheap". However, if ( intBPar7 & 3 ) < 3
  *   BundleSolver does take note that a linearization is there
  *   in order to avoid to touch it "unless strictly necessary". If the
  *   bit is 1 instead, then BundleSolver will right away add
  *   the linearization to its bundle (the master problem), which is
  *   likely the best strategy if producing linearizations is "costly" and
  *   therefore it makes sense to profit from the effort that the
  *   C05Function(s) has done on behalf of the other Solver(s).
  *
  *   The bit 3 ( intBPar7 & 8 ) has a similar role for the initialization
  *   phase: if it is == 1, then BundleSolver will also scan
  *   the global pool of each component when it is attached to the Block,
  *   and immediately add to the bundle every linearization it finds
  *   there.
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
  *   The first two bits control if heuristics are used compute a new value
  *   of t when increasing/decreasing it:
  *
  *    bit 0:    1 (+1) if the heuristic is used during a SS, 0 otherwise
  *
  *    bit 1:    1 (+2) if the heuristic is used during a NS, 0 otherwise
  *
  *   The following 2 bits (bit 2 and 3) of inttSPar1 tell which long-term
  *   t-strategy is used, with the following values:
  *
  *    00 (+ 0): none, only the heuristics (if any) are used
  *
  *    01 (+ 4): the "soft" long-term t-strategy is used: the value EpsU is
  *              maintained (cf. intBPar6) which estimates the current
  *              relative error, and decreases of t are inhibited whenever
  *              v < tSPar2 * EpsU * | FiVal |
  *
  *    10 (+ 8): the "hard" long-term t-strategy is used: the value EpsU is
  *              maintained as above, and t is increased whenever
  *              v < tSPar2 * EpsU * | FiVal |
  *
  *    11 (+12): the "balancing" long-term t-strategy is used, where the two
  *              terms D*_t( -z* ) and Sigma* are kept of "roughly the same
  *              size": if D*_1( -z* ) <= tSPar2 * Sigma* then t increases
  *              are inhibited (increasing t causes a decrease of D*_1( -z* )
  *              that is already small), if tSPar2 * D*_1( -z* ) >= Sigma*
  *              then t decreases are inhibited (decreasing t causes an
  *              increase of D*_1( -z* ) that is already big)
  *
  *   These three long-term t-strategies are mutually exclusive. The following
  *   2 bits (bit 4 and 5) of inttSPar1 can instead be used to specify
  *   t-strategies that can be activated in addition to these, with the
  *   following values:
  *
  *    bit 4: 1 (+16) if the "endgame" t-strategy is used, where if
  *           D*_1( -z* ) is "small" (~ 1/10 of the current absolute epsilon)
  *           t is decreased no matter what the other strategies dictated.
  *           The rationale is that we are "towards the end" of the
  *           optimization and here t needs decrease. However, note that
  *           having D*_1( -z* ) "small" is no guarantee that we actually
  *           are at the end, especially if the oracle dynamically
  *           generates its variables, so use with caution
  *
  *    bit 5: currently reserved
  *
  *   The following 2 bits (bits 6 and 7) dictate which of the available
  *   heuristics are used to select t during a SS, if this is activated (see
  *   the bit 0). There are four heuristics available. The first three are
  *   based on a quadratic interpolation of the function based on known data,
  *   differing about which three of the four data points available (function
  *   values and derivatives in 0 and t) are used, while the fourth one is
  *   the so-called reversal form of the poorman's quasi-Newton update.
  *   The setting of the bits is:
  *
  *     00 (+  0):  Heuristic1 is used
  *     01 (+ 64):  Heuristic2 is used
  *     10 (+128):  Heuristic3 is used
  *     11 (+192):  Heuristic4 is used
  *
  *   The following two bits (8 and 9) are analogous to the previous two ones
  *   in that they dictate which of the available heuristics are used to
  *   select t during a NS:
  *
  *     00 (+  0):  Heuristic1 is used
  *     01 (+256):  Heuristic2 is used
  *     10 (+512):  Heuristic3 is used
  *     11 (+768):  Heuristic4 is used
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
  * - intDoEasy [1]: this boolean parameter controls whether
  *                  BundleSolver uses the "easy components"
  *   approach on components that allow it (LagBFunction with linear
  *   constraints, objective and continuous variables only).
  *
  *   If it is 0, then all components are treated as "hard" even if
  *   they could be treated as "easy". If it is 1, then all "easy"
  *   components are treated as such.
  *
  *   NOTE: In a previous version intDoEasy controlled also which structures
  *   were allowed to change in an esy component. With the use of MPBlock
  *   this is not necessary anymore, as the latter allows for changes in
  *   any part of the structure of easy components.
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
  *   capability of BundleSolver to produce global valid lower
  *   bounds (for a minimization problem, upper bounds for a maximization
  *   one). Indeed, these can only be produced when z* is "0"; this is
  *   taken to mean "almost 0" in the specific sense dictated by this
  *   parameter together with dblZNEps.
  *
  * - intFrcLstSS [0]: bit-wise encoding that controls two "symmetric"
  *                    issues about "trusting" input/output state.
  *     bit 0: if 1, it ensures that all the non-easy components have been
  *            evaluated the last time on the point that is returned
  *            (first) by get_var_solution(). Some approaches using
  *            BundleSolver may require this because they use
  *            some other information provided by the compute()-tion
  *            process of the components that need be "current" with the
  *            optimal solution. This may happen automatically if the very
  *            last iteration that the algorithm performs before stopping
  *            is a "serious step", but in general this is not guaranteed,
  *            whence the need for this parameter. Note that setting it
  *            to 1 may be as expensive as computing all components is;
  *            in particular it cannot work if the maximum time limit has
  *            been exceeded already, and it may trigger a kStopTime
  *            return status where a kOK would have been produced exactly
  *            due to the cost of the extra compute()-tions.
  *     bit 1: if 1, it makes it so that the function values stored in a
  *            BundleSolverState [see] are not trusted when it
  *            is put() back in BundleSolver. This causes the
  *            re-computation of all function values at the first
  *            iteration before the algorithm can declare optimality.
  *            This is provided in case the state of the Block to which
  *            BundleSolver is registered is not the same as
  *            that when the BundleSolverState was get(),
  *            which may happen either if they are actually two different
  *            (similar, but not identical) Block, or if
  *            BundleSolver was detached, some changes were
  *            effected in the Block and then BundleSolver was
  *            re-attached. This way of using BundleSolverState
  *            is explicitly permitted by the definition of State, with
  *            the provision that the Solver must be able to identify
  *            somehow the inconsistencies that it can create. Since
  *            BundleSolver has no way to check what had
  *            happened to the Block "when it was not listening to
  *            Modification", this setting provides a (crude but
  *            functional) way to ensure consistency for this use case.
  *
  * - intTrgtMng [0]: bit-wise encoding of several details of the algorithm
  *                   pertaining to how the upper and lower model are used to
  *   set targets and accuracies for the compute()-tion of the components.
  *     bit 0: if 1, then a lower/upper target is set to each convex/concave
  *            component (the role of lower and upper target is exchanged
  *            for convex and concave functions)
  *     bit 1: if 1, then a upper/lower target is set to each convex/concave
  *            component (the role of lower and upper target is exchanged
  *            for convex and concave functions)
  *     bit 2: if 1, then the accuracy in compute() is set as the difference
  *            between the upper and lower target save for what bit 3 says
  *     bit 3: if 1, then the accuracy in compute() is set as (almost) the
  *            current achieved accuracy (EpsU as computed via tStar);
  *            note that if both bit 2 and bit 3 are 1, then the accuracy is
  *            set as the *maximum* between the two
  *     bit 4: if 1, then the Lipschitz constant is requested to all the
  *            components and used to compute the value of the upper model,
  *            which is useful to perform SS without having to necessarily
  *            compute() all the components; this is optional in that the
  *            computation of the Lipschitz constant may be costly
  *     bit 5- : the following bits encode different formulae for computing
  *              the "\beta_k" factors, i.e., a positive number for each
  *              component summing to one telling how the required global
  *              decrease or increase is partitioned between the different
  *              components. currently available formulae are:
  *
  *              = 0: all \beta_k = 1 / number of non-easy components
  *
  *              note that the values of these bits is only significant if
  *              at least one among the bits 0, 1 and 2 is 1.
  *
  * - intMPStbl [0]: type of stabilization to be used for the Master Problem.
  *                  Please see [MasterProblemBlock.h:207] for the currently
  *                  implemented stabilization type.
  *
  * - intMPPrimal [0]: tells which formulation should be used for the Master
  *                    Problem. If 1 then the primal version of the MP will be
  *                    initialized, otherwise MasterProblemBlock will use the
  *                    dual one. Note that if the parameter DoEasy is set to 1
  *                    and there are easy components, then the only possibility
  *                    is to solve the dual representation.
  *
  * - intRstAlg [2]: parameter to handle the reset of the algorithm when
  *                  a new Block is set, bit-wise coded:
  *                  0 bit == 1 -> don't reset algorithmic parameters
  *                  1 bit == 1 -> set current point to using current values
  *                                of the Variable (otherwise reset to all-0)
  *
  * - intMPV2Form [0]: storage frame used by the dual Master Problem:
  *                    0 selects the displacement form, while 1 selects the
  *                    iterate form. The parameter has no effect on the primal
  *                    Master Problem.
  *
  * - intMPHScaling [0]: bit-wise numerical scaling of the
  *                     PolyhedralFunctionBlock representing each hard
  *                     component in the Master Problem:
  *                     bit 0 enables local row scaling, bit 1 enables global
  *                     epigraph scaling. Hence 0 = none, 1 = local only,
  *                     2 = global only, 3 = both.
  * 
  * - intMaxLevelNR [5]: TBD (max steps of NR for level method allowed)
  */

 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameters of BundleSolver
 /** Set the double parameters specific of BundleSolver, together
  * with the parameters of CDASolver that BundleSolver actually
  * "listens to":
  *
  * - dblMaxTime [Inf< double >()]: maximum CPU time for the next call to
  *                               compute(), in seconds
  *
  * - dblRelAcc [1e-6]: relative accuracy for declaring a solution optimal
  *                     (the "easy part", see dbltStar below for the
  *                     "complicated part")
  *
  * - dblAbsAcc [Inf< double >()]: absolute accuracy for declaring a solution
  *                              optimal; if INF, it is disabled; see dbltStar
  *   below for a detailed recount on how this is used
  *
  * - dblEveryTTm [0]: periodicity of eEveryTTime events
  *
  * - dblNZEps [0]:    parameter controlling when the norm of the aggregated
  *                    subgradient z* is declared to be "almost 0". See
  *   intWZNorm for the details of how this is done in terms of which norm
  *   is used and how this constant is treated, as well as on the impact
  *   it has on the ability of BundleSolver to declare globally
  *   valid lower bounds (for a minimization problem, upper bounds for a
  *   maximization one). Choosing a very small value for dblZNEps may
  *   result in BundleSolver not being able to declare any
  *   global valid lower bound (especially if the alternative stopping
  *   criterion is used, see dbltStar), but on the other hand using a
  *   loose tolerance may result in declaring invalid global upper bound.
  *
  *   A relevant case where choosing a fair value for dblNZEps can be easier
  *   is that when the C05Function(s) is (are) Lagrangian function(s), since
  *   then z* is the violation of the relaxed constraints. Hence, the right
  *   value for dblNZEps (and its actual form, see intWZNorm) is the one
  *   corresponding to the tolerance required for declaring a solution of the
  *   original problem feasible.
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
  *   iteration, see intBPar6 for details  *
  *
  * - dblm1 [0.01]: m1 factor in the all-important NS/SS decision. This
  *                 factor sets the lower target for the objective function
  *   value at the new iterate Lambda1, with the following formula: if
  *   m1 > 0, then
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
  *   Fi( Lambda1 ) is found that is <= upper_target (which means that a
  *   "sizable decrease" has been achieved), a Serious Step can be safely
  *   declared. Note that it must be 0 < m1 <= m2 < 1; m1 = m2 = 0 is
  *   theoretically possible, but not practically advisable, for
  *   polyhedral functions provided that both function values and v* are
  *   computed without numerical errors (which is typically impossible).
  *   Also, note that whenever m1 < m2 both a SS and a NS may be possible
  *   at the same time, in which case the BundleSolver will
  *   typically favor the SS.
  *
  * - dblm3 [0.99]: factor governing the Noise Reduction for "unfaithful"
  *                 oracles that pretend to provide information with the
  *   required accuracy but in fact they do not. This results in negative
  *   linearization errors, and therefore possibly in detecting directions
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
  *                   after a SS the new value of t must be <= t * dblmxIncr
  *   (t is the previous value). It must be dblmxIncr >= dblmnIncr > 1.
  *
  * - dblmnIncr [1.5]: minimum increasing t-factor: each time t is increased
  *                    after a SS the new value of t must be >= t * dblmnIncr
  *   (t is the previous value). It must be dblmxIncr >= dblmnIncr > 1.
  *
  * - dblmxDecr [0.1]: maximum decreasing t-factor: each time t is decreased
  *                    after a NS the new value of t must be >= t * dblmxDecr
  *   (t is the previous value). It must be dblmxDecr <= dblmnDecr < 1.
  *
  * - dblmnDecr [0.66]: maximum decreasing t-factor: each time t is decreased
  *                     after a NS the new value of t must be <= t * dblmnDecr
  *   (t is the previous value). It must be dblmxDecr <= dblmnDecr < 1.
  *
  * - dbltMaior [1e+6]: maximum value of t
  *
  * - dbltMinor [1e-6]: minimum value of t
  *
  * - dbltInit [1]: initial value of t. Choosing the "right" initial value
  *                 of t can clearly help the BundleSolver to perform
  *   better in the initial iterations, although the t-strategies should see
  *   to the fact that blatantly wrong t values are rapidly corrected.
  *   Giving a reasonable value to this parameter (and, consequently, to
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
  * - dbltSPar3 [0]: numerical parameter for "small heuristic changes in t".
  *                  The general gist of t-management is that there are three
  *   hierarchical levels: 1) on top, the long-term t-strategies that try to
  *   consider the stage at which the optimization process is [see inttSPar1
  *   and dbltSPar2 for details]; 2) in the middle, multi-step strategies
  *   that take into account what has happened "in the last few iterations"
  *   [see intMnSSC, intMnNSC, dblmxIncr, dblmnIncr, dblmxDecr, dblmnDecr];
  *   3) at the bottom, heuristic t-strategies that only look at information
  *   from the very current iteration [see inttSPar1]. In particular, both
  *   the top and the middle layer would stipulate that t cannot change unless
  *   conditions are met, and when it does the change satisfies other
  *   conditions (must be larger/smaller than ...), and then the heuristics
  *   are only used to help in picking the actual value. However, it can be
  *   useful to allow the heuristics to do "small changes in t" whatever the
  *   other levels dictate. This is the meaning of this parameter: if
  *   abs( dbltSPar3 ) > 1, then the heuristics are allowed to change the
  *   value of t even if it should be fixed according to the other mechanisms.
  *   The rub is that the change should be "small", that is:
  *
  *    = if dbltSPar3 > 1, then t can only get values in the interval
  *      [ t / dbltSPar3 , t * dbltSPar3 ] (for instance, with dbltSPar3 = 2
  *      t can at most double or halve); note that this implies that t can
  *      decrease during a SS and increase during a NS, which is contrary to
  *      usual practices but not necessarily wrong (the heuristics may "know
  *      better");
  *
  *    = if dbltSPar3 < - 1, then t can only be increased if a SS is being
  *      performed and decreased in a NS is being performed; in particular,
  *      for a SS t can be chosen in [ t , t * ( - dbltSPar3 ) ], while for
  *      a NS t can be chosen in [ t / ( - dbltSPar3 ) , t ].
  *
  *   Any value of dbltSPar3 such that abs( dbltSPar3 ) <= 1 is equivalent
  *   to 0, which means "t cannot be changed by the heuristics only".
  *
  * Level stabilization keeps an expected decrease Delta and installs the
  * level
  *
  *        L = Fi( Lambda ) - Delta .
  *
  * There are two regimes. If a reliable lower bound LB is known, Delta is
  * driven by the certified gap Fi( Lambda ) - LB. If no reliable LB is known,
  * Delta is exogenous: it is initialized heuristically and adjusted until the
  * level master either produces useful steps or proves a first valid LB. Pure
  * level stabilization does not use the usual t-strategy machinery; the
  * consecutive SS/NS counters are reused to decide when a level update is
  * significant enough to be performed.
  *
  * If the level master is empty, the current L is a valid lower bound. The
  * solver records it, switches to the reliable-LB regime if necessary, and
  * recomputes Delta from the certified-gap formula below.
  *
  * - dblLStabM [0.5]: m_l parameter in (0,1) for level stabilization.
  *   When a reliable lower bound LB is available, the expected decrease is
  *
  *        Delta = (1 - m_l) * ( Fi( Lambda ) - LB )
  *
  *   and the level is L = Fi( Lambda ) - Delta. Under a Null Step, after
  *   the usual consecutive-NS gate controlled by intMnNSC allows a significant
  *   stabilization update, Delta is shortened as
  *
  *        Delta <- m_l * Delta .
  *
  *   Under a Serious Step, if LB is reliable and the usual consecutive-SS gate
  *   controlled by intMnSSC allows a significant update, Delta is capped at
  *   the new centre:
  *
  *        Delta <- min( Delta ,
  *                      (1 - m_l) * ( Fi( Lambda+ ) - LB ) ) .
  *
  * - dblLStabDlt [0.1]: fallback exogenous Delta fraction used while no
  *   reliable lower bound is known. If the one-shot initial level probe cannot
  *   provide a positive predicted decrease, the heuristic value is
  *
  *        Delta = dblLStabDlt * max( | Fi( Lambda ) | , 1 ) .
  *
  *   This exogenous initialization is abandoned as soon as the solver
  *   discovers a reliable lower bound, either from a global certificate or
  *   because the current level master is empty.
  *
  * - dblLStabIncr [2.0]: multiplicative factor used only by pure level
  *   stabilization while no reliable lower bound is known. If too many
  *   consecutive Serious Steps have been performed, i.e., the intMnSSC gate is
  *   open, and the model is accurate enough according to dblLStabSmall, the
  *   exogenous expected decrease is enlarged as
  *
  *        Delta <- dblLStabIncr * Delta .
  *
  *   The consecutive-SS counter is reset when this significant level update is
  *   performed.
  *
  * - dblLStabSmall [1e-2]: threshold used with dblLStabIncr while no reliable
  *   lower bound is known. The Delta increase above is performed only if the
  *   relative model-error ratio
  *
  *        ( Fi( Lambda1 ) - Fi_{B,Lambda}( d* ) )
  *        ------------------------------------------- <= dblLStabSmall
  *             max( | Fi_{B,Lambda}( d* ) | , 1 )
  *
  *   is small enough. Here Fi_{B,Lambda} is the cutting-plane lower model
  *   built from the current bundle B at the stability centre Lambda, and d*
  *   is the master displacement to Lambda1. Negative numerator values are
  *   clipped to zero in the implementation. Internally BundleSolver stores
  *   vStar as a predicted improvement, hence
  *
  *        Fi_{B,Lambda}( d* ) = Fi( Lambda ) + vStar .
  */

 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// move the string parameters of BundleSolver
 /** Move in the string parameters specific of BundleSolver,
  * together with the parameters of CDASolver that BundleSolver
  * actually "listens to":
  *
  * - strEasyCfg [empty]: filename from where the Configuration of the
  *                       "easy" components is taken. If not empty,
  *   strEasyCfg must be a filename (of either a text file or a netCDF
  *   one) out of which a ComputeConfiguration is loaded via a call to
  *   Configuration::deserialize( const std::string ); the ComputeConfig
  *   is then set to each of the "easy" components of the problem (via a
  *   call to set_ComputeConfig) at the time in which the
  *   BundleSolver is registered to the Block; if strEasyCfg
  *   is empty or deserialize() returns nullptr, then set_ComputeConfig()
  *   is not called.
  *
  * - strHardCfg [empty]: filename from where the Configuration of the
  *                       non-easy components is taken. If not empty,
  *   strHardCfg must be a filename (of either a text file or a netCDF one)
  *   out of which a ComputeConfiguration is loaded via a call to
  *   Configuration::deserialize( const std::string ); the ComputeConfig is
  *   then set to each of the non-easy components of the problem (via a call
  *   to set_ComputeConfig) at the time in which the BundleSolver is
  *   registered to the Block; if strHardCfg is empty or deserialize()
  *   returns nullptr, then set_ComputeConfig() is not called.
  *
  * - strMPBSolverCfg [empty]: filename from where the Configuration of
  *                            the Solver associated with the
  *   MasterProblemBlock is taken. If not empty, strMPBSolverCfg must be
  *   a filename from which a BlockSolverConfiguration is loaded via a
  *   call to Configuration::deserialize( const std::string ); the
  *   BlockSolverConfiguration is then used to identify which Solver
  *   should be attached to the MasterProblemBlock at the time in which
  *   the BundleSolver is registered to the Block. If
  *   strMPBSolverCfg is empty or deserialize() returns nullptr, the
  *   caller is required to attach a suitable Solver to the
  *   MasterProblemBlock by other means before compute() is invoked. */

 void set_par( idx_type par , std::string && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// move in the vector-of-int parameters of BundleSolver
 /** Move in the given vector-of-int parameters specific of
  * BundleSolver, together with the parameters of CDASolver that
  * BundleSolver actually "listens to":
  *
  * - vintNoEasy [empty]: the vector vintNoEasy is assumed to contain the
  *                       indices (numbers in 0, ..., total number of
  *   components - 1, ordered in increasing sense and therefore not repeated)
  *   of the components of the problem that must not be treated as "easy"
  *   even if they could. This clearly only applies if intDoEasy == 1,
  *   for otherwise no component is ever treated as "easy". The ordering of
  *   the components is as follows: if the Block only has a C05Function as
  *   Objective and no sub-Block then that it is component 0, otherwise the
  *   component i corresponds to the (the C05Function found in the
  *   FReal)Objective found in the i-th sub-Block of the Block
  *   (get_nested_Block( i )). */

 void set_par( idx_type par , std::vector< int > && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// move in the vector-of-string parameters of BundleSolver
 /** Move in the given vector-of-string parameters specific of
  * BundleSolver, together with the parameters of CDASolver that
  * BundleSolver actually "listens to":
  *
  * - vstrCmpCfg [empty]: the vector vstrCmpCfg is assumed to contain the
  *                       filenames of ComputeConfig to be passed to each
  *   of the components of the Block. The correspondence between vstrCmpCfg
  *   is positional: if vstrCmpCfg[ i ] is not empty then it is taken as a
  *   filename, used to load a ComputeConfig (with
  *   Configuration::deserialize( string ), see the comments to the
  *   method for details of the possible format of the string and of the
  *   supported file types) and used to configure the i-th component. The
  *   ordering of the components is as follows: if the Block only has a
  *   C05Function as Objective and no sub-Block then that it is component 0,
  *   otherwise the component i corresponds to the (the C05Function found
  *   in the FReal)Objective found in the i-th sub-Block of the Block
  *   (get_nested_Block( i )). This parameter overrides any other
  *   Configuration-related parameter, in particular strEasyCfg and
  *   strHardCfg: if a nonempty string is passed, then any other
  *   Configuration is ignored. If, instead, vstrCmpCfg[ i ].empty(), then
  *   the applicable one between strEasyCfg and strHardCfg, if nonempty,
  *   is applied. If everything is empty then the component is not
  *   configured (the existing configuration is not changed). Note that it
  *   is still possible to completely reset a component by passing it an
  *   "empty" ComputeConfig with diff() == true, which is different from
  *   no ComputeConfig at all.
  *
  * - vstr_C05_SPAR_Names [empty]: [vector-of-]string parameters names that
  *                                are set differently to each C05Function
  *   when BundleSolver is register()-ed to the Block. This
  *   parameter works in tandem with vstr_C05_SPAR_Vals [see].
  *
  * - vstr_C05_SPAR_Vals [empty]: baseline values for [vector-of-]string
  *                               parameters that are set differently to
  *   each C05Function when BundleSolver is register()-ed to
  *   the Block. This parameter works in tandem with vstr_C05_SPAR_Names
  *   as follows. They must have the same length. Then, for every
  *   h = 0, 1, ...,
  *   n_components() - 1, the parameter vstr_C05_SPAR_Names[ i ] is set to
  *   value <prefix>"_h"<suffix>, where <prefix> is the first part of
  *   vstr_C05_SPAR_Vals[ i ] up until the rightmost "." (if any) excluded,
  *   while <suffix> is the last part of vstr_C05_SPAR_Vals[ i ] from the
  *   rightmost "." included up untile the end (empty if there is no ".").
  *   The parameter is set as a vector-of-string one containing one single
  *   value if the first 4 characters of vstr_C05_SPAR_Names[ i ] are
  *   exactly "vstr", and as a single string parameter otherwise. This is
  *   geared towards setting different filenames (e.g., log files,
  *   Configuration files, instance files, ...) to each of the compute() of
  *   each C05Function.
  *
  * - vstr_C05_EI_SPAR_Names [empty]: [vector-of-]string parameters names
  *                                   that are set differently to each
  *   C05Function at every iteration (and call). This parameter works in
  *   tandem with vstr_C05_EI_SPAR_Vals [see].
  *
  * - vstr_C05_EI_SPAR_Vals [empty]: baseline values for [vector-of-]string
  *                                  parameters that are set differently to
  *   each C05Function at every iteration (and call). This parameter works in
  *   tandem with vstr_C05_EI_SPAR_Names as follows. They must have the same
  *   length. Then, for every h = 0, 1, ..., n_components() - 1, at every
  *   iteration k [see get_elapsed_iterations()] of every call z [see
  *   get_elapsed_calls()], the parameter vstr_C05_EI_SPAR_Names[ i ] is set
  *   to value <prefix>"_h_z_k"<suffix>, where <prefix> is the first part of
  *   vstr_C05_EI_SPAR_Vals[ i ] up until the rightmost "." (if any) excluded,
  *   while <suffix> is the last part of vstr_C05_EI_SPAR_Vals[ i ] from the
  *   rightmost "." included up untile the end (empty if there is no ".").
  *   The parameter is set as a vector-of-string one containing one single
  *   value if the first 4 characters of vstr_C05_EI_SPAR_Names[ i ] are
  *   exactly "vstr", and as a single string parameter otherwise. This is
  *   geared towards setting different filenames (e.g., log files,
  *   Configuration files, instance files, ...) to each of the compute() of
  *   each C05Function at every iteration (of every call). */

 void set_par( idx_type par , std::vector< std::string > && value ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the BundleSolver log

 void set_log( std::ostream *log_stream = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
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
  * the Block that the BundleSolver is solving. This method should
  * not be called if set_Block() has not been called, or has last been called
  * with nullptr argument. */

 Index n_components( void ) const { return( NrFi ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns a pointer to the i-th C05Function in the objective
 /** Returns a pointer to the i-th C05Function, i.e., the i-th term of the
  * sum of C05Function that (possibly together with one LinearFunction) makes up
  * the objective of the Block that the BundleSolver is solving. It
  * must ve 0 <= \p i <= n_components(). All returned pointers are not nullptr
  * provided that set_Block() has last been called with not nullptr argument
  * (otherwise this method should not be called). */

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
  * C05Function in the objective, if any. If the method returns nullptr, there
  * is no such LinearFunction (it is constantly 0, i.e., all its coefficients
  * are 0). This method should not be called if set_Block() has not been called,
  * or has last been called with nullptr argument. */

 LinearFunction * l_component( void ) const { return( f_lf ); }

/** @} ---------------------------------------------------------------------*/
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

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the MCF encoded by the current MCFBlock
 *  @{ */

 /// (try to) solve the Block

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any) and solution information
 *  @{ */

 /// returns the current number of variables in (all) the C0%Function(s)

 Index get_NumVar( void ) { return( NumVar ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the number of calls to compute() (the current one included)

 long get_elapsed_calls( void ) const override { return( SCalls ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the number of iterations in the current call to compute()
 /** Returns the number of iterations in the current call to compute(), the
  * last one included. Note that this is the number of master problem solutions,
  * as clearly the number of function evaluations may be rather different. */

 long get_elapsed_iterations( void ) const override { return( ParIter ); }

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

 double get_elapsed_time( void ) const override {
  auto end = std::chrono::system_clock::now();
  std::chrono::duration< double > elapsed = end - c_start;
  return( elapsed.count() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the best know global lower bound of the function.
 /** Returns the best know global lower bound of the function. If the
  *  objective function is convex, this is a true LB, otherwise the best
  *  solution found so far is returned. */

 VarValue get_lb( void ) override {
  if( f_convex ) {
   if( f_global_LB > - INFshift )
    return( f_global_LB + constant_value );

   // if LowerBound is not conditional (see Block.h:2728) then return it
   return( TrueLB ? LowerBound.back()  + constant_value : - INFshift );
   }
  else
   if( ( MaxSol > 1 ) && ( UpFiBest < UpFiLmb.back() ) )
    return( - UpFiBest + constant_value );
   else
    return( - UpFiLmb.back() + constant_value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the best know global upper bound of the function.
 /** Returns the best know global upper bound of the function. If the
  *  objective function is concave, this is a true UB, otherwise the best
  *  solution found so far is returned. */

 VarValue get_ub( void ) override {
  if( f_convex )
   if( ( MaxSol > 1 ) && ( UpFiBest < UpFiLmb.back() ) )
    return( UpFiBest + constant_value );
   else
    return( UpFiLmb.back() + constant_value );
  else {
   if( f_global_LB > - INFshift )
    return( - f_global_LB + constant_value );

   // if UpperBound is not conditional (see Block.h:2728) then return it
   return( TrueLB ? - LowerBound.back() + constant_value : INFshift );
   }
  }

/*--------------------------------------------------------------------------*/
 /// BundleSolver always returns a primal solution (maybe infeas.)

 bool has_var_solution( void ) override { return( true ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// BundleSolver always returns a dual solution, possibly unfeasible
 /** BundleSolver always returns a dual solution, possibly
  * infeasible. This
  *  in fact requires that the Master Problem has been solved at least once,
  *  but has_dual_solution() can only be called after compute() and therefore
  *  the Master Problem must have been solved at least once. */

 bool has_dual_solution( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/

 bool is_var_feasible( void ) override {
  return( ( Result != kInfeasible ) && ( UpFiLmb.back() < INFshift ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool is_dual_feasible( void ) override {
  return( f_global_LB > - INFshift );
  }

/*--------------------------------------------------------------------------*/
 /// write the "current" solution
 /** Write the  "current" optimal solution in the Block. There are two
  * different pieces of optimal solution:
  *
  * - The optimal values of the ColVariable in the Block;
  *
  * - If some of the C05Function are LagBFunction and are handled as "easy"
  *   components, the primal optimal solution of the dual problem that the
  *   BundleSolver is using is actually the dual optimal
  *   solution of the constraints in the master problem, and the reduced
  *   costs of the variables, that represent that component; this can be
  *   fished out of the master problem and written in the Block inside
  *   the LagBFunction (via the MILPSolver that is attached to it).
  *
  * If \p solc is nullptr, then only the optimal values of the ColVariable in
  * the Block are written. Otherwise, \p solc can be:
  *
  * - a pointer to a
  *   SimpleConfiguration< std::vector< std::pair< int , int > > >, assumed
  *   to contain pairs < i , h > such that:
  *
  *   = i is the index of an "easy" component of which the dual solution has
  *     to be written;
  *
  *   = h is coded bit-wise, with the first bit (+1) representing if the
  *     dual variables need be written, and the second bit (+2) if the
  *     reduced costs need be written;
  *
  *   note that indices i *must* be of "easy" components, otherwise an
  *   exception will be thrown;
  *
  * - a pointer to a SimpleConfiguration< int > containing the number h with
  *   the same meaning in the previous case (coded bit-wise, with the first
  *   bit (+1) representing if the dual variables need be written, and the
  *   second bit (+2) if the reduced costs need be written) which is applied
  *   to *all easy* components.
  *
  * IMPORTANT NOTE: getting access to dual variables / reduced costs of
  *                 any easy component requires access to that
  *   component's description, which is *deleted if not needed* as
  *   dictated by intDoEasy. Hence, if reduced costs are required then
  *   ( intDoEasy & 8 ) must be true; if dual variables are required then
  *   *both* ( intDoEasy & 4 ) and ( intDoEasy & 8 ) must be true, since
  *   the master pipeline reuses the same code path for reduced costs
  *   and dual variables. */

 void get_var_solution( Configuration *solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void get_var_solution_easy_pi( Index k );

 void get_var_solution_easy_rc( Index k );

/*--------------------------------------------------------------------------*/
 /// write the "current" dual solution
 /** Write the  "current" dual optimal solution in the Block. This is done
  * in two different ways for "easy" and "not easy" components:
  *
  * - for "not easy" ones, the optimal solution of the master problem for
  *   that component, which is a set of convex multipliers associated to
  *   the linearizations currently in the pool, is used for forming the
  *   "important linearization" of that component and adding them to the
  *   corresponding linearizations pool; this is unless the optimal
  *   aggregate linearization has been inserted in the bundle for other
  *   reasons (making space in a full bundle), in which case the
  *   coefficients of the "important linearization" of that component are
  *   just < 1 , index of the optimal aggregate linearization > (that is,
  *   this costs nothing);
  *
  * - for "easy" components, the dual solution for the dual problem that
  *   BundleSolver is solving is actually the primal optimal
  *   solution of the master problem for the variables that represent
  *   that component; this is fished out of the master problem and
  *   directly written in the Block inside the LagBFunction.
  *
  * If \p solc is nullptr, then the solution for all components is written.
  * Otherwise, \p solc must be a pointer to a SimpleConfiguration< std::vector<
  * int > >, assuming to contain the indices of the components of which the
  * primal solution has to be written.
  */

 void get_dual_solution( Configuration *solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void get_dual_solution_easy( Index k );

 void get_dual_solution_hard( Index k );

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

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR READING ALGORITHM PROGRESS DATA -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing data describing how the algorithm is progressing, useful
 *  e.g. inside events to track the algorithm performances or perform
 *  dynamic parameters adjustments
 *  @{ */

 /// returns the current value of the crucial t stabilization parameter

 VarValue get_t( void ) { return( t ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the current aggregated linearization error
 /** Returns the value of Sigma, the current aggregated linearization error.
  *  If the C05Function(s) are computed without errors (with a "faithful
  *  oracle", i.e., one that thoes not cheat on the lower bounds) then
  *  Sigma >= 0, and the current aggregated subgradient is a
  *  Sigma-subgradient at the current point. This justifies why the
  *  standard stopping condition of BundleSolver roughly
  *  speaking reads
  *
  *     Sigma is small and the norm of the aggregated subgradient is small
  *
  *  (see get_DSTS() and get_DST()) since the optimality condition would
  *  be "0 is a subgradient at the current point". A crucial issue is how
  *  to define "small" for the two terms, but this is easy at least for
  *  Sigma: in fact, if the aggregated subgradient is small, then
  *
  *    Sigma <= RelAcc * | current value of Fi |
  *
  *  means that the current point is RelAcc-optimal (relative). */

 VarValue get_Sigma( void ) { return( Sigma ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the scaled norm of the aggregated subgradient
 /** Returns D_{tStar}*( z* ), where z* is the aggregated subgradient. In
  *  the standard proximal stabilization, this is
  *
  *       ( tStar / 2 ) || z* ||_2^2
  *
  * (see get_DS()). Since z* is a Sigma-subgradient at the current point (see
  * get_Sigma()), justifies why one of the stopping condition of
  * BundleSolver is
  *
  *      Sigma + D_{tStar}*( z* ) <= RelAcc * | current value of Fi |
  *
  * assuming that tStar is a proper upper bound on the maximum possible step
  * that one could ever take in direction -z* in order to have the function
  * decrease (a rater big "if"). */

 VarValue get_DSTS( void ) { return( DSTS ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the unscaled norm of the aggregated subgradient
 /** Returns D_{1}*( z* ), where z* is the aggregated subgradient. In
  *  the standard proximal stabilization, this is || z* ||_2^2 / 2. See
  *  get_DSTS() and get_Sigma() for how this number enters in the stopping
  *  conditions of the method. */

 VarValue get_DS( void ) const { return( read_DStart( 1 ) ); }
 /**< D*_{1}(z*) = || z* ||^2 / 2 under standard proximal stabilization;
  *   see the doc on dbltStar. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the current estimate of the (un-)optimality of the current point
 /** Returns the current best estimate of how far from optimality the current
  *  point is, i.e.,
  *
  *    EpsU = Sigma + D_{tStar}*( z* ) / max( | current value of Fi | , 1 )
  *
  *  (see get_Sigma() and get_DSTS()). tStar is a proper upper bound on the
  *  maximum possible step that one could ever take in direction -z* in order
  *  to have the function decrease (a rater big "if"), then
  *
  *    EpsU <= RelAcc
  *
  *  should reliably indicate that the current point is RelAcc-optimal
  *  (relative). */

 VarValue get_EpsU( void ) { return( EpsU ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the number of consecutives Serious Steps
 /** Returns the current number of consecutives Serious Steps; this is > 0
  *  only if the last iteration has been a Serious Steps, as the counter is
  *  immediately reset when a Null Step is performed. Note that some
  *  t-strategies cause the master problem to be re-solved without the
  *  C05Function(s) being evaluated, these do not count as iterations. */

 Index get_CSSCntr( void ) { return( CSSCntr ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of consecutives Null Steps
 /** Returns the current number of consecutives Null Steps; this is > 0
  *  only if the last iteration has been a Null Steps, as the counter is
  *  immediately reset when a Serious Step is performed. Note that some
  *  t-strategies cause the master problem to be re-solved without the
  *  C05Function(s) being evaluated, these do not count as iterations. */

 Index get_CNSCntr( void ) { return( CNSCntr ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the predicted improvement of the model at the tentative point

 VarValue get_vStar( void ) { return( vStar.back() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the norm of the first inserted subgradient in the last iteration

 VarValue get_G1Norm( void ) { return( G1Norm ); }

/** @} ---------------------------------------------------------------------*/
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

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastBndSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastBndSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastBndSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_vint_par( void ) const override {
  return( idx_type( vintLastBndSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_vstr_par( void ) const override {
  return( idx_type( vstrLastBndSlvPar ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  static const std::array< int , 20 > dflt_int_par = {
    10 ,  // intBPar1
   100 ,  // intBPar2
     1 ,  // intBPar3
     1 ,  // intBPar4
     0 ,  // intBPar6
     3 ,  // intBPar7
     0 ,  // intMnSSC
     3 ,  // intMnNSC
    12 ,  // inttSPar1
     2 ,  // intMaxNrEvls
     1 ,  // intDoEasy
     2 ,  // intWZNorm
     0 ,  // intFrcLstSS
     0 ,  // intTrgtMng
     2 ,  // intMPStbl (default value is DoublyStabilized)
     0 ,  // intMPPrimal (default value is dual)
     2 ,  // intRstAlg, default value:
          // RstAlg = 0  -  reset algorithmic parameters
          // RstCrr = 1  -  set current point to using values of the Variable
     0 ,  // intMPV2Form (default value is displacement form)
     0 ,  // intMPHScaling (default value is no hard-component scaling)
     5    // intMaxLevelNR
     };

  if( ( par >= intLastParCDAS ) && ( par < intLastBndSlvPar ) )
   return( dflt_int_par[ par - intLastParCDAS ] );

  return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override {
  static const std::array< double , 20 > dflt_dbl_par = {
   0 ,      // dblNZEps
   1e+2 ,   // dbltStar
   0 ,      // dblMinNrEvls
   30 ,     // dblBPar5
   0.01 ,   // dblm1
   0.99 ,   // dblm2
   0.99 ,   // dblm3
   10 ,     // dblmxIncr
   1.5 ,    // dblmnIncr
   0.1 ,    // dblmxDecr
   0.66 ,   // dblmmDecr
   1e+6 ,   // dbltMaior
   1e-6,    // dbltMinor
   1 ,      // dbltInit
   1e-3 ,   // dbltSPar2
   0 ,      // dbltSPar3
   0.5 ,    // dblLStabM
   0.1 ,    // dblLStabDlt
   2.0 ,    // dblLStabIncr
   1e-2     // dblLStabSmall
   };

  if( ( par >= dblLastParCDAS ) && ( par < dblLastBndSlvPar ) )
   return( dflt_dbl_par[ par - dblLastParCDAS ] );

  return( CDASolver::get_dflt_dbl_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & get_dflt_str_par( idx_type par ) const override {
  static std::string __empty;
  if( ( par == strEasyCfg ) || ( par == strHardCfg ) ||
        ( par == strMPBSolverCfg ) )
   return( __empty );

  return( CDASolver::get_dflt_str_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* !! not necessary so far: CDASolver and Solver do not have vector-of-int
 *    parameters and the default is empty anyway

 const std::vector< int > & get_dflt_vint_par( idx_type par ) const override {
  static std::vector< int > __empty;
  if( par == vintNoEasy )
   return( __empty );

  return( CDASolver::get_dflt_vint_par( par ) );
  }
!!*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* !! not necessary so far: CDASolver and Solver do not have vector-of-string
 *    parameters and the default is empty anyway

 const std::vector< std::string > & get_dflt_vstr_par( idx_type par )
  const override {
  static std::vector< std::string > __empty;
  if( par == vstrCmpCfg )
   return( __empty );

  return( CDASolver::get_dflt_vstr_par( par ) );
  }
!!*/

/*--------------------------------------------------------------------------*/

 [[nodiscard]] int get_int_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::vector< int > & get_vint_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::vector< std::string > & get_vstr_par(
                                               idx_type par ) const override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  static const std::map< std::string , idx_type > int_pars_map = {
   { "intBPar1" , BundleSolver::intBPar1  } ,
   { "intBPar2" , BundleSolver::intBPar2  } ,
   { "intBPar3" , BundleSolver::intBPar3 } ,
   { "intBPar4" , BundleSolver::intBPar4 } ,
   { "intBPar6" , BundleSolver::intBPar6 } ,
   { "intBPar7" , BundleSolver::intBPar7 } ,
   { "intMnSSC" , BundleSolver::intMnSSC } ,
   { "intMnNSC" , BundleSolver::intMnNSC } ,
   { "inttSPar1" , BundleSolver::inttSPar1 } ,
   { "intMaxNrEvls" , BundleSolver::intMaxNrEvls } ,
   { "intDoEasy" , BundleSolver::intDoEasy } ,
   { "intWZNorm" , BundleSolver::intWZNorm } ,
   { "intFrcLstSS" , BundleSolver::intFrcLstSS } ,
   { "intTrgtMng" , BundleSolver::intTrgtMng } ,
   { "intMPStbl" , BundleSolver::intMPStbl } ,
   { "intMPPrimal" , BundleSolver::intMPPrimal } ,
   { "intRstAlg" , BundleSolver::intRstAlg } ,
   { "intMPV2Form" , BundleSolver::intMPV2Form } ,
   { "intMPHScaling" , BundleSolver::intMPHScaling } ,
   { "intMaxLevelNR" , BundleSolver::intMaxLevelNR }
   };

  const auto it = int_pars_map.find( name );
  if( it != int_pars_map.end() )
   return( it->second );
  else
   return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override {
  static const std::map< std::string , idx_type > dbl_pars_map = {
   { "dblNZEps" , BundleSolver::dblNZEps } ,
   { "dbltStar" , BundleSolver::dbltStar } ,
   { "dblMinNrEvls" , BundleSolver::dblMinNrEvls } ,
   { "dblBPar5" , BundleSolver::dblBPar5 } ,
   { "dblm1" , BundleSolver::dblm1 } ,
   { "dblm2" , BundleSolver::dblm2 } ,
   { "dblm3" , BundleSolver::dblm3 } ,
   { "dblmxIncr" , BundleSolver::dblmxIncr } ,
   { "dblmnIncr" , BundleSolver::dblmnIncr } ,
   { "dblmxDecr" , BundleSolver::dblmxDecr } ,
   { "dblmnDecr" , BundleSolver::dblmnDecr } ,
   { "dbltMaior" , BundleSolver::dbltMaior } ,
   { "dbltMinor" , BundleSolver::dbltMinor } ,
   { "dbltInit" , BundleSolver::dbltInit } ,
   { "dbltSPar2" , BundleSolver::dbltSPar2 } ,
   { "dbltSPar3" , BundleSolver::dbltSPar3 } ,
   { "dblLStabM" , BundleSolver::dblLStabM } ,
   { "dblLStabDlt" , BundleSolver::dblLStabDlt } ,
   { "dblLStabIncr" , BundleSolver::dblLStabIncr } ,
   { "dblLStabSmall" , BundleSolver::dblLStabSmall }
   };

  const auto it = dbl_pars_map.find( name );
  if( it != dbl_pars_map.end() )
   return( it->second );

  return( CDASolver::dbl_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override {
  if( name == "strEasyCfg" )
   return( strEasyCfg );
  if( name == "strHardCfg" )
   return( strHardCfg );
  if( name == "strMPBSolverCfg" )
   return( strMPBSolverCfg );

  return( CDASolver::str_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type vint_par_str2idx( const std::string & name )
  const override {
  if( name == "vintNoEasy" )
   return( vintNoEasy );

  return( CDASolver::vint_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type vstr_par_str2idx( const std::string & name )
  const override {
  static const std::map< std::string , idx_type > vstr_pars_map = {
   { "vstrCmpCfg" , BundleSolver::vstrCmpCfg } ,
   { "vstr_C05_SPAR_Names" , BundleSolver::vstr_C05_SPAR_Names } ,
   { "vstr_C05_SPAR_Vals" , BundleSolver::vstr_C05_SPAR_Vals } ,
   { "vstr_C05_EI_SPAR_Names" ,
                 BundleSolver::vstr_C05_EI_SPAR_Names } ,
   { "vstr_C05_EI_SPAR_Vals" , BundleSolver::vstr_C05_EI_SPAR_Vals }
   };

  const auto it = vstr_pars_map.find( name );
  if( it != vstr_pars_map.end() )
   return( it->second );

  return( CDASolver::vstr_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::array< std::string , 20 > int_pars_str = {
   "intBPar1" , "intBPar2" , "intBPar3" , "intBPar4" , "intBPar6" ,
   "intBPar7" , "intMnSSC" , "intMnNSC" , "inttSPar1" , "intMaxNrEvls" ,
   "intDoEasy" , "intWZNorm" , "intFrcLstSS" , "intTrgtMng" ,
   "intMPStbl" , "intMPPrimal" , "intRstAlg" , "intMPV2Form" ,
   "intMPHScaling" , "intMaxLevelNR" };

  if( ( idx >= intLastParCDAS ) && ( idx < intLastBndSlvPar ) )
   return( int_pars_str[ idx - intBPar1 ] );

  return( CDASolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override {
  static const std::array< std::string , 20 > dbl_pars_str = {
   "dblNZEps" , "dbltStar" , "dblMinNrEvls" , "dblBPar5" , "dblm1" ,
   "dblm2" , "dblm3" , "dblmxIncr" , "dblmnIncr" , "dblmxDecr" ,
   "dblmnDecr" , "dbltMaior" , "dbltMinor" , "dbltInit" , "dbltSPar2" ,
   "dbltSPar3" , "dblLStabM" , "dblLStabDlt" , "dblLStabIncr" ,
   "dblLStabSmall" };

 if( ( idx >= dblLastParCDAS ) && ( idx < dblLastBndSlvPar ) )
   return( dbl_pars_str[ idx - dblLastParCDAS ] );

  return( CDASolver::dbl_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override {
  static const std::array< std::string , 3 > str_pars_str = {
   "strEasyCfg" , "strHardCfg" , "strMPBSolverCfg" };

  if( ( idx >= strLastParCDAS ) && ( idx < strLastBndSlvPar ) )
   return( str_pars_str[ idx - strLastParCDAS ] );

  return( CDASolver::str_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & vint_par_idx2str( idx_type idx )
  const override {
  static const std::string __psname = "vintNoEasy";
  if( idx == vintNoEasy )
   return( __psname );

  return( CDASolver::vint_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & vstr_par_idx2str( idx_type idx )
  const override {
  static const std::array< std::string , 5 > vstr_pars_str = { "vstrCmpCfg" ,
   "vstr_C05_SPAR_Names" , "vstr_C05_SPAR_Vals " , "vstr_C05_EI_SPAR_Names" ,
   "vstr_C05_EI_SPAR_Vals" };

  if( ( idx >= vstrLastParCDAS ) && ( idx < vstrLastBndSlvPar ) )
   return( vstr_pars_str[ idx - vstrCmpCfg ] );

  return( CDASolver::vstr_par_idx2str( idx ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------ METHODS FOR HANDLING THE State OF THE BundleSolver -----*/
/*--------------------------------------------------------------------------*/
/** @name Handling the State of the BundleSolver
 *  @{ */

 State * get_State( void ) const override;

/*--------------------------------------------------------------------------*/
 /// sets the current "internal state" of the BundleSolver
 /** This method reads \p state, which must be a BundleSolverState
  * (otherwise exception is thrown) and uses it to set the "internal state" of
  * the BundleSolver. \p state is not changed (could not, it's const)
  * so that it can be re-used later.
  *
  *     IMPORTANT NOTE: IT IS NOT ALLOWED TO CALL THIS METHOD WHILE compute()
  *                     IS RUNNING, SAY FROM WITHIN AN EVENT. */

 void put_State( const State & state ) override;

/*--------------------------------------------------------------------------*/
 /// sets the current "internal state" of the BundleSolver
 /** This method reads \p state, which must be a BundleSolverState
  * (otherwise exception is thrown) and uses it to set the "internal state" of
  * the BundleSolver. \p state is "moved inside" the
  * BundleSolver, hence at the end of the call is in a consistent
  * state but "empty", and therefore it will probably be immediately destroyed
  * as it has little remaining use.
  *
  *     IMPORTANT NOTE: IT IS NOT ALLOWED TO CALL THIS METHOD WHILE compute()
  *                     IS RUNNING, SAY FROM WITHIN AN EVENT. */

 void put_State( State && state ) override;

/*--------------------------------------------------------------------------*/

 void serialize_State( netCDF::NcGroup & group ,
                       const std::string & sub_group_name = "" )
  const override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------------- FRIENDS ---------------------------------*/
/*--------------------------------------------------------------------------*/

 // make BundleSolverState friend
 friend class BundleSolverState;

/*--------------------------------------------------------------------------*/
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

 void guts_of_put_State( const BundleSolverState & state );

/*--------------------------------------------------------------------------*/
 /// read the linearization error of the item with global name \p name
 /** Routes to MasterPB->get_alpha mapped through ItemVcblr; returns 0
  * when MasterPB is not yet available. */

 VarValue read_alpha_global( Index name ) const;

/*--------------------------------------------------------------------------*/
 /// read the optimal multiplier theta of the item with global name \p name
 /** Routes to MasterPB->get_theta mapped through ItemVcblr; returns 0
  * when MasterPB is not yet available. */

 VarValue read_theta_global( Index name ) const;

/*--------------------------------------------------------------------------*/
 /// remove the cut with global name \p name from the master
 /** Routes the call to
  * MasterPB->remove_cut( ItemVcblr[ name ].first ,
  *                      ItemVcblr[ name ].second ). */

 void remove_cut_global( Index name );

/*--------------------------------------------------------------------------*/
 /// returns the upper bound of valid global names in the bundle
 /** Used as the loop upper bound (i.e.
  * `for( i = 0 ; i < get_max_name() ; ++i ) ...`); returns the highest
  * occupied global slot index plus one, always <= ItemVcblr.size() ==
  * vBPar2.back(). The loop body still has to filter
  * `ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ]` to skip
  * empty entries. */

 Index get_max_name( void ) const {
  return( f_max_name );
  }

/*--------------------------------------------------------------------------*/
 /// component index of the item with global name \p name
 /** Returns the 1-based component index of the item with global name
  * \p name: 1 .. NrFi for a regular component, NrFi + 1 if the item
  * belongs to the linear "0-th component", Inf< Index >() otherwise.
  * The information is read straight from ItemVcblr[ name ].first + 1. */

 Index wcomponent_global( Index name ) const {
  if( name < ItemVcblr.size() &&
      ItemVcblr[ name ].second < Inf< Index >() )
   return( ItemVcblr[ name ].first + 1 );
  return( Inf< Index >() );
  }

/*--------------------------------------------------------------------------*/
 /// is the item with global name \p name a "true" subgradient?
 /** Returns true iff the item is a diagonal linearization (subgradient),
  * false if it is vertical (feasibility cut) or empty. Routes through
  * MasterPB->is_subgradient( k , slot ) via ItemVcblr. */

 bool is_subgradient_global( Index name ) const {
  if( MasterPB && name < ItemVcblr.size() ) {
   const auto & loc = ItemVcblr[ name ];
   if( loc.second < Inf< Index >() )
    return( MasterPB->is_subgradient( int( loc.first ) ,
                                      int( loc.second ) ) );
   }
  return( false );
  }

/*--------------------------------------------------------------------------*/
 /// returns D*_{factor}( z* ) = ( factor / 2 ) * || z* ||^2
 /** Returns ( factor / 2 ) * MasterPB::get_dual_norm_squared(); returns
  * 0 if MasterPB is not yet available. */

 VarValue read_DStart( double factor ) const {
  if( MasterPB )
   return( factor / 2.0 * MasterPB->get_dual_norm_squared() );
  return( 0 );
  }

/*--------------------------------------------------------------------------*/
 /// total number of vertical (feasibility) cuts in the master problem
 /** Routes to MasterPB::get_vertical_count(), which scans every
  * HardCmps[ k ] and tallies the is_row_vertical hits. */

 Index get_bc_size( void ) const {
  if( MasterPB )
   return( Index( MasterPB->get_vertical_count() ) );
  return( 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns g_name . d for the cut with global name \p name
 /** Routes through MasterPB->get_Gid( k , slot ) via ItemVcblr. */

 VarValue read_Gid_global( Index name ) const {
  if( MasterPB && name < ItemVcblr.size() ) {
   const auto & loc = ItemVcblr[ name ];
   if( loc.second < Inf< Index >() )
    return( MasterPB->get_Gid( int( loc.first ) , int( loc.second ) ) );
   }
  return( 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the aggregated z* . d
 /** Routes through MasterPB->get_Gid_aggregate(). */

 VarValue read_Gid_aggregate( void ) const {
  if( MasterPB )
   return( MasterPB->get_Gid_aggregate() );
  return( 0 );
  }

/*--------------------------------------------------------------------------*/
 /* FormD() just calls SolveMP() once and calculates the direction d:
    however, it also implements some strategies to survive to "fatal"
    failures in the subproblem solver, typically eliminating some of the
    items in the bundle.

    Set the protected field Result to kOK if (eventually after some "fatal"
    failure) a tentative descent direction could be found, to kUnfsbl if the
    MP is dual unfeasible and to kError if this was returned by SolveMP(): in
    the latter cases, the whole algorithm must abort. */

 void FormD( void );

 /*--------------------------------------------------------------------------*/

 bool UsesLevelStabilization( void ) const {
  return( MPStbl == MasterProblemBlock::kLevel ||
          MPStbl == MasterProblemBlock::kDoublyStabilized );
  }

 /*--------------------------------------------------------------------------*/

 bool UsesPureLevelStabilization( void ) const {
  return( MPStbl == MasterProblemBlock::kLevel );
  }

 /*--------------------------------------------------------------------------*/

 bool UsesPrimalMaster( void ) const {
  return( IsMPPrimal && ! ( DoEasy && ( NrEasy > 0 ) ) );
  }

 /*--------------------------------------------------------------------------*/

 VarValue reliable_level_LB( void ) const;

 /*--------------------------------------------------------------------------*/

 void reset_level_stabilization( void );

 /*--------------------------------------------------------------------------*/

 void install_level_stabilization( void );

 /*--------------------------------------------------------------------------*/

 void refresh_level_after_master( void );

 /*--------------------------------------------------------------------------*/

 void update_level_after_step( bool serious_step , bool gated_update );

 /*--------------------------------------------------------------------------*/

 void record_level_lower_bound( VarValue lb );

/*--------------------------------------------------------------------------*/
 // Updates the out-of-base counters for all items in the Bundle.

 void UpdtCntrs( void );

/*--------------------------------------------------------------------------*/
 /* After a (successful) call to FormD(), sets the new tentative point Lambda1
  * as Lambda1 = Lambda + ( Tau / t ) * d. */

 void FormLambda1( double Tau );

/*--------------------------------------------------------------------------*/
 /* Performs the inner loop: repeatedly compute components up until the
  * conditions for either a SS or a NS (or both) are satisfied based on the
  * current value of UpTrgt and LwTrgt, or something very bad happens
  * (errors, out of time, ...). Sets MPchgs > 0 if the SS and/or NS
  * conditions are satisfied.
  *
  * It is virtual because this is precisely the point where a sequential and a
  * "basic" asynchronous implementation of the approach differ, and therefore
  * this is the obvious hook for a derived AsynchBundleSolver.
  *
  * The parameter extrastep, if true, means that InnerLoop() is not called from
  * within the normal main loop, but at the end to ensure that the components
  * are "current" with the returned optimal solution. This has two effects:
  *
  * - it disables all fancy early termination checks and forces all (non-easy)
  *   components to be evaluated;
  *
  * - it only compute function values (upper and lower estimates), but it does
  *   not collect any new linearization.
  *
  * The rationale for the latter choice is that the bundle contains enough
  * linearizations to stop already, hence new ones are not needed, and in fact
  * the Master Problem is not re-solved, hence even if they were collected they
  * would not be useful. This is potentially wasteful in that one does the
  * effort to compute() the functions but only collects "a small part" of the
  * ensuing information. However, producing linearizations may have a cost in
  * itself, so on the other hand it saves some time. Besides, if the C05Function
  * are re-compute()-d right after in Lambda, then if they are "complex and
  * costly" they will likely have checks to understand if the bulk of the
  * computation can be skipped because it has been done already.
  *
  * The method returns the number of different components evaluated. It may also
  * internally set Result == kError if an irrecoverable error happens during a
  * component compute()-tion or Result == kStopTime if the available running
  * time runs up. */

 virtual Index InnerLoop( bool extrastep = false );

/*--------------------------------------------------------------------------*/
 /* Computes Fi[ wFi ]( Lambda1 ). In the "default mode", where getgi == true,
  * it also inserts the obtained linearizations in the bundle. If getgi == false
  * instead this means that method is being called outside of the main loop,
  * with Lambda1 == Lambda, and that no linearizations need be collected.
  *
  * Returns true <=> at least one item was inserted. It also "sneakily" sets
  * MPchgs if appropriate. None of this clearly happens if getgi == false. */

 bool FiAndGi( Index wFi , bool getgi = true );

/*--------------------------------------------------------------------------*/
 // Set the component-specific string parameters, if any

 void SetupFiStrPar( Index wFi );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /* Prepares component wFi for computation on Lambda1 by setting the
  * thresholds and accuracy. */

 void SetupFiLambda1( Index wFi );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /* Prepares component wFi for computation on Lambda by setting the
  * thresholds and accuracy. */

 void SetupFiLambda( Index wFi );

/*--------------------------------------------------------------------------*/
 /* Gets the new linearizations out of the freshly computed component wFi,
  * returns true <=> at least one item was inserted. It also "sneakily" sets
  * MPchgs if appropriate. */

 bool GetGi( Index wFi );

/*--------------------------------------------------------------------------*/

 void update_UpFiLambd1( Index wFi , VarValue nval )
 {
  // we assume the previous upper estimate to still be valid, hence if the
  // new upper estimate is >= then the old one nothing needs be done
  if( UpFiLmb1[ wFi ] <= nval )
   return;

  if( UpFiLmb1[ wFi ] == INFshift )
   ++UpFiLmb1def;

   if( UpFiLmb1def == NrFi ) {
   ++UpFiLmb1def;  // all components + the sum computed
   UpFiLmb1[ wFi ] = nval;
   UpFiLmb1.back() = std::accumulate( UpFiLmb1.begin() , --(UpFiLmb1.end()) ,
                                      Fi0Lmb1 );
   // note that this is the point where the lower bound (if any) is
   // "baked in" the total function value
   if( TrueLB && ( UpFiLmb1.back() < LowerBound.back() ) )
    UpFiLmb1.back() = LowerBound.back();
   }
  else {
   if( UpFiLmb1def > NrFi ) {
    // the total value was already known, but since the value of this
    // component decreases, the total value decreases by the same amount
    UpFiLmb1.back() += nval - UpFiLmb1[ wFi ];

    // in fact the previous total value might have been == (or close to)
    // the global valid lower bound, so check that it remains >= that
    if( TrueLB && ( UpFiLmb1.back() < LowerBound.back() ) )
     UpFiLmb1.back() = LowerBound.back();
    }
   UpFiLmb1[ wFi ] = nval;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void update_LwFiLambd1( Index wFi , VarValue nval )
 {
  // we assume the previous lower estimate to still be valid, hence if the
  // new lower estimate is <= then the old one nothing needs be done
  if( LwFiLmb1[ wFi ] >= nval )
   return;

  if( LwFiLmb1[ wFi ] == -INFshift )
   ++LwFiLmb1def;

  // the value of the global lower bound is in principle "baked in" the total
  // lower estimate of the function value in Lambda1: if the contribution of
  // this component grows it may bring the total above the lower bound while
  // previously it was below, so ensure the total is recomputed
  if( TrueLB && ( LwFiLmb1def > NrFi ) )
   LwFiLmb1def = NrFi;

  if( LwFiLmb1def == NrFi ) {
   ++LwFiLmb1def;  // all components + the sum computed
   LwFiLmb1[ wFi ] = nval;
   LwFiLmb1.back() = std::accumulate( LwFiLmb1.begin() , --(LwFiLmb1.end()) ,
                                      Fi0Lmb1 );
   // note that this is the point where the lower bound (if any) is
   // "baked in" the total function value
   if( TrueLB && ( LwFiLmb1.back() < LowerBound.back() ) )
    LwFiLmb1.back() = LowerBound.back();
   }
  else {
   // the total value was already known, but since the value of this
   // component increases, the total value increases by the same amount;
   // note that this can only happen if TrueLB == false, which means that
   // there can be no global lower bound "baked in" the total value
   if( LwFiLmb1def > NrFi )
    LwFiLmb1[ NrFi ] += nval - LwFiLmb1[ wFi ];
   LwFiLmb1[ wFi ] = nval;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void update_UpFiLambd( Index wFi , VarValue nval )
 {
  // we assume the previous upper estimate to still be valid, hence if the
  // new upper estimate is >= then the old one nothing needs be done
  if( UpFiLmb[ wFi ] <= nval )
   return;

  if( UpFiLmb[ wFi ] == INFshift )
   ++UpFiLmbdef;

  if( UpFiLmbdef == NrFi ) {
   ++UpFiLmbdef;  // all components + the sum computed
   UpFiLmb[ wFi ] = nval;
   UpFiLmb.back() = std::accumulate( UpFiLmb.begin() , --(UpFiLmb.end()) ,
                                     Fi0Lmb );
   // note that this is the point where the lower bound (if any) is
   // "baked in" the total function value
   if( TrueLB && ( UpFiLmb.back() < LowerBound.back() ) )
    UpFiLmb.back() = LowerBound.back();
   }
  else {
   if( UpFiLmbdef > NrFi ) {
    // the total value was already known, but since the value of this
    // component decreases, the total value decreases by the same amount
    UpFiLmb.back() += nval - UpFiLmb[ wFi ];

    // in fact the previous total value might have been == (or close to)
    // the global valid lower bound, so check that it remains >= that
    if( TrueLB && ( UpFiLmb.back() < LowerBound.back() ) )
     UpFiLmb.back() = LowerBound.back();
    }
   UpFiLmb[ wFi ] = nval;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void update_LwFiLambd( Index wFi , VarValue nval )
 {
  // we assume the previous lower estimate to still be valid, hence if the
  // new lower estimate is <= then the old one nothing needs be done
  if( LwFiLmb[ wFi ] >= nval )
   return;

  if( LwFiLmb[ wFi ] == -INFshift )
   ++LwFiLmbdef;

  // the value of the global lower bound is in principle "baked in" the total
  // lower estimate of the function value in Lambda: if the contribution of
  // this component grows it may bring the total above the lower bound while
  // previously it was below, so ensure the total is recomputed
  if( TrueLB && ( LwFiLmbdef > NrFi ) )
   LwFiLmbdef = NrFi;

  if( LwFiLmbdef == NrFi ) {
   ++LwFiLmbdef;  // all components + the sum computed
   LwFiLmb[ wFi ] = nval;
   LwFiLmb.back() = std::accumulate( LwFiLmb.begin() , --(LwFiLmb.end()) ,
                                     Fi0Lmb );
   // note that this is the point where the lower bound (if any) is
   // "baked in" the total function value
   if( TrueLB && ( LwFiLmb.back() < LowerBound.back() ) )
    LwFiLmb.back() = LowerBound.back();
   }
  else {
   // the total value was already known, but since the value of this
   // component increases, the total value increases by the same amount;
   // note that this can only happen if TrueLB == false, which means that
   // there can be no global lower bound "baked in" the total value
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
 /* Compute the value of \beta_h for the given component, i.e., which
  * fraction of the total expected decrease is expected from that specific
  * component. */

 double BetaK( Index wFi );

/*--------------------------------------------------------------------------*/

 void Log1( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void Log2( double ft );

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
 /** Concave functions to be maximised are sneakily turned into convex
  * functions to be minimized inside by changing the sign of function values and
  * linearizations, but they have to be output with the right sign. */

 VarValue rs( const VarValue fv ) {
  return( f_convex ? fv : - fv );
  }

/*--------------------------------------------------------------------------*/
 /** Finds the next component to compute in the inner loop, writes is in
  * f_wFi; returns false if there is no other component to compute. */

 bool FindNext( void );

/*--------------------------------------------------------------------------*/
 /// method implementing the short-term t-strategies
 /** The method is called after a SS if inttSPar1 & 1 and after a NS if
  * tSPar1 & 2 and should return a proposed new value for t using only
  * information pertaining to the current iteration (typically, the aggregated
  * subgradient and its linearization error, the newly obtained subgradient and
  * its linearization error). The base class implements four simple heuristics
  * chosen on the basis of the bits 6 and 7 of inttSPar1; see the comments to
  * set_par( inttSPar1 ) and those inside Heuristic1() to Heuristic4() for
  * details. However, the method is virtual so that derived classes can
  * implement new ones, typically using the remaining bits of inttSPar1 (from 8
  * on) to select them. */

 virtual double Heuristic( Index whch );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 double Heuristic1( void );

 double Heuristic2( void );

 double Heuristic3( void );

 double Heuristic4( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -

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

 Index BPar1;       ///< parameter for removal of items (B-strategy)
 Index BPar2;       ///< max Bundle size
 Index BPar3;       ///< max number of items fetched from Fi() at each call
 Index BPar4;       ///< min number of items fetched from Fi() at each call
 double BPar5;      ///< control how the actual BPar3 changes over time
 int BPar6;         ///< control how the actual BPar3 changes over time
 int BPar7;         ///< if GBS "plays nice" with other Solver

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
 double tSPar3;     ///< double parameter for small heuristic t changes

 bool DoEasy;       ///< if "easy" components are managed

 char WZNorm;       ///< how to compute the norm of z*

 int FrcLstSS;      /**<  bit 0 = 1: all components must be computed in
                     *               the optimum
                     *    bit 1 = 1: the State is not to be trusted */

 Index TrgtMng;     ///< how targets on components are managed

 MasterProblemBlock::stabilization_type MPStbl;
                    ///< type of stabilization for the master problem

 bool IsMPPrimal;   ///< whether the master problem is solved in its primal
                    ///< form (true) or in its dual one (false)

 int RstAlgPrm;     ///< reset parameter, bit-wise coded

 int MPV2Form;      ///< dual MP storage frame: 0 = displacement, 1 = iterate

 int MPHScaling;    ///< bit-wise hard-component PFB scaling: local / global

 int MaxLevelNR;    ///< maximum number of consecutive NoiseReduction step allowed

 std::string EasyCfg;
 ///< filename for the Block[Solver]Config of easy components

 std::string HardCfg;
 ///< filename for the Block[Solver]Config of non-easy components

 std::string MPBSolverCfg;
 ///< filename for the SolverConfig of the Solver to be used in the MPBlock

 std::vector< int > NoEasy;  ///< which components never treat as "easy"

 std::vector< std::string > CmpCfg;  ///< individual Configurations

 std::vector< std::string > v_C05_SPAR_Names;
 ///< string parameters names that are set differently to each C05Function

 std::vector< std::string > v_C05_SPAR_Vals;
 ///< baseline values for different string parameters to each C05Function

 std::vector< std::string > v_C05_EI_SPAR_Names;
 /**< string parameters names that are changed at every iteration and call
  * differently to each C05Function */

 std::vector< std::string > v_C05_EI_SPAR_Vals;
 /**< baseline values for string parameters changed at every iteration and
  * call differently to each C05Function */

 // generic fields- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int Result;        ///< result of the latest call to Solve()

 Index NumVar;      ///< (current) number of variables
 Index NrFi;        ///< number of components of Fi()

 Index SCalls;      ///< number of calls to compute() (the current included)
 Index ParIter;     ///< number of iterations in this call to compute()

 std::vector< bool > IsEasy;
 ///< true if the component k has found to be easy

 Index NrEasy;      ///< number of "easy" component of Fi

 std::vector< FakeSolver * > v_FakeSolver;
 ///< FakeSolver used to handle Modification from the easy components

 Vec_VarValue Lambda;   ///< the current point

 Vec_VarValue Lambda1;  ///< the tentative point

 Vec_VarValue LmbdBst;  ///< the best point found so far

 bool LHasChgd;       /**< true if Lambda has changed since the latest call
                       * to FiAndGi(): allows repeated calls in the same Lambda,
                       * e.g. with increasing precision */
 bool tHasChgd;       ///< true if t has changed since the last MP

 char MPchgs;         ///< nonzero if we can prove no cycling will occur
                      /**< MPchgs == 1 means that the conditions for
                       * ensuring that no cycle will occur have been found
  * due to the function value (a SS can be done) or a diagonal linearization (a
  * NS can be done); MPchgs == 2 means that a vertical linearization (cutting
  * off Lambda1) has been found, and this by itself ensures no cycling. */

 Subset whisZ;     /**< the position in the bundle where the "aggregate
                    * subgradient" Z[ k ] of component k is kept in whisZ[ k ];
                    * Inf< Index >() == it is not in the bundle */
 std::vector< bool > Zvalid;  /**< Zvalid[ k ] == true if the item in position
                               * whisZ[ k ] is exactly Z[ k ] as computed by
 * the last master problem. Zvalid[ k ] == true ==> whisZ[ k ] < INF. if Zvalid[
 * k ] == false and whisZ[ k ] < INF, then Z[ k ] had been previously stored in
 * position whisZ[ k ], but the master problem has been re-solved since and
 * therefore Z[ k ] is no longer current. */

 Subset whisG1;    ///< "representative subgradient" for each component
                   /**< whisG1[ k ] is the first subgradient inserted in the
                    * bundle for component k in the latest call to InnerLoop
 * in which component k has actually been evaluated. This is used as the
 * "representative subgradient" for component k for the computation of the t
 * heuristics, which use information about an "aggregate representative
 * subgradient" obtained by summing all the individual "representative
 * subgradients". NOTE 1: the "representative subgradient" used to be selected
 * as the one
 *         with smallest Alfa1k and largest ScPr1k, which may have been a
 *         bit better but was heuristic anyway; choosing the first makes
 *         for a simpler logic, and ideally "the first should be the best"
 *         so it still makes sense.
 * NOTE 2: "easy" components do not have a "representative subgradient" and
 *         are therefore excluded from the "aggregate representative
 *         subgradient": THIS IS WRONG and probably means that most t
 *         heuristics are unreliable with "easy" components, but the solution
 *         is not obvious so this is it for now.
 * NOTE 3: since not all components may be evaluated during one InnerLoop
 *         (especially in the NS case, but also in the SS one), the
 *         "representative subgradient" for one component may in fact have
 *         been inserted several iterations ago or even be missing in the
 *         (unlikely) case it has been removed in the meantime: this may
 *         not be ideal but there is not much else we can do.
 * NOTE 4: A component k s.t. whisG1[ k ] < InINF && CurrNrEvls[ k ] == 0
 *         has not been compute()-d in the latest InnerLoop and therefore
 *         its contribution to ScPr1 and Alfa1 (see below) must be added
 *         separately, as it has not inside InnerLoop(). */

 Vec_VarValue G1;  ///< the aggregate "representative subgradient"
 VarValue G1Norm;  ///< norm( G1 , 2 )

 VarValue ScPr1;   ///< ideally, ScalarProduct( dir , G1 )
 VarValue Alfa1;   /**< ideally, the linearization error of G1 w.r.t. the
                    * current point Lambda. */

 Vec_VarValue LowerBound;  ///< Lower Bound over (each component of) Fi
 VarValue f_global_LB;     ///< an algorithmically discovered global LB

 VarValue f_level_Delta = 0;
 ///< expected decrease used by level stabilization

 VarValue f_level_value = INFshift;
 ///< current level target

 VarValue f_level_LB = -INFshift;
 ///< reliable lower bound that last drove f_level_Delta

 bool f_level_reliable_LB = false;
 ///< true once f_level_Delta is driven by a reliable global LB

 bool f_level_initialized = false;
 ///< true after the first level target has been installed

 VarValue LStabM;      ///< m_l parameter for level stabilization

 VarValue LStabDlt;    ///< initial exogenous Delta fraction for level stabilization

 VarValue LStabIncr;   ///< exogenous Delta increase factor for level stabilization

 VarValue LStabSmall;  ///< threshold for small model error in level stabilization

 VarValue t;           ///< the (tremendous) t parameter
 VarValue Prevt;       ///< what t were before being changed for funny reasons

 VarValue Sigma;       ///< Sigma*: convex combination of the Alfa's
 VarValue DSTS;        /**< D*_{t*}( -z* ), the other part of the dual
                        * objective */
 Vec_VarValue vStar;   ///< v*, the predicted improvement

 VarValue DeltaFi;     ///< FiLambda - FiLambda1
 VarValue EpsU;        ///< precision required by the long-term t-strategy

 Index CSSCntr;        ///< counter of consecutive SS

 Index CNSCntr;        ///< counter of consecutive NS

 Index LevelNRCntr;    ///< counter of consecutive level NR steps

 Subset vBPar2;        ///< size of the global pools of each component

 std::priority_queue< Index , std::vector< Index > ,
                      std::greater< Index > > FreList;
 ///< list of free positions in the bundle
 /**< FreList is the priority queue of free positions in the bundle, where
  * the position with higher priority is that with smaller name. This is why the
  * comparison parameter has to be explicitly set to std::greater< Index >,
  * since a priority queue with the default std::less< Index > spits out the
  * element with larger value. */

 /** NrItems[ k ] contains the number of items in the bundle (master problem)
  * for component k. If ( BPar7 & 3 ) < 3, this number may be strictly less than
  * the number of linearizations in the global pool of component k, since
  * removals from the bundle do not imply removals from the global pool
  */

 Subset NrItems;  ///< number of items in the bundle for each component

 /** FictLB[ k ] tracks whether a *fictitious* model lower bound is
  * currently installed on hard component k inside MasterProblemBlock,
  * to keep the dual master feasible while k's bundle is transiently
  * empty and k has no genuine individual lower bound. Set when k is
  * empty, cleared as soon as it receives a real cut. Only meaningful
  * when MasterPB != nullptr; sized NrFi. */

 std::vector< bool > FictLB;

 /** FrFItem[ k ] contains the index of the first position in the global
  * pool of component k where BundleSolver can put a new
  * linearization when the corresponding item is added to the bundle
  * (master problem). Note that whether a position is suitable to this
  * depends on BPar2: if ( BPar7 & 3 == 0 ), then the position must be
  * completely empty (InvItemVcblr[ k ][ i ] == INF), while if
  * ( BPar7 & 3 != 0 ) then another linearization can be in that position
  * already provided that there is no corresponding item in the bundle
  * (vBPar2[ k ] <= InvItemVcblr[ k ][ i ] < INF). */

 Subset FrFItem;  ///< the first free item in each global pool

 /** MaxItem[ k ] contains 1 + the maximum index of a position in the
  * global pool of component k where a linearization is stored; if MaxItem[ k ]
  * == 0, then the global pool is empty. Note that this ignores the fact that a
  * position in the global pool corresponds or not to an item in the bundle: all
  * linearizations count. As a consequence it must always be MaxItem[ k ] >=
  * FrFItem[ k ]. */

 Subset MaxItem;  ///< the first unused item in each global pool

 /** Vocabulary of items: ItemVcblr[ i ].first is the component name
  * and ItemVcblr[ i ].second is the name in the global pool of that component
  * for item in position i of the bundle (master problem). ItemVcblr[ i ].second
  * == INF means that position i in the bundle is not used. */

 std::vector< std::pair< Index , Index > > ItemVcblr;

 /** Inverse vocabulary of items. InvItemVcblr[ k ] is a Subset of size
  * vBPar2[ k ] and describes the global pool of component k. With p =
  * InvItemVcblr[ k ][ i ], if p < vBPar2[ NrFi ], then the linearization with
  * name i in the global pool of h is in the bundle at position p. If p == INF,
  * then there is no linearization with name i in the global pool of k. If
  * vBPar2[ NrFi ] <= p < INF, then there is a linearization with name i in the
  * global pool of h, but it is not in the bundle.
  *
  * NOTE: THE GLOBAL POOL OF SOME C05Function CAN BE LARGER THAN vBPar2[ k ],
  * BUT ALL ELEMENTS WITH NAME LARGER THAN vBPar2[ k ] ARE NEVER USED OR CHANGED
  * BY BundleSolver. */

 std::vector< Subset > InvItemVcblr;

 /** Highest global slot index occupied in ItemVcblr plus one, i.e. the
  * current upper bound of valid global names. Always in
  * [0, ItemVcblr.size()]; bumped by add_to_global_pool() and trimmed
  * by remove_from_global_pool(). */

 Index f_max_name{ 0 };

  /** Out-Of-Base counters: if OOBase[ i ]
   * = Inf< SIndex >() then there is no item in position i of the bundle = k > 0
   * means that the item in position i is out of base since k
   *   iterations
   * = 0 means in the current base but potentially removable = a *finite*
   * negative value - k means not removable for the next k
   *   iterations: note that some items in base may be such
   * = - Inf< SIndex >() means unremovable */

 Vec_SIndex OOBase;

 bool TrueLB;         /**< true if LowerBound is a "true" lower bound rather
                       * than a "conditional" one */
 bool SSDone;         ///< true if the last step was a SS

 Index f_wFi;         ///< which component was evaluated last

 Subset FiStatus;     ///< status of last computation of each component

 std::vector< C05Function * > v_c05f;
 ///< the vector of (pointers to) the components of the sum function

 OFValue constant_value{};
 ///< the summation of the constant terms of all "easy" components
 ///< collected outside the master MP, which cannot carry them on its own

 LinearFunction * f_lf;  ///< the 0-th component of the sum function

 bool f_convex;          ///< true if all objectives are convex

 MasterProblemBlock * MasterPB;  ///< the Master Problem Block

 std::vector< ColVariable * > LamVcblr;  ///< map Lambda -> ColVariable

 // per-component local→global Lambda index map for sparse Lambda mode.
 // v_local2global[ h ] has size get_num_active_var() + 1 for v_c05f[ h ];
 // entries [ 0 .. loc_NV - 1 ] are the indices in LamVcblr of h's active
 // Variables in the order get_linearization_coefficients writes them, and
 // the last slot is Inf< Index >() so v_local2global[ h ].data() is a
 // ready-to-use Inf-terminated SGBse for the gather logic. Empty
 // (size 0) when f_sparse_lambda is false (dense path).
 //
 // Owned by BundleSolver (not by the C05Function), in line
 // with the policy that Solver-specific bookkeeping never leaks into
 // the Function interface.
 std::vector< std::vector< Index > > v_local2global;

 // true iff at least one v_c05f[ h ] (or f_lf) exposes a strict subset
 // of LamVcblr as its active variables (or the same set in a different
 // order). Auto-detected in two places: at set_Block, by comparing per-
 // component active sets against the union LamVcblr; and at runtime in
 // process_outstanding_Modification's 4th loop, where a naked
 // FunctionModVars* (i.e. one that did NOT arrive as a lockstep
 // GroupModification covering all components) promotes a dense Solver
 // to sparse on the spot — materialising identity local-to-global maps
 // in v_local2global[ * ] and rebuilding Lambda2Idx / v_ref_count from
 // the dense invariant — before the sparse handlers below process the
 // Mod. When false, every gather site falls back to the dense
 // "all components see the same Lambda" code path.
 bool f_sparse_lambda = false;

 // pointer → global Lambda index map, the inverse of LamVcblr. Kept live
 // (and incrementally maintained) only when f_sparse_lambda == true;
 // used to resolve ColVariable * coming from a FunctionModVars* against
 // the global Lambda index space when v_c05f[ h ] is sparse (i.e. when
 // the dense invariant "all components have identical active vars in
 // the same order" does not hold and the Mod's first()/range()/subset()
 // entries cannot be interpreted globally). In dense mode this map is
 // cleared after set_Block to save memory.
 std::unordered_map< ColVariable * , Index > Lambda2Idx;

 // per-LamVcblr-slot refcount: v_ref_count[ i ] is the number of
 // v_c05f (+ f_lf if any) that have LamVcblr[ i ] as an active
 // variable. Built in set_Block alongside LamVcblr; decremented by the
 // sparse FunctionModVarsRngd / FunctionModVarsSbst handlers, and when
 // it reaches 0 the slot is queued for global removal — at the end of
 // the 4th Modification loop we compact LamVcblr / Lambda / Lambda2Idx
 // / v_local2global[ * ] and call MasterPB->remove_vars to reclaim the
 // master row. Kept live only when f_sparse_lambda == true.
 std::vector< Index > v_ref_count;

 VarValue UpTrgt;        ///< upper target
 VarValue LwTrgt;        ///< lower target

 VarValue UpFiBest;      ///< Fi best value vector

 Vec_VarValue UpRifFi;   /** The value of Fi[ k ]() where the zero of the
                          * translated Cutting Plane models are fixed */
 bool RifeqFi;           ///< true if UpRifFi == UpFiLmb
 bool f_linear_part_set = false;
 ///< guard that allows the constant linear gradient of f_lf to be handed
 /// to MasterPB exactly once during a solver instance lifetime: f_lf is
 /// the affine "0-th" component of the original sum-function and its
 /// LinearFunction coefficients are problem data, hence invariant across
 /// successive calls to compute()
 bool CmptdinL;          ///< true if all components are "current"
                         /**< true if the last point in which all non-easy
                          * components have been compute()-d is the current
                          * point Lambda. */

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

 double DST;             ///< D_t( z* ), used to compute the crucial Delta*
 double NrmD;            ///< Euclidean norm of the current direction d*
 double NrmZ;            ///< some norm of the aggregated subgradient z*
 double NrmZFctr;        ///< scaling factor to declare NrmZ "small"

 // fields for events - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::chrono::time_point< std::chrono::system_clock > c_start;
 ///< starting instant of last call to compute()

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void CreateMPB( void );

/*--------------------------------------------------------------------------*/

 static std::pair< double , double >
 effective_bounds( const ColVariable * var );

/*--------------------------------------------------------------------------*/

 void InitMPB( void );

/*--------------------------------------------------------------------------*/

 Index BStrategy( Index wFi );

/*--------------------------------------------------------------------------*/

 Index FindAPlace( Index wFi );

/*--------------------------------------------------------------------------*/

 bool NeedsAlfa1( void );

 bool NeedsScPr1( void );

 bool NeedsG1( void );

/*--------------------------------------------------------------------------*/

 void UpdateHeuristicInfo( void );

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

 void Delete( Index i , bool ModDelete = false );

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

 void add_to_global_pool( Index k , Index i , Index wh = Inf< Index >() );

 void add_to_bundle( Index k , Index i );

 void reset_bundle( void );

/*--------------------------------------------------------------------------*/

 Lst_sp_Mod::size_type num_outstanding_Modification( void );

 bool is_special_GroupMod( GroupModification & gmod );

 void flatten_Modification_list( Lst_sp_Mod & vmt , sp_Mod mod );

 void flatten_easy_Modification_list( Lst_sp_Mod & vmt , sp_Mod mod );

 void process_outstanding_easy_Modification( void );

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

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class BundleSolver )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS BundleSolverState ---------------------------*/
/*--------------------------------------------------------------------------*/
/// class to describe the "internal state" of a BundleSolver
/** Derived class from State to describe the "internal state" of a
 * BundleSolver: the current stability centre, the proximal
 * parameters, and the global pool of all non-easy components. In a better world
 * the State should comprise the State of the Solver used to solve the Master
 * Problem; this will hopefully happen one day. */

class BundleSolverState : public State {

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*----------------------------- PUBLIC TYPES -------------------------------*/

 using Index = BundleSolver::Index;
 using VarValue = BundleSolver::VarValue;
 using Vec_VarValue = BundleSolver::Vec_VarValue;

/*------- CONSTRUCTING AND DESTRUCTING BundleSolverState --------*/

 /// constructor, doing everything or nothing.
 /** Constructor of BundleSolverState. If provided with a pointer
  * to a BundleSolver it immediately copies its "internal state",
  * which is the only way in which the BundleSolverState can be
  * initialised out of an existing BundleSolverState. If nullptr is
  * passed (as by default), then an "empty" BundleSolverState is
  * constructed that can only be filled by calling deserialize(). */

 BundleSolverState( const BundleSolver * bs = nullptr )
   : State() {
  if( ! bs ) {
   NrFi = NumVar = 0;
   return;
   }
  t = bs->t;
  NrFi = bs->NrFi;
  NumVar = bs->NumVar;
  Lambda = bs->Lambda;
  UpFiLmbdef = bs->UpFiLmbdef;
  if( UpFiLmbdef )
   UpFiLmb = bs->UpFiLmb;
  LwFiLmbdef = bs->LwFiLmbdef;
  if( LwFiLmbdef )
   LwFiLmb = bs->LwFiLmb;
  Fi0Lmb = bs->Fi0Lmb;
  global_LB = bs->f_global_LB;
  v_comp_State.resize( NrFi , nullptr );
  for( BundleSolver::Index i = 0 ; i < NrFi ; ++i ) {
   if( bs->NrEasy && bs->IsEasy[ i ] )
    continue;
   v_comp_State[ i ] = bs->v_c05f[ i ]->get_State();
   }
  }

/*--------------------------------------------------------------------------*/
 /// de-serialize a BundleSolverState out of netCDF::NcGroup
 /** De-serialize a BundleSolverState out of netCDF::NcGroup; see
  * BundleSolverState::serialize() for a description of the format.
  */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// destructor

 virtual ~BundleSolverState() {
  for( auto el : v_comp_State )
   delete el;
  }

/*---- METHODS DESCRIBING THE BEHAVIOR OF A BundleSolverState ---*/

 /// serialize a BundleSolverState into a netCDF::NcGroup
 /** The method serializes the BundleSolverState into the provided
  * netCDF::NcGroup, so that it can later be read back by deserialize().
  *
  * After the call, \p group will contain:
  *
  * - The dimension "BundleSolver_NumVar" containing the (current) number
  *   of variables of all the components.
  *
  * - The variable "BundleSolver_Lambda", of type netCDF::NcDouble and
  *   indexed over the dimension BundleSolver_NumVar, which contains
  *   the current stability centre.
  *
  * - The scalar variable "BundleSolver_t", of type netCDF::NcDouble,
  *   which contains the current value of the (tremendous) t parameter.
  *
  * - The dimension "BundleSolver_NrFi" containing the total number of
  *   components (comprised the "easy" ones).
  *
  * - The dimension "BundleSolver_UpFiLmbdef" containing the number of
  *   components for which a non-+INF upper bound value in Lambda is known.
  *   The dimension is optional, if it is not present then 0 is assumed.
  *
  * - The dimension "BundleSolver_LwFiLmbdef" containing the number of
  *   components for which a non--INF lower bound value in Lambda is known.
  *   The dimension is optional, if it is not present then 0 is assumed.
  *
  * - The variable "BundleSolver_UpFiLmb", of type netCDF::NcDouble and
  *   indexed over the dimension BundleSolver_NrFi, which contains the
  *   known opper bound value on each component in Lambda. The variable
  *   must be there if "BundleSolver_UpFiLmbdef" (is defined and) has a
  *   value > 0, but it is ignored (and therefore is optional) otherwise.
  *
  * - The variable "BundleSolver_LwFiLmb", of type netCDF::NcDouble and
  *   indexed over the dimension BundleSolver_NrFi, which contains the
  *   lower opper bound value on each component in Lambda. The variable
  *   must be there if "BundleSolver_LwFiLmbdef" (is defined and) has a
  *   value > 0, but it is ignored (and therefore is optional) otherwise.
  *
  * - The variable "BundleSolver_Fi0Lmb", of type netCDF::NcDouble and
  *   indexed over no dimension (scalar), which contains the value of the
  *   linear part of the objective. The variable is optional, if it is not
  *   present then 0 (typically "no linear part") is assumed.
  *
  * - The variable "BundleSolver_global_LB", of type netCDF::NcDouble and
  *   indexed over no dimension (scalar), which contains the value of the
  *   algorithmic global LB. The variable is optional, if it is not present
  *   then - INF is assumed.
  *
  * - At most BundleSolver_NrFi netCDF::NcGroup with name
  *   "Component_State_X", with X an integer between 0 and
  *   BundleSolver_NrFi - 1, each one containing the State af the
  *   corresponding non-easy component. If a group is not there for some
  *   X then no State has been saved for that component, which surely
  *   happens for "easy" ones. */

 void serialize( netCDF::NcGroup & group ) const override;

/*-------------------------------- FRIENDS ---------------------------------*/

 friend class BundleSolver;  // make BundleSolver friend

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "BundleSolverState [" << this << "] with NrFi = " << NrFi
         << " and NumVar = " << NumVar;
  }

/*--------------------------- PROTECTED FIELDS -----------------------------*/

 Index NrFi;             ///< total number of components

 Index NumVar;           ///< number of variables

 Vec_VarValue Lambda;    ///< the current point

 VarValue t;             ///< the proximity parameter

 Index UpFiLmbdef;       ///< number of known upper bounds

 Index LwFiLmbdef;       ///< number of known lower bounds

 Vec_VarValue UpFiLmb;   ///< the upper bounds

 Vec_VarValue LwFiLmb;   ///< the lower bounds

 VarValue Fi0Lmb;        ///< value of the linear component

 VarValue global_LB;     ///< the global lower bound

 std::vector< State * > v_comp_State;  ///< the States of each component

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( BundleSolverState ) )

/** @} end( group( BundleSolver_CLASSES ) ) --------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BundleSolver.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File BundleSolver.h --------------------------*/
/*--------------------------------------------------------------------------*/
