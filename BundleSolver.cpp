/*--------------------------------------------------------------------------*/
/*------------------------ File BundleSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BunldeSolver class.
 *
 * \version 0.01
 *
 * \date 28 - 11 - 2019
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BundleSolver.h"
#include "LagBFunction.h"
#include "QPPnltMP.h"
#include "OSIMPSolver.h"
#include "ilcplex/cplex.h"
#include "OsiCpxSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define BLOG( l , x ) if( f_log && ( LogVerb > l ) ) *f_log << x

#define BLOG2( l , c , x ) if( f_log && ( LogVerb > l ) && c ) *f_log << x

#define BLOGb( l , x ) if( f_log && ( LogVerb & l ) ) *f_log << x

#define BLOG2b( l , c , x ) if( f_log && ( LogVerb & l ) && c ) *f_log << x

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BundleSolver to the Solver factory
SMSpp_insert_in_factory_cpp_0( BundleSolver );

/*--------------------------------------------------------------------------*/
// define and initialize here the vector of int parameters names
const std::vector< std::string > BundleSolver::int_pars_str = {
 "intBPar1" ,
 "intBPar2" ,
 "intBPar3" ,
 "intBPar4" ,
 "intBPar6" ,
 "intMnSSC" ,
 "intMnNSC" ,
 "inttSPar1" ,
 "intMaxNrEvls" ,
 "intMPName" ,
 "intMPlvl" ,
 "intQPmp1" ,
 "intQPmp2",
 "OSImp1" ,
 "OSImp2" ,
 "OSImp3" ,
 "intRstAlg"
 };

// define and initialize here the vector of double parameters names
const std::vector< std::string > BundleSolver::dbl_pars_str = {
 "dbltStar" ,
 "dblRelMPAcc" ,
 "dblRMPAccSol" ,
 "dblBPar5" ,
 "dblm1" ,
 "dblm2" ,
 "dblm3" ,
 "dblmxIncr" ,
 "dblmnIncr" ,
 "dblmxDecr" ,
 "dblmnDecr" ,
 "dbltMaior" ,
 "dbltMinor" ,
 "dbltInit" ,
 "dbltSPar2" ,
 "dblCtOff"
 };

// define and initialize here the map for int parameters names
const std::map< std::string , BundleSolver::idx_type >
 BundleSolver::int_pars_map = {
 { "intBPar1" , BundleSolver::intBPar1  } ,
 { "intBPar2" , BundleSolver::intBPar2  } ,
 { "intBPar3" , BundleSolver::intBPar3 } ,
 { "intBPar4" , BundleSolver::intBPar4 } ,
 { "intBPar6" , BundleSolver::intBPar6 } ,
 { "intMnSSC" , BundleSolver::intMnSSC } ,
 { "intMnNSC" , BundleSolver::intMnNSC } ,
 { "inttSPar1" , BundleSolver::inttSPar1 } ,
 { "intMaxNrEvls" , BundleSolver::intMaxNrEvls } ,
 { "intMPName" , BundleSolver::intMPName } ,
 { "intMPlvl" , BundleSolver::intMPlvl } ,
 { "intQPmp1" , BundleSolver::intQPmp1 } ,
 { "intQPmp2" , BundleSolver::intQPmp2 } ,
 { "intOSImp1" , BundleSolver::intOSImp1 } ,
 { "intOSImp2" , BundleSolver::intOSImp2 } ,
 { "intOSImp3" , BundleSolver::intOSImp3 } ,
 { "intRstAlg" , BundleSolver::intRstAlg } ,
 };

// define and initialize here the map for double parameters names
const std::map< std::string , BundleSolver::idx_type >
 BundleSolver::dbl_pars_map = {
 { "dbltStar" , BundleSolver::dbltStar } ,
 { "dblRelMPAcc" , BundleSolver::dblRelMPAcc } ,
 { "dblRMPAccSol" , BundleSolver::dblRMPAccSol } ,
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
 { "dblCtOff" , BundleSolver::dblCtOff }
 };

// define and initialize here the default int parameters
const std::vector< int > BundleSolver::dflt_int_par = {
 10 ,  // intBPar1
100 ,  // intBPar2
  1 ,  // intBPar3
  1 ,  // intBPar4
  0 ,  // intBPar6
  0 ,  // intMnSSC
  3 ,  // intMnNSC
 12 ,  // inttSPar1
  2 ,  // intMaxNrEvls
  1 ,  // intMPName
  0 ,  // intMPlvl
  0 ,  // intQPmp1
  0 ,  // intQPmp2
  4 ,  // intOSImp1
  0 ,  // intOSImp2
  1 ,  // intOSImp3
  62  // intRstAlg
 };

// define and initialize here the default double parameters
const std::vector<double> BundleSolver::dflt_dbl_par = {
 1e+2 ,   // dbltStar
 1e-8 ,   // dblRelMPAcc
 1e-8 ,   // dblRMPAccSol
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
 1e-1     // dblCtOff
 };

/*--------------------------------------------------------------------------*/

static const HpNum Nearly  = 1.01;
static const HpNum Nearly2 = 1.02;

static const char LogBnd = 16;        // log Bundle changes
static const char LogVar = 32;        // log variables changes

static cIndex tSP1Msk = ~ 3;          // mask for tSPar1
static cIndex kSLTTS =  4;            // "soft" long-term t-strategy
static cIndex kHLTTS =  8;            // "hard" long-term t-strategy
static cIndex kBLTTS = 12;            // "balancing" long-term t-strategy
static cIndex kEGTTS = 16;            // "endgame" long-term t-strategy

static const unsigned char RstAlg =  1;  // don't reset algorithmic parameters
static const unsigned char RstCrr =  2;  // don't reset current point
static const unsigned char RstSbg =  4;  // don't reset subgradients
static const unsigned char RstCnt =  8;  // don't reset constraints
static const unsigned char RstFiV = 16;  // don't reset FiVals
static const unsigned char NoStPt = 32;  // don't get an initial point

static cIndex InINF = SMSpp_di_unipi_it::Inf<Index>();

/*--------------------------------------------------------------------------*/

int BundleSolver::compute( bool changedvars )
{

 if( MaxIter == 0 )  // No iteration must be performed
  return kStopIter;

 // basic sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! Master )
  throw( std::logic_error( "Master is not set yet" ) );

 if( v_c05f.empty() )
  throw( std::logic_error( "C05Function is not set yet" ) );

 // first, process any outstanding Modification - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 process_outstanding_Modification();

 // initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ParIter = 0;
 Result = kStopIter;
 SCalls++;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 do {
  // construct the direction d- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormD();

  if( Result >= kError ) {  // problems in the Master Problem solver
   BLOG( 1 , " ~ error in the MPSolver" << std::endl );
   break;
   }

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( IsOptimal() ) {
   Result = kOK;
   break;
  }

  // check if "ex-ante" Noise Reduction is needed - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ensure that the \sigma* is "not too negative", if it is increase t (if
  // possible) and re-solve the MP; note that this kind of NR only happens if
  // the oracle is "unfaithful", i.e., it pretends to provide information with
  // the required accuracy but in fact it does not

  if( ( Sigma < 0 ) &&  // do not even call ReadDStart() if Sigma >= 0
      ( Sigma <= - m3 * Master->ReadDStart( t ) ) ) {
   if( t >= tMaior ) {
    BLOG( 1 , " ~ NR required but t maximum" << std::endl );
    Result = kError;
    break;
    }

   t = std::min( t * mxIncr , tMaior );
   BLOG( 2 , " ~ NR: t increased to " << t << std::endl );
   tHasChgd = true;
   continue;
   }

  // update out-of-base counters- - - - - - - - - - - - - - - - - - - - - - -

  UpdtCntrs();

  // some log - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log1();

  // Hard Long-Term t-strategy for quadratic stabilization- - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the hard long-term t-strategy requires t to increase if the step is too
  // small, and therefore has to be checked before the others
  // however, it is only viable under a quadratic stabilization

  if( ( ( ! MPName ) || ( MPName & 4 ) ) &&
      ( tStar > 0 ) && ( ( tSPar1 & tSP1Msk ) == kHLTTS ) &&
      ( UpFiLmb[ NrFi ] < Inf<double>() ) ) {

   double AFL = std::abs( UpFiLmb[ NrFi ] );
   if( AFL < 1 )
    AFL = 1;

   if( abs(vStar[ NrFi ]) <= tSPar2 * EpsU * AFL ) {
    BLOG( 1 , "small v => increase t" << std::endl << "           " );

    // collect two numbers vc and vl such that v( tNew ) >= vc + tNew * vl
    // we require that v( tNew ) >= vc + tNew * vl = tSPar2 * EpsU * AFL
    // ==> tNew = ( tSPar2 * EpsU * AFL - vc ) / vl

    double vl , vc;
    Master->SensitAnals( vl , vc );

    double tt;
    if( - vl < Eps<HpNum>() )  // v( t ) is [almost] constant ==> D*_t [~]= 0
     tt = tStar;               // ==> the CP model is ~bounded
    else
     tt = std::min( tStar , ( tSPar2 * EpsU * AFL * Nearly + vc ) /
		            ( - vl ) );

    if( ( tHasChgd = ( tt != t ) ) ) {
     t = tt;
     continue;         // loop only if t changes
     }
    }
   }  // end if( Hard t-strategy )  - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // a real iteration (iterations where Fi() is not evaluated do not count) -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ParIter++;

  // calculate Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  // calculate Fi( Lambda1 )- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index wFi = 0;
  CurrNrEvls.assign( NrFi , Index( 0 ) );

  bool MPchgs = false;  // true if no cycling will occur
  for( ; ; ) {   // ... possibly more than once due to precision issues

   MPchgs = FiAndGi( wFi );
   CurrNrEvls[ wFi ]++;

   if( MPchgs )  // if something changes
    break;       // all done

   if( ! FindNext( wFi ) )  // find next component
    break;                  // if none, end
   }

  if( ! MPchgs ) {  // noise reduction
   t = std::min( t * mxIncr , tMaior );
   BLOG( 1 , " ~ noise reduction: t increased to " << t << std::endl );
   tHasChgd = true;
   if( t >= tMaior )
    Result = kError;
   }

  // check whether the Lower Bounds have changed- - - - - - - - - - - - - - -

  UpdtLowerBound();

  // some log about the newly obtained information- - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log2();

  // check whether either any error has occurred or time has expired- - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1[NrFi] == - Inf<double>() ) {
   Result = kUnbounded;
   break;
   }

  if( ( Result == kError ) || ( Result == kStopTime ) )
   break;

  if( tHasChgd )  // "noise reduction": t has changed,
   continue;                      // so go solve the master problem again
                                  // (no NS/SS decision can be made)

  // check for the conditional lower bound- - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( ! TrueLB ) &&
      ( UpFiBest <= LowerBound *
	            ( 1 - ( LowerBound > 0 ? RelAcc : - RelAcc ) ) ) ) {
   Result = kUnbounded;
   break;
   }

  // avoid the t-changing phase if Lambda1 is unfeasible- - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // note: one possible alternative t-strategy would be to set t to the
  // largest value that would have produced a feasible point, i.e.
  // t := *Alfa1 / ( - *ScPr1 )

  if( UpFiLmb1[NrFi] == Inf<double>() )  // ???
   continue;
  else
   if( UpFiLmb[NrFi] == Inf<double>() ) {  // if reached feasibility  - - - - -
	GotoLambda1();             // go to the feasible point
    continue;                  // and start the actual minimization of Fi()
    }

  // the NS / SS decision - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SSDone = ( UpFiLmb1[ NrFi ] <= UpTrgt ) ? true : false;

  // compute the heuristic t- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  HpNum tt;
  if( ( SSDone && ( ! ( tSPar1 & 1 ) ) ) ||
       ( ( ! SSDone ) && ( tSPar1 & 2 ) ) )
   tt = Heuristic1();
  else
   tt = Heuristic2();

  if( SSDone ) {  // SS - - - - - - - - - - - - - - - - - - - - - - - - - - -

   BLOG( 1 , std::endl << " SS[" << CSSCntr << "]: DFi = " << DeltaFi
	     <<  " ~ Up1(" << UpFiLmb1[ NrFi ] << ") <= UpTrgt(" << UpTrgt
	     << ") ~ Ht = " << tt );

   tt = std::min( std::min( tMaior , t * mxIncr ) ,
 		  std::max( t * mnIncr , tt ) );

   if( CSSCntr < MnSSC )  // increasing t is inhibited
    tt = t;
   else
    if( ( tSPar1 & tSP1Msk ) == kBLTTS )  // "balancing" long-term t-strategy
     if( ( DSTS <= tSPar2 * Sigma ) && ( CSSCntr < 10 ) ) {  //!! 10!
      BLOG( 1 , " ~ small D*_t( 1 )" );
      tt = t;
      }

   BLOG( 1 , std::endl );

   GotoLambda1();
   CSSCntr++;
   CNSCntr = 0;
   }
  else {        // NS - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   BLOG( 1 , std::endl << " NS[" << CNSCntr << "]: DFi = " << DeltaFi
	     <<  " ~ Lw1(" << LwFiLmb1[ NrFi ] << ") >= LwTrgt(" << LwTrgt
	     << ") ~ Ht = " << tt );

   tt = std::max( std::max( tMinor , t * mxDecr ) ,
 		  std::min( t * mnDecr , tt ) );

   if( CNSCntr < MnNSC )  // decreasing t is inhibited
    tt = t;
   else
    if( Alfa1[ NrFi ] <= m1 * Sigma ) {
     BLOG( 1 , " ~ small Alfa1" );
     tt = t;
     }
    else
     switch( tSPar1 & tSP1Msk ) {
      case( kSLTTS ):
      case( kHLTTS ):
       if( abs( vStar[ NrFi ] ) <= tSPar2 * EpsU * max_error() ) {
        BLOG( 1 , " small v" );
	tt = t;
        }
       break;
      case( kBLTTS ):
       /*!! this version avoids problems which may occur with ill-set
 	    tStar or tSPar2, but it may give worse performances with
 	    "difficult" problems
        if( ( tSPar2 * DSTS >= Sigma ) && ( CNSCntr < 20 ) ) {
        !!*/

       if( tSPar2 * DSTS >= Sigma ) {
        BLOG( 1 , " ~ large D*_t( t* )" );
        tt = t;
        }
      }

   BLOG( 1 , std::endl );

   CNSCntr++;
   CSSCntr = 0;

   }   // end else( NS )- - - - - - - - - - - - - - - - - - - - - - - - - - -

  // actually update t- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( tSPar1 & kEGTTS ) && ( UpFiLmb[ NrFi ] < Inf< VarValue >() ) )
   // endgame t-strategy: note the "/ 10"!!
   if( DSTS < max_error() / 10 ) {
    tt = std::max( t * ( mxDecr + mnDecr ) / 2 , tMinor );
    BLOG( 1 , " ~ endgame, t = " << tt );
    }

  //!! the reverse should also be done: if sigma is small and D*( t* ) is
  //!! large, t should be increased --> but this would happen surely at
  //!! the beginning, it should be done only near the end

  if( ( tHasChgd = ( t != tt ) ) ) {
   CSSCntr = CNSCntr = 0;  // reset the counters as t changes
   t = tt;
   }

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  } while( ParIter < MaxIter );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

return( Result );

}  // end( BundleSolver::compute() ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::set_Block( Block * block )
{
 if( f_Block ) {  // changing from a previous oracle - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  guts_of_destructor();   // deallocate memory
  }

 Solver::set_Block( block );  // attach to the new Block

 if( ! f_Block )  // that was actually clearing the Block
  return;         // all done

 /* Two types of block can be handled by the BundleSolver:

     1. Only one single non-smooth function
     2. a sum of some non-smooth functions

    The algorithm here developed aims at solving non-constrained non-smooth
    optimization. The block can have box constraints at the most.
    It is expected the block to have in the first case a FRealObjective
    whose the function is a C05Function one and having no children, while in
    the second case a FRealObjective one whose the function is a LinearFunction
    and having as many sub-blocks as the number of components. In the latter case,
    each sub-block must not contain any Variable or Constraint.
    Variable may have a lower and upper bound. If the lower bound  has a finite
    value, it must be 0. */

 if( f_Block->get_nested_Blocks().empty() ) {

  // the objective function of the block must be a C05Function  - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
  if( obj == nullptr )
   throw( std::logic_error( "objective is not a FRealObjective" ) );

  auto c05f = dynamic_cast< C05Function * >( (obj)->get_function() );
  if( c05f == nullptr )
   throw( std::logic_error( "the objective is not a C05Function" ) );

  v_c05f.push_back( c05f );
  f_lf = nullptr;

  }
 else {

  // the objective function of each block must be a LinearFunction - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( f_Block->get_objective() )
   f_lf =  nullptr;
  else {
   auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
   if( obj == nullptr )
    throw( std::logic_error( "the objective is not a real function" ) );

   f_lf = dynamic_cast<LinearFunction *>( (obj)->get_function() );
   if( f_lf == nullptr )
    throw( std::logic_error( "the objective is not a LinearFunction" ) );
   }


  auto sb = f_Block->get_nested_Blocks();
  v_c05f.resize( sb.size() );

  for( Index i = 0 ; i < sb.size() ; ++i ) { // for each sub-block

   // the objective function of each sub-block must be a C05Function - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto obj = dynamic_cast< FRealObjective * >( sb[ i ]->get_objective() );
   if( obj == nullptr )
    throw( std::logic_error( "the objective is not a real function" ) );

   auto c05f = dynamic_cast<C05Function *>( (obj)->get_function() );
   if( c05f == nullptr )
	throw( std::logic_error( "the objective is not a C05Function" ) );
   v_c05f[ i ] = c05f;

   // nephew are not allowed - - - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( sb[ i ]->get_nested_Blocks().size() )
	throw( std::logic_error( "nephew are not allowed" ) );

   // Variable of Sub-Block are not expected, neither the Constraint - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( sb[ i ]->get_static_variables().size() )
	throw( std::logic_error( "static Variable are not allowed" ) );

   if( sb[ i ]->get_dynamic_variables().size() )
	throw( std::logic_error( "dynamic Variable are not allowed" ) );

   if( sb[ i ]->get_static_constraints().size() )
	throw( std::logic_error( "static Constraint are not allowed" ) );

   if( sb[ i ]->get_dynamic_constraints().size() )
	throw( std::logic_error( "dynamic Constraint are not allowed" ) );

   }
  } // end decomposed case - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if some Variable are present, they are of the ColVariable type - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar = 0;  // count the number of Variable
 auto v_s_Variable = f_Block->get_static_variables();
 for( auto & el : v_s_Variable ) {
  if( un_any_thing_0( ColVariable , el , ++NumVar ) )
   continue;
  if( un_any_thing_1( ColVariable , el , NumVar += var.size() ) )
   continue;
  if( un_any_thing_K( ColVariable , el , NumVar += var.num_elements() ) )
   continue;
  throw( std::logic_error( "some static Variable is not a ColVariable" ) );
  }

 // construct the vocabulary for Variable and sort it  - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LamVcblr.resize( NumVar );
 Index count = 0;
 for( auto & el : v_s_Variable )
  un_any_static( el , [ & ]( ColVariable & static_var ) {
                       LamVcblr[ count++ ] = & static_var;
                       } ,
		  un_any_type<ColVariable>() );

 assert( count == NumVar );
 std::sort( LamVcblr.begin() , LamVcblr.end() );

 // no dynaimic variables are allowed  - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_Block->get_dynamic_variables().size() )
  throw( std::logic_error( "dynamic Variable are not allowed" ) );

 // if some Constraint are present, their functions must be linear
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_Block->get_static_constraints().size() )
	throw( std::logic_error( "static Constraint are not allowed" ) );

 if( f_Block->get_dynamic_constraints().size() )
	throw( std::logic_error( "dynamic Constraint are not allowed" ) );

 // read information about the function  - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NrFi = v_c05f.size();

 MILP_s.resize( NrFi , nullptr );
 vStar.resize( NrFi + 1 , 0 );

 NrEasy = 0;
 if( NrFi > 1 ) {
  auto sb = f_Block->get_nested_Blocks();
  std::vector<Index> BNC( NrFi );

  IsEasy.resize( NrFi , false );
  for( Index k = 0 ; k < NrFi ; ++k ) {
   // ?? quando i solver distruggere dopo SetDim o
   // con il distruttore della classe??

   auto LagB = dynamic_cast<LagBFunction *>( v_c05f[ k ] );
   if( LagB  ) {
    MILP_s[ k ] = new MILPSolver();
    MILP_s[ k ]->set_Block( LagB->get_inner_block() );
    BNC[ k ] = MILP_s[ k ]->get_numcols();
    if( BNC[ k ] ) {
     IsEasy[ k ] = true;
     NrEasy++;
     }
    }
   }

  if( ! NrEasy )
   IsEasy.clear();
  }
 else
  MILP_s[ 0 ] = nullptr;

 // set the global pool size to all non-easy functions - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NrEasy )
  for( Index i = 0 ; i < NrFi ; ++i ) {
   if( IsEasy[ i ] )
    continue;
   auto gps = v_c05f[ i ]->get_int_par( C05Function::intGPMaxSz );
   if( gps < BPar2 )
    v_c05f[ i ]->set_par( C05Function::intGPMaxSz , BPar2 );
   }
 else
  for( auto fun : v_c05f ) {
   auto gps = fun->get_int_par( C05Function::intGPMaxSz );
   if( gps < BPar2 )
    fun->set_par( C05Function::intGPMaxSz , BPar2 );
   }

 // allocate memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 t = tInit;
 Prevt = Inf<double>();

 Lambda.resize( NumVar );    // the default starting point
 Lambda1.resize( NumVar );   // the tentative point

 if( MaxSol > 1 )  // best point found so far
  LmbdBst.resize( NumVar );

 OOBase.resize( BPar2 * ( NrFi - NrEasy ) , Inf<SIndex>() );
 // counter for eliminating outdated items: Inf<SIndex>() means empty

 FreList = priority_queue<Index>();     // list of free bundle slots
 whisZ.resize( NrFi );  // for each component, the name of its "Z" if it is
                        // in the bunlde

 FiStatus.resize( NrFi , kUnEval );
 LowerBound = -Inf<double>(); // lower bounds
 TrueLB = false;

 UpFiBest = Inf<VarValue>();     // best, ...
 UpRifFi.resize( NrFi + 1 , 0 ); // and reference Fi() values
 UpFiLmb1.resize( NrFi + 1 );    // upper and lower function value
 LwFiLmb1.resize( NrFi + 1 );    // ... at the tentative point

 UpFiLmb.resize( NrFi + 1 , Inf<VarValue>() );  // upper and lower function value
 LwFiLmb.resize( NrFi + 1 , -Inf<VarValue>() ); // ... at the current point

 whisG1.resize( NrFi , Inf<Index>() );  // no representative yet

 ScPr1.resize( NrFi + 1 , 0 );
 Alfa1.resize( NrFi + 1 , 0 );
 DeltaAlfa.resize( NrFi );

 Result = kError;
 SSDone = false;

 // warning: the following things can only be done *after* that
 // Oracle->SetMaxName() has been invoked, because they use methods of the
 // oracle which depends on knowledge of the MaxName to work properly
 // read b0- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // here one could initialize b0, if that was found to be of any use
 // b0 = Oracle->GetVal( BPar2 );

 // initialize the MP Solver, if any - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Master )        // a MP solver is set ??? dove metterlo ???
  Master->SetDim();  // clear all its internal state

 if( ! MPName ) {
  Master = new QPPenaltyMP( );
  QPPenaltyMP *qp = dynamic_cast<QPPenaltyMP*>( Master );
  qp->SetPricing( CtOff );
  qp->SetMaxVarAdd( MxAdd );
  qp->SetMaxVarRmv( MxRmv );
  }
 else {
  OSIMPSolver * osi_mps = new OSIMPSolver();
  Master = osi_mps;

  if( MPName & 2 ) {
   OsiCpxSolverInterface *osicpx = new OsiCpxSolverInterface();
   CPXENVptr env = osicpx->getEnvironmentPtr ();
   CPXsetintparam( env , CPX_PARAM_THREADS , threads );
   // 12.8
   // CPXsetlogfile( env , NULL );
   // 12.9
   // CPXsetlogfilename( env, "/dev/null" , "w" ) ;

   CPXsetintparam( env , CPXPARAM_ScreenOutput , CPX_OFF );
   CPXsetintparam( env , CPXPARAM_Barrier_Display , 0 );
   CPXsetintparam( env , CPXPARAM_Simplex_Display , 0 );
   CPXsetintparam( env , CPXPARAM_Sifting_Display , 0 );
   CPXsetintparam( env , CPXPARAM_Network_Display , 0 );
   CPXsetintparam( env , CPXPARAM_ParamDisplay  , CPX_OFF );
 
   osi_mps->SetOsi( osicpx );
   }
  else
   osi_mps->SetOsi( new OsiClpSolverInterface() );

  osi_mps->SetStabType( MPName & 4 ? OSIMPSolver::quadratic :
			             OSIMPSolver::boxstep );

  osi_mps->SetAlgo( OSIMPSolver::OsiAlg( algo ) ,
		    OSIMPSolver::OsiRed( reduction ) );
  }

 InitMP();

 // reset algorithm  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // default value  RstCrr | RstSbg | RstCnt | RstFiV | NoStPt = 62
 ReSetAlg( RstAlgPrm );  // Fi( Lambda ) is reset inside

 }  // end( BundleSolver::set_Block )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( const idx_type par , const int value ) {

 switch( par ) {
  case( intMaxIter ):
   if( value < 0 )
    throw( std::invalid_argument( "MaxIter must be >= 0" ) );
   MaxIter = value;
   break;
  case( intMaxSol ):
   if( value < 1 )
    throw( std::invalid_argument( "MaxSol must be >= 1" ) );
   MaxSol = value;
   break;
  case( intLogVerb ):
   LogVerb = value;
   break;
  case( intBPar1 ):
   if( value < 1 )
    throw( std::invalid_argument( "BPar1 must be >= 1" ) );
   BPar1 = value;
   break;
  case( intBPar2 ):
   if( value < 2 )
    throw( std::invalid_argument( "BPar2 must be >= 2" ) );
   BPar2 = value;
   break;
  case( intBPar3 ):
   if( value < BPar4 )
    throw( std::invalid_argument( "BPar3 must be >= BPar4" ) );
   BPar3 = value;
   break;
  case( intBPar4 ):
   if( value < 1 )
    throw( std::invalid_argument( "BPar4 must be >= 1" ) );
   BPar4 = value;
   if( BPar4 > BPar3 )
    BPar3 = BPar4;
   break;
  case( intBPar6 ):
   BPar6 = value;
   break;
  case( intMnSSC ):
   MnSSC = value;
   break;
  case( intMnNSC ):
   MnNSC = value;
   break;
  case( inttSPar1 ):
   tSPar1 = value;
   break;
  case( intMaxNrEvls ):
   MaxNrEvls = value;
   break;
  case( intMPName ):
   if( ( value < 0 ) || ( value > 7 ) )
    throw( std::invalid_argument( "MPName must be in [0, 7]" ) );
   MPName = value;
   break;
  case( intMPlvl ):
   MPlvl = value;
   break;
  case( intQPmp1 ):
   MxAdd = value;
   break;
  case( intQPmp2 ):
   MxRmv = value;
   break;
  case( intOSImp1 ):
   algo = value;
   break;
  case( intOSImp2 ):
   reduction = value;
   break;
  case( intOSImp3 ):
   threads = value;
   break;
  case( intRstAlg ):
   RstAlgPrm = value;
   break;
  default:
   CDASolver::set_par( par , value );
  }
 } // end (BundleSolver::set_par( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( const idx_type par , const double value )
{
 switch( par ) {
  case( dblMaxTime ):
   if( value <= 0 )
    throw( std::invalid_argument( "dblMaxTime must be > 0" ) );
   MaxTime = value;
   break;
  case( dblRelAcc ):
   if( ( value <= 0 ) || ( value >= Inf< VarValue >() ) )
    throw( std::invalid_argument( "RelAcc must be > 0 and finite" ) );
   RelAcc = value;
   break;
  case( dblAbsAcc ):
   if( value <= 0 )
    throw( std::invalid_argument( "AbsAcc must be > 0" ) );
   AbsAcc = value;
   break;
  case( dblRAccSol ):
   if( value <= 0 )
    throw( std::invalid_argument( "RAccSol must be > 0" ) );
   RAccSol = value;
   break;
  case( dblAAccSol ):
   if( value <= 0 )
    throw( std::invalid_argument( "AAccSol must be > 0" ) );
   AAccSol = value;
   break;
  case( dbltStar ):
   tStar = value;
   break;
  case( dblRelMPAcc ):
   if( value <= 0 )
    throw( std::invalid_argument( "RelMPAcc must be > 0" ) );
   RelMPAcc = value;
   break;
  case( dblRMPAccSol ):
   if( value <= 0 )
    throw( std::invalid_argument( "RMPAccSol must be > 0" ) );
   RMPAccSol = value;
   break;
  case( dblBPar5 ):
   BPar5 = value;
   break;
  case( dblm1 ):
   if( std::abs( value ) >= 1 )
    throw( std::invalid_argument( "| m1 | must be in (0, 1)" ) );
   m1 = value;
   break;
  case( dblm2 ):
   if( ( value < std::abs( m1 ) ) || ( value >= 1 ) )
    throw( std::invalid_argument( "m2 must be in [ | m1 |, 1)" ) );
   m2 = value;
   break;
  case( dblm3 ):
  if( ( value <= 0 ) || ( value >= 1 ) )
    throw( std::invalid_argument( "m3 must be in (0, 1)" ) );
   m3 = value;
   break;
  case( dblmxIncr ):
  if( value <= 1 )
    throw( std::invalid_argument( "mxIncr must be > 1" ) );
   mxIncr = value;
   break;
  case( dblmnIncr ):
   if( value <= 1 )
    throw( std::invalid_argument( "mnIncr must be > 1" ) );
   mnIncr = std::min( value , mxIncr );
   break;
  case( dblmxDecr ):
   if( ( value <= 0 ) || ( value > 1 ) )
    throw( std::invalid_argument( "mxDecr must be in (0, 1)" ) );
   mxDecr = value;
   break;
  case( dblmnDecr ):
   if( ( value <= 0 ) || ( value > 1 ) )
    throw( std::invalid_argument( "mnDecr must be in (0, 1)" ) );
   mnDecr = std::max( value , mxDecr );
   break;
  case( dbltMaior ):
   if( value <= 0 )
    throw( std::invalid_argument( "tMaior must be > 0" ) );
   tMaior = value;
   break;
  case( dbltMinor ):
   if( ( value <= 0 ) || ( value > tMaior ) )
    throw( std::invalid_argument( "tMinor must be in (0, tMaior]" ) );
   tMinor = value;
   break;
  case( dbltInit ):
   if( ( value < tMinor ) || ( value > tMaior ) )
    throw( std::invalid_argument( "tInit must be in [tMinor, tMaior]" ) );
   tInit = value;
   break;
  case( dbltSPar2 ):
   if( value <= 0 )
    throw( std::invalid_argument( "tSPar2 must be > 0" ) );
   tSPar2 = value;
   break;
  case( dblCtOff ):
   if( value < 0 )
    throw( std::invalid_argument( "CtOff must be >= 0" ) );
   CtOff = value;
   break;
  default:
   CDASolver::set_par( par , value );
  }

 } // end (BundleSolver::set_par( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::set_log( std::ostream * log_stream )
{
 f_log = log_stream;
 Master->SetMPLog( f_log , MPlvl );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::get_dual_solution( Configuration *solc ) {
 for( Index i = 0 ; i < zA.size() ; ++i ) {
  if( zA[i].second.empty() )
   throw( std::invalid_argument( "the combination is not present" ) );
  v_c05f[ i ]->set_important_linearization( std::move(zA[i].second) , zA[i].first );
  }
 } // end ( BundleSolver::get_dual_solution() )  - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

int BundleSolver::get_int_par( const idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):
   return( MaxIter );
   break;
  case( intMaxSol ):
   return( MaxSol );
   break;
  case( intLogVerb ):
   return( LogVerb );
   break;
  case( intBPar1 ):
   return( BPar1 );
   break;
  case( intBPar2 ):
   return( BPar2 );
   break;
  case( intBPar3 ):
   return( BPar3 );
   break;
  case( intBPar4 ):
   return( BPar4 );
   break;
  case( intBPar6 ):
   return( BPar6 );
   break;
  case( intMnSSC ):
   return( MnSSC );
   break;
  case( intMnNSC ):
   return( MnNSC );
   break;
  case( inttSPar1 ):
   return( tSPar1 );
   break;
  case( intMaxNrEvls ):
   return( MaxNrEvls );
   break;
  case( intMPName ):
   return( MPName );
   break;
  case( intMPlvl ):
   return( MPlvl );
   break;
  case( intQPmp1 ):
   return( CtOff );
   break;
  case( intQPmp2 ):
   return( MxRmv );
   break;
  case( intOSImp1 ):
   return( algo );
   break;
  case( intOSImp2  ):
   return( reduction );
   break;
  case( intOSImp3 ):
   return( threads  );
   break;
  case( intRstAlg ):
   return( RstAlgPrm  );
   break;
  default:
   return( CDASolver::get_dflt_int_par( par ) );
  }

 } // end( BundleSolver::get_int_par )  - - - - - - - - - - - - - - - - - - -

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double BundleSolver::get_dbl_par( const idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ):
   return( MaxTime );
   break;
  case( dblRelAcc ):
   return( RelAcc );
   break;
  case( dblAbsAcc ):
   return( AbsAcc );
   break;
  case( dblRAccSol ):
   return( RAccSol );
   break;
  case( dblAAccSol ):
   return( AAccSol );
   break;
  case( dbltStar ):
   return( tStar );
   break;
  case( dblRelMPAcc ):
   return( RelMPAcc );
   break;
  case( dblRMPAccSol ):
   return( RMPAccSol );
   break;
  case( dblBPar5 ):
   return( BPar5 );
   break;
  case( dblm1 ):
   return( m1 );
   break;
  case( dblm2 ):
   return( m2 );
   break;
  case( dblm3 ):
   return( m3 );
   break;
  case( dblmxIncr ):
   return( mxIncr );
   break;
  case( dblmnIncr ):
   return( mnIncr );
   break;
  case( dblmxDecr ):
   return( mxDecr );
   break;
  case( dblmnDecr ):
   return( mnDecr );
   break;
  case( dbltMaior ):
   return( tMaior );
   break;
  case( dbltMinor ):
   return( tMinor );
   break;
  case( dbltInit ):
   return( tInit );
   break;
  case( dbltSPar2 ):
   return( tSPar2 );
   break;
  case( dblCtOff ):
   return( CtOff );
   break;
  default:
   return( CDASolver::get_dflt_dbl_par( par ) );
  }

 } // end( BundleSolver::get_dbl_par ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*----------------------- OTHER PROTECTED METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::FormD( void )
{
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

 if( Master->BCSize() >= Master->BSize() ) {
  if( ( t > tMinor ) && ( Prevt == Inf<double>() ) ) {
   Prevt = t;
   t = tMinor;
   tHasChgd = true;
   }
  }
 else
  if( Prevt < Inf<double>() ) {
   if( t != Prevt ) {
    t = Prevt;
    tHasChgd = true;
    }
   Prevt = Inf<double>();
   }

 if( tHasChgd ) {
  Master->Sett( t );
  tHasChgd = false;
  }

 // collect and set individual and global lower bounds- - - - - - - - - - - -

 if( LBHasChgd && ( UpFiLmb[ NrFi ] < Inf<double>() ) ) {
  if( LowerBound > - Inf<double>() )
   Master->SetLowerBound( LowerBound - UpFiLmb[ NrFi ] );
  else
   Master->SetLowerBound( - Inf<double>() );

  for( Index k = 0 ; k < NrFi ; k++ ) {
   if( NrEasy && IsEasy[ k ] )  // skip easy components
    continue;

   if( LowerBound > - Inf<double>() )
    Master->SetLowerBound( LowerBound - UpFiLmb[ k ] , k + 1 );
   else
    Master->SetLowerBound( - Inf<double>() , k + 1 );
   }

  LBHasChgd = false;
  }

 // set termination criterion - - - - - - - - - - - - - - - - - - - - - - - -

 if( UpFiLmb[ NrFi ] < Inf<double>() )
  Master->SetPar( MPSolver::kZero ,
		  max_error() / std::max( tStar / t , HpNum( 1 ) ) );

 for(;;)  // error-handling loop - - - - - - - - - - - - - - - - - - - - - -
 {        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  MPSolver::MPStatus mps = Master->SolveMP();  // solve the MP

  if( mps == MPSolver::kOK )        // everything's alright
   break;

  if( mps == MPSolver::kUnfsbl ) {  // the feasible set is empty
   Result = kInfeasible;
   break;
   }

  if( mps == MPSolver::kUnbndd ) {  // the MP is unbounded: this can always
                                     // be mended by decreasing t ...
   if( ( t <= tMinor ) || ( Master->BCSize() >= Master->BSize() ) ) {
    // ... but t must always be >= tMinor, and it is already == tMinor in
    // the "empty" case of the initial iteration with empty bundle
    BLOG( 1 , std::endl << "Bundle::FormD: failure in MPSolver." );
    Result = kError;
    break;
    }

   BLOG( 1 , std::endl << "Bundle::FormD: MP unbounded, decreasing t" );
   Master->Sett( t = std::max( t / 2 , tMinor ) );
   continue;
   }

  if( mps == MPSolver::kStppd ) {  // stopped by time limit
   //!! so far, the time limit in the MPSolver is only due to the
   //!! global time limit in the NDOSolver, but one day we may want
   //!! to set it independently; then, some checks will have to be
   //!! done if the solution is feasible and it can still be used
   //!! and v is sufficiently < 0: in this case we can use the
   //!! solution as well, otherwise we have to give the MPSolver
   //!! more time
   Result = kStopTime;
   break;
   }

  // mps == MPSolver::kError, i.e., there has been a numerical problem in- -
  // the MP Solver; it's not yet time to despair, as by eliminating items- -
  // it may be possible to solve the problem - - - - - - - - - - - - - - - -

  Index MBDm;
  cIndex_Set MBse;
  cHpRow Mlt = Master->ReadMult( MBse , MBDm );
  Index i = InINF;

  // the last *removable* item in Base is eliminated - - - - - - - - - - - -

  if( MBse ) {
   for( ; MBDm-- ; )
    if( ( OOBase[ MBse[ MBDm ] ] >= 0 ) &&
        ( Mlt[ MBDm ] >= Eps<HpNum>() ) ) {
     i = MBse[ MBDm ];
     break;
     }
   }
  else
   for( ; MBDm-- ; )
    if( ( OOBase[ MBDm ] >= 0 ) && ( Mlt[ MBDm ] >= Eps<HpNum>() ) ) {
     i = MBDm;
     break;
     }

  if( i == InINF )  // there are no *removable* items in Base - - - -
   for( Index j = Master->MaxName() ; j-- ; )  // pick any removable item
    if( ( OOBase[ j ] >= 0 ) && ( OOBase[ j ] < Inf<SIndex>() ) ) {
     i = j;
     break;
     }

  if( i == InINF ) {  // there are no removable items at all- - - - -
   BLOG( 0 , std::endl << "Bundle::FormD: unrecoverable MP failure." );
   Result = kError;
   return;
   }

  Delete( i );      // just delete i

  }  // end ( error-handling loop )- - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Sigma = Master->ReadSigma();                  // read Sigma*
 vStar[ NrFi ] = Master->ReadFiBLambda();      // read v*

 if( IsEasy.size() ) {                         // there are easy components
  for( Index k = 0 ; k < NrFi ; k++ )          // read the *exact* Fi-value
   if( IsEasy[ k ] )                           // for all them
    UpFiLmb1[ k ] = Master->ReadFiBLambda( k + 1 );
   else
    vStar[ k ] = Master->ReadFiBLambda( k + 1 );

  if( UpFiLmb[NrFi] < Inf<double>() )
   for( Index k = 0 ; k < NrFi ; k++ )
    if( IsEasy[ k ] )
     vStar[ NrFi ] += UpRifFi[ k ];
  }
 else
  for( Index k = 0 ; k < NrFi ; k++ )
   vStar[ k ] = Master->ReadFiBLambda( k + 1 );

 if( tStar > 0 )
  DSTS = Master->ReadDStart( tStar );                  // D_{t*,\beta,x}
 else
  DSTS = std::abs( tStar ) * Master->ReadDStart( 1 );  // | t* | * || z* ||
  
 // Sigma* + D*_{t*}( -z* ) is the "maximum expected increase" used in
 // the stopping criterion, EpsU is that relative to Fi( Lambda )

 if( UpFiLmb[ NrFi ] < Inf<double>() )
  EpsU = ( DSTS + Sigma ) / std::max( std::abs( UpFiLmb[ NrFi ] ) ,
				      double( 1 ) );
 else
  EpsU = 1;  // ensure EpsU is initialized somehow

 // the z[ i ] have changed, so in principle they are no longer in the
 // bundle: it may be the case that they actually are, but this is
 // taken care of in UpdtCntrs()
 whisZ.assign( NrFi , Inf<Index>() );

 // the scalar products have changed
 ScPr1.assign( NrFi , Inf<double>() );

 // additional information not present in the Bundle implementation
 // for NDOSolver interface  - - - - - - - - - - - - - - - - - - - - - - - - -

 DeltaStar = Master->ReadDStart( t ) / 2.0 + Sigma;
 cLMRow tdir = Master->Readd( true );

 NrmD = 0;                                    // d-norm
 for( Index i = 0 ; i < NumVar ; i++ )
  NrmD += tdir[ i ] * tdir[ i ];
 NrmD = sqrt(  NrmD );

 }  // end( BundleSolver::FormD )  - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdtCntrs( void )
{
 // increase all the OOBase[] counters but those == +/-Inf<SIndex>() - - - - -
 // items whose OOBase[] becomes 0 (e.g. the newly entered items, which have
 // OOBase[] == -1) are set to +1, in such a way that only the items in the
 // optimal base have OOBase[] == 0; note that the converse is not true, as
 // items in the optimal base may have OOBase[] < 0 instead

 for( SIndex* tOO = OOBase.data() + Master->MaxName() ; tOO-- > OOBase.data() ; )
  if( ( *tOO < Inf<SIndex>() ) && ( *tOO > -Inf<SIndex>() ) ) {
   (*tOO)++;
   if( ! *tOO )
    (*tOO)++;
   }

 // set to 0 the OOBase[] counter for items in base (if not < 0)- - - - - - -
 // note that there is a case in which a component wFi has Z[ wFi ] "for free"
 // in the bundle: this is when wFi only has *one* subgradient in base (or, in
 // practice, a subgradient with multiplier very close to one). This is
 // checked here (it is basically for free), and in case whisZ[] is properly
 // set so as to avoid pointless aggregations and OOBase[] is set to -1,
 // because under no circumnstances such a subgradient can ever be removed
 // from the bundle

 const Index* MBse;
 const double* Mlt = Master->ReadMult( MBse , MBDim );
 if( MBse ) {
  for( Index i ; ( i = *(MBse++) ) < InINF ; Mlt++ )
   if( *Mlt >= Eps<HpNum>() ) {
    if( ( *Mlt >= 1 - RAccSol ) && Master->IsSubG( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ Master->WComponent( i ) - 1 ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }
  }
 else
  for( Index i = 0 ; i < MBDim ; i++ , Mlt++ )
   if( *Mlt >= Eps<double>() ) { // ?? come mai non da' errore ??
    if( ( *Mlt >= 1 - RAccSol ) && Master->IsSubG( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ Master->WComponent( i ) - 1 ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }

 }  // end( UpdtCntrs )  - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FormLambda1( HpNum Tau )
{
 Master->MakeLambda1( Lambda.data() , Lambda1.data() , Tau );

 if( Master->NumBxdVars() ) {
  // as the relative precision required to the MPSolver is not enough to
  // ensure that the bounds on the variables will be satisfied with the
  // precision required by the FiOracle, the (upper and lower) bounds are
  // strictly enforced here

  std::vector<VarValue> tL1 = Lambda1;

  if( Master->NumNNVars() )             // there are NN vars and UB vars
   if( Master->NumNNVars() == NumVar )  // actually, all variables are NN
    for( Index h = 0 ; h < NumVar ; h++ ) {
     if( tL1[ h ] < 0 )
      tL1[ h ] = 0;

     const double UBh = LamVcblr[h]->get_ub();
     if( tL1[ h ] > UBh )
      tL1[ h ] = UBh;
     }
   else                                 // not all variables are NN
    for( Index h = 0 ; h < NumVar ; h++ ) {
     if( Master->IsNN( h ) && ( tL1[ h ] < 0 ) )
      tL1[ h ] = 0;

     const double UBh = LamVcblr[h]->get_ub();
     if( tL1[ h ] > UBh )
      tL1[ h ] = UBh;
     }
  else  // there are only UB vars
   for( Index h = 0 ; h < NumVar ; h++ ) {
    const double UBh = LamVcblr[h]->get_ub();
    if( tL1[ h ] > UBh )
     tL1[ h ] = UBh;
    }

  Lambda1 = tL1;
  }  // end( if( the bounds have to be enforced ) )


 whisG1.assign( NrFi , Inf<Index>() );

 // Lambda has changed, pass the new one to the oracle - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiStatus.assign( NrFi ,  kUnEval );
 for( Index i = 0 ; i < NumVar ; i++ )
  LamVcblr[ i ]->set_value( Lambda1[i] );

 // compute the upper and lower model at the tentative point   - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_lf ) { // add the linear part to the "full function"
  f_lf->compute( true );
  UpFiLmb1[ NrFi ] = f_lf->get_upper_estimate();
  LwFiLmb1[ NrFi ] = f_lf->get_lower_estimate();
  }
 else
  UpFiLmb1[ NrFi ] = LwFiLmb1[ NrFi ] = 0;

 for( Index k = 0 ; k < NrFi ; k++ ) {
  if( IsEasy.size() && IsEasy[ k ] )  // if k is an easy component
   UpFiLmb1[ k ] =  LwFiLmb1[ k ] = Master->ReadFiBLambda( k );
  else {

   // initialize upper and lower bound for each component  - - - - - - - - - -

   if( v_c05f[ k ]->get_Lipschitz_constant() < Inf<VarValue>()
	   && UpFiLmb[ k ] < Inf<VarValue>() )
    UpFiLmb1[ k ] = UpRifFi[ k ] + v_c05f[ k ]->get_Lipschitz_constant() * NrmD;
   else
	UpFiLmb1[ k ] = Inf<VarValue>();

   if( LwFiLmb[ k ] > -Inf<VarValue>() )
    LwFiLmb1[ k ] = UpRifFi[ k ] + vStar[ k ];
   else
    LwFiLmb1[ k ] = -Inf<VarValue>();
   }

  // sum over the components, the zero-component is already there
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1[ NrFi ] < Inf<VarValue>() ) {
   if( UpFiLmb1[ k ] < Inf<VarValue>() )
    UpFiLmb1[ NrFi ] += UpFiLmb1[ k ];
   else
	UpFiLmb1[ NrFi ] = Inf<VarValue>();
   }

  if( LwFiLmb1[ NrFi ] > -Inf<VarValue>() ) {
   if( LwFiLmb1[ k ] < Inf<VarValue>() )
    LwFiLmb1[ NrFi ] += LwFiLmb1[ k ];
   else
    LwFiLmb1[ NrFi ] = -Inf<VarValue>();
   }
  }

 // update the upper and lower targets - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( UpFiLmb[ NrFi ] < Inf<VarValue>() )
  UpTrgt = UpRifFi[ NrFi ] + ( 1.0 - m2 ) * vStar[ NrFi ];
 else
  UpTrgt = Inf<VarValue>();

 if( LwFiLmb[ NrFi ] > -Inf<VarValue>() )
  if( m1 > 0 )
   LwTrgt = UpRifFi[ NrFi ] + vStar[ NrFi ] + m1 * DeltaStar;
  else
   LwTrgt = UpRifFi[ NrFi ] + ( 1.0 + m1 ) * vStar[ NrFi ];
 else
  LwTrgt = -Inf<VarValue>();

 }  // end( BundleSolver::FormLambda1 )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::FiAndGi( Index wFi )
{
 double UpCutOff, LwCutOff, LwFiK, EpsCurr;

 if( IsEasy.size() && IsEasy[ wFi ] )
  return( false );

 LwFiK =  UpRifFi[ wFi ] + vStar[ wFi ];

 if( UpFiLmb[ wFi ] < Inf<VarValue>() )
  if( UpTrgt < Inf<VarValue>() && UpFiLmb1[ NrFi ] < Inf<VarValue>() )
   UpCutOff = std::max( UpTrgt - ( UpFiLmb1[ NrFi ] - UpFiLmb1[ wFi ] ) ,
   		              LwFiK - m2 * BetaK( wFi ) * vStar[ NrFi ] );
  else
   UpCutOff = LwFiK - m2 * BetaK( wFi ) * vStar[ NrFi ];
 else
  UpCutOff = Inf<VarValue>();

 if( LwFiLmb[ wFi ] > -Inf<VarValue>() )
  if( LwTrgt > -Inf<VarValue>() && LwFiLmb1[ NrFi ] > -Inf<VarValue>() )
   LwCutOff = std::max( LwTrgt - ( LwFiLmb1[ NrFi ] - LwFiLmb1[ wFi ] ) ,
    		              LwFiK + m1 * BetaK( wFi ) * DeltaStar );
  else
   LwCutOff = LwFiK + m1 * BetaK( wFi ) * DeltaStar;
 else
  LwCutOff = -Inf<VarValue>();

 if( ( LwCutOff > -Inf<VarValue>() ) && ( UpCutOff < Inf<VarValue>() ) )
  EpsCurr = ( UpCutOff - LwCutOff ) / std::max( 1.0 , std::abs(UpRifFi[ wFi ] ) );
 else
  EpsCurr = RelAcc / Nearly;

 // assign the cutoff values to the c05Function - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_c05f[ wFi ]->set_par( dblUpCutOff , UpCutOff );
 v_c05f[ wFi ]->set_par( dblLwCutOff , LwCutOff );
 v_c05f[ wFi ]->set_par( dblRelAcc , EpsCurr );

 if( FiStatus[ wFi ] ==  kUnEval )
  FiStatus[ wFi ] = v_c05f[ wFi ]->compute( true );
 else
  FiStatus[ wFi ] = v_c05f[ wFi ]->compute( false );

 if( UpFiLmb1[ NrFi ] < Inf<VarValue>() )
  UpFiLmb1[ NrFi ] -= UpFiLmb1[ wFi ];

 if( LwFiLmb1[ NrFi ] > -Inf<VarValue>() )
  LwFiLmb1[ NrFi ] -= LwFiLmb1[ wFi ];

 UpFiLmb1[ wFi ] = std::min( v_c05f[ wFi ]->get_upper_estimate() , UpFiLmb1[ wFi ] );
 LwFiLmb1[ wFi ] = std::max( v_c05f[ wFi ]->get_lower_estimate() , LwFiLmb1[ wFi ] );

 if( UpFiLmb1[ NrFi ] < Inf<VarValue>() )
  UpFiLmb1[ NrFi ] += UpFiLmb1[ wFi ];
 else
  if( UpFiLmb1[ wFi ] < Inf<VarValue>() ) {
   if( f_lf )
	UpFiLmb1[ NrFi ] = f_lf->get_upper_estimate();
   else
    UpFiLmb1[ NrFi ] = 0;

   for( Index k = 0 ; k < NrFi ; k++ )
	if( UpFiLmb1[ k ] < Inf<VarValue>() )
     UpFiLmb1[ NrFi ] += UpFiLmb1[ wFi ];
	else {
     UpFiLmb1[ NrFi ] = Inf<VarValue>();
     break;
	 }

   }

 if( LwFiLmb1[ NrFi ] > -Inf<VarValue>() )
  LwFiLmb1[ NrFi ] += LwFiLmb1[ wFi ];
 else
  if( LwFiLmb1[ wFi ] > -Inf<VarValue>() ) {
   if( f_lf )
    LwFiLmb1[ NrFi ] = f_lf->get_lower_estimate();
   else
    LwFiLmb1[ NrFi ] = 0;

   for( Index k = 0 ; k < NrFi ; k++ )
    if( LwFiLmb1[ k ] > -Inf<VarValue>() )
     LwFiLmb1[ NrFi ] += LwFiLmb1[ wFi ];
    else {
     LwFiLmb1[ NrFi ] = -Inf<VarValue>();
     break;
     }
   }

 if( UpFiLmb1[ NrFi ] == Inf<VarValue>() )  // Fi() is not defined in Lambda1
  DeltaFi = Inf<VarValue>();
 else
  DeltaFi = UpFiLmb1[ NrFi ] - UpRifFi[ NrFi ];

 FiEvaltns++;

 // update FiBest, if necessary - - - - - - - - - - - - - - - - - - - - - - -

 if( UpFiLmb1[ NrFi ] < UpFiBest ) {
  UpFiBest = UpFiLmb1[ NrFi ];
  if( MaxSol > 1 )
   LmbdBst = Lambda1;
  }

 // get a new linearization - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool HasLinearization;
 bool diagonal;

 for( Index Ftchd = 0 ; Ftchd < aBP3 ; ) {

  diagonal = true;
  if( Ftchd ==  0 ) {

   // first look for a constraint then for a sub-gradient
   if( UpFiLmb1[ wFi ] == Inf<VarValue>() ) {
    HasLinearization = v_c05f[ wFi ]->has_linearization( diagonal = false );
    if( !HasLinearization )
     HasLinearization = v_c05f[ wFi ]->has_linearization( diagonal );
    }
   else
	HasLinearization = v_c05f[ wFi ]->has_linearization( diagonal );

   }
  else {

   if( UpFiLmb1[ wFi ] == Inf<VarValue>() ) {
    HasLinearization = v_c05f[ wFi ]->compute_new_linearization( diagonal = false );
    if( !HasLinearization )
     HasLinearization = v_c05f[ wFi ]->compute_new_linearization( diagonal );
    }
   else
	HasLinearization = v_c05f[ wFi ]->compute_new_linearization( diagonal );
   }

  if( !HasLinearization )
   break;

  // check if aggregation has to be performed - - - - - - - - - - - - - - - -
  // doing this now could occasionally result in useless aggregations, but it
  // avoids complications in the interface of MPSolver (inserting some
  // Z[ wFi ] while inserting the new item)

  Index wh = BStrategy( wFi );

  // get the space for the item from the MPSolver - - - - - - - - - - - - - -

  double* G1 = Master->GetItem( wFi + 1 );

  // fetch the item from the Oracle - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set SGBse = nullptr;
  v_c05f[ wFi ]->get_linearization_coefficients( G1 );
  HpNum eps = v_c05f[ wFi ]->get_linearization_constant();

  GiEvaltns++;

  // pass the base to the MP Solver - - - - - - - - - - - - - - - - - - - - -

  Master->SetItemBse( SGBse , NumVar );

  // calculate ScPr1k and Alfa1k- - - - - - - - - - - - - - - - - - - - - - -

  Index cp;
  HpNum ScPr1k;

  // update alpha value at Lambda1 point  - - - - - - - - - - - - - - - - - -

  eps = UpFiLmb1[ wFi ] - eps
		  - std::inner_product( Lambda1.begin() , Lambda1.end() , G1 , double(0) );

  HpNum Alfa1k = eps;

  if( !diagonal )  // it is a constraint
   cp = Master->CheckCnst( Alfa1k , ScPr1k , Lambda.data() );
  else                             // it is a subgradient
   cp = Master->CheckSubG( UpFiLmb1[ wFi ] - UpRifFi[ wFi ] ,
                          t , Alfa1k , ScPr1k );

  if( cp < InINF ) {  // the item is a copy- - - - - - - - - - - - - -
   BLOGb( LogBnd , std::endl << "New item is a copy of " << cp );

   cHpNum OrigA1k = (Master->ReadLinErr())[ cp ];

   if( OrigA1k > Alfa1k ) {        // if the copy has smaller Alfa than the
    Master->SubstItem( wh = cp );  // original, substitute it

    BLOGb( LogBnd , " with smaller Alfa" );
    }
   else
    wh = InINF;               // otherwise, nothing new has happened
   }
  else {              // insert the item, if there is space - - - - - - - - -
   if( wh < InINF )     // someone has been selected in BStrategy()
    Master->RmvItem( wh );     // remove it from the MP
   else
    wh = FindAPlace( wFi );    // find a spot in the bundle

   if( wh == InINF ) {  // no space found ...
    if( ! Ftchd ) {            // ... and this was the first item
     BLOG( 0 , std::endl << " ERROR: No space in the bundle" << std::endl );
     Result = kError;          // signal an error
                               // ensure that the outer Fi-cycle ends
     }
    else
     BLOG( 1 , std::endl << " WARNING: No space in the bundle" << std::endl );

     break;                     // the cycle ends
    }

   Master->SetItem( wh );      // insert the item in the MP Solver
   OOBase[ wh ] = -1;          // ensure it won't be touched again this round

   if( f_log && LogVerb ) {
    if( ! diagonal )
     *f_log << std::endl << "New constraint " << wh << ", rhs = " << Alfa1k;
    else
     *f_log << std::endl << "New eps-subgradient " << wh << " for Fi[ "
	    << wFi << " ] , eps = " << eps <<
	    " , Alfa1 = " << Alfa1k << ", gd = " << - ScPr1k;
    }
   }

  // if something was inserted, bookkeeping is needed - - - - - - - - - - - -

  if( wh < InINF ) {
   Ftchd++;                    // one more item
   v_c05f[ wFi ]->store_linearization( wh ); // tell the name of the item to the FiOracle

   if( UpFiLmb1[ wFi ] < Inf<double>() ) {  // it is a subgradient
    if( ( whisG1[ wFi ] == InINF ) || ( Alfa1k < Alfa1[ wFi ] ) ||
       ( ( Alfa1k == Alfa1[ wFi ] ) && ( ScPr1k > ScPr1[ wFi ] ) ) ) {
     whisG1[ wFi ] = wh;       // wh is the new representative of wFi
     Alfa1[ wFi ] = Alfa1k;
     ScPr1[ wFi ] = ScPr1k;
     }
    }
   else
    OOBase[ wh ] = - Inf<SIndex>();
     /* if the item is a constraint, mark it as permanently fixed: this may be
       a bad choice in practice, although it is required by the theory
       (we'll see ...) */
   }

  // compute *Alfa1 and *ScPr1 - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Alfa1[NrFi] = 0;
  ScPr1[NrFi] = Master->ReadGid();

  for( Index k = 0 ; k < NrFi ; k++ )
   if( whisG1[ k ] < InINF ) {
    if( Alfa1[ k ] == Inf<double>() )
     Alfa1[ k ] = (Master->ReadLinErr())[ whisG1[ k ] ];

    Alfa1[NrFi] += Alfa1[ k ];

    if( ScPr1[ k ] == Inf<double>() )
     ScPr1[ k ] = Master->ReadGid( whisG1[ k ] );

    ScPr1[NrFi] += ScPr1[ k ];
    }
   else
    Alfa1[ k ] = ScPr1[ k ] = 0;

  }

 // update lower and upper estimates  - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( LwFiLmb1[ NrFi ] > LwTrgt ) || ( UpFiLmb1[ NrFi ] < UpTrgt ) )
  return( true );
 else
  return( false );

 }  // end( BundleSolver::FiAndGi() )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::GotoLambda1( void )
{
 std::vector<VarValue>  DeltaFi( NumVar );
 std::transform( UpFiLmb1.begin(), UpFiLmb1.end(), UpRifFi.begin(),
		 DeltaFi.begin(), std::minus<double>() );

 // do the move - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Lambda.swap(Lambda1);
 UpFiLmb.swap(UpFiLmb1);
 LwFiLmb.swap(LwFiLmb1);
 UpRifFi = UpFiLmb;

 // change the current point in the MP Solver - - - - - - - - - - - - - - - -

 Master->ChangeCurrPoint( t , DeltaFi.data() );

 // signal that Alfa1[] is not reliable - - - - - - - - - - - - - - - - - - -

 Alfa1.assign( NrFi + 1 , Inf<double>() );

 }  // end( GotoLambda1 )

/*--------------------------------------------------------------------------*/

void BundleSolver::SimpleBStrat( void )
{
 for( SIndex* tOO = OOBase.data() + Master->MaxName() ; tOO-- > OOBase.data() ; )
  if( ( *tOO < Inf<SIndex>() ) && ( *tOO > SIndex( BPar1 ) ) ) {
   const Index h = tOO - OOBase.data();
   v_c05f[ Master->WComponent( h ) - 1 ]->delete_linearization( h );
   Delete( h );
   }
 }  // end( BundleSolver::SimpleBStrat ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdtLowerBound( void )
{
 // first of all, check if a "hard" lower bound is available
 double LwrBnd = f_Block->get_valid_lower_bound( false );
 if( LwrBnd > - Inf<double>() )
  TrueLB = true;
 else {
  // if not, check if at least a "conditional" one is available
  TrueLB = false;
  LwrBnd = f_Block->get_valid_lower_bound( true );
  }

 if( LwrBnd != LowerBound ) {
  LowerBound = LwrBnd;
  LBHasChgd = true;
  }
 } // end( BundleSolver::UpdtLowerBound )  - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

double BundleSolver::BetaK( Index wFi )
{
 return( 1.0 / double( NrEasy ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::Log1( void )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "{" << SCalls << "-" << ParIter << "-"
	<< Master->MaxName() - FreList.size() << "-" << MBDim << "} t = " << t
	<< " ~ D*_1( z* ) = " << Master->ReadDStart( 1 )
	<< " ~ Sigma = " << Sigma << std::endl << "           ";

 *f_log <<  " Fi = ";

 if( UpFiLmb[ NrFi ] == Inf<double>() )
  *f_log << " - INF";
 else
  *f_log << UpFiLmb[NrFi] << " ~ eU = " << EpsU;

 if( BPar6 )
  *f_log << " ~ BP3 = " << aBP3;

 } // end( BundleSolver::Log1 )  - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::Log2( void )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "           ";

 if( LowerBound > - Inf<double>() )
  *f_log << "LB = " << LowerBound << " ~ ";

 *f_log << "Fi1 = ";

 if( UpFiLmb1[ NrFi ] <= - Inf<double>() )
  *f_log << "+ INF => STOP." << std::endl;
 else
  if( UpFiLmb1[ NrFi ] >= Inf<double>() )
   *f_log << " - INF" << std::endl;
  else
   *f_log << UpFiLmb1[ NrFi ] << " ~ Alfa1 = " << Alfa1[ NrFi ]
	  << " ~ Gi1xd = " << - ScPr1[ NrFi ] << std::endl;

 } // end( BundleSolver::Log2 )  - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::InitMP( void )
{
 // this method is called only when *both* the FiOracle *and* the MPSolver
 // have been set, and it is re-called each time any one of the two changes
 // set the size- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FakeFi.initialize();
 Master->SetDim( BPar2 * ( NrFi - NrEasy ) , &FakeFi , false );

 Master->SetPar( MPSolver::kOptEps , RelMPAcc );
 Master->SetPar( MPSolver::kFsbEps , RMPAccSol );

 // insert the constant subgradient of the 0-th component - - - - - - - - - -
 // TODO: check if the 0-th component is sparse, if so pass a proper base
 //       to the MPSolver

 if( f_lf ) {
  f_lf->get_linearization_coefficients( Master->GetItem( 0 ) );
  Master->SetItemBse( nullptr , NumVar );
  Master->SetItem( InINF );
  }

 tHasChgd = LBHasChgd = true;

 }  // end( BundleSolver::InitMP( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::FindNext( Index & wFi )
{
 Index PrWFi = wFi;
 bool NextIsAccepted = false;
 do {
  wFi = ( wFi + 1 ) % NrFi;
  if( ( FiStatus[ wFi ] == kUnEval ) ||
      ( ( FiStatus[ wFi ] < kError ) && ( FiStatus[ wFi ] > kOK ) &&
	( CurrNrEvls[ wFi ] < MaxNrEvls ) ) )
   NextIsAccepted = true;
  } while( ( ! NextIsAccepted ) && ( wFi != PrWFi ) );

 return( NextIsAccepted );

 } // end( BundleSolver::FindNext( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::BStrategy( Index wFi )
{
 // this method implements the B-strategies of the code, i.e., which "old"
 // items are discarded if the bundle is full and a new item belonging to
 // component wFi has to be inserted
 // note: this is called *before* that we know if the place will actually be
 // required, so it may return InINF and nothing bad may happen

 // in particular, it returns InINF is there is plenty of space left
 // in the bundle so that no B-strategy (no removal or aggregation) is
 // required; picking a specific spot in the free space is the task of
 // FindAPlace(), which however is not called right away because the place
 // may end up not being needed

 if( ( ! FreList.empty() ) || ( Master->MaxName() <
				Index( BPar2 * ( NrFi - NrEasy ) ) ) )
  return( InINF );

 // there is not plenty of space- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // find the removable item with largest OOBase[]; among these with the
 // largest OOBase[], select that with largest Alfa[]
 // note: the Z[ wFi ] in the bundle (if any) are not removable and
 //       therefore cannot be selected, which in particular happens if
 //       wFi has only *one* subgradient in base

 Index wh = 0;
 cHpRow tA = Master->ReadLinErr();
 for( Index i = 0 ; ++i < Index( BPar2 * ( NrFi - NrEasy ) ) ; )
  if( ( OOBase[ i ] > OOBase[ wh ] ) ||
      ( ( OOBase[ i ] == OOBase[ wh ] ) && ( tA[ i ] > tA[ wh ] ) ) )
   wh = i;

 if( OOBase[ wh ] < 0 )    // all items are non-removable: nothing else to
  return( InINF );  // do (except maybe complaining very loudly)
 else
  if( OOBase[ wh ] > 0 )   // a place is found
   return( wh );           // there are no problems, all done

 // wh is a basic item- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this means that *all* items of *all* components are either in base or not
 // removable, for otherwise we would have selected an item with OOBase > 0;
 // we cannot discard anything before having performed aggregation, but in
 // order to do so we also need to free some space for the Z[]

 wh = InINF;
 for( Index wFi2 = wFi ; ; ) {
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is done component-wise, round-robin starting from wFi in order to
  // try to keep the balance between the space allocated to each component
  if( IsEasy.size() && IsEasy[ wFi2 ] ) {  // ... but skipping easy components
    wFi2 = ( wFi2 == NrFi ? 1 : wFi2 + 1 );
    if( wFi2 == wFi )
     break;
    else
     continue;
    }

  Index MBDm;
  cIndex_Set MBse;
  cHpRow Mlt = Master->ReadMult( MBse , MBDm , wFi2 , false );
  // note that since *all* items of *all* components are either in base or
  // not removable we can scan MBse[] for the items to be removed, possibly
  // ignoring items with Mlt[] == 0 -- but in fact not doing it because there
  // is not any

  if( whisZ[ wFi2 ] == InINF ) {
   // Z[ wFi2 ] is not already in: in principle aggregation might have to
   // be performed, so select the item to take Z[ wFi2 ] as the one with
   // smallest Mlt

   cHpRow tM = Mlt;
   Index whZ = InINF;
   HpNum tMin = Inf<double>();
   if( MBse ) {
    cIndex_Set tB = MBse;
    for( Index h ; ( h = *(tB++) ) < InINF ; tM++ )
     if( ( *tM < tMin ) && ( OOBase[ h ] >= 0 ) ) {
      whZ = h;
      tMin = *tM;
      }
    }
   else
    for( Index h = 0 ; h < MBDm ; h++ , tM++ )
     if( ( *tM < tMin ) && ( OOBase[ h ] >= 0 ) ) {
      whZ = h;
      tMin = *tM;
      }

   if( whZ < InINF ) {  // if there is space for Z[ wFi2 ]
    AggregateZ( Mlt , MBse , MBDm , wFi2 , whZ );  // put it in whZ
    Mlt = Master->ReadMult( MBse , MBDm , wFi2 , false );
    // and now ensure that Mlt and MBse are reliable again; this is
    // needed because AggregateZ() calls methods of the MPSolver which
    // may therefore "invalidate" temporary vectors like Mlt and MBse
    }
   else {
    // there is *no* removable item of this component: this in some sense is
    // no problem, because in this case aggregation is not needed, but on the
    // other hand it means that there is no way this component will provide
    // any place, so have to move to the next
    wFi2 = ( wFi2 == NrFi ? 1 : wFi2 + 1 );
    if( wFi2 == wFi )
     break;
    else
     continue;
    }
   }  // end if( Z[ wFi2 ] is not already in )

  // at this point, Z[ wFi2 ] is in the bundle: try to select the exiting
  // item as the one with the smallest Mult[] among these belonging to wFi2
  // this may fail, prompting to move to the next component, so one could
  // fear of having just done aggregation to no avoil; yet I don't think
  // this can happen, as the only case would be that of having only one
  // removable subgradient in the wFi2 base (that is taken by whZ), but the
  // only reasonable case in which this can happen is that there is only one
  // subgradient at all, in which case it is Z[ wFi2 ] and no aggregation is
  // done (in other words, you do aggregation only if you have at least two
  // subgradients in base, and there is no reason one of them should be non
  // removable)
  HpNum tMin = Inf<double>();
  if( MBse ) {
   for( Index h ; ( h = *(MBse++) ) < InINF ; Mlt++ )
    if( ( *Mlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
     wh = h;
     tMin = *Mlt;
     }
   }
  else
   for( Index h = 0 ; h < MBDm ; h++ , Mlt++ )
    if( ( *Mlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
     wh = h;
     tMin = *Mlt;
     }

  if( wh < InINF )  // a place has been found
   break;                  // all done
  else {                   // all the items of wFi2 are not removable
   wFi2 = ( wFi2 == NrFi ? 1 : wFi2 + 1 );  // try the next component
   if( wFi2 == wFi )                        // if any has remained
    break;                                  // otherwise give up
   }
  }  // end for( wFi2 ) - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( wh );

 }  // end( BundleSolver::BStrategy )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::FindAPlace( cIndex wFi )
{
 // this method is used to return the index of an available position in the
 // bundle where to store a new item belonging to "component" wFi; if there
 // are no possible positions left, then InINF is returned

 Index wh = InINF;

 if( FreList.size() ) {              // there are deleted items
  wh = FreList.top();   // pick one
  FreList.pop();
  }
 else                                       // there are no deleted items ...
  if( Master->MaxName() < Index( BPar2 * ( NrFi - NrEasy ) ) )
                                            // ... but there is still space
   wh = Master->MaxName();                  // next name

 assert( Master->MaxName() >= FreList.size() );

 return( wh );

 }  // end( BundleSolver::FindAPlace ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::AggregateZ( cHpRow Mlt , cIndex_Set MBse , Index MBDm ,
			  	cIndex wFi , cIndex whr )
{
 // note: this is *never* called in the "easy" case where aggregation is not
 // needed because there is only one subgradient in base (for the component
 // wFi) and therefore, its Mlt[] is == 1, since in this case whisZ[] has
 // been properly set in UpdtCntrs()

 // tell the C05Function what is going to happen  - - - - - - - - - - - - - - -

 LinearCombination coefficients( MBDm );
 for( Index i = 0 ; i < MBDm ; i++ ) {
  coefficients[i].first =  MBse[ i ];
  coefficients[i].second =  Mlt[ i ];
  }

 v_c05f[ wFi ]->store_combination_of_linearizations( coefficients , whr );

 // ask the MPSolver the memory for keeping Z[ wFi ]- - - - - - - - - - - - -
 // note: Mlt and MBse could very well be "temporary" memory belonging to the
 // MPSolver, and any call to a method of the MPSolver may invalidate it;
 // the calls start now, and in fact MBse and Mlt are no longer used

 SgRow tZ = Master->GetItem( wFi );

 // read Z[ wFi ] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index ZBDm;
 cIndex_Set ZBse;
 Master->ReadZ( tZ , ZBse , ZBDm , wFi );

 // now pass Z[ wFi ] back to the MP Solver - - - - - - - - - - - - - - - - -

 Master->SetItemBse( ZBse , ZBDm );

 HpNum ScPri;
 HpNum Ai = Master->ReadSigma( wFi );          // its alfa is Sigma[ wFi ]

 Master->CheckSubG( 0 , 0 , Ai , ScPri );      // DFi == Tau == 0

 Master->RmvItem( whr );  // remove the old item in position whr

 Master->SetItem( whr );  // set Z[ wFi ] in position whr

 whisZ[ wFi ] = whr;      // Z[ wFi ] is in the bundle ...
 OOBase[ whr ] = -1;      // ... and it won't be removed in this iteration

 BLOGb( LogBnd , std::endl << "Aggregation performed into " << whr );

 }  // end( BundleSolver::AggregateZ() ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::Heuristic1( void )
{
 if( Alfa1[0] < Eps<double>() )
  return( DeltaFi > Eps<double>() ? tMaior : tMinor );
 else
  return( t * ( ( DeltaFi + Alfa1[NrFi] ) / ( 2 * Alfa1[NrFi] ) ) );
 } // end( BundleSolver::Heuristic1() ) - -  - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::Heuristic2( void )
{
 if( std::abs( vStar[ NrFi ] + DeltaFi ) < Eps<double>() )
  return( tMaior );
 else
  return( t * abs( vStar[ NrFi ] / ( 2 * ( vStar[ NrFi ] + DeltaFi ) ) ) );
 } // end( BundleSolver::Heuristic2() ) - -  - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::RemoveItems( void )
{
 if( Master )
  Master->RmvItems();  // remove all items from the MPSolver (if any)

 if( v_c05f.size() )
  for( auto & fun : v_c05f ) { // tell the c05Function (if any) about it
   int GPMaxSz = fun->get_int_par( C05Function::intGPMaxSz );
   for( Index i = 0 ; i < GPMaxSz ; i++ )
	fun->delete_linearization( i );
    }

 if( f_lf ) {
  int GPMaxSz = f_lf->get_int_par( C05Function::intGPMaxSz );
  for( Index i = 0 ; i < GPMaxSz ; i++ )
   f_lf->delete_linearization( i );
  }

 FreList = priority_queue<Index>();

 }  // end( BundleSolver::RemoveItems )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::guts_of_destructor( void )
{
 LmbdBst.clear(); // does it make sense??
 Lambda1.clear();
 Lambda.clear();

 Alfa1.clear();
 ScPr1.clear();
 whisG1.clear();

 FiStatus.clear();

 whisZ.clear();
 FreList = priority_queue<Index>();
 OOBase.clear();

 if( !IsEasy.empty() )
  IsEasy.clear();

 }  // end( BundleSolver::MemDealloc( ) )  - - - - - - - - - - - - - - - - - -

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

 //!! check if MPSolver != NULL !!

 // if( ! ( RstLvl & RstCrr ) )  // reset the current point to all-0- - - - - - -
 // SetLambda();

 if( ! ( RstLvl & ( RstSbg | RstCnt ) ) )  // reset everything- - - - - - - -
  RemoveItems();
 else
  if( ! ( RstLvl & RstSbg ) ) {  // reset subgrads (but not constrs)- - - - -
   if( Master->BSize() ) {       // if the bundle is nonempty
    if( ! Master->BCSize() )     // and it contains only subgradients
     RemoveItems();              // remove everything
    else
     for( Index i = Master->MaxName() ; i-- ; )
      if( Master->IsSubG( i ) ) {
       v_c05f[ Master->WComponent( i ) - 1 ]->delete_linearization( i );
       Delete( i );
       }
    }
   }
  else
   if( ! ( RstLvl & RstCnt ) ) // reset constrs (but not subgrads) - - - - - -
    if( Master->BSize() ) {    // if the bundle is nonempty
     if( Master->BSize() == Master->BCSize() )
                               // and it contains only constrs
      RemoveItems();           // remove everything
     else
      for( Index i = Master->MaxName() ; i-- ; )
       if( ! Master->IsSubG( i ) ) {
    	v_c05f[ Master->WComponent( i ) - 1 ]->delete_linearization( i );
	Delete( i );
        }
     }

 if( ! ( RstLvl & RstFiV ) )  // reset the current value of Fi( Lambda ) - - -
  UpFiLmb[ NrFi ] = Inf< VarValue >();

 if( !( RstLvl & NoStPt ) ) {  // get an initial point
  UpFiLmb[ NrFi ] = Inf< VarValue >();
  for( Index i = 0 ; i < NumVar ; i++ )
   Lambda[i] = LamVcblr[ i ]->get_value( );
  }

 }  // end( BundleSolver::ReSetAlg ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::Delete( cIndex i )
{
 if( Master ) {
  // check if this item was the "representative" for its component- - - - - -

  cIndex wFi = Master->WComponent( i ) - 1;

  if( Master->IsSubG( i ) )  // it is a subgradient
   if( whisG1[ wFi ] == i )  // it is the representative of wFi
    whisG1[ wFi ] = InINF;   // a new representative is needed

  // delete the item with name `i' from the MP- - - - - - - - - - - - - - - -

  Master->RmvItem( i );
  }

 BLOGb( LogBnd , std::endl << "Item " << i << " removed" );

 // bookkeeping of internal data structures - - - - - - - - - - - - - - - - -

 FreList.push( i );
 OOBase[ i ] = Inf<SIndex>();

 // compacting FreList[] if it's too big- - - - - - - - - - - - - - - - - - -
 // remove from FreList[] every name >= Master->MaxName(); note that every
 // ordered set *is* a Heap. apart from efficiency reasons, this is
 // needed because Master->MaxName() - FreDim is the only way in which the
 // Bundle can compute the number of "live" items

 cIndex MxNm = Master->MaxName();
 if( FreList.size() > MxNm ) {
  FreList = priority_queue<Index>();
  for( Index i = 0 ; i < MxNm ; i++ )
   if( OOBase[ i ] == Inf<SIndex>() )
    FreList.push( i );
  }

 assert( Master->MaxName() >= FreList.size() );

 }  // end( BundleSolver::Delete() ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdtaBP3( void )
{
 if( BPar5 == 0 )
  return;

 switch( BPar6 ) {
  case( 4 ):
   if( UpFiLmb[ NrFi ] > -Inf<double>() )
    aBP3 = ( BPar5 > 0 ? BPar4 : BPar3 ) +
           Index( BPar5 / std::log10( EpsU / RelAcc ) );
   break;
  case( 3 ):
   if( UpFiLmb[NrFi] > -Inf<double>() )
    aBP3 = ( BPar5 > 0 ? BPar4 : BPar3 ) +
           Index( BPar5 / std::sqrt( EpsU / RelAcc ) );
   break;
  case( 2 ):
   if( UpFiLmb[NrFi] > -Inf<double>() )
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

 if( aBP3 > BPar3 )
  aBP3 = BPar3;
 else
  if( aBP3 < BPar4 )
   aBP3 = BPar4;

 }  // end( BundleSolver::UpdtaBP3 ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::IsOptimal( double eps ) const
{
 if( eps <= 0 )
  eps = RelAcc;

 c_VarValue err = max_error( eps );
 if( err >= Inf< VarValue >() )
  return( false );

 if( tStar > 0 )
  return( DSTS + Sigma <= err );
 else
  return( ( Sigma <= err ) &&
	  ( DSTS <= std::min( AAccSol , RAccSol * std::abs( tStar ) ) ) );

 } // end( BundleSolver::IsOptimal() ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::CheckAlfa( const bool All )
{
 return( Sigma >= - t * m3 * Master->ReadDStart( t ) );
 }  // end( CheckAlfa )  - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::process_outstanding_Modification( void )
{
 // no-frills loop: do them in order, with no attempt at optimizing
 while( ! v_mod.empty() ) {
  auto mod = v_mod.front();  // pick (a reference to) the first Modification


  /* Use a Lambda to define a "guts" of the method that can be called
     recursively. Note the trick of defining the std::function object and
     "passing" it to the lambda, which allows recursive calls. Note the need
     to explicitly capture "this" to use fields/methods of the class. */

  // auto MCFB = static_cast< MCFBlock * >( f_Block );

  std::function< void( sp_Mod )> guts_of_poM;
  guts_of_poM = [ this , & guts_of_poM ]( sp_Mod mod ) {

   // process Modification - - - - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   /* This requires to patiently sift through the possible Modification types
       to find what this Modification exactly is, and call the appropriate
       method. */

   // GroupModification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   {
    const auto tmod = std::dynamic_pointer_cast<GroupModification>( mod );
    if( tmod ) {
     for( const auto & submod : tmod->v_sub_Modifications )
      guts_of_poM( submod );

     return;
     }
    } // end GroupModification - - - - - - - - - - - - - - - - - - - - - - - -

   // C05FunctionMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   {
    const auto tmod = std::dynamic_pointer_cast<C05FunctionMod>( mod );
    if( tmod ) {
     Index wFi = get_index_of_component(tmod->function());
     std::vector<double> Alfa1;
     switch( tmod->type() ) {
      case( C05FunctionMod::AllLinearizationChanged ):
	   Master->ChgSubG( 0 , NumVar , wFi+1 );
      case( C05FunctionMod::AlphaChanged ):
	   // a finite f_shift should be treated in a different way but
	   // as of now the finite shifts are ignored by MPSolver
	   Alfa1.resize( Master->MaxName(wFi+1) );
       for( Index i = 0 ; i < Master->MaxName(wFi+1) ; i++ )
        if( Master->WComponent( i ) == wFi+1 ) {
         Alfa1[ i ] = v_c05f[wFi]->get_linearization_constant( i );
         // alpha has to be referred to \Lambda
         std::vector<double> G1(NumVar);
         Range range = make_pair( 0, NumVar );
         v_c05f[ wFi ]->get_linearization_coefficients( G1.data() , range, i );
         Alfa1[ i ] = UpRifFi[ wFi ] - Alfa1[ i ]
          		  - std::inner_product( Lambda.begin() , Lambda.end() , G1.data() , double(0) );
         Master->ChgAlfa( Alfa1.data() , wFi );
        }
       break;
      case( C05FunctionMod::AllEntriesChanged ):
       Master->ChgSubG( 0 , NumVar , wFi+1 );
       break;
      } // switch( tmod->f_type )
     }  // end  if( tmod )
    } // end C05FunctionMod  - - - - - - - - - - - - - - - - - - - - - - - - -

   // C05FunctionModRngd - - - - - - - - - - - - - - - - - - - - - - - - - - -
   {
    const auto tmod = std::dynamic_pointer_cast<C05FunctionModRngd>( mod );
    if( tmod ) {
     Index wFi = get_index_of_component(tmod->function());
     std::vector<double> Alfa1;
     switch( tmod->type() ) {
      case( C05FunctionMod::AllLinearizationChanged ):
        Master->ChgSubG( tmod->range().first , tmod->range().second , wFi+1 );
      case( C05FunctionMod::AlphaChanged ):
       // a finite f_shift should be treated in a different way but
       // as of now the finite shifts are ignored by MPSolver
       Alfa1.resize( Master->MaxName(wFi+1) );
       for( Index i = 0 ; i < Master->MaxName(wFi+1) ; i++ )
        if( Master->WComponent( i ) == wFi+1 ) {
         Alfa1[ i ] = v_c05f[wFi]->get_linearization_constant( i );
         // alpha has to be referred to \Lambda
         std::vector<double> G1(NumVar);
         Range range = make_pair( 0, NumVar );
         v_c05f[ wFi ]->get_linearization_coefficients( G1.data() , range , i );
         Alfa1[ i ] = UpRifFi[ wFi ] - Alfa1[ i ]
		          		  - std::inner_product( Lambda.begin() , Lambda.end() , G1.data() , double(0) );
         Master->ChgAlfa( Alfa1.data() , wFi );
         }
       break;
      case( C05FunctionMod::AllEntriesChanged ):
       Master->ChgSubG( tmod->range().first , tmod->range().second , wFi+1 );
       break;
      } // switch( tmod->f_type )
     }  // end  if( tmod )
    } // end C05FunctionModRngd  - - - - - - - - - - - - - - - - - - - - - - -

   // C05FunctionModSbst - - - - - - - - - - - - - - - - - - - - - - - - - - -
   {
    const auto tmod = std::dynamic_pointer_cast<C05FunctionModSbst>( mod );
    if( tmod ) {
     Index wFi = get_index_of_component(tmod->function());
     std::vector<double> Alfa1;
     switch( tmod->type() ) {
      case( C05FunctionMod::AllLinearizationChanged ):
       Master->ChgSubG( v_c05f[wFi]->is_active(tmod->vars()[0]) ,
               v_c05f[wFi]->is_active(tmod->vars()[tmod->vars().size()-1]) , wFi+1 );
      case( C05FunctionMod::AlphaChanged ):
       // a finite f_shift should be treated in a different way but
       // as of now the finite shifts are ignored by MPSolver
       Alfa1.resize( Master->MaxName(wFi+1) );
       for( Index i = 0 ; i < Master->MaxName(wFi+1) ; i++ )
        if( Master->WComponent( i ) == wFi+1 ) {
         Alfa1[ i ] = v_c05f[wFi]->get_linearization_constant( i );
         // alpha has to be referred to \Lambda
         std::vector<double> G1(NumVar);
         Range range = make_pair( 0, NumVar );
         v_c05f[ wFi ]->get_linearization_coefficients( G1.data() , range , i );
         Alfa1[ i ] = UpRifFi[ wFi ] - Alfa1[ i ]
		          		  - std::inner_product( Lambda.begin() , Lambda.end() , G1.data() , double(0) );
         Master->ChgAlfa( Alfa1.data() , wFi );
         }
       break;
      case( C05FunctionMod::AllEntriesChanged ):
       Master->ChgSubG( v_c05f[wFi]->is_active(tmod->vars()[0]) ,
               v_c05f[wFi]->is_active(tmod->vars()[tmod->vars().size()-1]) , wFi+1 );
       break;
      } // switch( tmod->f_type )
     }  // end  if( tmod )
    } // end C05FunctionModSbst  - - - - - - - - - - - - - - - - - - - - - - -

   // C05FunctionModLin  - - - - - - - - - - - - - - - - - - - - - - - - - - -
   {
    const auto tmod = std::dynamic_pointer_cast<C05FunctionModLin>( mod );
    if( tmod ) {
     Index wFi = get_index_of_component(tmod->function());
     if( wFi == Inf<Index>() ) { // 0th component
      double * G1 = Master->GetItem( 0 );
      f_lf->get_linearization_coefficients( G1 );
      for( Index i = 0 ; i < tmod->vars().size() ; ++i )
       G1[ f_lf->is_active(tmod->vars()[i]) ] +=
        tmod->delta()[ i ];
      const Index* SGBse = nullptr;
      Master->SetItemBse( SGBse , NumVar );
      Master->SetItem( InINF );
      } // end if( f_lf )
     else
      for( Index i = 0 ; i < Master->MaxName(wFi+1) ; i++ )
       if( Master->WComponent( i ) == wFi+1 ) {
        std::vector<double> G1(NumVar);
        Range range = make_pair( 0, NumVar );
        v_c05f[ wFi ]->get_linearization_coefficients( G1.data() , range , i );
        for( Index i = 0 ; i < tmod->vars().size() ; ++i )
          G1[ v_c05f[ wFi ]->is_active(tmod->vars()[i]) ] +=
           tmod->delta()[ i ];
        throw( std::logic_error( "expected to be completed" ) );
        }
     } // end  if( tmod )
    } // end C05FunctionModLin   - - - - - - - - - - - - - - - - - - - - - - -

   };  // end( guts_of_poM ) - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // finally, call the "guts of" - - - - - - - - - - - - - - - - - - - - - - -

  guts_of_poM( mod );  // now the actual call

  v_mod.pop_front();   // now the Modification is processed: remove it

  }  // end( while( there are Modification ) )
 }  // end( BundleSolver::process_outstanding_Modification ) - - - - - - - - -

/*--------------------------------------------------------------------------*/

BundleSolver::FakeFiOracle::FakeFiOracle( BundleSolver *solver ) : FiOracle()
{
 bslv = solver;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetNDOSolver( NDOSolver * NwSlvr )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetFiLog( ostream *outs , const char lvl )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetFiTime( const bool TimeIt )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetMaxName( cIndex MxNme )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetNumVar( void ) const
{
 return( bslv->NumVar );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetNrFi( void ) const
{
 return( bslv->v_c05f.size() );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetMaxName( void ) const
{
 return( bslv->BPar2 * ( bslv->NrFi - bslv->NrEasy ) );
 }

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetMinusInfinity( void )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetMaxNZ( cIndex wFi ) const
{
 if( wFi != Inf<Index>() )
  throw( std::logic_error( "GetMaxNZ can be called with wFi = Inf only" ) );

 return( bslv->NumVar );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetMaxCNZ( cIndex wFi ) const
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

bool BundleSolver::FakeFiOracle::GetUC( cIndex i )
{
 double lb_value = bslv->LamVcblr[ i ]->get_lb();
 if( lb_value == -Inf<ColVariable::VarValue>() )
  return( true );

 if( lb_value != ColVariable::VarValue( 0 ) )
  throw( std::logic_error( "any value different from zero is not allowed" ) );

 return( false );
 }

/*--------------------------------------------------------------------------*/

LMNum BundleSolver::FakeFiOracle::GetUB( cIndex i )
{
 return( bslv->LamVcblr[ i ]->get_ub() );
 }

/*--------------------------------------------------------------------------*/

LMNum BundleSolver::FakeFiOracle::GetBndEps( void )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetGlobalLipschitz( cIndex wFi )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetBNC( cIndex wFi )
{
 if( bslv->NrEasy && bslv->IsEasy[ wFi - 1 ] )
  return( bslv->MILP_s[ wFi - 1 ]->get_numcols() );
 else
  return( 0 );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetBNR( cIndex wFi )
{
 return( bslv->MILP_s[ wFi -1 ]->get_numrows() );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetBNZ( cIndex wFi )
{
 return( bslv->MILP_s[ wFi -1 ]->get_nzelements() );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::GetBDesc( cIndex wFi , int *Bbeg ,
					   int *Bind , double *Bval ,
					   double *lhs , double *rhs ,
					   double *cst , double *lbd ,
					   double *ubd )
{
 auto MILPSlv = bslv->MILP_s[wFi-1];

 int num_col = MILPSlv->get_numcols();

 std::copy( MILPSlv->get_matbeg().begin() , MILPSlv->get_matbeg().end() ,
	    Bbeg );
 std::copy( MILPSlv->get_matind().begin() , MILPSlv->get_matind().end() ,
	    Bind );
 std::copy( MILPSlv->get_matval().begin() , MILPSlv->get_matval().end() ,
	    Bval );

 std::copy( MILPSlv->get_lb().begin() , MILPSlv->get_lb().end() , lbd );
 std::copy( MILPSlv->get_ub().begin() , MILPSlv->get_ub().end() , ubd );

 for( Index i = 0 ; i < num_col ; i++ )
  if( MILPSlv->get_sense()[ i ] == 'L' ) {
   rhs[ i ] = MILPSlv->get_rhs()[ i ];
   lhs[ i ] = -Inf<double>();
   }
  else
   if( MILPSlv->get_sense()[ i ] == 'E' ) {
    rhs[ i ] = MILPSlv->get_rhs()[ i ];
    lhs[ i ] = MILPSlv->get_rhs()[ i ];
    }
   else
    if( bslv->MILP_s[ wFi -1 ]->get_sense()[ i ] == 'G' ) {
     rhs[ i ] = Inf<double>();
     lhs[ i ] = MILPSlv->get_rhs()[ i ];
     }
    else {
     double rngval = MILPSlv->get_rngval()[ i ];
     if( rngval > double(0) ) {
      rhs[ i ] = MILPSlv->get_rhs()[ i ] + rngval;
      lhs[ i ] = MILPSlv->get_rhs()[ i ];
      }
     else {
      rhs[ i ] = MILPSlv->get_rhs()[ i ];
      lhs[ i ] = MILPSlv->get_rhs()[ i ] + rngval;
      }
     }

 } // end ( BundleSolver::FakeFiOracle::GetBDesc() ) - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetANZ( cIndex wFi ,
					  cIndex strt , Index stp )
{
 if( ! bslv->IsEasy[ wFi - 1 ] )
  throw( std::logic_error( "the Function is not a Lagrangian one" ) );

 auto LagB = static_cast<LagBFunction *>( bslv->v_c05f[ wFi - 1 ] );
 //?? DA FARE return( LagB->get_Amat_nzelements() );

 } // end ( BundleSolver::FakeFiOracle::GetANZ() ) - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::GetADesc( cIndex wFi , int *Abeg , int *Aind ,
					   double *Aval , cIndex strt ,
					   Index stp )
{
 if( ! bslv->IsEasy[ wFi - 1 ] )
  throw( std::logic_error( "the Function is not a Lagrangian one" ) );

 auto LagB = static_cast<LagBFunction *>( bslv->v_c05f[ wFi - 1 ] );
 //?? DA FARE LagB->get_Amat_desc( Abeg , Aind , Aval , strt , stp );

 } // end ( BundleSolver::FakeFiOracle::GetANZ() ) - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

NDOSolver * BundleSolver::FakeFiOracle::GetNDOSolver( void )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SETTING LAMBDA ------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetLambda( cLMRow Lmbd )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

void BundleSolver::FakeFiOracle::SetLamBase( cIndex_Set LmbdB  ,
					     cIndex LmbdBD )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

bool BundleSolver::FakeFiOracle::SetPrecision( HpNum Eps )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR COMPUTING Fi() ----------------------*/
/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::Fi( cIndex wFi )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/

bool BundleSolver::FakeFiOracle::NewGi( cIndex wFi )
{
 if( wFi == 0 )
  throw( std::invalid_argument( "asking for the 0th component" ) );
 last_c05 =  wFi - 1;
 return( true );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetGi( SgRow SubG , cIndex_Set &SGBse ,
					 cIndex Name , cIndex strt , Index stp
					 )
{
 auto range = make_pair( strt , stp );
 bslv->v_c05f[ std::get< 1 >( GiNameVcblr[ Name ] ) ]->
  get_linearization_coefficients( SubG , range ,
				  std::get< 0 >( GiNameVcblr[ Name ] ) );
 SGBse = nullptr;
 return( stp - strt );
 }

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetVal( cIndex Name )
{
 return( bslv->v_c05f[ std::get< 1 >( GiNameVcblr[ Name ] ) ]->
	 get_linearization_constant( std::get< 0 >( GiNameVcblr[ Name ] ) ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetGiName( cIndex Name )
{
 auto it = GiNameVcblr.begin();
 for( ; it != GiNameVcblr.end() ; ++it  )
  if( std::get<1>( *it ) == last_c05 && std::get<2>( *it ) == true ) {
   std::get<2>( *it ) = false;
   break;
   }

 if( it == GiNameVcblr.end() )
  throw( std::invalid_argument( "the global pool is full" ) );

 bslv->v_c05f[ last_c05 ]->store_linearization( std::get< 0 >( *it ) );

 } // end ( FakeFiOracle::SetGiName( ) ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetLowerBound( cIndex wFi )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

FiOracle::FiStatus BundleSolver::FakeFiOracle::GetFiStatus( Index wFi )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::Deleted( cIndex i )
{
 bslv->v_c05f[ std::get< 1 >( GiNameVcblr[ i ] ) ]->
  delete_linearization( std::get< 0 >( GiNameVcblr[ i ] ) );

 std::get< 2 >( GiNameVcblr[ i ] ) = true;
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::Aggregate( cHpRow Mlt , cIndex_Set NmSt ,
					    cIndex Dm , cIndex NwNm )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::initialize( void )
{
 GiNameVcblr.resize( GetMaxName() );
 auto it =  GiNameVcblr.begin();
 for( Index i = 0 ; i < bslv->v_c05f.size() ; ++i )
  for( Index j = 0 ; j < bslv->BPar2 ; ++j )
   *it = std::make_tuple( j , i , true );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File BundleSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
