/*--------------------------------------------------------------------------*/
/*------------------------ File BundleSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BunldeSolver class.
 *
 * \version 0.13
 *
 * \date 25 - 04 - 2020
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
 * Copyright &copy 2019 - 2020 by Antonio Frangioni, Enrico Gorgone
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
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define BLOG( l , x ) if( f_log && ( LogVerb > l ) ) *f_log << x

#define BLOG2( l , c , x ) if( f_log && ( LogVerb > l ) && c ) *f_log << x

#define USE_MPTESTER 0

// if USE_MPTESTER is nonzero, the MPSolver is a MPTester whose master is
// an OSIMPSolver and whose slave is a QPPenaltyMP

#if USE_MPTESTER
 #include "MPTester.h"
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

void Compact( BundleSolver::Vec_VarValue & g , BundleSolver::Subset & B )
{
 // takes a "dense" n-vector g and "compacts" it deleting the elements whose
 // indices are in B; all elements of B must be in the range 0 .. n, B must
 // be ordered in increasing sense
 // the remaining entries in g are shifted left of the minimum possible
 // amount in order to fill the holes left by the deleted ones
 // g is *not* resized in here

 auto Bit = B.begin();
 auto i = *(Bit++);
 auto git = g.begin() + (i++);

 for( ; Bit != B.end() ; ++i ) {
  auto h = *(Bit++);
  while( i < h )
   *(git++) = g[ i++ ];
  }

 std::copy( g.begin() + i , g.end() , git );

 }  // end( Compact )

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
 "intBPar7" ,
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
 { "intBPar7" , BundleSolver::intBPar7 } ,
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
  1 ,  // intBPar7
  0 ,  // intMnSSC
  3 ,  // intMnNSC
 12 ,  // inttSPar1
  2 ,  // intMaxNrEvls
  0 ,  // intMPName
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
/*----------------------- METHODS OF BundleSolver --------------------------*/
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
 // v_mod is atomically copied in a temporary data structure to be processed,
 // but while the latter happens new Modification may come in; hence,
 // process_outstanding_Modification() may be called more than once

 while( ! v_mod.empty() )
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

  // some log - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log1();

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( IsOptimal() ) {
   BLOG( 1 , " ~ stop (optimal)" << std::endl );
   Result = kOK;
   break;
   }

  // check if "ex-ante" Noise Reduction is needed - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ensure that the \sigma* is "not too negative", if it is increase t (if
  // possible) and re-solve the MP; note that this kind of NR only happens if
  // the oracle is "unfaithful", i.e., it pretends to provide information with
  // the required accuracy but in fact it does not
  //
  // do the check only if Sigma is "negative enough to matter", otherwise do
  // not even call ReadDStart()
  //
  // however, avoid doing any of this if Fi( Lambda ) is defined, because if
  // not then negative linearization errors are "normal" (the reference value
  // is "random" and there is no reason to believe it's >= than the true value)

  if( ( UpFiLmb[ NrFi ] < Inf<double>() ) &&
      ( Sigma < - max_error( UpRifFi[ NrFi ] , RelAcc ) ) &&
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

  // Hard Long-Term t-strategy for quadratic stabilization- - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the hard long-term t-strategy requires t to increase if the step is too
  // small, and therefore has to be checked before the others
  // however, it is only viable under a quadratic stabilization

  if( ( ( ! ( MPName & 1 ) ) || ( MPName & 4 ) ) &&
      ( tStar > 0 ) && ( ( tSPar1 & tSP1Msk ) == kHLTTS ) &&
      ( UpFiLmb[ NrFi ] < Inf<double>() ) ) {

   double AFL = std::abs( UpFiLmb[ NrFi ] );
   if( AFL < 1 )
    AFL = 1;

   if( abs( vStar[ NrFi ] ) <= tSPar2 * EpsU * AFL ) {
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
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

  // compute the "aggregated" Alfa1 and ScPr1 - - - - - - - - - - - - - - - -
  // ... using the "representatives" of all components: these are used in
  // some "global" formulae, such as the t heuristics

  Alfa1[ NrFi ] = 0;
  ScPr1[ NrFi ] = Master->ReadGid();

  for( Index k = 0 ; k < NrFi ; k++ )
   if( whisG1[ k ] < InINF ) {
    if( Alfa1[ k ] == Inf<double>() )
     Alfa1[ k ] = (Master->ReadLinErr())[ whisG1[ k ] ];

    Alfa1[ NrFi ] += Alfa1[ k ];

    if( ScPr1[ k ] == Inf<double>() )
     ScPr1[ k ] = Master->ReadGid( whisG1[ k ] );

    ScPr1[ NrFi ] += ScPr1[ k ];
    }
   else
    Alfa1[ k ] = ScPr1[ k ] = 0;

  // check whether the Lower Bounds have changed- - - - - - - - - - - - - - -

  UpdtLowerBound();

  // some log about the newly obtained information- - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log2();

  // check whether either any error has occurred or time has expired- - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1[NrFi] == - Inf<double>() ) {
   BLOG( 1 , " ~ stop (unbounded)" << std::endl );
   Result = kUnbounded;
   break;
   }

  if( ( Result == kError ) || ( Result == kStopTime ) ) {
   BLOG2( 1 , ( Result == kError ) , " ~ stop (error)" << std::endl );
   BLOG2( 1 , ( Result == kStopTime ) , " ~ stop (time)" << std::endl );
   break;
   }

  if( tHasChgd ) {  // "noise reduction": t has changed
   BLOG( 1 , " ~ NR" << std::endl );  // so go solve the master problem again
   continue;                          // (no NS/SS decision can be made)
   }

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

  if( UpFiLmb1[ NrFi ] == Inf<double>() )  // ???
   continue;
  else
   if( UpFiLmb[ NrFi ] == Inf<double>() ) {  // if reached feasibility- - - -
    BLOG( 1 , "           Fi1 < INF ==> SS " << std::endl );
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

 //!! PrintBundle();

 return( Result );

 }  // end( BundleSolver::compute() ) - - - - - - - - - - - - - - - - - - - -

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
    the second case a FRealObjective one whose the function is a
    LinearFunction and having as many sub-blocks as the number of components.
    In the latter case, each sub-block must not contain any Variable or
    Constraint.cVariable may have a lower and upper bound. If the lower
    bound  has a finite value, it must be 0. */

 if( f_Block->get_nested_Blocks().empty() ) {
  // the objective function of the block must be a C05Function  - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
  if( ! obj )
   throw( std::logic_error( "objective is not a FRealObjective" ) );

  auto c05f = dynamic_cast< C05Function * >( (obj)->get_function() );
  if( ! c05f )
   throw( std::logic_error( "the objective is not a C05Function" ) );

  v_c05f.push_back( c05f );
  f_lf = nullptr;
  }
 else {
  // the objective function of the block must be a LinearFunction- - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ! f_Block->get_objective() )  // there is no Objective
   f_lf =  nullptr;
  else {
   auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
   if( ! obj )
    throw( std::logic_error( "the objective is not a real function" ) );

   if( ! obj->get_function() ) // the FRealObjective has no Function
    f_lf = nullptr;
   else {
    f_lf = dynamic_cast<LinearFunction *>( obj->get_function() );
    if( ! f_lf )
     throw( std::logic_error( "the objective is not a LinearFunction" ) );

    if( ! f_lf->get_num_active_var() )  // the LinarFunction has no Variable
     f_lf = nullptr;
    }
   }

  auto sb = f_Block->get_nested_Blocks();
  v_c05f.resize( sb.size() );

  for( Index i = 0 ; i < sb.size() ; ++i ) {  // for each sub-block

   // the objective function of each sub-block must be a C05Function - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto obj = dynamic_cast< FRealObjective * >( sb[ i ]->get_objective() );
   if( ! obj )
    throw( std::logic_error( "the objective is not a real function" ) );

   auto c05f = dynamic_cast<C05Function *>( (obj)->get_function() );
   if( ! c05f )
    throw( std::logic_error( "the objective is not a C05Function" ) );
   v_c05f[ i ] = c05f;

   // nephews are not allowed- - - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( sb[ i ]->get_nested_Blocks().size() )
    throw( std::logic_error( "nephew are not allowed" ) );

   // Variable of sub-Block are not expected, neither are Constraint - - - - -
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

  // the set of "active" Variable in all Function must be the same - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( f_lf ) {
   if( f_lf->get_num_active_var() != v_c05f[ 0 ]->get_num_active_var() )
    throw( std::logic_error( "the list of active Variable do not match" ) );

   auto vi = v_c05f[ 0 ]->begin();
   for( auto v : *f_lf )
    if( v != *(vi++) ) 
     throw( std::logic_error( "the list of active Variable do not match" ) );
   }

  for( Index i = 1 ; i < sb.size() ; ++i ) {
   if( v_c05f[ i - 1 ]->get_num_active_var() !=
       v_c05f[ i ]->get_num_active_var() )
    throw( std::logic_error( "the list of active Variable do not match" ) );

   auto vi = v_c05f[ i ]->begin();
   for( auto v : *v_c05f[ i - 1 ] )
    if( v != *(vi++) ) 
     throw( std::logic_error( "the list of active Variable do not match" ) );
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

 // set the global pool size to all non-easy functions - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 vBPar2.resize( NrFi + 1, 0 );
 if( NrEasy )
  for( Index i = 0 ; i < NrFi ; ++i ) {
   if( IsEasy[ i ] )
    continue;
   auto gps = v_c05f[ i ]->get_int_par( C05Function::intGPMaxSz );
   if( BPar2 == 0 ) {
    vBPar2[ NrFi ] += gps;
    vBPar2[ i ] = gps;
    }
   else {
    if( gps < BPar2 )
     v_c05f[ i ]->set_par( C05Function::intGPMaxSz , BPar2 );
    vBPar2[ NrFi ] += gps;
    vBPar2[ i ] = BPar2;
    }
   }
 else
  for( Index i = 0 ; i < NrFi ; ++i ) {
   auto gps = v_c05f[ i ]->get_int_par( C05Function::intGPMaxSz );
   if( BPar2 == 0 ) {
    vBPar2[ NrFi ] += gps;
    vBPar2[ i ] = gps;
    }
   else {
    if( gps < BPar2 )
     v_c05f[ i ]->set_par( C05Function::intGPMaxSz , BPar2 );
    vBPar2[ NrFi ] += BPar2;
    vBPar2[ i ] = BPar2;
    }
   }

 // allocate memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 t = tInit;
 Prevt = Inf<double>();

 Lambda.resize( NumVar );    // the default starting point
 Lambda1.resize( NumVar );   // the tentative point

 if( MaxSol > 1 )  // best point found so far
  LmbdBst.resize( NumVar );

 OOBase.resize( vBPar2[ NrFi ] , Inf<SIndex>() );
 // counter for eliminating outdated items: Inf<SIndex>() means empty

 ItemVcblr.resize( vBPar2[ NrFi ] );
 for( Index i = 0 ; i < vBPar2[ NrFi ] ; i++ )
  ItemVcblr[ i ] = make_pair( InINF , Inf<SIndex>() );

 NrItems.resize( NrFi , 0 );
 DFItems.resize( NrFi , 0 );
 NFItems.resize( NrFi , 0 );

 FreList = priority_queue<Index>();     // list of free bundle slots
 whisZ.resize( NrFi );  // for each component, the name of its "Z" if it is
                        // in the bunlde

 FiStatus.resize( NrFi , kUnEval );
 LowerBound = -Inf<double>();     // globaò lower bounds
 TrueLB = false;

 UpFiBest = Inf<VarValue>();      // best, ...
 UpRifFi.resize( NrFi + 1 , 0 );  // and reference Fi() values
 UpFiLmb1.resize( NrFi + 1 );     // upper and lower function value
 LwFiLmb1.resize( NrFi + 1 );     // ... at the tentative point

 UpFiLmb.resize( NrFi + 1 , Inf<VarValue>() );   // upper 
 LwFiLmb.resize( NrFi + 1 , -Inf<VarValue>() );  // ... and lower Fi-value
                                                 // ... at the current point

 whisG1.resize( NrFi , InINF );  // no representative yet

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

 #if( ! USE_MPTESTER )
  if( MPName & 1 ) {  // the MPSolver is a OSIMPSolver
 #endif
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
 #if( ! USE_MPTESTER )
   }
  else {  // the MPSolver is a QPPenaltyMP
 #endif
   QPPenaltyMP *qp = new QPPenaltyMP();
   qp->SetPricing( CtOff );
   qp->SetMaxVarAdd( MxAdd );
   qp->SetMaxVarRmv( MxRmv );
 #if( USE_MPTESTER )
   Master = new MPTester( Master , qp );
 #else
   Master = qp;
   }
 #endif

 InitMP();

 // reset algorithm  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // default value  RstCrr | RstSbg | RstCnt | RstFiV | NoStPt = 62
 ReSetAlg( RstAlgPrm );  // Fi( Lambda ) is reset inside

 }  // end( BundleSolver::set_Block )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( const idx_type par , const int value )
{
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
  case( intBPar7 ):
   BPar7 = value;
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
   if( ( value < 0 ) || ( value > 15 ) )
    throw( std::invalid_argument( "MPName must be in [0, 15]" ) );
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
  case( intBPar7 ):
   return( BPar7 );
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
  if( TrueLB && ( LowerBound > - Inf<double>() ) )
   Master->SetLowerBound( LowerBound - UpFiLmb[ NrFi ] );
  else
   Master->SetLowerBound( - Inf<double>() );

  if( MPName & 1 )  // QPPenaltyMP does not allow individual lower bounds
   for( Index k = 0 ; k < NrFi ; k++ ) {
    if( NrEasy && IsEasy[ k ] )  // skip easy components
     continue;

    if( TrueLB && ( LowerBound > - Inf<double>() ) )
     Master->SetLowerBound( LowerBound - UpFiLmb[ k ] , k + 1 );
    else
     Master->SetLowerBound( - Inf<double>() , k + 1 );
    }

  LBHasChgd = false;
  }

 /* set termination criterion - - - - - - - - - - - - - - - - - - - - - - - -

 if( UpFiLmb[ NrFi ] < Inf<double>() )
  Master->SetPar( MPSolver::kZero ,
		  max_error() / std::max( tStar / t , HpNum( 1 ) ) );
 */
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

  if( UpFiLmb[ NrFi ] < Inf<double>() )
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
 whisZ.assign( NrFi , InINF );

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

 // Lambda has changed, pass the new one to the oracle - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiStatus.assign( NrFi , kUnEval );
 whisG1.assign( NrFi , InINF );

 for( Index i = 0 ; i < NumVar ; i++ )
  LamVcblr[ i ]->set_value( Lambda1[ i ] );

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

   c_VarValue Lk = v_c05f[ k ]->get_Lipschitz_constant();
   if( ( Lk < Inf<VarValue>() ) && ( UpFiLmb[ k ] < Inf<VarValue>() ) )
    UpFiLmb1[ k ] = UpRifFi[ k ] + Lk * NrmD;
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

 LwFiK = UpRifFi[ wFi ] + vStar[ wFi ];

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
  EpsCurr = ( UpCutOff - LwCutOff ) / std::max( 1.0 ,
						std::abs( UpRifFi[ wFi ] ) );
 else
  EpsCurr = RelAcc / Nearly;

 // assign the cutoff values to the C05Function - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto fwFi = v_c05f[ wFi ];

 fwFi->set_par( dblUpCutOff , UpCutOff );
 fwFi->set_par( dblLwCutOff , LwCutOff );
 fwFi->set_par( dblRelAcc , EpsCurr );

 // now compute the C05Function and retrieve upper and lower estimates- - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiStatus[ wFi ] = fwFi->compute( ( FiStatus[ wFi ] == kUnEval ) );

 if( UpFiLmb1[ NrFi ] < Inf<VarValue>() )
  UpFiLmb1[ NrFi ] -= UpFiLmb1[ wFi ];

 if( LwFiLmb1[ NrFi ] > -Inf<VarValue>() )
  LwFiLmb1[ NrFi ] -= LwFiLmb1[ wFi ];

 UpFiLmb1[ wFi ] = std::min( fwFi->get_upper_estimate() , UpFiLmb1[ wFi ] );
 LwFiLmb1[ wFi ] = std::max( fwFi->get_lower_estimate() , LwFiLmb1[ wFi ] );

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

 // get new linearizations- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool HasLinearization;
 bool diagonal;

 for( Index Ftchd = 0 ; Ftchd < aBP3 ; ) {
  diagonal = true;

  if( Ftchd == 0 ) {
   // first look for a constraint then for a sub-gradient

   if( UpFiLmb1[ wFi ] == Inf<VarValue>() ) {
    HasLinearization = fwFi->has_linearization( diagonal = false );
    if( ! HasLinearization )
     HasLinearization = fwFi->has_linearization( diagonal );
    }
   else
    HasLinearization = fwFi->has_linearization( diagonal );
   }
  else {
   if( UpFiLmb1[ wFi ] == Inf<VarValue>() ) {
    HasLinearization = fwFi->compute_new_linearization( diagonal = false );
    if( ! HasLinearization )
     HasLinearization = fwFi->compute_new_linearization( diagonal );
    }
   else
    HasLinearization = fwFi->compute_new_linearization( diagonal );
   }

  if( ! HasLinearization )
   break;

  // check if aggregation has to be performed - - - - - - - - - - - - - - - -
  // doing this now could occasionally result in useless aggregations, but it
  // avoids complications in the interface of MPSolver (inserting some
  // Z[ wFi ] while inserting the new item)

  Index wh = BStrategy( wFi );

  // get the space for the item from the MPSolver - - - - - - - - - - - - - -

  double *G1 = Master->GetItem( wFi + 1 );

  // fetch the item from the Oracle - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set SGBse = nullptr;
  fwFi->get_linearization_coefficients( G1 );
  auto eps = fwFi->get_linearization_constant();

  GiEvaltns++;

  // pass the base to the MP Solver - - - - - - - - - - - - - - - - - - - - -

  Master->SetItemBse( SGBse , NumVar );

  // calculate ScPr1k and Alfa1k- - - - - - - - - - - - - - - - - - - - - - -

  Index cp;
  HpNum ScPr1k;

  // update alpha value at Lambda1 point- - - - - - - - - - - - - - - - - - -

  eps = UpFiLmb1[ wFi ] - eps -
   std::inner_product( Lambda1.begin() , Lambda1.end() , G1 , double( 0 ) );

  HpNum Alfa1k = eps;

  if( ! diagonal )                 // it is a constraint
   cp = Master->CheckCnst( Alfa1k , ScPr1k , Lambda.data() );
  else                             // it is a subgradient
   cp = Master->CheckSubG( UpFiLmb1[ wFi ] - UpRifFi[ wFi ] ,
			   t , Alfa1k , ScPr1k );

  if( cp < InINF ) {  // the item is a copy- - - - - - - - - - - - - -
   BLOG( 2 , std::endl << "            New " );
   BLOG2( 2 , diagonal , "subgradient" );
   BLOG2( 2 , ! diagonal , "constraint" );
   BLOG( 2 , " for Fi[ " << wFi << " ] is a copy of " << cp );

   cHpNum OrigA1k = (Master->ReadLinErr())[ cp ];

   assert( ItemVcblr[cp].first != InINF &&
           ItemVcblr[cp].second < vBPar2[ItemVcblr[cp].first] &&
  		   ItemVcblr[cp].second >= 0 );

   if( OrigA1k > Alfa1k ) {   // the copy has smaller Alfa than the original
    BLOG( 2 , " with smaller Alfa" );
    Master->SubstItem( wh = cp );  // substitute it
    if( BPar7 ) {
     ItemVcblr[ cp ].second += vBPar2[ wFi ]; // the linearization ItemVcblr[ cp ].second
     DFItems[ wFi ]++;                        // is free in the global pool
     }
    else {
     ItemVcblr[ cp ].second -= vBPar2[ wFi ]; // the linearization ItemVcblr[ cp ].second
     NFItems[ wFi ]++;  // could be assigned if there is no more free items
     }
    NrItems[ wFi ]--;
    }
   else {
	wh = InINF;    // otherwise, nothing new has happened
	Ftchd++;       // anyhow, this counts as a new item
    }
   }
  else {               // insert the item, if there is space - - - - - - - - -
   if( wh < InINF ) {  // someone has been selected in BStrategy()
    Master->RmvItem( wh ); // remove it from the MP
    if( BPar7 ) {
     ItemVcblr[ wh ].second += vBPar2[ wFi ]; // the linearization ItemVcblr[ cp ].second
     DFItems[ wFi ]++;                        // is free in the global pool
     }
    else {
     ItemVcblr[ wh ].second -= vBPar2[ wFi ]; // the linearization ItemVcblr[ cp ].second
     NFItems[ wFi ]++;  // could be assigned if there is no more free items
     }
    NrItems[ wFi ]--;
    }
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

   if( f_log && ( LogVerb > 2 ) ) {
    *f_log << std::endl << "            New ";
    if( ! diagonal )
     *f_log << "constraint " << wh << " ~ rhs = " << Alfa1k;
    else
     *f_log << "eps-subgradient " << wh << " for Fi[ " << wFi << " ] ~ eps = "
	    << eps << " ~ Alfa1 = " << Alfa1k << " ~ gd = " << - ScPr1k;
    }
   }

  // if something was inserted, bookkeeping is needed - - - - - - - - - - - -

  if( wh < InINF ) {
   Ftchd++;              // one more item

   SetItemName( wFi , wh );
   inhibit_Modification( true );
   v_c05f[ wFi ]->store_linearization( ItemVcblr[wh].second );
   inhibit_Modification( false );

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
  }  // end( items-collecting loop )- - - - - - - - - - - - - - - - - - - - -

 // update lower and upper estimates  - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( ( LwFiLmb1[ NrFi ] > LwTrgt ) || ( UpFiLmb1[ NrFi ] < UpTrgt ) );

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
   Delete( h );
   }
 }  // end( BundleSolver::SimpleBStrat ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdtLowerBound( void )
{
 // note: set LBHasChgd only if the new lower bound has to be set into the
 //       MPSolver, which only happens if a new non-conditional lower bound
 //       if found or if a previously non-conditional lower bound disappears
 //       and only leaves a conditional one (in the latter case, the lower
 //       bound in the MPSolver has to be set to -INF)

 // first of all, check if a "hard" lower bound is available
 double LwrBnd = f_Block->get_valid_lower_bound( false );
 if( LwrBnd > - Inf<double>() ) {
  LBHasChgd = ( ! TrueLB ) || ( LwrBnd != LowerBound );
  LowerBound = LwrBnd;
  TrueLB = true;
  }
 else {
  // if not, check if at least a "conditional" one is available
  LwrBnd = f_Block->get_valid_lower_bound( true );
  LBHasChgd = TrueLB;
  LowerBound = LwrBnd;
  TrueLB = false;
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

 *f_log << std::endl << "            ";

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

 Master->SetDim( vBPar2[ NrFi ] , &FakeFi , false );

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

 tHasChgd = true;
 LBHasChgd = false;

 if( MPName & 8 )
  Master->CheckIdentical();
 
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

 if( NrItems[ wFi ] < vBPar2[ wFi ] )
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
 for( Index i = 0 ; ++i < vBPar2[ NrFi ] ; )
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

  // throw( std::logic_error( "Fre List" ) );

  }
 else                                       // there are no deleted items ...
  if( Master->MaxName() < Index( vBPar2[ NrFi ] ) )
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

 SetItemName( wFi , whr );
 inhibit_Modification( true );
 v_c05f[ wFi ]->store_combination_of_linearizations( coefficients , ItemVcblr[whr].second );
 inhibit_Modification( false );

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

 BLOG( 2 , std::endl << "Aggregation performed into " << whr );

 }  // end( BundleSolver::AggregateZ() ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::Heuristic1( void )
{
 if( Alfa1[ NrFi ] < Eps<double>() )
  return( DeltaFi > Eps<double>() ? tMaior : tMinor );
 else
  return( t * ( ( DeltaFi + Alfa1[ NrFi ] ) / ( 2 * Alfa1[ NrFi ] ) ) );

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

 FreList = priority_queue<Index>();

 ItemVcblr.resize( vBPar2[ NrFi ] );
  for( Index i = 0 ; i < vBPar2[ NrFi ] ; i++ )
   ItemVcblr[ i ] = make_pair( InINF , Inf<SIndex>() );

 NrItems.resize( NrFi , 0 );
 DFItems.resize( NrFi , 0 );
 NFItems.resize( NrFi , 0 );

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
 NrItems.clear();
 DFItems.clear();
 NFItems.clear();

 ItemVcblr.clear();
 LamVcblr.clear();

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
      if( Master->IsSubG( i ) )
       Delete( i );
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
       if( ! Master->IsSubG( i ) )
    	Delete( i );
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

void BundleSolver::Delete( cIndex i , bool ModDelete )
{
 cIndex wFi = ItemVcblr[ i ].first;

 if( Master ) {
  // check if this item was the "representative" for its component- - - - - -

  if( Master->IsSubG( i ) )  // it is a subgradient
   if( whisG1[ wFi ] == i )  // it is the representative of wFi
    whisG1[ wFi ] = InINF;   // a new representative is needed

  // delete the item with name `i' from the MP- - - - - - - - - - - - - - - -

  Master->RmvItem( i );
  }

 BLOG( 2 , std::endl << "Item " << i << " removed" );

 // bookkeeping of internal data structures - - - - - - - - - - - - - - - - -

 FreList.push( i );
 OOBase[ i ] = Inf<SIndex>();

 // because of the deletion, if the item is the last one
 // in the pool the updating of max item must be performed
 // in any case update the max item of the pool


 if( ModDelete ) {
  ItemVcblr[ i ].second += vBPar2[ wFi ];   // slot can be rewritten
  DFItems[ wFi ]++;
  }
 else {
  if( BPar7 ) {
   ItemVcblr[ i ].second += vBPar2[ wFi ];
   DFItems[ wFi ]++;
   }
  else {
   ItemVcblr[ i ].second -= vBPar2[ wFi ];
   NFItems[ wFi ]++;
   }
  }

 NrItems[ wFi ]--;

 // compacting FreList[] if it's too big- - - - - - - - - - - - - - - - - - -
 // remove from FreList[] every name >= Master->MaxName(); note that every
 // ordered set *is* a Heap. apart from efficiency reasons, this is
 // needed because Master->MaxName() - FreDim is the only way in which the
 // Bundle can compute the number of "live" items

 cIndex MxNm = Master->MaxName();
 if( FreList.size() > MxNm ) {
  FreList = priority_queue<Index>();
  for( Index j = 0 ; j < MxNm ; j++ )
   if( OOBase[ j ] == Inf<SIndex>() )
    FreList.push( j );
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

 if( vStar[ NrFi ] >= Inf< VarValue >() )  // there are no subgradients
  return( false );

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

void BundleSolver::FModChg( VarValue f_shift , Index wFi )
{
 if( f_shift == INFshift ) {           // function changed monotonically up
  UpFiLmb[ wFi ] = Inf<VarValue>();    // reset upper function values
  UpFiLmb[ NrFi ] = Inf<VarValue>();
  UpFiBest = Inf<VarValue>();          // comprised best one
  return;
  }

 if( f_shift == -INFshift ) {          // function changed monotonically dn
  LwFiLmb[ wFi ] = -Inf<VarValue>();   // reset lower function values
  LwFiLmb[ NrFi ] = -Inf<VarValue>();
  return;
  }

 if( std::isnan( f_shift ) ) {         // function changed unpredictably
  UpFiLmb[ wFi ] = Inf<VarValue>();    // reset both upper ...
  LwFiLmb[ wFi ] = -Inf<VarValue>();   // ... and lower function values
  UpFiLmb[ NrFi ] = Inf<VarValue>();
  LwFiLmb[ NrFi ] = -Inf<VarValue>();
  UpFiBest = Inf<VarValue>();          // and of course best one
  return;
  }

 // function changed by shift():  just update everything

 if( UpFiLmb[ wFi ] < Inf<VarValue>() )
  UpFiLmb[ wFi ] += f_shift;

 if( UpFiLmb[ NrFi ] < Inf<VarValue>() )
  UpFiLmb[ NrFi ] += f_shift;

 if( UpFiBest < Inf<VarValue>() )
  UpFiBest += f_shift;

 if( LwFiLmb[ wFi ] > -Inf<VarValue>() )
  LwFiLmb[ wFi ] += f_shift;

 if( LwFiLmb[ NrFi ] > -Inf<VarValue>() )
  LwFiLmb[ NrFi ] += f_shift;

 UpRifFi[ wFi ] += f_shift;
 UpRifFi[ NrFi ] += f_shift;

 } // end ( BundleSolver::FModChg )  - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::SetItemName( Index wFi , Index wh ) {

 // we distinguish three different type of *free* slots of the
 // global pool:
 //	- the free one (a), never assigned, roughly speaking it
 //   is one of the last free positions;
 // - used slot (b), its linearization has been removed
 //   or could be removed; that slot is in the set "DFItems"
 // - almost busy slot (c), its linearization could be removed
 //   but *only if* the slots (a) and (b) are not available;
 //   it belongs to the set "NFItems"

 // two cases can occur:
 // 1. BPar7 == 1, the search is in order performed among the
 //  slots of type b) and the ones of type a)
 // 2. BPar7 == 0, the search is in sequel performed among the
 //  slots of type b), the slots of type a) and finally
 //  the ones of type c)

 // the item wh cannot be currently used
 assert( ! (ItemVcblr[wh].first != InINF &&
         ItemVcblr[wh].second < vBPar2[ItemVcblr[wh].first] &&
		 ItemVcblr[wh].second >= 0 ) );

 if( DFItems[ wFi] ) {
  if( !( ItemVcblr[wh].first == wFi
	  && ItemVcblr[wh].second >= vBPar2[wFi]) ) {
   Index i = 0;
   for( ; i < vBPar2[ NrFi ] ; i++  )
    if( ItemVcblr[i].first == wFi
		&& ItemVcblr[i].second >= vBPar2[wFi] ) {
     swap( ItemVcblr[wh] , ItemVcblr[i] );

     break;
	 }
   assert( i != vBPar2[NrFi] );
   }
  ItemVcblr[wh].second = ItemVcblr[wh].second - vBPar2[wFi];
  DFItems[ wFi ]--;
  }
 else
  if( NrItems[ wFi ] + NFItems[ wFi] < vBPar2[wFi] ) {
   ItemVcblr[wh].first = wFi; // set component name
   ItemVcblr[wh].second = NrItems[ wFi ] + NFItems[ wFi];
   }
  else {
   // if BPar7 == 1 --> (NFItems==0)
   assert( NFItems[ wFi] != 0 );
   if( !( ItemVcblr[wh].first == wFi
	  && ItemVcblr[wh].second < 0 ) ) {
    Index i = 0;
    for( ; i < vBPar2[ NrFi ] ; i++  )
	 if( ItemVcblr[i].first == wFi
			&& ItemVcblr[i].second < 0 ) {
	  swap( ItemVcblr[wh] , ItemVcblr[i] );
      break;
      }
    assert( i != vBPar2[NrFi] );
    }
   ItemVcblr[wh].second = ItemVcblr[wh].second + vBPar2[wFi];
   NFItems[ wFi ]--;
   }

 NrItems[ wFi ]++; // update the number of item of wFi

 } // end ( BundleSolver::SetItemName )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::is_special_GroupMod( GroupModification & gmod )
{
 // recognise "special" GroupModification for changing the set of "active"
 // Variable of all the Objective at the same time; note that these
 // contain FunctionModVars* not necessarily C05FunctionModVars* because
 // the Modification may not be strongly quasi-additive

 Index nsm = gmod.v_sub_Modifications.size();
 if( nsm != NrFi + ( f_lf ? 1 : 0 ) )
  return( false );

 auto sm0 = gmod.v_sub_Modifications[ 0 ]
 for( Index i = 1 ; i < nsm ; ++i )
  if( typeid( sm0 ) != typeid( tmod.v_sub_Modifications[ i ] ) )
   return( false );

 // check FunctionModVarsAddd
 {
  const auto mod0 = std::dynamic_pointer_cast<FunctionModVarsAddd>( sm0 );
  if( mod0 ) {
   for( Index i = 1 ; i < nsm ; ++i ) {
    auto modi = std::static_pointer_cast<FunctionModVarsAddd>(
					    tmod.v_sub_Modifications[ i ] );
    if( ( mod0->first() != modi->first() ) ||
	( mod0->vars() != modi->vars() ) )
     throw( std::logic_error( "different Variable change in components" ) );
    }

   return( true );
   }
  }

 // check FunctionModVarsRngd
 {
  const auto mod0 = std::dynamic_pointer_cast<FunctionModVarsRngd>( sm0 );
  if( mod0 ) {
   for( Index i = 1 ; i < nsm ; ++i ) {
    auto modi = std::static_pointer_cast<FunctionModVarsRngd>(
					    tmod.v_sub_Modifications[ i ] );
    if( mod0->range() != modi->range() )
     throw( std::logic_error( "different Variable change in components" ) );
    }

   return( true );
   }
  }

 // check FunctionModVarsSbst
 {
  const auto mod0 = std::dynamic_pointer_cast<FunctionModVarsSbst>( sm0 );
  if( mod0 ) {
   for( Index i = 1 ; i < nsm ; ++i ) {
    auto modi = std::static_pointer_cast<FunctionModVarsSbst>(
					    tmod.v_sub_Modifications[ i ] );
    if( mod0->subset() != modi->subset() )
     throw( std::logic_error( "different Variable change in components" ) );
    }

   return( true );
   }
  }

 return( false );

 }  // end( BundleSolver::is_special_GroupMod )- - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::flatten_Modification_list( Lst_sp_Mod & vmt , sp_Mod mod )
{
 const auto tmod = std::dynamic_pointer_cast<GroupModification>( mod );
 if( tmod && ( ! is_special_GroupMod( *tmod ) ) )
  for( auto submod : tmod->v_sub_Modifications )
   flatten_Modification_list( vmt , submod );
 else
  vmt.push_back( mod );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::compute_inverse_dictionary( Dctnry & id )
{
 if( ! id.empty() )  // inverse dictionary computed already
  return;            // nothing to do

 // the inverse dictionary is a std::vector with one entry per component; the
 // entry has size NrItems[ h ], and NrItems[ h ][ i ] is the position in the
 // bundle of the i-th item in the h-th global pool, or Inf< Index > if the
 // i-th item in the global pool is not in the bundle

 // allocate memory
 id.resize( NrFi );
 for( Index h = 0 ; h < NrFi : ++h )
  id[ h ].resize( NrItems[ h ] , Inf< Index >() );

 // now construct the inverse vocabulary
 for( Index i = 0 ; i < Master->MaxName() ; ++i )
  if( ( ItemVcblr[ i ].second >= 0 ) &&
      ( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) )
   id[ ItemVcblr[ i ].first ][ ItemVcblr[ i ].second ] = i;

 }  // end( BundleSolver::compute_inverse_dictionary ) - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::process_outstanding_Modification( void )
{
 // multiple-loop version, where several passes are done in order to gather
 // which kind of Modification have occurred and avoid doing costly work
 // more than once
 //
 // All loops use a Lambda to define a "guts" of the method that can be
 // called recursively. Note the trick of defining the std::function object
 // and "passing" it to the lambda, which allows recursive calls. Note the
 // need to explicitly capture "this" to use fields/methods of the class.

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 0-th loop: "atomically flatten" v_mod into a temporary list to better
 // handle it, then clear it

 Lst_sp_Mod v_mod_tmp;

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;  // try to acquire lock, spin on failure

 for( auto mod : v_mod )
  flatten_Modification_list( v_mod_tmp , mod );

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the 1st loop is made in reverse, from the latest Modification to the
 // earlies, and does the following:
 // - reset/change the upper/lower bounds that need to (check all shift() of
 //   all *FunctionMod*
 // - check if some global pool has been "hard" reset, i.e., all the
 //   linearization in there have been deleted; this is brought about by
 //   a C05FunctionMod with type() == GlobalPoolRemoved and which().empty()
 //   or by any FunctionMod that does *not* imply strong quasi-additivity
 //   (i.e., that it is not a *C05*FunctionMod*)
 // - if a *FunctionMod* changing linearizations happens before a reset of the
 //   global pool (meaning it is found afterwards in the reverse order) it
 //   is deleted since it is useless (after having checked if it also impacts
 //   the upper/lower bounds)
 // - if there is more than one component, check that no "naked"
 //   *FunctionModVars* is there
 // - check that no ConstraintMod or VariableMod are there, since they are
 //   not handled (yet)
 //
 // TODO: during the 1st loop we could compute the set of components that have
 //       been modified anyhow and use this information to avoid constructing
 //       the numerous data structures like reset[] that are indexed over NrFi.
 //       this might be important if, say, NrFi is 10000 but only a smattering
 //       of the components (say, one) change

 std::vector<bool> reset( NrFi , false );

 for( auto rimod = v_mod_tmp.rbegin() ; rimod != v_mod_tmp.rend() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to rimod
      [ & to_delete , & v_mod_tmp ]( decltype( rimod ) & ri ) {
       if( to_delete )
	ri = decltype( ri )( v_mod_tmp.erase( std::next( ri ).base() ) );
       else
	++ri;
       }( rimod ) ) {
  bool to_delete = false;
  auto mod = *rimod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // first check if it is any kinf of FunctionMod, since this gives immediate
  // access to the component, and any FunctionMod pertaining to an already
  // reset component can be almost immediately deleted
  {
   const auto tmod = std::dynamic_pointer_cast<FunctionMod>( mod );
   if( tmod ) {
    auto wFi = get_index_of_component( tmod->function() );

    // adjust or reset upper/lower values as needed
    // note that the list is scanned in reverse, hence these changes are
    // applied in reverse order. however, if the upper/lower values are
    // reset at any point in the list they stay reset forever. indeed,
    // even if a function has a finite shift after a reset, this says
    // nothing because there are no known values to shift. if, rather, the
    // values are only shifted by finite amounts, the total shift is the
    // sum of the shift, and the order of additions do not change the result
    if( wFi < NrFi )
     FModChg( tmod->shift() , wFi );
    else
     // this is a FunctionMod coming from the linear 0-th component,
     // it surely does not reset any component
     continue;

    if( reset[ wFi ] ) {
     // any kind of FunctionMod after (before) one that completely reset the
     // component is useless, delete it and move forward (backward)
     to_delete = true;
     continue;
     }

    // if the component is not reset (yet), one must look in details what
    // exact type the *FunctionMod* is and react accordingly

    {
     // a C05FunctionModRngd only changes existing linearizations, and
     // therefore is never a "hard" reset; the only easy case is
     // NothingChanged, which by definition does nothing save for the
     // shift(), that has been dealt with already
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionModRngd>( tmod );
     if( ttmod ) {
      switch( ttmod->type() ) {
       case( C05FunctionMod::NothingChanged ):
	to_delete = true;
       case( C05FunctionMod::AllLinearizationChanged ):
       case( C05FunctionMod::AllEntriesChanged ):
        continue;
       default:
	throw( std::invalid_argument( "wrong type in C05FunctionModRngd" ) );
       }  // end( switch( ttmod->f_type ) )
      }  // end( if( ttmod ) )
     }  // end C05FunctionModRngd

    {
     // a C05FunctionModSbst only changes existing linearizations, and
     // therefore is never a "hard" reset; the only easy case is
     // NothingChanged, which by definition does nothing save for the
     // shift(), that has been dealt with already
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionModSbst>( tmod );
     if( ttmod ) {
      switch( ttmod->type() ) {
       case( C05FunctionMod::NothingChanged ):
	to_delete = true;
       case( C05FunctionMod::AllLinearizationChanged ):
       case( C05FunctionMod::AllEntriesChanged ):
        continue;
       default:
	throw( std::invalid_argument( "wrong type in C05FunctionModSbst" ) );
       }  // end( switch( ttmod->f_type ) )
      }  // end( if( ttmod ) )
     }  // end C05FunctionModSbst

    {
     // a C05FunctionMod of type GlobalPoolRemoved with which.empty() resets
     // all the component. NothingChanged by definition does nothing (save
     // for the shift(), that has been dealt with already). all other cases
     // will have to be dealt with later
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionMod>( tmod );
     if( ttmod ) {
      switch( ttmod->type() ) {
       case( C05FunctionMod::GlobalPoolRemoved ):
	if( ttmod->which().empty() )
	 reset[ wFi ] = true;
       case( C05FunctionMod::NothingChanged ):
	to_delete = true;
       case( C05FunctionMod::AllLinearizationChanged ):
       case( C05FunctionMod::AllEntriesChanged ):
       case( C05FunctionMod::AlphaChanged ):
       case( C05FunctionMod::GlobalPoolAdded ):
	continue;
       default:
	throw( std::invalid_argument( "wrong type in C05FunctionMod" ) );
       }  // end( switch( ttmod->f_type ) )
      }  // end( if( ttmod ) )
     }  // end C05FunctionMod

    {
     // a C05FunctionModLin* only changes existing linearizations, and
     // therefore is never a "hard" reset
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionModLinRngd>(
								       tmod );
     if( ttmod )
      continue;
     }

    // if control reaches this point, mod is (indistinguishable from) a base
    // FunctionMod, and in particular it is not a C05FunctionMod*. hence the
    // change in the Function is not quasi-additive, and therefore a fortiori
    // not strongly quasi-additive. as a result, this is a "hard" reset

    reset[ wFi ] = true;
    to_delete = true;

    }  // end( if( tmod ) )
   }  // end FunctionMod

  {
   // a "naked" FunctionModVars is only allowed if there is only one
   // component (comprised the linear one). if it is allowed, it is
   // of no consequence here, except for the possible effect on the
   // function values, if it is a C05FunctionModVars*, meaning that it
   // represents a strongly quasi-additive variable change. if not, the
   // variable change also implies a reset
   // in no case, however, the Modification is removed from the list

   const auto tmod = std::dynamic_pointer_cast<FunctionModVars>( mod );
   if( tmod ) {
    if( ( NrFi > 1 ) || f_lf )
     throw( std::invalid_argument( "naked FunctionModVars not allowed" ) );

    auto wFi = get_index_of_component( tmod->function() );

    FModChg( tmod->shift() , wFi );  // change/reset upper/lower values

    {
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionModVarsAddd>(
								       tmod );
     if( ttmod )
      continue;
     }

    {
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionModVarsRngd>(
								       tmod );
     if( ttmod )
      continue;
     }

    {
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionModVarsSbst>(
								       tmod );
     if( ttmod )
      continue;
     }

    // if control reaches here, this is a FunctionModVars* that is not a
    // C05FunctionModVars*, i.e., a non strongly quasi-additive variable
    // change, which implies a "hard" reset for the component
    reset[ wFi ] = true;
    continue;

    }  // end( if( tmod ) )
   }  // end FunctionModVars

  {
   // a GroupModification here can only be a bunch of identical
   // *FunctionModVar*: pick the first one and act on it
   const auto tmod = std::dynamic_pointer_cast<GroupModification>( mod );
   if( tmod ) {
    auto fmod = tmod->v_sub_Modifications.front();

    {
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionModVarsAddd>(
								       fmod );
     if( ttmod )
      continue;
     }

    {
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionModVarsRngd>(
								       fmod );
     if( ttmod )
      continue;
     }

    {
     const auto ttmod = std::dynamic_pointer_cast<C05FunctionModVarsSbst>(
								       fmod );
     if( ttmod )
      continue;
     }

    // if control reaches here, this is a FunctionModVars* that is not a
    // C05FunctionModVars*, i.e., a non strongly quasi-additive variable
    // change, which implies a "hard" reset for *all* components
    reset.assign( NrFi , true );
    continue;
    }
   }

  {
   const auto tmod = std::dynamic_pointer_cast<ConstraintMod>( mod );
   if( tmod )
    throw( std::invalid_argument( "ConstraintMod not handled (yet)" ) );
   }

  {
   const auto tmod = std::dynamic_pointer_cast<VariableMod>( mod );
   if( tmod )
    throw( std::invalid_argument( "VariableMod not handled (yet)" ) );
   }

  {
   const auto tmod = std::dynamic_pointer_cast<BlockMod>( mod );
   if( tmod )
    throw( std::invalid_argument( "BlockMod not handled (yet)" ) );
   }

  {
   const auto tmod = std::dynamic_pointer_cast<BlockModAD>( mod );
   if( tmod )
    throw( std::invalid_argument( "BlockModAD not handled (yet)" ) );
   }

  // if control reaches here, the Modification is "unknown", probably a
  // "physical" Modification that BundleSolver does not care about

  to_delete = true;

  }  // end( 1st loop, in reverse )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if any global pool has been reset (or, even better, *all* of them have
 // reset) then delete all items in the (<corresponding) bundle

 Dctnry inv_dict;
 
 if( reset.find( reset.begin() , reset.end() , false ) == reset.end() ) {
  // all components have been reset

  NrItems.assign( NrFi , 0 );
  OOBase.assign( vBPar2[ NrFi ] , Inf<SIndex>() );
  ItemVcblr.assign( vBPar2[ NrFi ] , make_pair( InINF , Inf<SIndex>() ) );
  whisG1.assign( NrFi , InINF );
  FreList.clear();

  Master->RmvItems();
  }
 else
  if( reset.find( reset.begin() , reset.end() , true ) != reset.end() ) {
   // at least a component has been reset: need to construct the inverse
   // dictionary < component , global pool position > --> bundle position
   // (in linear time) to do removals efficiently

   compute_inverse_dictionary( inv_dict );

   for( Index k = 0 ; k < NrFi ; ++k )
    if( reset[ k ] )
     for( auto i : inv_dict[ k ] )
      Delete( i , true );
   }

 // After this point, all the Modification adding, deleting or modifying
 // linearizations are significant: they either pertain to components that
 // have never been reset, or are the remaining ones after the (last) one
 // resetting the component

 if( v_mod_tmp.empty() )  // no more Modification to process
  return;                 // all done

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 2nd loop, again in reverse: check for "soft" reset of components, i.e.,
 // when all existing linearization changes. any Modification that changes
 // the linearizations happening before a "soft" reset of the global pool
 // (meaning it is found afterwards in the reverse order) is deleted since it
 // is useless. Actually, for this to be true we need to keep track already if
 // all the constants change (so that we can handle AllLinearizationChanged)
 //
 // note that Modification changing the linearizations happening *after* a
 // "soft" reset of the global pool (meaning it is found *before* in the
 // reverse order is also useless, since a reset forces the re-reading of all
 // linearizations, which by definition happens at their current (final) state.
 // yet, this is not done immediately
 //
 // note that we make no serious attempt at keeping track of the combined
 // effect of all changes, in order to detect if a large set of small
 // changes actually imples a reset. this is complicated for "horizontal"
 // changes (for all linearizations, a range/subset of entries) because the
 // names of the changed Variable may not be current (additions/deletions may
 // happen in the meantime), and keeping track is too burdensome. similarly
 // for "vertical" changes (a set of specific linearizations). some steps
 // in this direction will be done in later stages
 //
 // another note is that BundleSolver (due to limitations in the interface of
 // MPSolver) has an all-or-nothing approach to changing the constants. thus,
 // AlphaC[ wFi ] is set to true whenever a change in any constant happens,
 // even if it is on a subset of the linearizations

 reset.assign( NrFi , false );  // reset reset (couldn't resist)
 std::vector<bool> AlphaC( NrFi , false );

 for( auto rimod = v_mod_tmp.rbegin() ; rimod != v_mod_tmp.rend() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to rimod
      [ & to_delete , & v_mod_tmp ]( decltype( rimod ) & ri ) {
       if( to_delete )
	ri = decltype( ri )( v_mod_tmp.erase( std::next( ri ).base() ) );
       else
	++ri;
       }( rimod ) ) {
  bool to_delete = false;
  auto mod = *rimod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  {
   // a C05FunctionModRngd only changes a range of the linearizations,
   // and therefore is not considered a "soft" reset even if which().empty()
   // in fact the range could be so large as to be (almost) all the
   // variables, which would count as a reset, but so far we don't attempt
   // at detecting this. the only easy case would be NothingChanged, but
   // any such C05FunctionModRngd has been deleted already. however, if the
   // component is "soft" reset already, it can be deleted

   const auto tmod = std::dynamic_pointer_cast<C05FunctionModRngd>( mod );
   if( tmod ) {
    auto wFi = get_index_of_component( tmod->function() );
    if( reset[ wFi ] )
     to_delete = true;
    continue;
    }  // end( if( ttmod ) )
   }  // end C05FunctionModRngd

  {
   // a C05FunctionModSbst only changes a subet of the linearizations,
   // and therefore is not considered a "soft" reset even if which().empty()
   // in fact the subet could be so large as to be (almost) all the
   // variables, which would count as a reset, but so far we don't attempt
   // at detecting this. the only easy case would be NothingChanged, but
   // any such C05FunctionModSbst has been deleted already. however, if the
   // component is "soft" reset already, it can be deleted

   const auto tmod = std::dynamic_pointer_cast<C05FunctionModSbst>( mod );
   if( tmod ) {
    auto wFi = get_index_of_component( tmod->function() );
    if( reset[ wFi ] )
     to_delete = true;
    continue;
    }  // end( if( ttmod ) )
   }  // end C05FunctionModSbst

  {
   // a C05FunctionMod of type AllLinearizationChanged or AllEntriesChanged
   // with which.empty() "soft" resets all the component, and in the former
   // case also Alpha; AlphaChanged only changes the constants (obviously)

   const auto tmod = std::dynamic_pointer_cast<C05FunctionMod>( mod );
   if( tmod ) {
    auto wFi = get_index_of_component( tmod->function() );

    // AlphaChanged is considered by default applied to all the
    // linearizations, irrespectively of which()
    if( tmod->type() == C05FunctionMod::AlphaChanged ) {
     AlphaC[ wFi ] = true;
     to_delete = true;
     continue;
     }

    // in all other cases we only react to which().empty() 
    if( ttmod->type() == C05FunctionMod::AllLinearizationChanged ) {
     AlphaC[ wFi ] = true;
     if( tmod->which().empty() ) {
      reset[ wFi ] = true;
      to_delete = true;
      }
     continue;
     }

    if( tmod->type() == C05FunctionMod::AllEntriesChanged )
     if( tmod->which().empty() ) {
      reset[ wFi ] = true;
      to_delete = true;
      }

    continue;
    }  // end( if( tmod ) )
   }  // end C05FunctionMod

  {
   // a C05FunctionModLinRngd implies that a specific range in all the
   // linearizations must be changed by adding; this is never considered a
   // "soft" reset even if in fact the range could be so large as to be
   // (almost) all the variables, which would count as a reset, but so far
   // we don't attempt at detecting this. however, if the component is "soft"
   // reset already, it can be deleted
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModLinRngd>( mod );
   if( tmod ) {
    auto wFi = get_index_of_component( tmod->function() );
    if( reset[ wFi ] )
     to_delete = true;
    continue;
    }  // end( if( ttmod ) )
   }  // end C05FunctionModLinRngd

  {
   // a C05FunctionModLinSbst implies that a specific subset in all the
   // linearizations must be changed by adding; this is never considered a
   // "soft" reset even if in fact the subset could be so large as to be
   // (almost) all the variables, which would count as a reset, but so far
   // we don't attempt at detecting this. however, if the component is "soft"
   // reset already, it can be deleted
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModLinSbst>( mod );
   if( tmod ) {
    auto wFi = get_index_of_component( tmod->function() );
    if( reset[ wFi ] )
     to_delete = true;
    continue;
    }  // end( if( ttmod ) )
   }  // end C05FunctionModLinSbst

  {
   // a C05FunctionModLin implies that *all* the linearizations must be
   // changed by adding them \delta; this may in principle be handled in
   // a specialised way by BundleSolver, but is currently not, and
   // therefore it is a "full" reset
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModLin>( mod );
   if( tmod ) {
    auto wFi = get_index_of_component( tmod->function() );
    reset[ wFi ] = true;
    to_delete = true;
    }  // end( if( ttmod ) )
   }  // end C05FunctionModLin
  }  // end( 2nd loop, again in reverse )

 // note that even if there were no more Modification to process we could not
 // stop because this means that reset[ wFi ] == true for some wFi. in fact
 // v_mod_tmp as not empty(), and elements can be remove from it only if
 // some component is "soft" reset. 

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 3rd loop, forward: prepare for addition/removal/changes of individual
 // linearization for each component by computing the three sets
 // - linearizations that need be removed
 // - linearizations that need be added
 // - linearizations that need be changed
 // Note that:
 // - if a linearization that is added/changed is later removed, it is no
 //   longer added/changed
 // - if a linearization that is removed/changed is later added it is no
 //   longer removed/changed (adding over an existing linearization changes
 //   it anyway, no reason to remove it)
 // - changes to a reset component can be ignored
 // - due to limitations in the MPSolver interface, "horizontal" changes (of
 //   a given range/subset of entries) to a subset of linearizations are
 //   not supported. these must be either mapped in "horizontal" changes to
 //   *all* linearizations (of a given component), or to "vertical" changes
 //   of all components to a subset of linearization. somewhat arbitrarily,
 //   the second option is chosen here. as a consequence, C05FunctionModRngd
 //   and C05FunctionModSbst are considered C05FunctionMod. note that
 //   C05FunctionMod* with which().empty() still remain untreated, as well
 //   as C05FunctionModLinRngd and C05FunctionModLinSbst (that by definition
 //   always concern all the linearizations), while C05FunctionModLin have
 //   been dealt with already

 Dctnry Addd( NrFi );
 Dctnry Rmvd( NrFi );
 Dctnry Chgd( NrFi );

 for( auto imod = v_mod_tmp.begin() ; imod != v_mod_tmp.end() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to imod
      [ & to_delete , & v_mod_tmp ]( decltype( imod ) & it ) {
       if( to_delete )
	it =  v_mod_tmp.erase( it );
       else
	++it;
       }( imod ) ) {
  bool to_delete = false;
  auto mod = *imod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly
  //
  // actually, only C05FunctionMod need be treated here. we do not need to
  // distinguish between the base and the derived classes because we take
  // the change of a range/subset of the entries for a change of the whole
  // linearization. however, we do in fact distinguish them because we
  // ignore the Modification if ttmod->which().empty() and type() ==
  // AllEntriesChanged or == AllLinearizationChanged. any such C05FunctionMod
  // has been deleted (at worst) in the 2nd loop, and therefore here it must
  // be a C05FunctionModRngd or C05FunctionModSbst
  //
  // also, note that NothingChanged and AlphaChanged cannot happen, since
  // again these cases have been dealt with already

  const auto tmod = std::dynamic_pointer_cast<C05FunctionMod>( mod );
  if( tmod ) {
   auto wFi = get_index_of_component( tmod->function() );

   switch( ttmod->type() ) {
    case( C05FunctionMod::AllLinearizationChanged ):
    case( C05FunctionMod::AllEntriesChanged ):
     // if tmod->which().empty(), this must actually be either a
     // C05FunctionModRngd or a C05FunctionModSbst: save it, for it
     // will be dealt with in the next loop
     if( tmod->which().empty() )
      continue;
     to_delete = true;   // otherwise, delete it
     if( reset[ wFi ] )  // changes in reset components
      continue;          // are ignored
     // add to Chgd[ wFi ] the names in tmod->which(), save those that are
     // in either Addd[ wFi ] or Rmvd[ wFi ]
     if( Addd[ wFi ].empty() && Rmvd[ wFi ].empty() ) {
      // no items are added/removed, all changed are changed
      if( Chgd[ wFi ].empty() )
       Chgd[ wFi ] = tmod->which();
      else {
       Subset tmp( std::min( Chgd[ wFi ].size() + tmod->which().size() ,
			     vBPar2[ wFi ] ) );
       std::set_union( Chgd[ wFi ].begin() , Chgd[ wFi ].end() ,
		       tmod->which().begin() , tmod->which().end() ,
		       tmp.begin() );
       Chgd[ wFi ] = std::move( tmp );
       }
      }
     else {
      // only those items that are not added and/or removed need be changed
      Subset tmp( tmod->which() );
      if( ! Addd[ wFi ].empty() ) {
       Subset tmp2( tmp.size() );
       std::set_difference( tmp.begin() , tmp.end() ,
			    Addd[ wFi ].begin() , Addd[ wFi ].end() ,
			    tmp2.begin() );
       tmp = std::move( tmp2 );
       }
      if( ( ! Rmvd[ wFi ].empty() ) && ( ! tmp.empty() ) ) {
       Subset tmp2( tmp.size() );
       std::set_difference( tmp.begin() , tmp.end() ,
			    Rmvd[ wFi ].begin() , Rmvd[ wFi ].end() ,
			    tmp2.begin() );
       tmp = std::move( tmp2 );
       }
      if( Chgd[ wFi ].empty() )
       Chgd[ wFi ] = std::move( tmp );
      else {
       Subset tmp2( std::min( Chgd[ wFi ].size() + tmp.size() ,
			      vBPar2[ wFi ] ) );
       std::set_union( Chgd[ wFi ].begin() , Chgd[ wFi ].end() ,
		       tmp.begin() , tmp.end() , tmp2.begin() );
       Chgd[ wFi ] = std::move( tmp2 );
       }
      }
    case( C05FunctionMod::GlobalPoolAdded ):
     // add to Addd[ wFi ] the names in tmod->which(), and remove them
     // from Rmvd[ wFi ], and Chgd[ wFi ] if the component is not reset
     if( Addd[ wFi ].empty() )
      Addd[ wFi ] = tmod->which();
     else {
      Subset tmp( std::min( Addd[ wFi ].size() + tmod->which().size() ,
			    vBPar2[ wFi ] ) );
      std::set_union( Addd[ wFi ].begin() , Addd[ wFi ].end() ,
		      tmod->which().begin() , tmod->which().end() ,
		      tmp.begin() );
      Addd[ wFi ] = std::move( tmp );
      }
     if( ! Rmvd[ wFi ].empty() ) {
      Subset tmp( Rmvd[ wFi ].size() );
      std::set_difference( Rmvd[ wFi ].begin() , Rmvd[ wFi ].end() ,
			   tmod->which().begin() , tmod->which().end() ,
			   tmp.begin() );
      Rmvd[ wFi ] = tmp;
      }
     if( ( ! reset[ wFi ] ) && ( ! Chgd[ wFi ].empty() ) ) {
      Subset tmp( Chgd[ wFi ].size() );
      std::set_difference( Chgd[ wFi ].begin() , Chgd[ wFi ].end() ,
			   tmod->which().begin() , tmod->which().end() ,
			   tmp.begin() );
      Chgd[ wFi ] = std::move( tmp );
      }
     to_delete = true;
     continue;
    case( C05FunctionMod::GlobalPoolRemoved ):
     // add to Rmvd[ wFi ] the names in tmod->which(), and remove them
     // from Addd[ wFi ], and Chgd[ wFi ] if the component is not reset
     if( Rmvd[ wFi ].empty() )
      Rmvd[ wFi ] = tmod->which();
     else {
      Subset tmp( std::min( Rmvd[ wFi ].size() + tmod->which().size() ,
			    vBPar2[ wFi ] ) );
      std::set_union( Rmvd[ wFi ].begin() , Rmvd[ wFi ].end() ,
		      tmod->which().begin() , tmod->which().end() ,
		      tmp.begin() );
      Rmvd[ wFi ] = std::move( tmp );
      }
     if( ! Addd[ wFi ].empty() ) {
      Subset tmp( Addd[ wFi ].size() );
      std::set_difference( Addd[ wFi ].begin() , Addd[ wFi ].end() ,
			   tmod->which().begin() , tmod->which().end() ,
			   tmp.begin() );
      Addd[ wFi ] = std::move( tmp );
      }
     if( ( ! reset[ wFi ] ) && ( ! Chgd[ wFi ].empty() ) ) {
      Subset tmp( Chgd[ wFi ].size() );
      std::set_difference( Chgd[ wFi ].begin() , Chgd[ wFi ].end() ,
			   tmod->which().begin() , tmod->which().end() ,
			   tmp.begin() );
      Chgd[ wFi ] = std::move( tmp );
      }
     to_delete = true;
    }  // end( switch( tmod->type() ) )
   }  // end( if( ttmod ) )
  }  // end( 3rd loop, forward )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now act on the just gathered information, i.e., delete all linearization
 // that need to, if any

 if( Rmvd.find_if( Rmvd.begin() , Rmvd.end() ,
		   []( Subset & Rk ) { return( ! Rk.empty() ); }
		   ) != Rmvd.end() ) {
  // at least a component has had lnearizations removed: need to construct
  // the inverse dictionary < component , global pool position > --> bundle
  // position (if not constructed already) to do removals efficiently

  compute_inverse_dictionary( inv_dict );

  for( Index k = 0 ; k < NrFi ; ++k )
   for( auto i : Rmvd[ k ] )
    Delete( inv_dict[ k ][ i ] , true );
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 4th loop: handle exclusively Variable addition and removal, in forward
 // order. all removals of existing Variable are immediately performed,
 // diminishing NumVar. additions of Variable are "cached", and only performed
 // once (if ever) after the end of the cycle. removals of Variable that have
 // not been added yet just decreases the number of new Variable to be added
 // at the end. all corresponding Modification are removed from the list

 Index to_add = 0;
 bool addd_vars = false;  // if any Variable has ever been added
 bool rmvd_vars = false;  // if any Variable has ever been removed

 for( auto imod = v_mod_tmp.begin() ; imod != v_mod_tmp.end() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to imod
      [ & to_delete , & v_mod_tmp ]( decltype( imod ) & it ) {
       if( to_delete )
	it =  v_mod_tmp.erase( it );
       else
	++it;
       }( imod ) ) {
  bool to_delete = false;
  auto mod = *imod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // not that we do not distinguish C05FunctionModVars* from "plain"
  // FunctionModVars*, since the only difference is whether or not the
  // operation is strongly quasi-additive, i.e., it implies or not a "hard"
  // reset, but this has already been acted upon
  
  {
   // a "naked" FunctionModVars
   const auto tmod = std::dynamic_pointer_cast<FunctionModVars>( mod );
   if( ! tmod ) {
    // if it is not a "naked" FunctionModVars, it can still be a group of
    // identical *FunctionModVars* "dressed" into a GroupModification
    const auto gmod = std::dynamic_pointer_cast<GroupModification>( mod );
    if( gmod )  // if so, pick the first one and act on it
     tmod = std::static_pointer_cast<FunctionModVars>(
				        tmod->v_sub_Modifications.front() );
     }

   if( tmod ) {
    // if we have a *FunctionModVars*, we have to distinguish its exact type
    // and add/delete Variable accordingly; in all cases, however, the
    // Modification is processed and can be deleted
    to_delete = true;

    {
     const auto ttmod = std::dynamic_pointer_cast<FunctionModVarsAddd>(
								       tmod );
     if( ttmod ) {
      addd_vars = true;
      if( ! to_add ) {
       // the first time, check that the Modification data agrees with what
       // we expect
       if( ttmod->first() != NumVar )
	throw( std::logic_error( "wrong Variable names in FunctionModVars" )
	       );
       }
      
      to_add += vars.size();
      continue;
      }
     }

    {
     const auto ttmod = std::dynamic_pointer_cast<FunctionModVarsRngd>(
								       tmod );
     if( ttmod ) {
      rmvd_vars = true;
      Range rng = ttmod->range();
      if( rng.first >= NumVar ) {  // all the Variable are deleted already
       auto nr = rng.second - rng.first;
       if( nr > to_add )
	throw( std::logic_error( "removing non-existing Variable" ) );
       to_add -= nr;               // "virtually" remove them
       continue;                   // nothing else to do
       }
      if( rng.second >= NumVar ) {  // some of the Variable are deleted already
       auto nr = rng.second - NumVar;
       if( nr > to_add )
	throw( std::logic_error( "removing non-existing Variable" ) );
       to_add -= nr;               // "virtually" remove them
       rng.second = NumVar;
       }
      if( rng.second < NumVar ) {
       // if deleting the last range of Variable nothing has to be done,
       // but deleting Variable "in the middle" rather requires moving
       // down the remaining range of values in Lambda
       std::copy( Lambda.begin() + rng.second ,
		  Lambda.begin() + Lambda.end() ,
		  Lambda.begin() + rng.first );
       }
      Subset tdlt( rng.second - rng.first );
      NumVar -= tdlt.size();
      Lambda.resize( NumVar );  // adjust Lambda
      Lambda1.resize( NumVar );
      if( MaxSol > 1 )
       LmbdBst.resize( NumVar );
      std::iota( tdlt.begin() , tdlt.end() , rng.first );
      Master->RmvVars( tdlt.data() , tdlt.size() );  // remove from MP
      continue;
      }
     }

    {
     const auto ttmod = std::dynamic_pointer_cast<FunctionModVarsSbst>(
								       tmod );
     if( ttmod ) {
      rmvd_vars = true;
      if( ttmod->subset().front() >= NumVar ) {
       // all the Variable are deleted already
       if( ttmod->subset().size() > to_add )
	throw( std::logic_error( "removing non-existing Variable" ) );
       to_add -= ttmod->subset();  // "virtually" remove them
       continue;                   // nothing else to do
       }

      Subset & sbst;
      Subset tsbst;
      if( ttmod->subset().back() < NumVar )  // no Variable deleted already
       sbst = & ttmod->subset();             // delete them all
      else {                                 // construct the subset to delete
       auto sbstit = ttmod->subset().end();
       while( *(--it) >= NumVar );
       tsbst = Subset( ttmod->subset().begin() , ++it );
       auto nr = ttmod->subset().size() - tsbst-size();
       if( nr > to_add )
	throw( std::logic_error( "removing non-existing Variable" ) );
       to_add -= nr;               // "virtually" remove them
       }

      Compact( Lambda , sbst );  // adjust Lambda
      NumVar -= sbst.size();
      Lambda.resize( NumVar );
      Lambda1.resize( NumVar );
      if( MaxSol > 1 )
       LmbdBst.resize( NumVar );
      Master->RmvVars( sbst.data() , sbst.size() );  // remove from MP
      continue;
      }
     }

    // if control reaches here, this is an unknown *FunctionModVars* (??)
    throw( std:.logic_error( "unknown FunctionModVars" ) );

    }  // end( if( tmod ) )
   }  // end FunctionModVars
  }  // end( 4th loop, forward )

 // at this point, the set of Variable in the BundleSolver/Master Problem
 // coincides with the set of Variable in the C05Function(s), save for the
 // Variable to be added: in other words, the positions from 0 no NumVar - 1
 // in the linearizations corresponds to what BundleSolver expects

 // if there are no more Modification to process, no Variable to add, and
 // no component in need of a reset, all done
  
 if( v_mod_tmp.empty() && ( ! to_add ) &&
     ( reset.find( reset.begin() , reset.end() , true ) == reset.end() ) )
  return;

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 5th loop: handle "horizontal" changes, i.e., changes of a given range
 // (subset) of entries in all the linearizations, i.e., C05FunctionMod* and
 // C05FunctionModLin* with which().empty(). note that due to limitations
 // of the MPSolver interface, subsets are anyway translated to a range,
 // thereby possibly requiring also entries that have not changed. if
 // Variable have been removed the original indices in the Modification need
 // be translated, and in fact if Variable in the range have been removed
 // the range can be shrank, up to disappearing altogether. the mapping is
 // more difficult if Variable have also been added. however we can exploit
 // the property that if a Variable ever has a name that is larger than its
 // index in the Modification, this can only mean that the Variable has been
 // removed and then re-added after that the Modification has been issued.
 // these Variable can therefore be ignored, since (unless they have been
 // re-removed) they will be added in the end
 //
 // this is the final loop, so the list must be empty at the end

 for( ; ! v_mod.empty() ; v_mod.pop_front() ) {
  auto mod = v_mod.front();  // pick (a reference to) the first Modification

  Range range( NumVar , 0 );    // an empty range
  c_Subset * subset = nullptr;  // an empty subset
  c_Vec_p_Var & vars;           // the affected Variable

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  {
   // a C05FunctionModRngd, that at this point can only have which().empty()
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModRngd>( mod );
   if( tmod ) {
    if( ! tmod->which().empty() )
     throw( std::logic_error( "unexpected nonempty C05FunctionModRngd" ) );

    vars = tmod->vars();
    range = tmod->range();
    }
   }

  {
   // a C05FunctionModSbst, that at this point can only have which().empty()
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModSbst>( mod );
   if( tmod ) {
    if( ! tmod->which().empty() )
     throw( std::logic_error( "unexpected nonempty C05FunctionModSbst" ) );

    vars = tmod->vars();
    subset = & tmod->subset();
    }
   }

  {
   // a C05FunctionModLinRngd implies that a specific range in all the
   // linearizations must be changed (by adding something)
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModLinRngd>( mod );
   if( tmod ) {
    vars = tmod->vars();
    range = tmod->range();
    }
   }

  {
   // a C05FunctionModLinSbst implies that a specific subset in all the
   // linearizations must be changed (by adding something)
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModLinSbst>( mod );
   if( tmod ) {
    vars = tmod->vars();
    subset = & tmod->subset();
    }
   }

  if( ( range.first >= range.second ) && ( ! subset ) )
   // it is neither of the above: this should not happen
   throw( std::logic_error( "unexpected Modification slipped in" ) );

  if( ! rmvd_vars ) {
   // Variable have never been removed, hence the names can be used directly
   if( subset ) ) {  // turn the subset into a range
    range.first = subset->front();
    range.second = subset()->back() + 1;
    }
   }
  else {
   // Variable have been removed, and hence names need be actualised
   // this is done by directly checking vars() against the "active"
   // Variable of v_c05f[ 0 ], which is fairly taken as a representative
   // since all the C05Function have the same "active" Variable
   if( ! add_vars ) {
    // ... but never added: names can have only decreased, but even more
    // importantly must have remained ordered, i.e., the first "active"
    // Variable in vars() is the first variable of the range, the last
    // "active" Variable vars() is the last variable of the range
    // note that we do not use subset and range here, as the range is
    // reconstructed from scratch using vars
    auto lit = vars.begin();
    for( ; lit != vars.end() ; ++lit ) {
     range.first = v_c05f[ 0 ]->is_active( *lit );
     if( range.first < v_c05f[ 0 ]->get_num_active_var() )
      break;
     }
    if( lit == vars.end() )  // no Variable in vars is still "active"
     continue;               // nothing else to do
    // since we know that here are some "active" Variable in vars(), this
    // second loop will necessarily end
     for( auto rit = vars.rbegin() ; ; ++rit ) {
      range.second = v_c05f[ 0 ]->is_active( *lit );
      if( range.second < v_c05f[ 0 ]->get_num_active_var() )
       break;
      }
     ++range.second;  // the range is [ first , second )
     }
   else {
    // the complicated case: Variable have both been removed and added
    // names can have changed in an almost arbitrary way, except that if
    // a name has increased then the Variable has been deleted and re-added
    // and therefore need not be included
    Subset newnames( vars.size() );
    auto lit = vars.begin();
    auto nni = newnames.begin();
    if( subset ) {
     auto sit = subset->begin();
     for( ; lit != vars.end() ; ++lit , ++sit ) {
      auto i = v_c05f[ 0 ]->is_active( *lit );
      if( ( i <= *sit ) && ( i < v_c05f[ 0 ]->get_num_active_var() ) )
       *(nni++) = i;
      }
     }
    else {
     for( ; lit != vars.end() ; ++lit , ++range.first ) {
      auto i = v_c05f[ 0 ]->is_active( *lit );
      if( ( i <= range.first ) && ( i < v_c05f[ 0 ]->get_num_active_var() ) )
       *(nni++) = i;
      }
     }
    if( nni == newnames.begin() )  // no Variable in vars is still "active"
     continue;                     // nothing else to do
    newnames.resize( std::distance( newnames.begin() , nni ) );
    range.first = newnames.front();
    range.second = newnames.back();
    ++range.second;  // the range is [ first , second )
    }
   }

  // now actually do it
  Master->ChgSubG( range.first , range.second , NrFi + 1 );

  }  // end( 5th loop, forward )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are Variable to add, do it now in one blow
 // there is a trade-off here: doing this now causes Master->AddVars() to
 // (indirectly) call get_linearization_coefficients() on a smaller set of
 // linearizations, if additions are done, but on the other hand increases
 // NumVar and therefore the work done in later stages. hence, this is done
 // here only if no additions are done

 bool toadd = Addd.find_if( Addd.begin() , Addd.end() ,
			    []( Subset & Ak ) { return( ! Ak.empty() ); }
			    ) != Addd.end();
 if( to_add && ( ! toadd ) ) {
  NumVar += to_add;
  Lambda.resize( NumVar , 0 );
  Lambda1.resize( NumVar , 0 );
  if( MaxSol > 1 )
   LmbdBst.resize( NumVar , 0 );
  Master->AddVars( to_add );
  to_add = 0;  // done already
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are linearization to add/change, do it now in one blow
 // note that due to limitations in the MPSolver interface, changing a
 // linearization is identical to adding one, even if the change was limited
 // to a range/subset of the entries

 if( toadd ||
     ( Chgd.find_if( Chgd.begin() , Chgd.end() ,
		     []( Subset & Ck ) { return( ! Ck.empty() ); }
		     ) != Chgd.end() ) ) {
  // at least a component has had lnearizations added or changed: need to
  // construct the inverse dictionary < component , global pool position >
  // --> bundle position (if not constructed already) to do additions/changes
  // efficiently

  compute_inverse_dictionary( inv_dict );

  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( Chgd[ k ].size() >= NrItems[ k ] )  // all items change
    AlphaC[ k ] = reset[ k ] = false;      // this component is served

   // compute the union between Addd[ k ] and Chgd[ k ] into Addd[ k ]
   if( Addd[ k ].empty() )
    if( Chgd[ k ].empty() )
     continue;
    else
     Addd[ k ] = std::move( Chgd[ k ] );
   else
    if( ! Chgd[ k ].empty() ) {
     Subset tmp( Addd[ k ].size() + Chgd[ k ].size() );
     std::set_union( Addd[ k ].begin() , Addd[ k ].end() ,
		     Chgd[ k ].begin() , Chgd[ k ].end() , tmp.begin() );
     Addd[ wFi ] = std::move( tmp );
     }

   for( auto i : Addd[ k ] ) {
    double *G1 = Master->GetItem( k + 1 );
    v_c05f[ k ]->get_linearization_coefficients( G1 ,
						 make_pair( 0 , NumVar ) , i );
    auto Ai = v_c05f[ wFi ]->get_linearization_constant( i );
    Master->SetItemBse( nullptr , NumVar );
    double ScPri;
    if( v_c05f[ wFi ]->is_linearization_vertical( i ) )
     Master->CheckCnst( Ai , ScPri , Lambda.data() );
    else {
     Ai = UpRifFi[ k ] - Ai -
          std::inner_product( Lambda.begin() , Lambda.end() ,
			      G1.data() , VarValue( 0 ) );
     Master->CheckSubG( 0 , 0 , Ai , ScPri );
     }
    Master->SetItem( inv_dict[ k ][ i ] );
    }
   }
  }  // end( if( additions or changes ) )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if some component need be reset

 if( reset.find( reset.begin() , reset.end() , true ) != reset.end() )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( reset[ k ] )
    Master->ChgSubG( 0 , NumVar , k + 1 );

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are Alphas to change, do it now in one blow

 if( AlphaC.find( AlphaC.begin() , AlphaC.end() , true ) != AlphaC.end() ) {
  std::vector< VarValue > Gi( NumVar );

  for( Index k = 0 ; k < NrFi ; ++k )
   if( AlphaC[ k ] ) {
    std::vector< VarValue > Alfa( Master->MaxName( k + 1 ) );

    for( Index i = 0 ; i < Master->MaxName( k + 1 ) ; ++i )
     if( ( ItemVcblr[ i ].first == k ) &&
	 ( ItemVcblr[ i ].second < vBPar2[ k ] ) &&
	 ( ItemVcblr[ i ].second >= 0 )  ) {
      auto Ai = v_c05f[ k ]->get_linearization_constant(
						     ItemVcblr[ i ].second );
      if( std::isnan( Ai ) )  // linearization no longer valid
       throw( std::logic_error( "inconsistent ItemVcblr" ) );

      // compute the linearization error in Lambda
      v_c05f[ k ]->get_linearization_coefficients( G1.data() ,
						   Range( 0 , NumVar ) ,
						   ItemVcblr[ i ].second );
      Alfa[ i ] = UpRifFi[ wFi ] - Ai -
                  std::inner_product( Lambda.begin() , Lambda.end() ,
				      Gi.data() , VarValue( 0 ) );
      }

    Master->ChgAlfa( Alfa.data() , k + 1 );
    }
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are (still) Variable to add, do it now in one blow

 if( to_add ) {
  NumVar += to_add;
  Lambda.resize( NumVar , 0 );
  Lambda1.resize( NumVar , 0 );
  if( MaxSol > 1 )
   LmbdBst.resize( NumVar , 0 );
  Master->AddVars( to_add );
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // and this, finally, is all!!

 }  // end( BundleSolver::process_outstanding_Modification ) - - - - - - - - -

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

void BundleSolver::PrintBundle( void )
{
 if( ! f_log )
  return;

 auto Alfa = Master->ReadLinErr();
 std::vector< VarValue > G( NumVar );
 // 
 *f_log << std::endl;
 for( Index i = 0 ; i < Master->MaxName() ; ++i ) {
  *f_log << i << "\t";
  if( ItemVcblr[ i ].second >= vBPar2[ ItemVcblr[ i ].first ]
	  || ItemVcblr[ i ].second < 0 ) {
   *f_log << "[empty]" << std::endl;
   continue;
   }

  auto wFi = ItemVcblr[ i ].first;
  auto j = ItemVcblr[ i ].second;
  *f_log << wFi << "\t" << j << "\t[ ";

  v_c05f[ wFi ]->get_linearization_coefficients( G.data() ,
						 Range( 0 , NumVar ) , j );
  for( Index h = 0 ; h < NumVar - 1 ; ++h )
   *f_log << G[ h ] << ", ";
   
  *f_log << G.back() << " ]\t"
	 << v_c05f[ wFi ]->get_linearization_constant( j )
	 << "\t" << Alfa[ i ] << std::endl;
  }
 }
   
#endif

/*--------------------------------------------------------------------------*/
/*--------------- METHODS OF BundleSolver::FakeFiOracle --------------------*/
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
 return( bslv->vBPar2[ bslv->NrFi ] );
 }

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetMinusInfinity( void )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetMaxNZ( cIndex wFi ) const
{
 if( wFi != InINF )
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
 return( LagB->get_NzMat() );

 } // end ( BundleSolver::FakeFiOracle::GetANZ() ) - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::GetADesc( cIndex wFi , int *Abeg , int *Aind ,
					   double *Aval , cIndex strt ,
					   Index stp )
{
 if( ! bslv->IsEasy[ wFi - 1 ] )
  throw( std::logic_error( "the Function is not a Lagrangian one" ) );

 auto LagB = static_cast<LagBFunction *>( bslv->v_c05f[ wFi - 1 ] );
 LagB->get_MatDesc( Abeg , Aind , Aval , strt , stp );

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

 if( Name == bslv->vBPar2[ bslv->NrFi ] ) // get the zero-component subgradient
  bslv->f_lf->get_linearization_coefficients( SubG , range );
 else
  bslv->v_c05f[ bslv->ItemVcblr[ Name ].first ]->
   get_linearization_coefficients( SubG , range ,
				   bslv->ItemVcblr[ Name ].second );
 SGBse = nullptr;
 return( stp - strt );
 }

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetVal( cIndex Name )
{
 if( Name == bslv->vBPar2[ bslv->NrFi ] ) // get the zero-component subgradient
  return( bslv->f_lf->get_linearization_constant( ) );
 else
  return( bslv->v_c05f[ bslv->ItemVcblr[ Name ].first ]->
		 get_linearization_constant( bslv->ItemVcblr[ Name ].second ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetGiName( cIndex Name )
{
 throw( std::logic_error( "this method cannot be called" ) );
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
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::Aggregate( cHpRow Mlt , cIndex_Set NmSt ,
					    cIndex Dm , cIndex NwNm )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File BundleSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
