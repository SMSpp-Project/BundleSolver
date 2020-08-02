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

/*--------------------------------------------------------------------------*/

#define USE_MPTESTER 0

// if USE_MPTESTER is nonzero, the MPSolver is a MPTester. in particular, if
// USE_MPTESTER == 1 then the master of the MPTester is an OSIMPSolver and
// the slave is a QPPenaltyMP, while if USE_MPTESTER != 1 then the master is
// a QPPenaltyMP and the slave is an OSIMPSolver

#if USE_MPTESTER
 #include "MPTester.h"
#endif

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 1
 /* Perform long and costly checks on the data structures, coded bit-wise:
  *
  * - CHECK_DS & 1 == checks the data structures representing the bundle and
  *                   the global pools them against the MPSolver and the
  *                   C05Function(s)
  * - CHECK_DS & 2 == checks that the aggregated linearization produced by
  *                   the C05Function agrees with that produced by the
  *                   MPSolver
  */
#else
 #define CHECK_DS 0
 // never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

static void Compact( BundleSolver::Vec_VarValue & g ,
		     BundleSolver::c_Subset & B )
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

static void set_difference_in_place( BundleSolver::Subset & S1 ,
				     BundleSolver::c_Subset & S2 )
{
 // removes from S1 all elements in S2, resizing it accordingly
 // both S1 and S2 are assumed to be ordered and with unique elements

 if( S1.empty() )  // nothing to delete from
  return;          // nothing to do

 auto S1it = S1.begin();
 auto S2it = S2.begin();

 // first phase: find the first element present in both S1 and S2

 for( ; ; ) {
  while( ( S1it != S1.end() ) && ( *S1it < *S2it ) )
   ++S1it;
  if( S1it == S1.end() )
   break;
  while( ( S2it != S2.end() ) && ( *S1it > *S2it ) )
   ++S2it;
  if( S2it == S2.end() )
   break;
  if( *S1it == *S2it )
   break;
  }

 if( ( S1it == S1.end() ) || ( S2it == S2.end() ) ) // if there are none
  return;                                           // all done

 // now S1it points to the first element in S1 == than the first in S2
 // elements in S1 after the common one(s) will have to be moved
 auto S1wit = S1it++;  // skip the first equal element
 S2it++;

 for( ; ( S1it != S1.end() ) && ( S2it != S2.end() ) ; ) {
  while( ( S1it != S1.end() ) && ( *S1it < *S2it ) )
   *(S1wit++) = *(S1it++);
  if( S1it == S1.end() )
   break;
  while( ( S2it != S2.end() ) && ( *S1it > *S2it ) )
   ++S2it;
  if( S2it == S2.end() )
   break;
  if( *S1it == *S2it ) { ++S1it; ++S2it; }
  }

 while( S1it != S1.end() )  // copy the part remaining after the end of S2
  *(S1wit++) = *(S1it++);
 
 S1.resize( std::distance( S1.begin() , S1wit ) );

 }  // end( set_difference_in_place )

/*--------------------------------------------------------------------------*/

static void set_union_in_place( BundleSolver::Subset & S1 ,
				BundleSolver::c_Subset & S2 )
{
 // make S1 to be the union of S1 and S2
 if( S2.empty() )
  return;

 if( S1.empty() )
  S1 = S2;
 else {
  BundleSolver::Subset tmp;
  std::set_union( S1.begin() , S1.end() , S2.begin() , S2.end() ,
		  std::back_inserter( tmp ) );
  S1 = std::move( tmp );
  }
 }  // end( set_union_in_place )

/*--------------------------------------------------------------------------*/

static void set_union_in_place( BundleSolver::Subset & S1 ,
				BundleSolver::Subset && S2 )
{
 // make S1 to be the union of S1 and S2, if useful destroy S2 in the process
 if( S2.empty() )
  return;

 if( S1.empty() )
  S1 = std::move( S2 );
 else {
  BundleSolver::Subset tmp;
  std::set_union( S1.begin() , S1.end() , S2.begin() , S2.end() ,
		  std::back_inserter( tmp ) );
  S1 = std::move( tmp );
  }
 }  // end( set_union_in_place )

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
  3 ,  // intBPar7
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
  2    // intRstAlg, default value:
       // RstAlg = 0  -  reset algorithmic parameters
       // RstCrr = 1  -  set current point to using values of the Variable 
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

static const unsigned char RstAlg = 1;  // don't reset algorithmic parameters
static const unsigned char RstCrr = 2;  // don't reset current point to all-0
                                        // (use value from the Variable)

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

 if( ! f_Block )
  return( kBlockLocked );

 // first, process any outstanding Modification - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // v_mod is atomically copied in a temporary data structure to be processed,
 // but while the latter happens new Modification may come in; hence,
 // process_outstanding_Modification() may be called more than once

 while( ! v_mod.empty() ) {
  bool owned = f_Block->is_owned_by( f_id );       // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) )  // if not try to read_lock
   return( kBlockLocked );                         // return error on failure

  process_outstanding_Modification();

  if( ! owned )             // if the Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it
  }

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

  //!! PrintBundle();
 
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
  // is "random" and there is no reason to believe it's >= than the true
  // value)

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

  // some log about the newly obtained information- - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log2();

  // check whether either any error has occurred or time has expired- - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1[ NrFi ] == - Inf<double>() ) {
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
      ( UpFiBest <= LowerBound[ NrFi ] *
	            ( 1 - ( LowerBound[ NrFi ] > 0 ? RelAcc : - RelAcc ) ) )
      ) {
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

 const auto & sb = f_Block->get_nested_Blocks();

 if( sb.empty() ) {
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
  }

 // the set of "active" Variable in all Function must be the same- - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar = v_c05f[ 0 ]->get_num_active_var();
 LamVcblr.resize( NumVar );
 auto vi = std::as_const( v_c05f[ 0 ] )->begin();
 for( Index i = 0 ; i < NumVar ; ++vi )
  LamVcblr[ i++ ] = static_cast< ColVariable * >( & (*vi) );

 if( f_lf ) {
  if( f_lf->get_num_active_var() != NumVar )
   throw( std::logic_error( "the list of active Variable do not match" ) );

  auto v = f_lf->begin();
  for( auto vi = LamVcblr.begin() ; vi != LamVcblr.end() ; ++v , ++vi )
   if( static_cast< ColVariable * >( & (*v) ) != *vi )
    throw( std::logic_error( "the list of active Variable do not match" ) );
  }

 for( Index i = 1 ; i < sb.size() ; ++i ) {
  if( v_c05f[ i ]->get_num_active_var() != NumVar )
   throw( std::logic_error( "the list of active Variable do not match" ) );

  auto v = v_c05f[ i ]->begin();
  for( auto vi = LamVcblr.begin() ; vi != LamVcblr.end() ; ++v , ++vi )
   if( static_cast< ColVariable * >( & (*v) ) != *vi ) 
    throw( std::logic_error( "the list of active Variable do not match" ) );
  }

 // if some Variable are present, they are of the ColVariable type - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index NumBVar = 0;  // count the number of Variable in the Block
 auto v_s_Variable = f_Block->get_static_variables();
 for( auto & el : v_s_Variable ) {
  if( un_any_thing_0( ColVariable , el , ++NumBVar ) )
   continue;
  if( un_any_thing_1( ColVariable , el , NumBVar += var.size() ) )
   continue;
  if( un_any_thing_K( ColVariable , el , NumBVar += var.num_elements() ) )
   continue;
  throw( std::logic_error( "some static Variable is not a ColVariable" ) );
  }

 if( NumBVar < NumVar )
  throw( std::logic_error( "too few ColVariable in the Block" ) );

 // check that the Variable in the Block agree with that in the C05Function- -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector<ColVariable *> LamBVcblr( NumBVar );

 Index cnt = 0;
 for( auto & el : v_s_Variable )
  un_any_static( el , [ & ]( ColVariable & sv ) { LamBVcblr[ cnt++ ] = & sv;
                  } , un_any_type<ColVariable>() );

 std::sort( LamBVcblr.begin() , LamBVcblr.end() );

 std::vector<ColVariable *> LamVcblrO( LamVcblr );
 std::sort( LamVcblrO.begin() , LamVcblrO.end() );

 if( ! std::includes( LamBVcblr.begin() , LamBVcblr.end() ,
		      LamVcblrO.begin() , LamVcblrO.end() ) )
 throw( std::logic_error(
		   "some ColVariable in C05Function are not in the Block" ) );

 LamVcblrO.clear();
 LamBVcblr.clear();

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

 NrEasy = 0;
 NrFi = v_c05f.size();

 if( NrFi > 1 ) {
  IsEasy.resize( NrFi , false );
  MILP_s.resize( NrFi , nullptr );
  for( Index k = 0 ; k < NrFi ; ++k ) {
   auto LagB = dynamic_cast< LagBFunction * >( v_c05f[ k ] );
   if( LagB ) {
    MILP_s[ k ] = new MILPSolver();
    MILP_s[ k ]->set_Block( LagB->get_inner_block() );
    IsEasy[ k ] = true;
    ++NrEasy;
    }
   }

  if( ! NrEasy ) {
   IsEasy.clear();
   MILP_s.clear();
   }
  }

 // set the global pool size to all non-easy functions - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // with BPar2 == 0, BundleSolver just takes whatever size of the global pool
 // it finds already in the C05Function. with BPar2 > 0, BundleSolver ensures
 // that the size of the global pool is *at least* BPar2 by increasing it if
 // it is below. this means that:
 // - BundleSolver never *decreases* the size of the global pool
 // - BundleSolver only uses the first BPar2 linearizations in each global
 //   pool; if there are more, the other ones are ignored
 // meanwhile, also set the accuracy of multipliers to RMPAccSol

 vBPar2.resize( NrFi + 1, 0 );
 for( Index k = 0 ; k < NrFi ; ++k ) {
  if( NrEasy && IsEasy[ k ] )
   continue;
  // note: we don't really trust the accuracy of MPSolver, so we give the
  // C05Function more slack
  v_c05f[ k ]->set_par( C05Function::dblAAccMlt , 100 * RMPAccSol );
  auto gps = v_c05f[ k ]->get_int_par( C05Function::intGPMaxSz );
  if( BPar2 == 0 ) {  // use the current global pool size
   if( gps < 2 )
    throw( std::logic_error( "BPar2 == 0 but too small global pool" ) );
   vBPar2[ NrFi ] += gps;
   vBPar2[ k ] = gps;
   }
  else {              // force the global pool size to be *at least* BPar2
   if( gps < BPar2 )
    v_c05f[ k ]->set_par( C05Function::intGPMaxSz , BPar2 );
   vBPar2[ NrFi ] += BPar2;
   vBPar2[ k ] = BPar2;
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

 ItemVcblr.resize( vBPar2[ NrFi ] , make_pair( InINF , InINF ) );

 InvItemVcblr.resize( NrFi );
 for( Index k = 0 ; k < NrFi ; ++k )
  InvItemVcblr[ k ].resize( vBPar2[ k ] , InINF );

 NrItems.resize( NrFi + 1 , 0 );
 FrFItem.resize( NrFi , 0 );
 MaxItem.resize( NrFi , 0 );

 FreList = {};
 whisZ.resize( NrFi , InINF );
 Zvalid.resize( NrFi , false );

 CurrNrEvls.resize( NrFi , 0 );

 FiStatus.resize( NrFi , kUnEval );
 TrueLB = false;

 UpFiBest = Inf<VarValue>();      // best, ...
 UpRifFi.resize( NrFi + 1 , 0 );  // and reference Fi() values
 UpFiLmb1.resize( NrFi + 1 );     // upper and lower function value
 LwFiLmb1.resize( NrFi + 1 );     // ... at the tentative point
 UpFiLmb.resize( NrFi + 1 ,  Inf<VarValue>() );  // upper 
 LwFiLmb.resize( NrFi + 1 , -Inf<VarValue>() );  // ... and lower Fi-value
                                                 // ... at the current point
 LowerBound.resize( NrFi + 1 , -Inf<VarValue>() );  // global lower bounds
 vStar.resize( NrFi + 1 , 0 );
 whisG1.resize( NrFi , InINF );  // no representative yet

 ScPr1.resize( NrFi + 1 , 0 );
 Alfa1.resize( NrFi + 1 , 0 );

 Result = kError;
 SSDone = false;

 // initialize the MP Solver - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Master )        // a MPSolver is set already
  Master->SetDim();  // clear all its internal state
 else {              // a MPSolver is not set yet, create it now
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
     #if( USE_MPTESTER == 1 )
      Master = new MPTester( Master , qp );
     #else
      Master = new MPTester( qp , Master );
     #endif
    #else
     Master = qp;
    }
    #endif
  }

 InitMP();

 // reset algorithm  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this has to be done before the next step since it sets Lambda, and
 //       if linearizations are added to the bundle their linearization error
 //       depends on Lambda

 ReSetAlg( RstAlgPrm );  // Lambda is reset inside

 // deal with existing linearizations- - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if BPar7 & 8, read from all the C05Function and immediately add to the
 // master problem each and every linearization found in their global pools
 //
 // otherwise read from all the C05Function and mark into InvItemVcblr each
 // and every linearization found in their global pools
 //
 // note that if ( BPar7 & 3 ) >= 2, BundleSolver will happily delete from
 // the global pool any linearization it deletes from the bundle; yet we
 // do not immediately delete existing linearizations from the global pools
 // here. these will likely be overwritten during the optimization, and if
 // memory is a problem they can be cleaned up by the user before set_Block()
 // is called. besides, in many scenarios there will be no linearizations
 // anyway

 if( BPar7 & 8 ) {
  for( Index k = 0 ; k < NrFi ; ++k )
   for( Index i = 0 ; i < vBPar2[ k ] ; ++i )
    if( v_c05f[ k ]->is_linearization_there( i ) )
     add_to_bundle( k , i );    
  }
 else
  for( Index k = 0 ; k < NrFi ; ++k )
   for( Index i = 0 ; i < vBPar2[ k ] ; ++i )
    if( v_c05f[ k ]->is_linearization_there( i ) )
     add_to_global_pool( k , i );    


 //!! PrintBundle();

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
   if( BPar2 == value )
    break;
   if( f_Block )
    throw( std::invalid_argument( "changing BPar2 not supported yet" ) );
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
 }  // end( BundleSolver::set_par( int ) )- - - - - - - - - - - - - - - - - -

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
 }  // end( BundleSolver::set_par( double ) ) - - - - - - - - - - - - - - - -

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

void BundleSolver::get_dual_solution( Configuration *solc )
{
 for( Index i = 0 ; i < zA.size() ; ++i ) {
  if( zA[ i ].second.empty() )
   throw( std::invalid_argument( "the combination is not present" ) );
  v_c05f[ i ]->set_important_linearization( std::move( zA[ i ].second ) ,
					    zA[ i ].first );
  }
 }  // end( BundleSolver::get_dual_solution() )  - - - - - - - - - - - - - - -

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
 }  // end( BundleSolver::get_int_par ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

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
 }  // end( BundleSolver::get_dbl_par ) - - - - - - - - - - - - - - - - - - -

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
 // first of all, check if a "hard" lower bound is available
 // note: this used to be done elsewhere, in particular after each function
 //       evaluation. the rationale was that in the Lagrangian case one would
 //       do heuristics as a part of the computation, and these could produce
 //       a better lower bound that one may immediately want to check. but
 //       the truth is that one does not: the real way in which a lower bound
 //       is useful is when it enters the master problem and helps in getting
 //       the optimality conditions. thus, the right place to check and
 //       update the lower bound(s) is right before the master problem is
 //       solved. the conditional lower bound ( TrueLB == false ) is indeed
 //       useful right after the function computation to prove unboundedness,
 //       but it is available then, and anyway it typically does not change
 //       when the function is computed, unlike the "hard" one

 auto LwrBnd = f_Block->get_valid_lower_bound( false );
 if( LwrBnd > - Inf<double>() ) {
  TrueLB = true;
  if( LwrBnd != LowerBound[ NrFi ] ) {
   LowerBound[ NrFi ] = LwrBnd;
   Master->SetLowerBound( LowerBound[ NrFi ] - UpRifFi[ NrFi ] );
   }
  }
 else {
  TrueLB = false;
  if( LwrBnd != LowerBound[ NrFi ] ) {
   LowerBound[ NrFi ] = - Inf<double>();
   Master->SetLowerBound( - Inf<double>() );
   }
  }

 if( ! TrueLB )  // if not, at least pick the a "conditional" one (if any)
  LowerBound[ NrFi ] = f_Block->get_valid_lower_bound( true );

 // now, if the MPSolver accepts them, collect and if necessary set the
 // individual lower bounds. note that if all of them are finite and the
 // 0-th component is not there their sum would give an alterntive valid
 // global lower bound. however, the same information is already encoded
 // in the individual bounds, hence it's of no use. the exception is that
 // QPPenaltyMP does not allow individual lower bounds, but this is not a
 // permanent issue and it'll go away when we'll get rid of MPSolver;
 // besides, it very unlikely to ever really happen

 #if( ! USE_MPTESTER )
  // QPPenaltyMP does not allow individual lower bounds, and if a MPTester
  // is used then a QPPenaltyMP is involved anyway

  if( MPName & 1 ) 
   for( Index k = 0 ; k < NrFi ; ++k ) {
    if( NrEasy && IsEasy[ k ] )  // skip easy components
     continue;

    auto LwrBndk = v_c05f[ k ]->get_global_lower_bound();
    if( LwrBndk != LowerBound[ k ] ) {
     LowerBound[ k ] = LwrBndk;
     Master->SetLowerBound( LowerBound[ k ] - UpRifFi[ k ] , k + 1 );
     }
    }
 #endif

 /* set termination criterion - - - - - - - - - - - - - - - - - - - - - - - -
  * leftover code for a previous version of MPSolver having a MPSolver::kZero
  * parameter, now removed; to be deleted

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

 Sigma = Master->ReadSigma();              // read Sigma*
 vStar[ NrFi ] = Master->ReadFiBLambda();  // read v*

 if( IsEasy.empty() )                      // there are no easy components
  for( Index k = 0 ; k < NrFi ; ++k )
   vStar[ k ] = Master->ReadFiBLambda( k + 1 );  // read model value
 else {                                    // there are easy components
  for( Index k = 0 ; k < NrFi ; ++k )
   if( IsEasy[ k ] )                                 // for easy components
    UpFiLmb1[ k ] = Master->ReadFiBLambda( k + 1 );  // read *exact* Fi-value
   else                                              // for hard components
    vStar[ k ] = Master->ReadFiBLambda( k + 1 );     // read model value

  // add the contribution of easy components to the total function value
  if( UpFiLmb[ NrFi ] < Inf<double>() )
   for( Index k = 0 ; k < NrFi ; k++ )
    if( IsEasy[ k ] )
     vStar[ NrFi ] += UpRifFi[ k ];
  }


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

 // the z[ i ] are no longer valid
 Zvalid.assign( NrFi , false );

 // the scalar products have changed
 ScPr1.assign( NrFi , Inf<double>() );

 // additional information not present in the Bundle implementation
 // for NDOSolver interface  - - - - - - - - - - - - - - - - - - - - - - - - -

 DeltaStar = Master->ReadDStart( t ) / 2.0 + Sigma;
 cLMRow tdir = Master->Readd( true );

 NrmD = 0;                                    // d-norm
 for( Index i = 0 ; i < NumVar ; ++i )
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

 for( auto OOit = OOBase.begin() ;
      OOit != OOBase.begin() + Master->MaxName() ; ++OOit )
  if( ( *OOit < Inf<SIndex>() ) && ( *OOit > -Inf<SIndex>() ) ) {
   ++(*OOit);
   if( ! *OOit )
    ++(*OOit);
   }

 // set to 0 the OOBase[] counter for items in base (if not < 0)- - - - - - -
 // note that chechking if the multiplier is strictly positive should be
 // redundant, if one was trusting the MPSolver
 
 const Index* MBse;
 const double* Mlt = Master->ReadMult( MBse , MBDim );
 if( MBse ) {
  for( Index i ; ( i = *(MBse++) ) < InINF ; ++Mlt )
   if( ( *Mlt >= RMPAccSol ) && ( OOBase[ i ] > 0 ) )
    OOBase[ i ] = 0;
  }
 else
  for( Index i = 0 ; i < MBDim ; ++i , ++Mlt )
   if( ( *Mlt >= RMPAccSol ) && ( OOBase[ i ] > 0 ) )
    OOBase[ i ] = 0;

 /*!!
 // note that there is a case in which a component wFi has Z[ wFi ] "for free"
 // in the bundle: this is when wFi only has *one* subgradient in base (or, in
 // practice, a subgradient with multiplier very close to one). This is
 // checked here (it is basically for free), and in case whisZ[] is properly
 // set so as to avoid pointless aggregations and OOBase[] is set to -1,
 // because under no circumnstances such a subgradient can ever be removed
 // from the bundle
 //
 // it now seems to me that this is stupid, since if there is only one
 // subgradient in base no aggregation is ever performed; the only issue
 // is if the base is all (but one) taken by constraints, but then even
 // aggregating does not help

 if( MBse ) {
  for( Index i ; ( i = *(MBse++) ) < InINF ; Mlt++ )
   if( *Mlt >= Eps<HpNum>() ) {
    if( ( *Mlt >= 1 - RMPAccSol ) && Master->IsSubG( i ) ) {
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
   if( *Mlt >= Eps<double>() ) {
    if( ( *Mlt >= 1 - RAccSol ) && Master->IsSubG( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ Master->WComponent( i ) - 1 ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }
    !!*/

 }  // end( UpdtCntrs ) - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
    for( Index i = 0 ; i < NumVar ; ++i ) {
     if( tL1[ i ] < 0 )
      tL1[ i ] = 0;

     const double UBh = LamVcblr[ i ]->get_ub();
     if( tL1[ i ] > UBh )
      tL1[ i ] = UBh;
     }
   else                                 // not all variables are NN
    for( Index i = 0 ; i < NumVar ; ++i ) {
     if( Master->IsNN( i ) && ( tL1[ i ] < 0 ) )
      tL1[ i ] = 0;

     const double UBh = LamVcblr[ i ]->get_ub();
     if( tL1[ i ] > UBh )
      tL1[ i ] = UBh;
     }
  else  // there are only UB vars
   for( Index i = 0 ; i < NumVar ; ++i ) {
    const double UBh = LamVcblr[ i ]->get_ub();
    if( tL1[ i ] > UBh )
     tL1[ i ] = UBh;
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

 if( f_lf ) {  // add the linear part to the "full function"
  f_lf->compute( true );
  UpFiLmb1[ NrFi ] = f_lf->get_upper_estimate();
  LwFiLmb1[ NrFi ] = f_lf->get_lower_estimate();
  }
 else
  UpFiLmb1[ NrFi ] = LwFiLmb1[ NrFi ] = 0;

 for( Index k = 0 ; k < NrFi ; ++k ) {
  if( ( ! IsEasy.empty() ) && IsEasy[ k ] )  // if k is an easy component
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

  // sum over the components, the 0th-component is already there - - - - - - -

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

 }  // end( BundleSolver::FormLambda1 ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::FiAndGi( Index wFi )
{
 double UpCutOff, LwCutOff, LwFiK, EpsCurr;

 if( ( ! IsEasy.empty() ) && IsEasy[ wFi ] )
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

 for( Index Ftchd = 0 ; Ftchd < aBP3 ; ++Ftchd ) {
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
  // is necessary due to limitations in the MPSolver interface (there can be
  // only one "un-named item being inserted", so inserting Z[ wFi ] while
  // inserting the new item is complicated

  auto wh = BStrategy( wFi );

  // get the space for the item from the MPSolver - - - - - - - - - - - - - -

  auto G1 = Master->GetItem( wFi + 1 );

  // fetch the item from the Oracle - - - - - - - - - - - - - - - - - - - - -

  fwFi->get_linearization_coefficients( G1 );
  auto Alfa1k = fwFi->get_linearization_constant();
  HpNum eps;

  // pass the base to the MP Solver - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set SGBse = nullptr;
  Master->SetItemBse( SGBse , NumVar );

  // compute ScPr1k and Alfa1k- - - - - - - - - - - - - - - - - - - - - - - -

  Index cp;
  HpNum ScPr1k;

  if( diagonal ) {  // it is a subgradient
   // update alpha value at Lambda1 point
   Alfa1k = UpFiLmb1[ wFi ] - Alfa1k -
    std::inner_product( Lambda1.begin() , Lambda1.end() , G1 , double( 0 ) );
   eps = Alfa1k;  // this is how much G1 is an eps-subgradent in Lambda1
   // CheckSubG changes Alfa1k so that G1 is an Alfa1k-subgradent in Lambda
   cp = Master->CheckSubG( UpFiLmb1[ wFi ] - UpRifFi[ wFi ] , t ,
			   Alfa1k , ScPr1k );
   }
  else              // it is a constraint
   cp = Master->CheckCnst( Alfa1k , ScPr1k , Lambda.data() );

  Index gpp = Inf<Index>();  // position in the global pool where to put it

  if( f_log && ( LogVerb > 2 ) ) {
   *f_log << std::endl << "            New ";
   if( diagonal ) {
    if( eps >= std::max( std::abs( UpRifFi[ wFi ] ) , double( 1 ) )
	       * RelAcc / 10 )
     *f_log << "eps-subgradient with eps = " << eps;
    else
     *f_log << "subgradient";
    *f_log << " for Fi[ " << wFi << " ] ~ Alfa1 = " << Alfa1k
	   << " ~ gd = " << - ScPr1k;
    }
   else
    *f_log << "constraint " << wh << " ~ rhs = " << Alfa1k;
   }

  bool to_insert = true;  // if it has to be inserted

  if( cp < InINF ) {  // the item is a copy - - - - - - - - - - - - - - - - -
   BLOG( 2 , " is copy of " << cp << " (" << ItemVcblr[ cp ].second << ")" );

   wh = cp;  // we have it already

   auto OldA1k = (Master->ReadLinErr())[ cp ];

   assert( ( ItemVcblr[ cp ].first == wFi ) &&
           ( ItemVcblr[ cp ].second < vBPar2[ wFi ] ) &&
	   ( InvItemVcblr[ wFi ][ ItemVcblr[ cp ].second ] == cp ) );

   if( OldA1k >= Alfa1k + std::max( std::abs( Alfa1k ) , double( 1 ) )
                          * RelAcc / 10 ) {
    // if the copy has a *substantially* smaller Alfa than the original,
    // replace the original with the copy; in principle relative differences
    // smaller than RelAcc could be ignored, but we use RelAcc / 10 for safety

    BLOG( 2 , " with smaller Alfa" );

    gpp = ItemVcblr[ cp ].second;
    if( ( BPar7 & 3 ) < 3 ) {
     // BundleSolver does not immediately replace the copy unless necessary,
     // but clearly if one linearization in the global pool has to be
     // sacrificed, it'll be the copy
     auto ngpp = find_place_in_global_pool( wFi );
     if( ngpp < Inf<Index>() ) {       // a free place has been found
      // although the old linearization is kept, it is removed from the
      // bundle: the position cp is now associated with ngpp, which
      // means that position gpp is now free
      remove_from_global_pool( wFi , gpp , false );
      gpp = ngpp;                      // store the copy there
      BLOG( 2 , " (" << gpp << ")" );  // print the chosen place
      }
     }

    Master->SubstItem( cp );  // substitute it in the master problem
    // note that the number of items of component wFi in the master problem
    // is unchanged
    }
   else                 // the item is a copy, not better than the original
    to_insert = false;  // do nothing

   BLOG( 2 , std::endl );
   }
  else {           // the item is not a copy- - - - - - - - - - - - - - - - -
   // insert the item, if there is space

   if( wh == InINF )  // the position has not been selected in BStrategy()
    wh = FindAPlace( wFi );  // find a free spot in the bundle

   if( wh == InINF ) {  // no space found ...
    if( ! Ftchd ) {     // ... and this was the first item
     BLOG( 0 , std::endl << " ERROR: No space in the bundle" << std::endl );
     Result = kError;   // signal an error to end the outer Fi-cycle
     }
    else
     BLOG( 1 , std::endl << " WARNING: No space in the bundle" << std::endl );
    break;              // the cycle ends
    }

   if( ItemVcblr[ wh ].second < vBPar2[ wFi ] )
    // the place is occupied already: this happens if the bundle was full
    // (and, possibly aggregation has been performed for safety)
    Master->RmvItem( wh );  // the old item has to be removed first
   else {                   // the place is unoccupied
    ++NrItems[ wFi ];       // one more item in the bundle (otherwise the
    ++NrItems[ NrFi ];      // number remains the same as one is replaced)
    }

   Master->SetItem( wh );   // insert the new item in the MP Solver

   // now find a position in the glonal pool of component wFi where to store
   // the new linearization
   gpp = find_place_in_global_pool( wFi );

   if( gpp == Inf<Index>() ) {  // there is none
    // this means that not only the global pool is full, but also the bundle
    // also full: one can therefore put it in the very place of the item it
    // replaces, which must be an item of the same component because
    // BStrategy() ensures this
    assert( ItemVcblr[ wh ].first == wFi );    
    gpp = ItemVcblr[ wh ].second;
    assert( gpp < vBPar2[ NrFi ] );  
    }

   BLOG( 2 , " stored in " << wh << " (" << gpp << ")" << std::endl  );
   }

  // in all (subgradient) cases, check and update whisG1- - - - - - - - - - -

  if( diagonal ) {     // it is a subgradient
   if( ( whisG1[ wFi ] == InINF ) || ( Alfa1k < Alfa1[ wFi ] ) ||
       ( ( Alfa1k == Alfa1[ wFi ] ) && ( ScPr1k > ScPr1[ wFi ] ) ) ) {
    whisG1[ wFi ] = wh;  // wh is the new representative of wFi
    Alfa1[ wFi ] = Alfa1k;
    ScPr1[ wFi ] = ScPr1k;
    }
   }
 
  // if something was inserted, bookkeeping is needed - - - - - - - - - - - -

  if( to_insert ) {

   inhibit_Modification( true );
   v_c05f[ wFi ]->store_linearization( gpp );
   inhibit_Modification( false );

   add_to_global_pool( wFi , gpp , wh );

   if( diagonal )       // it is a subgradient
    OOBase[ wh ] = -1;  // ensure it won't be touched again this round
   else                 // it is a constraint
    // mark it as permanently fixed: this may be a bad choice in practice,
    // although it is required by the theory (we'll see ...)
    OOBase[ wh ] = - Inf<SIndex>();
   }

  #if CHECK_DS & 1
   CheckBundle();
  #endif

  }  // end( items-collecting loop )- - - - - - - - - - - - - - - - - - - - -

 // update lower and upper estimates  - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( ( LwFiLmb1[ NrFi ] > LwTrgt ) || ( UpFiLmb1[ NrFi ] < UpTrgt ) );

 }  // end( BundleSolver::FiAndGi() )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::GotoLambda1( void )
{
 std::vector<VarValue> DeltaFi( NrFi + 1 );  // DeltaFi = UpFiLmb1 - UpRifFi
 std::transform( UpFiLmb1.begin() , UpFiLmb1.end() , UpRifFi.begin() ,
		 DeltaFi.begin() , std::minus<double>() );

 // do the move - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // Lambda = Lambda1

 Lambda.swap( Lambda1 );
 UpFiLmb.swap( UpFiLmb1 );
 LwFiLmb.swap( LwFiLmb1 );
 UpRifFi = UpFiLmb;

 // change the current point in the MP Solver - - - - - - - - - - - - - - - -

 Master->ChangeCurrPoint( t , DeltaFi.data() );

 // signal that Alfa1[] is not reliable - - - - - - - - - - - - - - - - - - -

 Alfa1.assign( NrFi + 1 , Inf<double>() );

 }  // end( GotoLambda1 )

/*--------------------------------------------------------------------------*/

void BundleSolver::SimpleBStrat( void )
{
 if( ( BPar7 & 3 ) == 3 ) {  // "eager" deletion
  std::vector<Subset> tbdltd( NrFi );
  for( Index i = 0 ; i < Master->MaxName() ; ++i )
   if( ( OOBase[ i ] < Inf<SIndex>() ) &&
       ( OOBase[ i ] > SIndex( BPar1 ) ) ) {
    tbdltd[ ItemVcblr[ i ].first ].push_back( ItemVcblr[ i ].second );
    Delete( i );
    }

  inhibit_Modification( true );
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ! tbdltd[ k ].empty() )
    v_c05f[ k ]->delete_linearizations( std::move( tbdltd[ k ] ) , false );
  inhibit_Modification( false );
  }
 else                        // "lazy" deletion
  for( Index i = 0 ; i < Master->MaxName() ; ++i )
   if( ( OOBase[ i ] < Inf<SIndex>() ) && ( OOBase[ i ] > SIndex( BPar1 ) ) )
    Delete( i );

 #if CHECK_DS & 1
  CheckBundle();
 #endif

 }  // end( BundleSolver::SimpleBStrat )- - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

double BundleSolver::BetaK( Index wFi ) { return( 1.0 / double( NrEasy ) ); }

/*--------------------------------------------------------------------------*/

void BundleSolver::Log1( void )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "{" << SCalls << "-" << ParIter << "-"
	<< NrItems[ NrFi ] << "-" << MBDim << "} t = " << t
	<< " ~ D*_1( z* ) = " << Master->ReadDStart( 1 )
	<< " ~ Sigma = " << Sigma << std::endl << "           ";

 *f_log <<  " Fi = ";

 if( UpFiLmb[ NrFi ] == Inf<double>() )
  *f_log << " - INF";
 else
  *f_log << UpFiLmb[ NrFi ] << " ~ eU = " << EpsU;

 if( BPar6 )
  *f_log << " ~ BP3 = " << aBP3;

 } // end( BundleSolver::Log1 )  - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::Log2( void )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "            ";

 if( LowerBound[ NrFi ] > - Inf<double>() )
  *f_log << "LB = " << LowerBound[ NrFi ] << " ~ ";

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
 //
 // this is called *before* that we know if the place will actually be
 // required, so it has a "loose attitude"; in particular, it will return
 // InINF under two opposite set of conditions:
 //
 // - there is plenty of space left in the bundle, so that no B-strategy (no
 //   removal or aggregation) is required;
 //
 // - there is no way in which space can be found, i.e., the bundle (for this
 //   component) is full, and removal/aggregation are not successful (a very
 //   strange occurrence due to an extremely small bundle or it being chock
 //   full of constraints)
 //
 // Picking a specific spot in the free space is the task of FindAPlace(),
 // which however is not called right away because the place may end up not
 // being needed. If BStrategy() returne InINF because there is plenty of
 // space then FindAPlace() will suceed, if BStrategy() returne InINF because
 // there is no way space can be found then FindAPlace() will fail and it
 // will be clear that disaster looms
 //
 // important note: if the returned wFi is not InINF, *and* it is the name
 // of an item still in the bundle, then:
 //
 // - the item belongs to the component wFi, so as to keep the number of
 //   items of that component constant
 //
 // - the item is *not* deleted, since it is not 100% sure this will need
 //   to be done, so deletion will be responsibility of the caller

 // there are "free" items in the global pool: there is "plenty of space"
 if( FrFItem[ wFi ] < vBPar2[ wFi ] )
  return( InINF );

 // there is not plenty of space, take 1- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // there are no "free" items in the gloal pool, but there are items in the
 // global pool that are not in the bundle; these will go first

 if( NrItems[ wFi ] < vBPar2[ wFi ] )
  return( InINF );

 // there is not plenty of space, take 2- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the bundle for component wFi is full, which implies that the global pool
 // is also full: among the items in the bundle for component wFi, find the
 // removable one with largest OOBase[], snd among these with the largest
 // OOBase[] select that with largest Alfa[]
 // note: the Z[ wFi ] in the bundle (if any) is not removable and
 //       therefore cannot be selected, which in particular happens if
 //       wFi has only *one* subgradient in base

 Index wh;
 SIndex OOwh = - Inf<SIndex>();
 HpNum Awh = -Inf<HpNum>();
 cHpRow tA = Master->ReadLinErr();
 for( auto i : InvItemVcblr[ wFi ] ) {
  assert( i < vBPar2[ NrFi ] );
  if( ( OOBase[ i ] > OOwh ) ||
      ( ( OOBase[ i ] == OOwh ) && ( tA[ i ] > Awh ) ) ) {
   wh = i;
   OOwh = OOBase[ i ];
   Awh = tA[ i ];
   }
  }

 if( OOBase[ wh ] < 0 )    // all items are non-removable: nothing else to
  return( InINF );         // do (except maybe complaining very loudly)
 else
  if( OOBase[ wh ] > 0 )   // a place is found
   return( wh );           // there are no problems, all done

 // wh is a basic item- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this means that *all* items of component wFi are either in base or not
 // removable, for otherwise we would have selected an item with OOBase > 0;
 // we cannot discard anything before having performed aggregation, but in
 // order to do so we also need to free some space for the Z[ wFi ]
 //
 // note: there is an easy case where z is "naturally" in the base without
 //       any aggregation: there is only one subgradient in base (for the
 //       component wFi), and therefore, its Mlt[] is == 1. however this
 //       can never happen here because it would mean vBPar2[ wFi ] == 1
 //       which is not permitted. the other case in which this could happen
 //       is that all the items in the bundle for component wF are
 //       constraints, but then aggregation would not be useful (in fact
 //       this is checked and reported as a failure)
 //
 // note: this also means that wh is the item in base with largest Alpha
 //
 // note: that since *all* items of this component are either in base or
 //       not removable, we can scan MBse[] for the items to be removed,
 //       possibly ignoring items with Mlt[] == 0 -- but in fact not doing
 //       it because there is not any

 Index MBDm;
 cIndex_Set MBse;
 cHpRow Mlt = Master->ReadMult( MBse , MBDm , wFi + 1 , false );

 if( ! Zvalid[ wFi ] ) {
  // a valid Z[ wFi ] is not already in: aggregation has to be performed

  Index whZ = InINF;  // the position where Z[ wFi ] has to go

  if( ( whisZ[ wFi ] < InINF ) && Master->IsSubG( whisZ[ wFi ] ) )
   whZ = whisZ[ wFi ];  // preferably re-use the last position
  else {
   // there is no last position for Z[ wFi ], choose the one with min Mult
   // among all the removable ones different from wh

   cHpRow tMlt = Mlt;
   HpNum tMin = Inf<HpNum>();
   if( MBse ) {
    cIndex_Set tMBse = MBse;
    for( Index h ; ( h = *(tMBse++) ) < InINF ; ++tMlt )
     if( ( h != wh ) && ( *tMlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
      whZ = h;
      tMin = *tMlt;
      }
    }
   else
    for( Index h = 0 ; h < MBDm ; ++h , ++tMlt )
     if( ( h != wh ) && ( *tMlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
      whZ = h;
      tMin = *tMlt;
      }
   }

  if( whZ == InINF )  // there is no removable item apart from wh
   return( InINF );   // nothing else to do except complaining very loudly

  // tell the C05Function what is going to happen - - - - - - - - - - - - - -
  // note that this only happens when the bundle (for component wFi) is "very
  // full", and therefore also the global pool (for component wFi) is such.
  // hence, whZ is an item already in the bundle, and therefore in the global
  // pool. the natural choice is to put the new aggregate linearization in the
  // same position in the global pool where whZ was, i.e.,
  // ItemVcblr[ whZ ].second. hence, ItemVcblr, InvItemVcblr and so on need
  // not be changed
  
  LinearCombination coeff( MBDm );
  if( MBse )
   for( Index i = 0 ; i < MBDm ; ++i ) {
    coeff[ i ].first = ItemVcblr[ MBse[ i ] ].second;
    coeff[ i ].second = Mlt[ i ];
    }
  else
   for( Index i = 0 ; i < MBDm ; ++i ) {
    coeff[ i ].first = ItemVcblr[ i ].second;
    coeff[ i ].second = Mlt[ i ];
    }

  inhibit_Modification( true );
  v_c05f[ wFi ]->store_combination_of_linearizations( coeff ,
						   ItemVcblr[ whZ ].second );
  inhibit_Modification( false );

  Master->RmvItem( whZ );  // remove the old item in position whZ

  // ask the MPSolver the memory for keeping Z[ wFi ] - - - - - - - - - - - -
  // note: Mlt and MBse could very well be "temporary" memory belonging to the
  // MPSolver, and any call to a method of the MPSolver may invalidate it;
  // the calls start now, and in fact MBse and Mlt are no longer used

  SgRow tZ = Master->GetItem( wFi + 1 );

  // read Z[ wFi ]- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index ZBDm;
  cIndex_Set ZBse;
  Master->ReadZ( tZ , ZBse , ZBDm , wFi + 1 );

  #if CHECK_DS & 2
   std::vector<VarValue> Z( NumVar );
   v_c05f[ wFi ]->get_linearization_coefficients( Z.data() ,
						  Range( 0 , NumVar ) ,
						  ItemVcblr[ whZ ].second );
   if( ZBse ) {
    Index j = 0;
    for( Index i = 0 ; i < NumVar ; ++i )
     if( ( j < ZBDm ) && ( ZBse[ j ] == i ) ) {
      if( std::abs( Z[ i ] - tZ[ j ] ) >=
	  RMPAccSol * std::max( Z[ i ] , double( 1 ) ) )
       std::cerr << "CZ[ " << i << " ] = " << Z[ i ]
		 << " ~ MZ[ " << i << " ] = " << tZ[ j ] << std::endl;
      ++j;
      }
     else
      if( std::abs( Z[ i ] ) >= RMPAccSol )
       std::cerr << "CZ[ " << i << " ] = " << Z[ i ]
		 << " ~ MZ[ " << i << " ] = 0" << std::endl;
    }
   else
    for( Index i = 0 ; i < NumVar ; ++i )
     if( std::abs( Z[ i ] - tZ[ i ] ) >=
	 RMPAccSol * std::max( Z[ i ] , double( 1 ) ) )
      std::cerr << "CZ[ " << i << " ] = " << Z[ i ]
		<< " ~ MZ[ " << i << " ] = " << tZ[ i ] << std::endl;
  #endif
 
  // now pass Z[ wFi ] back to the MP Solver- - - - - - - - - - - - - - - - -

  Master->SetItemBse( ZBse , ZBDm );

  HpNum ScPri;
  HpNum Ai = Master->ReadSigma( wFi + 1 );      // its alfa is Sigma[ wFi ]

  #if CHECK_DS & 2
   HpNum tAi = v_c05f[ wFi ]->get_linearization_constant(
						    ItemVcblr[ whZ ].second );
   tAi = UpRifFi[ wFi ] - tAi -
                          std::inner_product( Lambda.begin() , Lambda.end() ,
					      Z.begin() , double( 0 ) );

   if( std::abs( tAi - Ai ) >=
       RMPAccSol * std::max( std::max( Ai , UpRifFi[ wFi ] ) , double( 1 ) ) )
    std::cerr << "Csigma = " << tAi << " ~ Msigma = " << Ai << std::endl;
  #endif

  // note that Tau == -1, meaning that Ai need not be changed since
  // Ai is already the linearization error in Lambda, but still the
  // ScPri need be computed
  Master->CheckSubG( 0 , -1 , Ai , ScPri );
  
  Master->SetItem( whZ );  // set Z[ wFi ] in position whZ

  whisZ[ wFi ] = whZ;      // Z[ wFi ] is in the bundle in position whZ
  Zvalid[ wFi ] = true;    // ... and it is valid
  OOBase[ whZ ] = -1;      // ... and it won't be removed in this iteration

  BLOG( 2 , std::endl << "Aggregation performed into " << whZ );
  }

 // at this point, Z[ wFi ] is in the bundle, hence it is safe to replace wh
 // with anything the oracle provides us

 return( wh );

 }  // end( BundleSolver::BStrategy ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::FindAPlace( Index wFi )
{
 // this method is used to return the index of an available position in the
 // bundle where to store a new item belonging to "component" wFi; if there
 // are no possible positions left, then InINF is returned
 //
 // note that positions in the bundle are not "reserved to components", so
 // the task of FindAPlace() is easy. it is the business of other parts of
 // the code to ensure that a component does not use more than its fair
 // share of positions in the bundle

 Index wh = InINF;

 if( FreList.size() ) {     // there are deleted items
  wh = FreList.top();       // pick the first one (smaller name)
  FreList.pop();
  }
 else                       // there are no deleted items ...
  if( Master->MaxName() < vBPar2[ NrFi ] )
                            // ... but there is still space
   wh = Master->MaxName();  // next name

 return( wh );

 }  // end( BundleSolver::FindAPlace )- - - - - - - - - - - - - - - - - - - -

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

void BundleSolver::guts_of_destructor( void )
{
 if( Master ) {
  Master->SetDim();
  delete Master;
  Master = nullptr;
  }

 Alfa1.clear();
 ScPr1.clear();
 whisG1.clear();
 vStar.clear();

 LowerBound.clear();
 LwFiLmb.clear();
 UpFiLmb.clear();
 LwFiLmb1.clear();
 UpFiLmb1.clear();
 UpRifFi.clear();
 
 FiStatus.clear();

 CurrNrEvls.clear();

 Zvalid.clear();
 whisZ.clear();
 FreList = {};

 MaxItem.clear();
 FrFItem.clear();
 NrItems.clear();

 InvItemVcblr.clear();
 ItemVcblr.clear();

 OOBase.clear();

 LmbdBst.clear();
 Lambda1.clear();
 Lambda.clear();

 vBPar2.clear();

 for( auto milpp : MILP_s )
  delete( milpp );
 MILP_s.clear();
 IsEasy.clear();

 LamVcblr.clear();

 v_c05f.clear();

 }  // end( BundleSolver:guts_of_destructor ) - - - - - - - - - - - - - - - -

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

 if( RstLvl & RstCrr )    // get an initial point - - - - - - - - - - - - - -
  for( Index i = 0 ; i < NumVar ; ++i )
   Lambda[ i ] = LamVcblr[ i ]->get_value();
 else                     // reset the current point to all-0 - - - - - - - - -
  Lambda.assign( NumVar , 0 );

 }  // end( BundleSolver::ReSetAlg ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::Delete( cIndex i , bool ModDelete )
{
 // deletes from the bundle the item in position i
 //
 // ModDelete == true means that this is called in response of a Modification
 // where the linearization has been removed from the global pool
 //
 // whatever BPar7 says, no linearization is physically deleted here inside,
 // this has to be done by the caller (if needed)

 cIndex k = ItemVcblr[ i ].first;

 // check if this item was the "representative" for its component - - - - - -

 if( whisG1[ k ] == i )  // it is the representative of k
  whisG1[ k ] = InINF;   // a new representative is needed

 // check if this item was the z* for its component - - - - - - - - - - - - -

 if( whisZ[ k ] == i ) {  // it is the aggregate subgradient of k
  whisZ[ k ] = InINF;     // no aggregate subgradient is in the bundle
  Zvalid[ k ] = false;    // a fortiori, no valid one
  }

 // delete the item from the MP - - - - - - - - - - - - - - - - - - - - - - -

 Master->RmvItem( i );

 BLOG( 2 , std::endl << "Item " << i << " removed" );

 // bookkeeping of internal data structures - - - - - - - - - - - - - - - - -
 // note that any item whose name is >= Master->MaxName() is surely not in
 // the bundle (master problem), and therefore it need not be in FreList

 cIndex MxNm = Master->MaxName();
 if( i < MxNm )
  FreList.push( i );

 OOBase[ i ] = Inf<SIndex>();
 --NrItems[ k ];
 --NrItems[ NrFi ];

 // remove from the global pool: the removal is "hard" if either BPar7 says
 // so, or the linearization had been deleted anyway

 remove_from_global_pool( k , ItemVcblr[ i ].second ,
			  ( ( BPar7 & 3 ) == 3 ) || ModDelete );
 ItemVcblr[ i ].second = Inf<Index>();

 // check if compacting FreList is appropriate- - - - - - - - - - - - - - - -
 // the issue with having indices of "free" position in the bundle stored in
 // a priority_queue is the following: if the bundle gets "full", but then is
 // "emptied", FreList may end up containing "many" elements, and in
 // particular elements that are >= Master->MaxName(), which therefore are
 // useless since they are obviously not in the bundle. the check above tries
 // to avoid that, but it may clearly fail (say, if small items are deleted
 // before large ones). checking if there are items with name >=
 // Master->MaxName() in FreList and deleting them is not cheap. the only
 // easy-to-check case is the one where FreList.size() > Master->MaxName():
 // if this happens, FreList is cleared and re-initialized

 if( FreList.size() > MxNm ) {
  FreList = {};
  for( Index h = 0 ; h < MxNm ; ++h )
   if( ItemVcblr[ h ].second == Inf<Index>() )
    FreList.push( h );
  }
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
   if( UpFiLmb[ NrFi ] > -Inf<double>() )
    aBP3 = ( BPar5 > 0 ? BPar4 : BPar3 ) +
           Index( BPar5 / std::sqrt( EpsU / RelAcc ) );
   break;
  case( 2 ):
   if( UpFiLmb[ NrFi ] > -Inf<double>() )
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

void BundleSolver::remove_from_global_pool( Index k , Index i , bool hard )
{
 // update InvItemVcblr and all the associated fields to the fact that the
 // linearization currently in position i of the global pool of component k
 // is removed; actually, there may as well not be any linearization in
 // position i already
 //
 // the removal is "hard" (meaning the linearization is actually deleted) if
 // hard == true, and "soft" (the linearization is kept in the global pool
 // but deleted from the bundle) otherwise
 //
 // the actual removal from the global pool is not handled here

 InvItemVcblr[ k ][ i ] = hard ? Inf<Index>() : vBPar2[ NrFi ];
 while( MaxItem[ k ] &&
	( InvItemVcblr[ k ][ MaxItem[ k ] - 1 ] == Inf<Index>() ) )
  --MaxItem[ k ];
 if( i < FrFItem[ k ] )   // creating a new "hole" before the FrFItem
  FrFItem[ k ] = i;       // this is the new FrFItem
 else                     // deleting something that may be FrFItem
  while( FrFItem[ k ] &&
	 ( InvItemVcblr[ k ][ FrFItem[ k ] - 1 ] >=
	   ( ( BPar7 & 3 ) ? vBPar2[ NrFi ] : Inf<Index>() ) ) )
   --FrFItem[ k ];

 }  // end( remove_from_global_pool )

/*--------------------------------------------------------------------------*/

Index BundleSolver::find_place_in_global_pool( Index k )
{
 // returns a suitable position in the global pool of component k, or InINF
 // if there is no free space

 if( FrFItem[ k ] < vBPar2[ k ] )  // there are free positions
  return( FrFItem[ k ] );          // return the first of them

 // if there are no free positions but there are less items in the bundle
 // (for component k) than the size of the global pool, there must be a
 // some linearization in the global pool that is not in the bundle: find
 // and return the first one (smallest position)
 Index gpp = Inf<Index>();
 if( NrItems[ k ] < vBPar2[ k ] )
  for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
   if( InvItemVcblr[ k ][ i ] >= vBPar2[ NrFi ] ) {
    gpp = i;
    break;
    }

 return( gpp );

 }  // end( find_place_in_global_pool )

/*--------------------------------------------------------------------------*/

void BundleSolver::add_to_global_pool( Index k , Index i , Index wh )
{
 // update ItemVcblr, InvItemVcblr and all the associated fields the to fact
 // that the linearization in position i in the global global_pool of
 // component k will be kept in the bundle at position wh
 //
 // the actual addition to the global pool is not handled here

 ItemVcblr[ wh ].first = k;
 ItemVcblr[ wh ].second = i;

 InvItemVcblr[ k ][ i ] = wh;
 while( ( MaxItem[ k ] < vBPar2[ k ] ) &&
	( InvItemVcblr[ k ][ MaxItem[ k ] ] < Inf<Index>() ) )
  ++MaxItem[ k ];
 while( ( FrFItem[ k ] < MaxItem[ k ] ) &&
	( InvItemVcblr[ k ][ FrFItem[ k ] ] <
	  ( ( BPar7 & 3 ) ? vBPar2[ NrFi ] : Inf<Index>() ) ) )
  ++FrFItem[ k ];

 }  // end( BundleSolver::add_to_global_pool( k , i , wh ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::add_to_global_pool( Index k , Index i )
{
 // update InvItemVcblr and all the associated fields the to fact that there
 // is a linerization in position i in the global global_pool of component k,
 // although there is no corresponding item in the bundle

 InvItemVcblr[ k ][ i ] = vBPar2[ NrFi ];
 while( ( MaxItem[ k ] < vBPar2[ k ] ) &&
	( InvItemVcblr[ k ][ MaxItem[ k ] ] < Inf<Index>() ) )
  ++MaxItem[ k ];
 if( ! ( BPar7 & 3 ) )  // if items not in the bundle are not "free"
  return;               // nothing else to do
 while( ( FrFItem[ k ] < MaxItem[ k ] ) &&
	( InvItemVcblr[ k ][ FrFItem[ k ] ] < vBPar2[ NrFi ] ) )
  ++FrFItem[ k ];

 }  // end( BundleSolver::add_to_global_pool( k , i ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::add_to_bundle( Index k , Index i )
{
 // add to the bundle (master problem) the item corresponding to the
 // linearization to be found at position i in the global pool of component
 // k; this assumes that the linearization is already there in the global
 // pool. if InvItemVcblr[ k ][ i ] < vBPar2[ NrFi ], i.e., the item is
 // already in the bundle, then it is replaced, otherwise it is added
 //
 // note that CheckSubG() or CheckCnst() need be called, but even if the
 // item is identical to some in the bundle already this information is
 // ignored and the item is inserted anyway; hence, if the check is active,
 // it is temporarily deactivated (and then re-activated)

 auto wh = InvItemVcblr[ k ][ i ];
 if( wh >= vBPar2[ NrFi ] ) {  // the item is not there already
  wh = FindAPlace( k );        // find a "free" spot in the bundle
  if( wh == InINF )            // one must be there
   throw( std::logic_error( "no space found in the bundle" ) );

  ++NrItems[ k ];              // keep count
  ++NrItems[ NrFi ];
  add_to_global_pool( k , i , wh );  // update dictionaries
  }
 else                          // the item is there already
  Master->RmvItem( wh );       // remove it so that it can be replaced

 // ask the MPSolver for the space to write the item to
 auto G1 = Master->GetItem( k + 1 );

 // recover the linearization from the C05Function
 v_c05f[ k ]->get_linearization_coefficients( G1 , Range( 0 , NumVar ) , i );

 // recover the constant and "translate" it w.r.t. Lambda
 auto Ai = v_c05f[ k ]->get_linearization_constant( i );

 Master->SetItemBse( nullptr , NumVar );

 if( MPName & 8 )                   // if checking for copies is active
  Master->CheckIdentical( false );  // temporarily de-activate it now

 double ScPri;
 if( v_c05f[ k ]->is_linearization_vertical( i ) )
  Master->CheckCnst( Ai , ScPri , Lambda.data() );
 else {
  Ai = UpRifFi[ k ] - Ai -
   std::inner_product( Lambda.begin() , Lambda.end() , G1 , double( 0 ) );
  Master->CheckSubG( 0 , 0 , Ai , ScPri );
  }

 if( MPName & 8 )                   // if checking for copies is active
  Master->CheckIdentical();         // re-activate it now

 Master->SetItem( wh );  // add the item to the master problem

 }  // end( BundleSolver::add_to_bundle )

/*--------------------------------------------------------------------------*/

void BundleSolver::reset_bundle( void )
{
 // completely resets the bundle, because a (bunch of) Modification(s) saying
 // so has(ve) been received. this only affects the BundleSolver data
 // structures and the MPSolver, not the C05Function(s)

 OOBase.assign( vBPar2[ NrFi ] , Inf<SIndex>() );

 ItemVcblr.assign( vBPar2[ NrFi ] , make_pair( InINF , InINF ) );

 for( Index k = 0 ; k < NrFi ; ++k )
  InvItemVcblr[ k ].assign( vBPar2[ k ] , InINF );

 NrItems.assign( NrFi + 1 , 0 );
 FrFItem.assign( NrFi , 0 );
 MaxItem.assign( NrFi , 0 );

 FreList = {};
 whisZ.assign( NrFi , InINF );
 Zvalid.assign( NrFi , false );

 whisG1.assign( NrFi , InINF );

 Master->RmvItems();
 }

/*--------------------------------------------------------------------------*/

bool BundleSolver::is_special_GroupMod( GroupModification & gmod )
{
 // recognise "special" GroupModification for changing the set of "active"
 // Variable of all the Objective at the same time; note that these
 // contain FunctionModVars* not necessarily C05FunctionModVars* because
 // the Modification may not be strongly quasi-additive

 if( gmod.v_sub_Modifications.size() != NrFi + ( f_lf ? 1 : 0 ) )
  return( false );

 auto smi = gmod.v_sub_Modifications.begin();
 auto sm0 = *(smi++);
 for( ; smi !=  gmod.v_sub_Modifications.end() ; ++smi )
  if( typeid( sm0 ) != typeid( *smi ) )
   return( false );

 smi = gmod.v_sub_Modifications.begin();
 ++smi;

 // check FunctionModVarsAddd
 {
  const auto mod0 = std::dynamic_pointer_cast<FunctionModVarsAddd>( sm0 );
  if( mod0 ) {
   for( ; smi != gmod.v_sub_Modifications.end() ; ++smi ) {
    auto modi = std::static_pointer_cast<FunctionModVarsAddd>( *smi );
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
   for( ; smi != gmod.v_sub_Modifications.end() ; ++smi ) {
    auto modi = std::static_pointer_cast<FunctionModVarsRngd>( *smi );
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
   for( ; smi != gmod.v_sub_Modifications.end() ; ++smi ) {
    auto modi = std::static_pointer_cast<FunctionModVarsSbst>( *smi );
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

void BundleSolver::process_outstanding_Modification( void )
{
 // multiple-loop version, where several passes are done in order to gather
 // which kind of Modification have occurred and avoid doing costly work
 // more than once

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
 //   all *FunctionMod*)
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

 bool to_delete;  // should have been defined inside, but there is not
                  // visible by the lambda

 for( auto rimod = v_mod_tmp.rbegin() ; rimod != v_mod_tmp.rend() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to rimod
      [ & to_delete , & v_mod_tmp ]( decltype( rimod ) & ri ) {
       if( to_delete )
	ri = std::reverse_iterator( v_mod_tmp.erase( std::next( ri ).base()
						     ) );
       else
	++ri;
       }( rimod ) ) {
  to_delete = false;
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
	if( ttmod->which().empty() ) {
	 reset[ wFi ] = true;
	 to_delete = true;
	 }
	continue;
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
     const auto ttmod =
                    std::dynamic_pointer_cast<C05FunctionModVarsAddd>( tmod );
     if( ttmod )
      continue;
     }

    {
     const auto ttmod =
                    std::dynamic_pointer_cast<C05FunctionModVarsRngd>( tmod );
     if( ttmod )
      continue;
     }

    {
     const auto ttmod =
                    std::dynamic_pointer_cast<C05FunctionModVarsSbst>( tmod );
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
     const auto ttmod =
                    std::dynamic_pointer_cast<C05FunctionModVarsAddd>( fmod );
     if( ttmod )
      continue;
     }

    {
     const auto ttmod =
                    std::dynamic_pointer_cast<C05FunctionModVarsRngd>( fmod );
     if( ttmod )
      continue;
     }

    {
     const auto ttmod =
                    std::dynamic_pointer_cast<C05FunctionModVarsSbst>( fmod );
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
 // reset) then delete all corresponding items in the bundle (all the bundle)

 if( std::find( reset.begin() , reset.end() , false ) == reset.end() )
  reset_bundle();  // all components have been reset
 else
  if( std::find( reset.begin() , reset.end() , true ) != reset.end() ) {
   // at least a component has been reset

   for( Index k = 0 ; k < NrFi ; ++k )
    if( reset[ k ] ) {
     for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
      if( InvItemVcblr[ k ][ i ] < vBPar2[ NrFi ] )
       Delete( InvItemVcblr[ k ][ i ] , true );
     }
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
 // is useless.
 //
 // note that, due to limitations in the MPSolver interface, changing a
 // linearization implies changing its constant; therefore, reset[ k ] == true
 // means that everything is changing. when reset[ k ] == true we take
 // AlphaC[ k ] == false since in this case che change of constants is not
 // needed because it is implied by the soft reset. thus, AlphaC[ k ] == true
 // means that only the constants need be changed, but not all the rest
 //
 // note that Modification changing the linearizations happening *after* a
 // "soft" reset of the global pool (meaning it is found *before* in the
 // reverse order) is also useless, since a reset forces the re-reading of all
 // linearizations, which by definition happens at their current (final)
 // state. yet, this is not done immediately
 //
 // note that we make no serious attempt at keeping track of the combined
 // effect of all changes, in order to detect if a large set of small
 // changes actually imples a reset. this is complicated for "horizontal"
 // changes (for all linearizations, a range/subset of entries) because the
 // names of the changed Variable may not be current (additions/deletions may
 // happen in the meantime), and keeping track is too burdensome. similarly
 // for "vertical" changes (a set of specific linearizations). some steps
 // in this direction will perhaps be done in later stages of develpment

 reset.assign( NrFi , false );  // reset reset (couldn't resist)
 std::vector<bool> AlphaC( NrFi , false );

 for( auto rimod = v_mod_tmp.rbegin() ; rimod != v_mod_tmp.rend() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to rimod
      [ & to_delete , & v_mod_tmp ]( decltype( rimod ) & ri ) {
       if( to_delete )
	ri = std::reverse_iterator( v_mod_tmp.erase( std::next( ri ).base()
						     ) );
       else
	++ri;
       }( rimod ) ) {
  to_delete = false;
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
   //
   // the Modification can also change constants, so this has to be checked.
   // if the Modification changes all the constants and the component is
   // reset already, it can be deleted

   const auto tmod = std::dynamic_pointer_cast<C05FunctionModRngd>( mod );
   if( tmod ) {
    auto wFi = get_index_of_component( tmod->function() );
    switch( tmod->type() ) {
     case( C05FunctionMod::AllEntriesChanged ):
      if( reset[ wFi ] )            // component reset already
       to_delete = true;            // nothing else to do
      continue;
     case( C05FunctionMod::AllLinearizationChanged ):
      if( tmod->which().empty() ) {  // reset of the constants
       if( reset[ wFi ] )            // component reset already
	to_delete = true;            // nothing else to do
       else                          // component not reset
	AlphaC[ wFi ] = true;        // reset the constants
       }
      continue;
     default:  // this must not happen
      throw( std::invalid_argument( "wrong type() in C05FunctionModRngd" ) );
     }
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
    switch( tmod->type() ) {
     case( C05FunctionMod::AllEntriesChanged ):
      if( reset[ wFi ] )            // component reset already
       to_delete = true;            // nothing else to do
      continue;
     case( C05FunctionMod::AllLinearizationChanged ):
      if( tmod->which().empty() ) {  // reset of the constants
       if( reset[ wFi ] )            // component reset already
	to_delete = true;            // nothing else to do
       else                          // component not reset
	AlphaC[ wFi ] = true;        // reset the constants
       }
      continue;
     default:  // this must not happen
      throw( std::invalid_argument( "wrong type() in C05FunctionModSbst" ) );
     }
    }  // end( if( ttmod ) )
   }  // end C05FunctionModSbst

  {
   // a C05FunctionMod of type AllLinearizationChanged or AllEntriesChanged
   // with which.empty() "soft" resets all the component, and in the former
   // case also Alpha; AlphaChanged only changes the constants (obviously)

   const auto tmod = std::dynamic_pointer_cast<C05FunctionMod>( mod );
   if( tmod ) {
    // we only react to which().empty() 
    if( ! tmod->which().empty() )
     continue;

    auto wFi = get_index_of_component( tmod->function() );

    switch( tmod->type() ) {
     case( C05FunctionMod::AlphaChanged ):
      AlphaC[ wFi ] = true;
      break;
     case( C05FunctionMod::AllEntriesChanged ):
     case( C05FunctionMod::AllLinearizationChanged ):
      reset[ wFi ] = AlphaC[ wFi ] = true;
      break;
     default:  // this must not happen, as GlobalPoolRemoved with
               // which.empty() has been dealt with and deleted before
      throw( std::invalid_argument(
		        "wrong type in C05FunctionMod with empty which()" ) );
     }

    to_delete = true;
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
    reset[ wFi ] = AlphaC[ wFi ] = true;
    to_delete = true;
    }  // end( if( ttmod ) )
   }  // end C05FunctionModLin
  }  // end( 2nd loop, again in reverse )

 // note that even if there were no more Modification to process we could not
 // stop because this means that reset[ k ] == true and/or AlphaC[ k ] == true
 // for some k. in fact v_mod_tmp was not empty(), and elements can be removed
 // from it only if some component is "soft" reset. 

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 3rd loop, forward: prepare for addition/removal/changes of individual
 // linearization for each component by computing the three sets
 // - linearizations that need be removed
 // - linearizations that need be added
 // - linearizations that need be changed
 // - constants that need be changed
 // Note that:
 // - due to limitations in the MPSolver interface, changing a linearization
 //   implies changing its constant; therefore, Cchg[] contains changes in
 //   the constants only and nothing else, which means that has empty
 //   intersection with all other three sets
 // - if a linearization that is added/changed is later removed, it is no
 //   longer added/changed
 // - if a linearization that is removed/changed is later added it is no
 //   longer removed/changed (adding over an existing linearization changes
 //   it anyway, no reason to remove it)
 // - thus, Addd[], Rmvd[], Chgd[] and Cchg[] all have empty intersection
 // - linearization changes to a reset component can be ignored
 // - constant changes when all constants change can be ignored
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

 std::vector< Subset > Addd( NrFi );  // added linearizations
 std::vector< Subset > Rmvd( NrFi );  // removed linearizations
 std::vector< Subset > Chgd( NrFi );  // changed linearizations
 std::vector< Subset > Cchg( NrFi );  // changed constants

 for( auto imod = v_mod_tmp.begin() ; imod != v_mod_tmp.end() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to imod
      [ & to_delete , & v_mod_tmp ]( decltype( imod ) & it ) {
       if( to_delete )
	it =  v_mod_tmp.erase( it );
       else
	++it;
       }( imod ) ) {
  to_delete = false;
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

   switch( tmod->type() ) {
    case( C05FunctionMod::AlphaChanged ):
     assert( ! tmod->which().empty() );
     // this cannot happen, the case has been dealt with and deleted
     if( AlphaC[ wFi ] )  // constants are reset already
      break;              // nothing else to do

     // note that reset[ wFi ] ==> AlphaC[ wFi ], hence surely reset[ wFi ]
     // is false here
     // add to Cchg[ wFi ] the names in tmod->which(), save those that are
     // in either Addd[ wFi ], or Rmvd[ wFi ], or Chgd[ wFi ]
     if( Addd[ wFi ].empty() && Rmvd[ wFi ].empty() && Chgd[ wFi ].empty() )
      set_union_in_place( Cchg[ wFi ] , tmod->which() );
     else {
      // only those constants corresponding to linearizations that are not
      // added or removed or changed need be changed
      Subset tmp( tmod->which() );

      if( ! Addd[ wFi ].empty() )
       set_difference_in_place( tmp , Addd[ wFi ] );

      if( ( ! Rmvd[ wFi ].empty() ) && ( ! tmp.empty() ) )
       set_difference_in_place( tmp , Rmvd[ wFi ] );

      if( ( ! Chgd[ wFi ].empty() ) && ( ! tmp.empty() ) )
       set_difference_in_place( tmp , Chgd[ wFi ] );

      set_union_in_place( Cchg[ wFi ] , std::move( tmp ) );
      }
     break;
    case( C05FunctionMod::AllLinearizationChanged ):
    case( C05FunctionMod::AllEntriesChanged ):
     // AllEntriesChanged is completely equivalent to AllLinearizationChanged
     // since the change in the linearization implies that of the constant
     // if tmod->which().empty(), this must actually be either a
     // C05FunctionModRngd or a C05FunctionModSbst: save it, for it
     // will be dealt with in the next loop
     if( tmod->which().empty() )
      continue;

     if( reset[ wFi ] )  // changes in reset components
      break;             // are ignored

     // add to Chgd[ wFi ] the names in tmod->which(), save those that are
     // in either Addd[ wFi ] or Rmvd[ wFi ]
     if( Addd[ wFi ].empty() && Rmvd[ wFi ].empty() )
      set_union_in_place( Chgd[ wFi ] , tmod->which() );
     else {
      // only those items that are not added and/or removed need be changed
      Subset tmp( tmod->which() );

      if( ! Addd[ wFi ].empty() )
       set_difference_in_place( tmp , Addd[ wFi ] );

      if( ( ! Rmvd[ wFi ].empty() ) && ( ! tmp.empty() ) )
       set_difference_in_place( tmp , Rmvd[ wFi ] );

      set_union_in_place( Chgd[ wFi ] , std::move( tmp ) );
      }

     // changed linearizations imply changed constants
     if( ! AlphaC[ wFi ] )
      set_difference_in_place( Cchg[ wFi ] , Chgd[ wFi ] );

     break;

    case( C05FunctionMod::GlobalPoolAdded ):
     // add to Addd[ wFi ] the names in tmod->which(), and remove them
     // from Rmvd[ wFi ], and from Chgd[ wFi ] if the component is not
     // reset, and from Cchg[ wFi ] is constants are not reset

     set_union_in_place( Addd[ wFi ] , tmod->which() );

     set_difference_in_place( Rmvd[ wFi ] , tmod->which() );

     if( ! reset[ wFi ] )
      set_difference_in_place( Chgd[ wFi ] , tmod->which() );

     if( ! AlphaC[ wFi ] )
      set_difference_in_place( Cchg[ wFi ] , tmod->which() );
     
     break;

    case( C05FunctionMod::GlobalPoolRemoved ):
     // add to Rmvd[ wFi ] the names in tmod->which(), and remove them
     // from Addd[ wFi ], and from Chgd[ wFi ] if the component is not
     // reset, and from Cchg[ wFi ] is constants are not reset

     set_union_in_place( Rmvd[ wFi ] , tmod->which() );

     set_difference_in_place( Addd[ wFi ] , tmod->which() );

     if( ! reset[ wFi ] )
      set_difference_in_place( Chgd[ wFi ] , tmod->which() );

     if( ! AlphaC[ wFi ] )
      set_difference_in_place( Cchg[ wFi ] , tmod->which() );
 
    }  // end( switch( tmod->type() ) )

   to_delete = true;   // delete it

   }  // end( if( ttmod ) )
  }  // end( 3rd loop, forward )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now delete all linearization that need to, if any

 if( std::find_if( Rmvd.begin() , Rmvd.end() ,
		   []( Subset & Rk ) { return( ! Rk.empty() ); }
		   ) != Rmvd.end() ) {
  // at least a component has had lnearizations removed, but note that not
  // all lnearizations need be in the bundle; if they are not they are
  // just removed from the global pool (if they are there)

  for( Index k = 0 ; k < NrFi ; ++k )
   for( auto i : Rmvd[ k ] ) {
    auto h = InvItemVcblr[ k ][ i ];
    if( h < vBPar2[ NrFi ] )
     Delete( h , true );
    else
     remove_from_global_pool( k , i , true );
    }
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
  to_delete = false;
  auto mod = *imod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // not that we do not distinguish C05FunctionModVars* from "plain"
  // FunctionModVars*, since the only difference is whether or not the
  // operation is strongly quasi-additive, i.e., it implies or not a "hard"
  // reset, but this has already been acted upon
  
  {
   // a "naked" FunctionModVars
   auto tmod = std::dynamic_pointer_cast<FunctionModVars>( mod );
   if( ! tmod ) {
    // if it is not a "naked" FunctionModVars, it can still be a group of
    // identical *FunctionModVars* "dressed" into a GroupModification
    const auto gmod = std::dynamic_pointer_cast<GroupModification>( mod );
    if( gmod )  // if so, pick the first one and act on it
     tmod = std::static_pointer_cast<FunctionModVars>(
				        gmod->v_sub_Modifications.front() );
    }

   if( tmod ) {
    // if we have a *FunctionModVars*, we have to distinguish its exact type
    // and add/delete Variable accordingly; in all cases, however, the
    // Modification is processed and can be deleted
    to_delete = true;

    {
     const auto ttmod =
                      std::dynamic_pointer_cast<FunctionModVarsAddd>( tmod );
     if( ttmod ) {
      addd_vars = true;
      if( ! to_add ) {
       // the first time, check that the Modification data agrees with what
       // we expect
       if( ttmod->first() != NumVar )
	throw( std::logic_error( "wrong Variable names in FunctionModVars" )
	       );
       }
      
      to_add += ttmod->vars().size();
      continue;
      }
     }

    {
     const auto ttmod =
                      std::dynamic_pointer_cast<FunctionModVarsRngd>( tmod );
     if( ttmod ) {
      rmvd_vars = true;
      Range rng = ttmod->range();

      if( ( rng.first == 0 ) && ( rng.second >= NumVar ) ) {
       NumVar = 0;                      // deleting *all* Variable
       Lambda.clear();
       Lambda1.clear();
       LmbdBst.clear();
       Master->RmvVars( nullptr , 0 );  // remove from MP
       continue;                        // nothing else to do
       }
      
      if( rng.first >= NumVar ) {  // all the Variable are deleted already
       auto nr = rng.second - rng.first;
       if( nr > to_add )
	throw( std::logic_error( "removing non-existing Variable" ) );
       to_add -= nr;               // "virtually" remove them
       continue;                   // nothing else to do
       }
      if( rng.second >= NumVar ) {  // some of the Variable deleted already
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
       std::copy( Lambda.begin() + rng.second , Lambda.end() ,
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
     const auto ttmod =
                      std::dynamic_pointer_cast<FunctionModVarsSbst>( tmod );
     if( ttmod ) {
      rmvd_vars = true;

      if( ttmod->subset().empty() ) {  // deleting *all* Variable
       NumVar = 0;
       Lambda.clear();
       Lambda1.clear();
       LmbdBst.clear();
       Master->RmvVars( nullptr , 0 );  // remove from MP
       continue;                        // nothing else to do
       }

      if( ttmod->subset().front() >= NumVar ) {
       // all the Variable are deleted already
       if( ttmod->subset().size() > to_add )
	throw( std::logic_error( "removing non-existing Variable" ) );
       to_add -= ttmod->subset().size();  // "virtually" remove them
       continue;                          // nothing else to do
       }

      Subset tsbst;
      c_Subset * sbst = & tsbst;
      if( ttmod->subset().back() < NumVar )  // no Variable deleted already
       sbst = & ttmod->subset();             // delete them all
      else {                                 // construct the subset to delete
       auto sbstit = ttmod->subset().end();
       while( *(--sbstit) >= NumVar );
       tsbst = Subset( ttmod->subset().begin() , ++sbstit );
       auto nr = ttmod->subset().size() - tsbst.size();
       if( nr > to_add )
	throw( std::logic_error( "removing non-existing Variable" ) );
       to_add -= nr;               // "virtually" remove them
       }

      Compact( Lambda , *sbst );  // adjust Lambda
      NumVar -= sbst->size();
      Lambda.resize( NumVar );
      Lambda1.resize( NumVar );
      if( MaxSol > 1 )
       LmbdBst.resize( NumVar );
      Master->RmvVars( sbst->data() , sbst->size() );  // remove from MP
      continue;
      }
     }

    // if control reaches here, this is an unknown *FunctionModVars* (??)
    throw( std::logic_error( "unknown FunctionModVars" ) );

    }  // end( if( tmod ) )
   }  // end FunctionModVars
  }  // end( 4th loop, forward )

 // at this point, the set of Variable in the BundleSolver/Master Problem
 // coincides with the set of Variable in the C05Function(s), save for the
 // Variable to be added: in other words, the positions from 0 no NumVar - 1
 // in the linearizations corresponds to what BundleSolver expects

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
  c_Vec_p_Var * vars;           // the affected Variable

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  {
   // a C05FunctionModRngd, that at this point can only have which().empty()
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModRngd>( mod );
   if( tmod ) {
    if( ! tmod->which().empty() )
     throw( std::logic_error( "unexpected nonempty C05FunctionModRngd" ) );

    vars = & tmod->vars();
    range = tmod->range();
    }
   }

  {
   // a C05FunctionModSbst, that at this point can only have which().empty()
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModSbst>( mod );
   if( tmod ) {
    if( ! tmod->which().empty() )
     throw( std::logic_error( "unexpected nonempty C05FunctionModSbst" ) );

    vars = & tmod->vars();
    subset = & tmod->subset();
    }
   }

  {
   // a C05FunctionModLinRngd implies that a specific range in all the
   // linearizations must be changed (by adding something)
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModLinRngd>( mod );
   if( tmod ) {
    vars = & tmod->vars();
    range = tmod->range();
    }
   }

  {
   // a C05FunctionModLinSbst implies that a specific subset in all the
   // linearizations must be changed (by adding something)
   const auto tmod = std::dynamic_pointer_cast<C05FunctionModLinSbst>( mod );
   if( tmod ) {
    vars = & tmod->vars();
    subset = & tmod->subset();
    }
   }

  if( ( range.first >= range.second ) && ( ! subset ) )
   // it is neither of the above: this should not happen
   throw( std::logic_error( "unexpected Modification slipped in" ) );

  if( ! rmvd_vars ) {
   // Variable have never been removed, hence the names can be used directly
   if( subset ) {  // turn the subset into a range
    range.first = subset->front();
    range.second = subset->back() + 1;
    }
   }
  else {
   // Variable have been removed, and hence names need be actualised
   // this is done by directly checking vars() against the "active"
   // Variable of v_c05f[ 0 ], which is fairly taken as a representative
   // since all the C05Function have the same "active" Variable
   if( ! addd_vars ) {
    // ... but never added: names can have only decreased, but even more
    // importantly must have remained ordered, i.e., the first "active"
    // Variable in vars() is the first variable of the range, the last
    // "active" Variable vars() is the last variable of the range
    // note that we do not use subset and range here, as the range is
    // reconstructed from scratch using vars
    auto lit = vars->begin();
    for( ; lit != vars->end() ; ++lit ) {
     range.first = v_c05f[ 0 ]->is_active( *lit );
     if( range.first < v_c05f[ 0 ]->get_num_active_var() )
      break;
     }
    if( lit == vars->end() )  // no Variable in vars is still "active"
     continue;                // nothing else to do
    // since we know that here are some "active" Variable in vars(), this
    // second loop will necessarily end
     for( auto rit = vars->rbegin() ; ; ++rit ) {
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
    Subset newnames( vars->size() );
    auto lit = vars->begin();
    auto nni = newnames.begin();
    if( subset ) {
     auto sit = subset->begin();
     for( ; lit != vars->end() ; ++lit , ++sit ) {
      auto i = v_c05f[ 0 ]->is_active( *lit );
      if( ( i <= *sit ) && ( i < v_c05f[ 0 ]->get_num_active_var() ) )
       *(nni++) = i;
      }
     }
    else {
     for( ; lit != vars->end() ; ++lit , ++range.first ) {
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

 bool toadd = std::find_if( Addd.begin() , Addd.end() ,
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
 //
 // note that due to limitations in the MPSolver interface, changing a
 // linearization is identical to adding one, even if the change was limited
 // to a range/subset of the entries
 //
 // yet, handling of additions and changes differs depending on BPar7.
 // in particular, if ( BPar7 & 4 ), then additions to the global pool
 // imply additions to the bundle, otherwise just marking that the
 // linearization exists
 //
 // note that if ( BPar7 & 3 ) >= 2, BundleSolver will happily delete from
 // the global pool any linearization it deletes from the bundle; if also
 // ! ( BPar7 & 4 ), one could therefore think it appropriate to delete
 // from the global pool any linearization that is added. however, we do not
 // do that, just refraining to add them to the bundle. these will likely be
 // overwritten during the optimization, and if memory is a problem then the
 // global pools should just be sized accordingly

 if( toadd ||
     ( std::find_if( Chgd.begin() , Chgd.end() ,
		     []( Subset & Ck ) { return( ! Ck.empty() ); }
		     ) != Chgd.end() ) ) {
  // at least a component has had lnearizations added or changed

  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( Addd[ k ].empty() && Chgd[ k ].empty() )  // nothing to see here
    continue;                                    // move off

   if( Chgd[ k ].size() >= NrItems[ k ] ) {  // all items change
    reset[ k ] = AlphaC[ k ] = true;         // this is a reset
    continue;
    }

   // first, cleanup Chgd[ k ] from linearizations not in the bundle
   if( ! Chgd[ k ].empty() ) {
    Subset tmp;
    for( auto i : Chgd[ k ] )
     if( InvItemVcblr[ k ][ i ] >= vBPar2[ NrFi ] )
      tmp.push_back( i );

    if( ! tmp.empty() ) {
     if( tmp.size() >= Chgd[ k ].size() ) {
      Chgd[ k ].clear();
      if( Addd[ k ].empty() )  // nothing more to see here
       continue;               // move off
      }
     else
      set_difference_in_place( Chgd[ k ] , tmp );     
     }
    }

   // now, if ! ( BPar7 & 4 ), also cleanup Addd[ k ] from linearizations
   // not in the bundle, but in doing so mark them into InvItemVcblr
   // note that for linearizations in the bundle, being in Addd[ k ] is the
   // same as being in Chgd[ k ]: the linearization has changed, and the
   // master problem must be changed to reflect this
   if( ( ! ( BPar7 & 4 ) ) && ( ! Addd[ k ].empty() ) ) {
    Subset tmp;
    for( auto i : Addd[ k ] )
     if( InvItemVcblr[ k ][ i ] >= vBPar2[ NrFi ] ) {
      InvItemVcblr[ k ][ i ] = vBPar2[ NrFi ];
      tmp.push_back( i );
      }

    if( ! tmp.empty() ) {
     if( tmp.size() >= Addd[ k ].size() ) {
      Addd[ k ].clear();
      if( Chgd[ k ].empty() )  // nothing more to see here
       continue;               // move off
      }
     else
      set_difference_in_place( Addd[ k ] , tmp );     
     }
    }

   // compute the union between Addd[ k ] and Chgd[ k ] into Addd[ k ]
   set_union_in_place( Addd[ k ] , std::move( Chgd[ k ] ) );

   // finally, add the resulting stuff to the bundle
   for( auto i : Addd[ k ] )
    add_to_bundle( k , i );
   }
  }  // end( if( additions or changes ) )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if some component need be reset, reset the linearizations: since
 // reset[ k ] ==> AlphaC[ k ], later on also the constants will be reset

 if( std::find( reset.begin() , reset.end() , false ) == reset.end() )
  Master->ChgSubG( 0 , NumVar , NrFi + 1 );  // all components have been reset
 else
  if( std::find( reset.begin() , reset.end() , true ) != reset.end() )
   // some components have been reset
   for( Index k = 0 ; k < NrFi ; ++k )
    if( reset[ k ] )
     Master->ChgSubG( 0 , NumVar , k + 1 );

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are constants to change entirely, do it now in one blow
 // note: in case of a full reset, get_linearization_coefficients() is
 //       called twice, once in ChgSubG() (via GetGi()) and once in the
 //       loop below. this has the potential to be horribly inefficient,
 //       but the only clean way out is to do away with MPSolver entirely

 if( std::find( AlphaC.begin() , AlphaC.end() , false ) == AlphaC.end() ) {
  // all components have been reset
  std::vector< VarValue > Gi( NumVar );
  std::vector< VarValue > Alfa( Master->MaxName() );

  for( auto & cchk : Cchg )
   cchk.clear();

  for( Index i = 0 ; i < Master->MaxName() ; ++i )
   if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) {
    auto Ai = v_c05f[ ItemVcblr[ i ].first ]->get_linearization_constant(
						     ItemVcblr[ i ].second );
    if( std::isnan( Ai ) )  // linearization no longer valid
     throw( std::logic_error( "inconsistent ItemVcblr" ) );

    // compute the linearization error in Lambda
    v_c05f[ ItemVcblr[ i ].first ]->get_linearization_coefficients( Gi.data() ,
				 Range( 0 , NumVar ) , ItemVcblr[ i ].second );
    Alfa[ i ] = UpRifFi[ ItemVcblr[ i ].first ] - Ai -
                std::inner_product( Lambda.begin() , Lambda.end() , Gi.data() ,
				    double( 0 ) );
    }

  Master->ChgAlfa( Alfa.data() , NrFi + 1 );  
  }
 else
  if( std::find( AlphaC.begin() , AlphaC.end() , true ) != AlphaC.end() ) {
   // some components have been reset
   std::vector< VarValue > Gi( NumVar );

   for( Index k = 0 ; k < NrFi ; ++k )
    if( AlphaC[ k ] ) {  // all the constants of this component are reset
     Cchg[ k ].clear();  // no need to change them individually

     std::vector< VarValue > Alfa( Master->MaxName( k + 1 ) );

     for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
      if( InvItemVcblr[ k ][ i ] < vBPar2[ NrFi ] ) {
       auto Ai = v_c05f[ k ]->get_linearization_constant( i );
       if( std::isnan( Ai ) )  // linearization no longer valid
	throw( std::logic_error( "inconsistent InvItemVcblr" ) );

       // compute the linearization error in Lambda
       v_c05f[ k ]->get_linearization_coefficients( Gi.data() ,
						    Range( 0 , NumVar ) , i );
       Alfa[ InvItemVcblr[ k ][ i ] ] = UpRifFi[ k ] - Ai -
	       std::inner_product( Lambda.begin() , Lambda.end() , Gi.data() ,
				   double( 0 ) );
       }

     Master->ChgAlfa( Alfa.data() , k + 1 );
     }
   }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there subsets of Alphas to change, do it now

 if( std::find_if( Cchg.begin() , Cchg.end() ,
		   []( Subset & Ck ) { return( ! Ck.empty() ); }
		   ) != Cchg.end() ) {
  std::vector< VarValue > Gi( NumVar );

  for( Index k = 0 ; k < NrFi ; ++k )
   for( auto i : Cchg[ k ] )
    if( InvItemVcblr[ k ][ i ] < vBPar2[ NrFi ] ) {
     auto Ai = v_c05f[ k ]->get_linearization_constant( i );
     if( std::isnan( Ai ) )  // linearization no longer valid
      throw( std::logic_error( "inexistent linearization" ) );

     // compute the linearization error in Lambda
     v_c05f[ k ]->get_linearization_coefficients( Gi.data() ,
						  Range( 0 , NumVar ) , i );
     Ai = UpRifFi[ k ] - Ai - std::inner_product( Lambda.begin() ,
						  Lambda.end() ,
						  Gi.begin() , double( 0 ) );
     Master->ChgAlfa( InvItemVcblr[ k ][ i ] , Ai );
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

 #if CHECK_DS & 1
  CheckBundle();  // save possibly some checks
  //!! PrintBundle();
 #endif
 
 }  // end( BundleSolver::process_outstanding_Modification ) - - - - - - - - -

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

void BundleSolver::CheckBundle( void )
{
 // check vBPar2
 for( Index k = 0 ; k < NrFi ; ++k )
  if( v_c05f[ k ]->get_int_par( C05Function::intGPMaxSz ) != vBPar2[ k ] )
   std::cerr << "size of global pool " << k << " does not match" << std::endl;

 // check ItemVcblr against InvItemVcblr and Master
 Subset tmp( NrFi , 0 );
 for( Index i = 0 ; i < Master->MaxName() ; ++i )
  if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) {
   ++tmp[ ItemVcblr[ i ].first ];
   if( Master->WComponent( i ) != ItemVcblr[ i ].first + 1 ) {
    std::cerr << "position " << i << " in the bundle should be of component "
	      << ItemVcblr[ i ].first << " but Master says ";
    if( Master->WComponent( i ) == Inf<Index>() )
     std::cerr << "empty" << std::endl;
    else
     std::cerr << Master->WComponent( i ) - 1 << std::endl;
    }

   if( InvItemVcblr[ ItemVcblr[ i ].first ][ ItemVcblr[ i ].second ] != i )
    std::cerr << "position " << i << " in the bundle shoud be linearization "
	      << ItemVcblr[ i ].second << " of component "
	      << ItemVcblr[ i ].first << " but InvItemVcblr disagrees"
	      << std::endl;
   }
  else
   if( Master->WComponent( i ) < Inf<Index>() )
    std::cerr << "position " << i
	      << " in the bundle should be empty but Master says "
	      << Master->WComponent( i ) << std::endl;

 // check NrItems
 if( int diff = ( NrItems[ NrFi ] - std::accumulate( NrItems.begin() ,
						     NrItems.begin() + NrFi ,
						     Index( 0 ) ) ) != 0 )
  std::cerr << " NrItems[ NrFi ] = " << NrItems[ NrFi ] << " but the sum is "
	    <<  NrItems[ NrFi ] + diff  << std::endl;

 for( Index k = 0 ; k < NrFi ; ++k )
  if( tmp[ k ] != NrItems[ k ] )
   std::cerr << "counted " << tmp[ k ]
	     << " items in the bundle for component " << k
	     << " but NrItems says " << NrItems[ k ] << std::endl;

 // check InvItemVcblr against ItemVcblr and C05Function
 for( Index k = 0 ; k < NrFi ; ++k )
  for( Index i = 0 ; i < vBPar2[ k ] ; ++i ) {
   if( ( InvItemVcblr[ k ][ i ] < Inf<Index>() ) &&
       ( ! v_c05f[ k ]->is_linearization_there( i ) ) )
    std::cerr << "linearization " << i << " in pool " << k
	      << " does not exist" << std::endl;

   if( ( InvItemVcblr[ k ][ i ] == Inf<Index>() ) &&
       v_c05f[ k ]->is_linearization_there( i ) )
    std::cerr << "linearization " << i << " in pool " << k
	      << " unaccounted for" << std::endl;

   if( ( InvItemVcblr[ k ][ i ] < vBPar2[ NrFi ] ) &&
       ( ( ItemVcblr[ InvItemVcblr[ k ][ i ] ].first != k ) ||
	 ( ItemVcblr[ InvItemVcblr[ k ][ i ] ].second != i ) ) )
    std::cerr << "linearization " << i << " in pool " << k
	      << " should be in bundle in position "
	      << InvItemVcblr[ k ][ i ] << " but ItemVcblr disagrees"
	      << std::endl;

   if( ( InvItemVcblr[ k ][ i ] < Inf<Index>() ) &&
       ( i >= MaxItem[ k ] ) )
    std::cerr << "free item in position " << i << " of pool " << k
	      << " but MaxItem says " << MaxItem[ k ] << std::endl;

   if( ( InvItemVcblr[ k ][ i ] >=
	 ( ( BPar7 & 3 ) ? vBPar2[ NrFi ] : Inf<Index>() ) ) &&
       ( i < FrFItem[ k ] ) )
     std::cerr << "free item in position " << i << " of pool " << k
	       << " but FrFItem says " << FrFItem[ k ] << std::endl;

   }  // end( for( i ) )

 }  // end( BundleSolver::CheckBundle )
 
/*--------------------------------------------------------------------------*/

void BundleSolver::PrintBundle( void )
{
 if( ! f_log )
  return;

 *f_log << std::endl << "Lambda = [ ";
 for( Index h = 0 ; h < NumVar - 1 ; ++h )
  *f_log << Lambda[ h ] << ", ";
 *f_log << Lambda.back() << " ]";

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
 throw( std::logic_error( "this method cannot be called" ) );
 // the implementation below is wrong, it gives the linearization constant
 // rather than the linearization error in the point where the component is
 // computed. this would require translating the value w.r.t. Lambda1,
 // which is possible but ugly, and anyway not needed

 if( Name == bslv->vBPar2[ bslv->NrFi ] ) // get the 0th-component subgradient
  return( bslv->f_lf->get_linearization_constant() );
 else
  return( bslv->v_c05f[ bslv->ItemVcblr[ Name ].first ]->
	       get_linearization_constant( bslv->ItemVcblr[ Name ].second ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetGiName( cIndex Name )
{
 throw( std::logic_error( "this method cannot be called" ) );
 }

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
