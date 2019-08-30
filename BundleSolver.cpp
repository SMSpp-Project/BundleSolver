/*--------------------------------------------------------------------------*/
/*------------------------ File BundleSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BunldeSolver class.
 *
 * \version 0.01
 *
 * \date 28 - 08 - 2019
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
#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( LOG_BND )
 #define BLOG( l , x ) if( LogVerb > l ) *f_log << x

 #define BLOG2( l , c , x ) if( ( LogVerb > l ) && c ) *f_log << x

 #define BLOGb( l , x ) if( LogVerb & l ) *f_log << x

 #define BLOG2b( l , c , x ) if( ( LogVerb & l ) && c ) *f_log << x
#else
 #define BLOG( l , x )

 #define BLOG2( l , c , x )

 #define BLOGb( l , x )

 #define BLOG2b( l , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BundleSolver to the Solver factory
SMSpp_insert_in_factory_cpp_0( BundleSolver );

/*--------------------------------------------------------------------------*/
// define and initialize here the vector of int parameters names
const std::vector< std::string > BundleSolver::int_pars_str =
             { "intBPar1" , "intBPar6" ,
               "intMnSSC" , "intMnNSC" ,
               "inttSPar1" , "intMaxNrEvls" , "intKpBstL" , "intMPName" ,
			   "intQPmp1" ,  "intQPmp2",
			   "OSImp1" , "OSImp2" , "OSImp3" , "OSImp4" , "OSImp5"  };

// define and initialize here the vector of double parameters names
const std::vector< std::string > BundleSolver::dbl_pars_str =
		     { "dbltStar"  , "dblEInit" ,
		       "dblBPar3" ,  "dblBPar4"  ,
		       "dblBPar5"  , "dblm1" , "dblm3" ,
			   "dblmxIncr" ,  "dblmnIncr" ,  "dblmxDecr" ,
			   "dblmnDecr" ,  "dbltMaior" ,  "dbltMinor" ,
			   "dbltInit" ,  "dbltSPar2" , "dblQPmp1" };

// define and initialize here the map for int parameters names
const std::map< std::string , BundleSolver::idx_type > BundleSolver::int_pars_map =
                   { { "intBPar1"  , BundleSolver::intBPar1  } ,
		     { "intBPar6" , BundleSolver::intBPar6 } ,
		     { "intMnSSC" , BundleSolver::intMnSSC } ,
		     { "intMnNSC" , BundleSolver::intMnNSC } ,
		     { "inttSPar1" , BundleSolver::inttSPar1 } ,
			 { "intMaxNrEvls" , BundleSolver::intMaxNrEvls } ,
		     { "intKpBstL" , BundleSolver::intKpBstL } ,
			 { "intMPName" , BundleSolver::intMPName } ,
			 { "intQPmp2" , BundleSolver::intQPmp1 } ,
			 { "intQPmp3" , BundleSolver::intQPmp2 } ,
			 { "intOSImp1" , BundleSolver::intOSImp1 } ,
			 { "intOSImp2" , BundleSolver::intOSImp2 } ,
			 { "intOSImp3" , BundleSolver::intOSImp3 } ,
			 { "intOSImp4" , BundleSolver::intOSImp4 } ,
			 { "intOSImp5" , BundleSolver::intOSImp5 } ,
              };

// define and initialize here the map for double parameters names
const std::map< std::string , BundleSolver::idx_type > BundleSolver::dbl_pars_map =
                   { { "dbltStar" , BundleSolver::dbltStar } ,
		     { "dblEInit" , BundleSolver::dblEInit } ,
			 { "dblBPar3" , BundleSolver::dblBPar3 } ,
			 { "dblBPar4" , BundleSolver::dblBPar4 } ,
			 { "dblBPar5" , BundleSolver::dblBPar5 } ,
			 { "dblm1" , BundleSolver::dblm1 } ,
			 { "dblm3" , BundleSolver::dblm3 } ,
			 { "dblmxIncr" , BundleSolver::dblmxIncr } ,
			 { "dblmnIncr" , BundleSolver::dblmnIncr } ,
			 { "dblmxDecr" , BundleSolver::dblmxDecr } ,
			 { "dblmnDecr" , BundleSolver::dblmnDecr } ,
			 { "dbltMaior" , BundleSolver::dbltMaior } ,
			 { "dbltMinor" , BundleSolver::dbltMinor } ,
			 { "dbltInit" , BundleSolver::dbltInit } ,
			 { "dbltSPar2" , BundleSolver::dbltSPar2 } ,
			 { "dblQPmp1" , BundleSolver::dblQPmp1 } };

// define and initialize here the default int parameters
const std::vector<int> BundleSolver::dflt_int_par =
        {    10 ,  // intBPar1
			  0 ,  // intBPar6
			  0 ,  // intMnSSC
			  3 ,  // intMnNSC
			 12 ,  // inttSPar1
			  2 ,  // intMaxNrEvls
			  0 ,  // intKpBstL
			  1 ,  // intMPName
			  0 ,  // intQPmp1
			  0 ,  // intQPmp2
			  4 ,  // intOSImp1
			  0 ,  // intOSImp2
			  1 ,  // intOSImp3
			  3 ,  // intOSImp4
			  1    // intOSImp5
		};

// define and initialize here the default double parameters
const std::vector<double> BundleSolver::dflt_dbl_par =
           { 1e2 ,    // dbltStar
			 1e-2 ,   // dblEInit
			 - 1 ,    // dblBPar3
			 - 1 ,    // dblBPar4
			 30 ,     // dblBPar5
			 0.1 ,    // dblm1
			 3 ,      // dblm3
			 10 ,     // dblmxIncr
			 1.5 ,    // dblmnIncr
			 0.1 ,    // dblmxDecr
			 0.66 ,   // dblmmDecr
			 1e6 ,    // dbltMaior
			 1e-6,    // dbltMinor
			  1 ,     // dbltInit
			0.001 ,   // dbltSPar2
			0.1       // dblQPmp1
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

static cIndex InINF = SMSpp_di_unipi_it::Inf<Index>();

/*--------------------------------------------------------------------------*/

int BundleSolver::compute( bool changedvars )
{
 // basic sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! Master )
  throw( std::logic_error( "Master is not set yet" ) );

 if( v_c05f.empty() )
  throw( std::logic_error( "C05Function is not set yet" ) );


 bool MPisQuad = false;

 auto OsiMP = dynamic_cast<OSIMPSolver*>( Master );
 auto QppMP = dynamic_cast<OSIMPSolver*>( Master );

 if( ( OsiMP &&  ( stblztn == OSIMPSolver::quadratic ) ) || QppMP )
  MPisQuad = true;

 // initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Result = kOK;
 SCalls++;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 do {
  // construct the direction d- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormD();

  if( !CheckAlfa() ) {
   t = std::min( t * mxIncr , tMaior );
   BLOG( 1 , " ~ noise reduction: t increased to " << t << std::endl );
   tHasChgd = true;
   if( t >= tMaior )
	Result = kError;
   else
    continue;
   }

  if( Result )  // problems in the Master Problem solver
   break;

  // update out-of-base counters- - - - - - - - - - - - - - - - - - - - - - -

  UpdtCntrs();

  // some log - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log1();

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( IsOptimal() )
   break;

  // Hard Long-Term t-strategy for quadratic stabilization- - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the hard long-term t-strategy requires t to increase if the step is too
  // small, and therefore has to be checked before the others

  if( MPisQuad && ( tStar > 0 ) && ( ( tSPar1 & tSP1Msk ) == kHLTTS )
	   && ( UpFiLmb[NrFi] < Inf<double>() ) ) {

   double AFL = std::abs( UpFiLmb[NrFi] );
   if( AFL < 1 )
    AFL = 1;

   if( vStar[ NrFi ] <= tSPar2 * EpsU * AFL ) {
    BLOG( 1 , "small v => increase t" << std::endl << "           " );

    // collect two numbers vc and vl such that v( tNew ) >= vc + tNew * vl
    // we require that v( tNew ) >= vc + tNew * vl = tSPar2 * EpsU * AFL
    // ==> tNew = ( tSPar2 * EpsU * AFL - vc ) / vl

    double vl , vc;
    Master->SensitAnals( vl , vc );

    double tt;
    if( - vl < Eps<HpNum>() )  // v( t ) is [almost] constant ==> D*_t [~]= 0
     tt = tStar;       // ==> the CP model is ~bounded
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
  // This is done *after* the call to Master->SensitAnals() in the Hard Long-Term t-strategy and
  // to FormLambda1(), because elimination of items from the bundle may make
  // the current solution of the master problem invalid, and therefore all
  // solution information may be lost. In theory this should not happen, since
  // only items "out of base" are eliminated, and therefore the solution
  // remains optimal; however, not all MPSolvers may behave in this respect.

  SimpleBStrat();

  // calculate Fi( Lambda1 )- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index wFi = 0;
  CurrNrEvls.assign( NrFi , Index(0) );

  bool MPchgs = false;  // true if no cycling will occur
  for( ; ; ) {   // ... possibly more than once due to precision issues

   MPchgs = FiAndGi( wFi ); // the component wFi needs more time to be solved
   CurrNrEvls[ wFi ]++;

   if( MPchgs ) // if something changes
	break;                                          // all done

   // find next component   - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( !FindNext( wFi ) )
	break;

   } // end Fi and Gi computation - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  if( !MPchgs ) { // noise reduction
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

  // check the Lower Bound- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb[NrFi] < Inf<double>() ) {  // .. but only if Fi( Lambda ) is defined
   if( TrueLB )
    if( UpFiBest[NrFi] - RelAcc * std::abs( UpFiBest[NrFi] ) <= LowerBound )
     break;

   if( UpFiBest[NrFi] <= LowerBound *
                  ( 1 - ( LowerBound > 0 ? RelAcc : - RelAcc ) ) ) {
    Result = kUnbounded;
    break;
    }
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
  SSDone = ( UpFiLmb1[ NrFi ] < UpTrgt )? true : false;

  // compute the heuristic t- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  HpNum tt;
  if( ( SSDone && ( ! ( tSPar1 & 1 ) ) ) ||
       ( ( ! SSDone ) && ( tSPar1 & 2 ) ) )
   tt = Heuristic1();
  else
   tt = Heuristic2();

  if( SSDone ) {  // SS - - - - - - - - - - - - - - - - - - - - - - - - - - -
   BLOG( 1 , std::endl << " SS[" << CSSCntr << "]: DFi (" << DeltaFi
 	          << ") >= m1 * Dv (" << std::abs( m1 ) * Deltav << ") ~ Ht = "
 	          << tt );

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
   ParSS++;
   CSSCntr++;
   CNSCntr = 0;
   }
  else {        // NS - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   BLOG( 1 , std::endl << " NS[" << CNSCntr << "]: DFi (" << DeltaFi
 	          << ") < m1 * Dv (" << std::abs( m1 ) * Deltav
 	          << ") ~ Ht = " << tt );

   tt = std::max( std::max( tMinor , t * mxDecr ) ,
 		  std::min( t * mnDecr , tt ) );

   if( CNSCntr < MnNSC )  // decreasing t is inhibited
    tt = t;
   else
    if( Alfa1[ NrFi ] <= m3 * Sigma ) {
     BLOG( 1 , " ~ small Alfa1" );
     tt = t;
     }
    else
     switch( tSPar1 & tSP1Msk ) {
      case( kSLTTS ):
      case( kHLTTS ):
       if( vStar[ NrFi ] <= tSPar2 * EpsU * std::max( std::abs( UpFiLmb[NrFi] ) , HpNum( 1 ) ) ) {
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

  if( tSPar1 & kEGTTS )  // endgame t-strategy: note the "/ 10"!!
   if( DSTS < RelAcc * std::max( std::abs( UpFiLmb[NrFi] ) , HpNum( 1 ) ) / 10 ) {
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

  } while( ( ! MaxIter ) || ( ParIter < MaxIter ) );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if( MaxIter && ( ParIter >= MaxIter ) && ( ! Result ) )
 Result = kStopIter;

return( Result );

}  // end( BundleSolver::compute() ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::set_Block( Block * block )
{
 if( f_Block ) {  // changing from a previous oracle - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RemoveItems();  // clear the bundle
  guts_of_destructor();   // deallocate memory
  }

 Solver::set_Block( block );  // attach to the new Block

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
  linear_function = nullptr;

  }
 else {

  // the objective function of each block must be a LinearFunction - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( f_Block->get_objective() )
   linear_function =  nullptr;
  else {
   auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
   if( obj == nullptr )
    throw( std::logic_error( "the objective is not a real function" ) );

   linear_function = dynamic_cast<LinearFunction *>( (obj)->get_function() );
   if( linear_function == nullptr )
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
   break;
  if( un_any_thing_1( ColVariable , el , NumVar += var.size() ) )
   break;
  if( un_any_thing_K( ColVariable , el , NumVar += var.size() ) )
   break;
  throw( std::logic_error( "some static Variable is not a ColVariable" ) );
  }

 // construct the vocabulary for Variable and sort it  - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LamVcblr.resize( NumVar );
 int count = 0;
 for( auto & el : v_s_Variable )
  un_any_static( el , [ & ]( ColVariable & static_var ){ LamVcblr[count++]= &static_var; } ,
		  un_any_type<ColVariable>() );

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

 // set Bpar2 as the sum of the global pool of the components  - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 BPar2 = 0;
 for( auto fun : v_c05f )
  BPar2 += fun->get_int_par( C05Function::intGPMaxSz );

 // read information about the function  - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NrFi = v_c05f.size();

 NrEasy = 0;
 if( NrFi > 1 ) {

  auto sb = f_Block->get_nested_Blocks();
  std::vector<Index> BNC( NrFi );

  bool HasEasy = false;
  IsEasy.resize( NrFi );
  for( Index k = 0 ; k < NrFi ; ++k ) {

   // ?? quando i solver distruggere dopo SetDim o
   // con il distruttore della classe??

   auto LagB = dynamic_cast<LagBFunction *>( v_c05f[ k ] );
   if( LagB  ) {

    MILP_s[ k ] = new MILPSolver();
    MILP_s[ k ]->set_Block( LagB->get_inner_block() );

    BNC[ k ] = MILP_s[ k ]->get_numcols();

    if( BNC[ k ] ) {
     IsEasy[ k ] = HasEasy = true;
     NrEasy++;
     }
    else
     IsEasy[ k ] = false;

    }
   }

   if( ! HasEasy )
    IsEasy.clear();

   }

 // allocate memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 t = tInit;
 Prevt = Inf<double>();

 Lambda.resize( NumVar );    // the default starting point
 Lambda1.resize( NumVar );   // the tentative point

 if( KpBstL )  // best point found so far
  LmbdBst.resize( NumVar);

 OOBase.resize( BPar2 , Inf<SIndex>() );  // counter for eliminating outdated
                                          // items: Inf<SIndex>() means empty

 FreList.resize( BPar2 );       // list of free bundle slots
 whisZ.resize( NrFi );          // for each component, the name of its "Z" if it is
                                // in the bunlde

 FiStatus.resize( NrFi , kUnEval );
 LowerBound = -Inf<double>(); // lower bounds
 TrueLB = false;

 UpFiBest.resize( NrFi + 1 , Inf<OFValue>() ); // best, ...
 UpRifFi.resize( NrFi + 1 , 0 ); // and reference Fi() values
 UpFiLmb1.resize( NrFi + 1 );    // upper and lower function value
 LwFiLmb1.resize( NrFi + 1 );    // ... at the tentative point

 UpFiLmb.resize( NrFi + 1 , Inf<OFValue>() );  // upper and lower function value
 LwFiLmb.resize( NrFi + 1 , -Inf<OFValue>() ); // ... at the current point

 whisG1.resize( NrFi , Inf<Index>() );  // no representative yet

 ScPr1.resize( NrFi + 1 , 0 );
 Alfa1.resize( NrFi + 1 , 0 );
 DeltaAlfa.resize( NrFi );

 FreDim = 0;
 Result = kError;
 SSDone = false;

 ReSetAlg( RstCrr | RstSbg | RstCnt );  // Fi( Lambda ) is reset inside

 // warning: the following things can only be done *after* that
 // Oracle->SetMaxName() has been invoked, because they use methods of the
 // oracle which depends on knowledge of the MaxName to work properly
 // read b0- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // here one could initialize b0, if that was found to be of any use
 // b0 = Oracle->GetVal( BPar2 );

 // initialize the MP Solver, if any - - - - - - - - - - - - - - - - - - - -

 if( Master )        // a MP solver is set ??? dove metterlo ???
  Master->SetDim();  // clear all its internal state

 if( MPName ) {
  ofstream qp_ofs ("param.qp");
  if( qp_ofs.is_open() ) {
   qp_ofs << CtOff << endl;
   qp_ofs << MxAdd << endl;
   qp_ofs << MxRmv << endl;
   qp_ofs.close();
   }
  else
   throw( std::logic_error( "errors in loading parameters" ) );
  ifstream qp_ifs ("param.qp");
  if( !qp_ifs.is_open() )
   throw( std::logic_error( "errors in loading parameters" ) );
  Master = new QPPenaltyMP( &qp_ifs );
  }
 else {
  Master = new OSIMPSolver( );
  OSIMPSolver *osi_mps = dynamic_cast<OSIMPSolver*>( Master );
  if( osi_type ) {
   OsiCpxSolverInterface *osicpx = new OsiCpxSolverInterface();
   osi_mps->SetOsi( osicpx );
   CPXENVptr env = osicpx->getEnvironmentPtr ();
   CPXsetintparam( env , CPX_PARAM_THREADS , threads );
   }
  else {
   OsiClpSolverInterface *osiclp = new OsiClpSolverInterface();
   osi_mps->SetOsi( osiclp );
   }
  osi_mps->SetStabType( OSIMPSolver::StabFun( stblztn ) );
  osi_mps->SetAlgo( OSIMPSolver::OsiAlg( algo ) , OSIMPSolver::OsiRed( reduction ) );
  }

 InitMP();

 }  // end( BundleSolver::set_Block )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( const idx_type par , const int value ) {

 switch( par ) {
  case( intMaxIter ):
   MaxIter = value;
   break;
  case( intMaxSol ):
   MaxSol = value;
   break;
  case( intLogVerb ):
   LogVerb = value;
   break;
  case( intBPar1 ):
   BPar1 = value;
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
  case( intKpBstL ):
   KpBstL = value;
   break;
  case( intMPName ):
   MPName = bool(value);
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
  case( intOSImp4 ):
   stblztn = value;
   break;
  case( intOSImp5 ):
   osi_type = bool(value);
   break;
  default:
   CDASolver::set_par( par , value );
  }

 } // end (BundleSolver::set_par( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( const idx_type par , const double value ) {

 switch( par ) {
  case( dblMaxTime ):
   MaxTime = value;
   break;
  case( dblRelAcc ):
   RelAcc = value;
   break;
  case( dblAbsAcc ):
   AbsAcc = value;
   break;
  case( dblRAccSol ):
   RAccSol = value;
   break;
  case( dblAAccSol ):
   AAccSol = value;
   break;
  case( dblFAccSol  ):
   FAccSol = value;
   break;
  case( dbltStar ):
   tStar = value;
   break;
  case( dblEInit ):
   EInit = value;
   break;
  case( dblBPar3 ):
   BPar3 = value;
   break;
  case( dblBPar4 ):
   BPar4 = value;
   break;
  case( dblBPar5 ):
   BPar5 = value;
   break;
  case( dblm1 ):
   m1 = value;
   break;
  case( dblm3 ):
   m3 = value;
   break;
  case( dblmxIncr ):
   mxIncr = value;
   break;
  case( dblmnIncr ):
   mnIncr = value;
   break;
  case( dblmxDecr ):
   mxDecr = value;
   break;
  case( dblmnDecr ):
   mnDecr = value;
   break;
  case( dbltMaior ):
   tMaior = value;
   break;
  case( dbltMinor ):
   tMinor = value;
   break;
  case( dbltInit ):
   tInit = value;
   break;
  case( dbltSPar2 ):
   tSPar2 = value;
   break;
  case( dblQPmp1 ):
   CtOff = value;
   break;
  default:
   CDASolver::set_par( par , value );
  }

 } // end (BundleSolver::set_par( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

// end( BundleSolver::compute ) - - - - - - - - - - - - - - - - - - - - - -

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

int BundleSolver::get_dflt_int_par( const idx_type par ) const
{
 if( ( par >= intBPar1 ) && ( par < intLastBndSlvPar ) )
  return( dflt_int_par[ par - intBPar1 ] );
 else
  return( CDASolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double BundleSolver::get_dflt_dbl_par( const idx_type par ) const
{
 if( ( par >= dbltStar ) && ( par < dblLastBndSlvPar ) )
  return( dflt_dbl_par[ par - dbltStar ] );
 else
  return( CDASolver::get_dflt_dbl_par( par ) );
 }

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
  case( intKpBstL ):
   return( KpBstL );
   break;
  case( intMPName ):
   return( MPName );
   break;
  case( intQPmp1 ):
   return( MxAdd );
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
  case( intOSImp4 ):
   return( threads );
   break;
  case( intOSImp5 ):
   return( stblztn );
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
  case( dblFAccSol  ):
   return( FAccSol );
   break;
  case( dbltStar ):
   return( tStar );
   break;
  case( dblEInit ):
   return( EInit );
   break;
  case( dblBPar3 ):
   return( BPar3 );
   break;
  case( dblBPar4 ):
   return( dblBPar4 );
   break;
  case( dblBPar5 ):
   return( BPar5 );
   break;
  case( dblm1 ):
   return( m1 );
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
  case( dblQPmp1 ):
   return( CtOff );
   break;
  default:
   return( CDASolver::get_dflt_dbl_par( par ) );
  }

 } // end( BundleSolver::get_dbl_par ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type BundleSolver::int_par_str2idx(
		const std::string & name ) const
{
 // these may be many enough as to warrant using a map
 const auto it = int_pars_map.find( name );
 if( it != int_pars_map.end() )
  return( it->second );
 else
 return( CDASolver::int_par_str2idx( name ) );

 } // end( BundleSolver::int_par_str2idx ) - - - - - - - - - - - - - - - - - -

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

ThinComputeInterface::idx_type BundleSolver::dbl_par_str2idx( const std::string & name ) const
{
 // these may be many enough as to warrant using a map
 const auto it = dbl_pars_map.find( name );
 if( it != dbl_pars_map.end() )
  return( it->second );
 else
  return( CDASolver::dbl_par_str2idx( name ) );
 } // end( BundleSolver::dbl_par_str2idx ) - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

const std::string & BundleSolver::int_par_idx2str(
	const ThinComputeInterface::idx_type idx ) const
{
 if( ( idx >= intBPar1 ) && ( idx < intLastBndSlvPar ) )
  return( int_pars_str[ idx - intBPar1 ] );
 else
  return( CDASolver::int_par_idx2str( idx ) );
 } // end( BundleSolver::int_par_idx2str ) - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

const std::string & BundleSolver::dbl_par_idx2str( const idx_type idx )
 const
{
 if( ( idx >= dbltStar ) && ( idx < dblLastBndSlvPar ) )
  return( dbl_pars_str[ idx - intBPar1 ] );
 else
  return( CDASolver::dbl_par_idx2str( idx ) );
 } // end( BundleSolver::dbl_par_idx2str ) - - - - - - - - - - - - - - - - - -

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

 if( LBHasChgd && ( UpFiLmb[NrFi] < Inf<double>() ) ) {
  if( LowerBound > - Inf<double>() )
   Master->SetLowerBound( LowerBound - UpFiLmb[NrFi] );
  else
   Master->SetLowerBound( - Inf<double>() );

  for( Index k = 0 ; k < NrFi ; k++ ) {
   if( IsEasy.size() && IsEasy[ k ] )  // skip easy components
    continue;

  if( LowerBound > - Inf<double>() )
    Master->SetLowerBound( LowerBound - UpFiLmb[ k ] , k + 1 );
   else
    Master->SetLowerBound( - Inf<double>() , k + 1 );
   }

  LBHasChgd = false;
  }

 // set termination criterion - - - - - - - - - - - - - - - - - - - - - - - -

 if( UpFiLmb[NrFi] < Inf<double>() )
  Master->SetPar( MPSolver::kZero ,
		  RelAcc * std::max( std::abs( UpFiLmb[NrFi] ) , double( 1 ) )
		         / std::max( tStar / t , HpNum( 1 ) ) );

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
 vStar[ NrFi ] = - Master->ReadFiBLambda();            // read v*

 if( IsEasy.size() ) {                         // there are easy components
  for( Index k = 0 ; k < NrFi ; k++ )          // read the *exact* Fi-value
   if( IsEasy[ k ] )                           // for all them
    UpFiLmb1[ k ] = Master->ReadFiBLambda( k );
   else
    vStar[ k ] = - Master->ReadFiBLambda( k );

  if( UpFiLmb[NrFi] < Inf<double>() )
   for( Index k = 0 ; k < NrFi ; k++ )
    if( IsEasy[ k ] )
     vStar[ NrFi ] += UpRifFi[ k ];
  }
 else
  for( Index k = 0 ; k < NrFi ; k++ )
   if( !IsEasy[ k ] )
    vStar[ k ] = - Master->ReadFiBLambda( k );

 DSTS = Master->ReadDStart( tStar );           // D_{t*,\beta,x}
 Deltav = vStar[ NrFi ];
 if( m1 < 0 )                                  // use - z( P_{t,\beta,x} )
  Deltav -= Master->ReadDt( t );

 // Sigma* + D*_{t*}( -z* ) is the "maximum expected increase" used in
 // the stopping criterion, EpsU is that relative to Fi( Lambda )

 if( UpFiLmb[NrFi] < Inf<double>() && tStar > 0 )
  EpsU = ( DSTS + Sigma ) / std::max( std::abs( UpFiLmb[NrFi] ) , double( 1 ) );
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

 if( linear_function ) { // add the linear part to the "full function"
  linear_function->compute( true );
  UpFiLmb1[ NrFi ] = linear_function->get_upper_estimate();
  LwFiLmb1[ NrFi ] = linear_function->get_lower_estimate();
  }
 else
  UpFiLmb1[ NrFi ] = LwFiLmb1[ NrFi ] = 0;

 for( Index k = 0 ; k < NrFi ; k++ ) {
  if( IsEasy.size() && IsEasy[ k ] )  // if k is an easy component
   UpFiLmb1[ k ] =  LwFiLmb1[ k ] = Master->ReadFiBLambda( k );
  else {

   // initialize upper and lower bound for each component  - - - - - - - - - -

   if( v_c05f[ k ]->get_Lipschitz_constant() < Inf<FunctionValue>()
	   && UpFiLmb[ k ] < Inf<FunctionValue>() )
    UpFiLmb1[ k ] = UpRifFi[ k ] + v_c05f[ k ]->get_Lipschitz_constant() * NrmD;
   else
	UpFiLmb1[ k ] = Inf<FunctionValue>();

   if( LwFiLmb[ k ] > -Inf<FunctionValue>() )
    LwFiLmb1[ k ] = UpRifFi[ k ] + vStar[ k ];
   else
    LwFiLmb1[ k ] = -Inf<FunctionValue>();
   }

  // sum over the components, the zero-component is already there
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1[ NrFi ] < Inf<FunctionValue>() ) {
   if( UpFiLmb1[ k ] < Inf<FunctionValue>() )
    UpFiLmb1[ NrFi ] += UpFiLmb1[ k ];
   else
	UpFiLmb1[ NrFi ] = Inf<FunctionValue>();
   }

  if( LwFiLmb1[ NrFi ] > -Inf<FunctionValue>() ) {
   if( LwFiLmb1[ k ] < Inf<FunctionValue>() )
    LwFiLmb1[ NrFi ] += LwFiLmb1[ k ];
   else
    LwFiLmb1[ NrFi ] = -Inf<FunctionValue>();
   }
  }

 // update the upper and lower targets - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( UpFiLmb[ NrFi ] < Inf<FunctionValue>() )
  UpTrgt = UpRifFi[ NrFi ] + (1.0 - m2) * vStar[ NrFi ];
 else
  UpTrgt = Inf<FunctionValue>();

 if( LwFiLmb[ NrFi ] > -Inf<FunctionValue>() )
  if( m1 > 0 )
   LwTrgt = UpRifFi[ NrFi ] + vStar[ NrFi ] + m1 * DeltaStar;
  else
   LwTrgt = UpRifFi[ NrFi ] + ( 1.0 - m1 ) * vStar[ NrFi ];
 else
  LwTrgt = -Inf<FunctionValue>();

 }  // end( BundleSolver::FormLambda1 )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::FiAndGi( Index wFi )
{

 double UpCutOff, LwCutOff, LwFiK, EpsCurr;

 if( IsEasy.size() && IsEasy[ wFi ] )
  return( false );

 LwFiK =  UpRifFi[ wFi ] + vStar[ wFi ];

 if( UpFiLmb[ wFi ] < Inf<OFValue>() )
  if( UpTrgt < Inf<OFValue>() && UpFiLmb1[ NrFi ] < Inf<OFValue>() )
   UpCutOff = std::max( UpTrgt - ( UpFiLmb1[ NrFi ] - UpFiLmb1[ wFi ] ) ,
   		              LwFiK - m2 * BetaK( wFi ) * vStar[ NrFi ] );
  else
   UpCutOff = LwFiK - m2 * BetaK( wFi ) * vStar[ NrFi ];
 else
  UpCutOff = Inf<OFValue>();

 if( LwFiLmb[ wFi ] > -Inf<OFValue>() )
  if( LwTrgt > -Inf<OFValue>() && LwFiLmb1[ NrFi ] > -Inf<OFValue>() )
   LwCutOff = std::max( LwTrgt - ( LwFiLmb1[ NrFi ] - LwFiLmb1[ wFi ] ) ,
    		              LwFiK + m1 * BetaK( wFi ) * DeltaStar );
  else
   LwCutOff = LwFiK + m1 * BetaK( wFi ) * DeltaStar;
 else
  LwCutOff = -Inf<OFValue>();

 if( LwCutOff > -Inf<OFValue>() && UpCutOff < Inf<OFValue>() )
  EpsCurr = ( UpCutOff - LwCutOff ) / std::max( 1.0 , std::abs(UpRifFi[ wFi ] ) );
 else
  EpsCurr = EInit;

 // assign the cutoff values to the c05Function - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_c05f[ wFi ]->set_par( dblUpCutOff , UpCutOff );
 v_c05f[ wFi ]->set_par( dblLwCutOff , LwCutOff );
 v_c05f[ wFi ]->set_par( dblRelAcc , EpsCurr );

 if( FiStatus[ wFi ] ==  kUnEval )
  FiStatus[ wFi ] = v_c05f[ wFi ]->compute( true );
 else {
  FiStatus[ wFi ] = v_c05f[ wFi ]->compute( false );

  if( UpFiLmb1[ NrFi ] < Inf<OFValue>() )
   UpFiLmb1[ NrFi ] -= UpFiLmb1[ wFi ];

  if( LwFiLmb1[ NrFi ] > -Inf<OFValue>() )
   LwFiLmb1[ NrFi ] -= LwFiLmb1[ wFi ];

  }

 UpFiLmb1[ wFi ] = std::min( v_c05f[ wFi ]->get_upper_estimate() , UpFiLmb1[ wFi ] );
 LwFiLmb1[ wFi ] = std::max( v_c05f[ wFi ]->get_lower_estimate() , LwFiLmb1[ wFi ] );

 if( UpFiLmb1[ NrFi ] < Inf<OFValue>() )
  UpFiLmb1[ NrFi ] += UpFiLmb1[ wFi ];
 else
  if( UpFiLmb1[ wFi ] < Inf<OFValue>() ) {
   if( linear_function )
	UpFiLmb1[ NrFi ] = linear_function->get_upper_estimate();
   else
    UpFiLmb1[ NrFi ] = 0;
   for( Index k ; k < NrFi ; k++ )
	if( UpFiLmb1[ k ] < Inf<OFValue>() )
     UpFiLmb1[ NrFi ] += UpFiLmb1[ wFi ];
	else {
     UpFiLmb1[ NrFi ] = Inf<OFValue>();
     break;
	 }
   }

 if( LwFiLmb1[ NrFi ] > -Inf<OFValue>() )
  LwFiLmb1[ NrFi ] += LwFiLmb1[ wFi ];
 else
  if( LwFiLmb1[ wFi ] > -Inf<OFValue>() ) {
   if( linear_function )
    LwFiLmb1[ NrFi ] = linear_function->get_lower_estimate();
   else
    LwFiLmb1[ NrFi ] = 0;
   for( Index k ; k < NrFi ; k++ )
    if( LwFiLmb1[ k ] > -Inf<OFValue>() )
     LwFiLmb1[ NrFi ] += LwFiLmb1[ wFi ];
    else {
     LwFiLmb1[ NrFi ] = -Inf<OFValue>();
     break;
     }
   }

 if( UpFiLmb1[ NrFi ] == Inf<OFValue>() )  // Fi() is not defined in Lambda1
  DeltaFi = Inf<OFValue>();
 else
  DeltaFi = UpRifFi[ NrFi ] - UpFiLmb1[ NrFi ];

 FiEvaltns++;

 // update FiBest, if necessary - - - - - - - - - - - - - - - - - - - - - - -

 if( UpFiLmb1[ NrFi ] < UpFiBest[ NrFi ] ) {
  UpFiBest = UpFiLmb1;
  if( KpBstL )
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
   if( UpFiLmb1[ wFi ] == Inf<OFValue>() ) {
    HasLinearization = v_c05f[ wFi ]->has_linearization( diagonal = false );
    if( !HasLinearization )
     HasLinearization = v_c05f[ wFi ]->has_linearization( diagonal );
    }
   else
	HasLinearization = v_c05f[ wFi ]->has_linearization( diagonal );

   }
  else {

   if( UpFiLmb1[ wFi ] == Inf<OFValue>() ) {
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

  double* G1 = Master->GetItem( wFi );

  // fetch the item from the Oracle - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set SGBse = nullptr;
  v_c05f[ wFi ]->get_linearization_coefficients( G1 );
  HpNum Alfa1k = v_c05f[ wFi ]->get_linearization_constant();

  GiEvaltns++;

  // pass the base to the MP Solver - - - - - - - - - - - - - - - - - - - - -

  Master->SetItemBse( SGBse , NumVar );

  // calculate ScPr1k and Alfa1k- - - - - - - - - - - - - - - - - - - - - - -

  Index cp;
  HpNum ScPr1k;

  // update alpha value at Lambda1 point  - - - - - - - - - - - - - - - - - -

  Alfa1k = UpFiLmb1[ wFi ] - Alfa1k
		  - std::inner_product( Lambda1.begin() , Lambda1.end() , G1 , 0 ); //??

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

   #if( LOG_BND )
    if( LogVerb ) {
     if( !diagonal )
      *f_log << std::endl << "New constraint " << wh << ", rhs = " << Alfa1k;
     else
      *f_log << std::endl << "New eps-subgradient " << wh << " for Fi[ "
	     << wFi << " ], eps = " << Alfa1k << ", gd = " << - ScPr1k;
     }
   #endif
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

 // do the move - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Lambda.swap(Lambda1);
 UpFiLmb.swap(UpFiLmb1);
 UpRifFi = UpFiLmb;

 // change the current point in the MP Solver - - - - - - - - - - - - - - - -

 Master->ChangeCurrPoint( t , UpFiLmb1.data() );

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
 double LwrBnd = f_Block->get_valid_lower_bound( false ); // ??
 if( LwrBnd > - Inf<double>() )
  TrueLB = true;
 else {
  TrueLB = false;
  LwrBnd = f_Block->get_valid_lower_bound( true );
  }

if( LwrBnd != LowerBound ) {
  LowerBound = LwrBnd;
  LBHasChgd = true;
  }

 } // end( BundleSolver::UpdtLowerBound )  - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

double BundleSolver::BetaK( Index wFi ) {

return( 1.0 / double( NrEasy ) );

} // end( BundleSolver::BetaK )  - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::Log1( void )
{
 #if( LOG_BND )
  if( LogVerb > 1 ) {
   *f_log << std::endl << "{" << SCalls << "-" << ParIter << "-"
	   << Master->MaxName() - FreDim << "-" << MBDim << "} t = " << t
	   << " ~ D*_1( z* ) = " << Master->ReadDStart( 1 )
	   << " ~ Sigma = " << Sigma << std::endl << "           ";

   *f_log <<  " Fi = ";

   if( UpFiLmb[NrFi] == Inf<double>() )
    *f_log << " - INF";
   else
    *f_log <<  -UpFiLmb[NrFi] << " ~ eU = " << EpsU;

   if( BPar6 )
    *f_log << " ~ BP3 = " << aBP3;
   }
 #endif
 } // end( BundleSolver::Log1 )  - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::Log2( void )
{
 #if( LOG_BND )
  if( LogVerb > 1 ) {
   *f_log << std::endl << "           ";

   if( LowerBound > - Inf<double>() )
    *f_log << "UB = " << - LowerBound << " ~ ";

   *f_log << "Fi1 = ";

   if( UpFiLmb1[NrFi] == - Inf<double>() )
    *f_log << "+ INF => STOP." << std::endl;
   else
    if( UpFiLmb1[NrFi] == Inf<double>() )
     *f_log << " - INF" << std::endl;
    else
     *f_log << - UpFiLmb1[NrFi] << " ~ Alfa1 = " << Alfa1[NrFi]
	     << " ~ Gi1xd = " << - ScPr1[NrFi] << std::endl;
   }
 #endif
 } // end( BundleSolver::Log2 )  - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::InitMP( void )
{
 // this method is called only when *both* the FiOracle *and* the MPSolver
 // have been set, and it is re-called each time any one of the two changes
 // set the size- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Master->SetDim( BPar2 , &FakeFi , false );

 Master->SetPar( MPSolver::kOptEps , RAccSol );
 Master->SetPar( MPSolver::kFsbEps , FAccSol );

 // insert the constant subgradient of the 0-th component - - - - - - - - - -

 if( linear_function ) {
  linear_function->get_linearization_coefficients( Master->GetItem( 0 ) );
  const Index* SGBse = nullptr;
  Master->SetItemBse( SGBse , NumVar );
  Master->SetItem( InINF );
  }

 tHasChgd = LBHasChgd = true;

 }  // end( BundleSolver::InitMP( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::FindNext( Index &wFi ) {

 Index PrWFi = wFi;
 bool NextIsAccepted = false;
 do {
  wFi = ( wFi == NrFi - 1 )? 0 : wFi + 1;
  if( ( FiStatus[ wFi ] == kUnEval ) ||
	 ( FiStatus[ wFi ] < kError && FiStatus[ wFi ] > kOK
	         && CurrNrEvls[ wFi ] != MaxNrEvls ) )
   NextIsAccepted = true;
  } while( !NextIsAccepted && ( wFi != PrWFi ) );

 if( NextIsAccepted )
  return( true );
 else
  return( false );

 } // end( BundleSolver::FindNext( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::BStrategy( cIndex wFi )
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
 if( FreDim || ( Master->MaxName() < Index( BPar2 ) ) )
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
 for( Index i = 0 ; ++i < Index( BPar2 ) ; )
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

 if( FreDim ) {              // there are deleted items
  wh = FreList.front();      // pick one
  std::pop_heap (FreList.begin(),FreList.end());
  FreList.pop_back();
  FreDim = FreList.size();
  }
 else                                       // there are no deleted items ...
  if( Master->MaxName() < Index( BPar2 ) )  // ... but there is still space
   wh = Master->MaxName();                  // next name

 assert( Master->MaxName() >= FreDim );

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
 if( std::abs( vStar[ NrFi ] - DeltaFi ) < Eps<double>() )
  return( tMaior );
 else
  return( t * ( vStar[ NrFi ] / ( 2 * ( vStar[ NrFi ] - DeltaFi ) ) ) );
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

 if( linear_function ) {
  int GPMaxSz = linear_function->get_int_par( C05Function::intGPMaxSz );
  for( Index i = 0 ; i < GPMaxSz ; i++ )
   linear_function->delete_linearization( i );
  }

 FreDim = 0;

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
 FreList.clear();
 OOBase.clear();

 if( !IsEasy.empty() )
  IsEasy.clear();

 }  // end( BundleSolver::MemDealloc( ) )  - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::ReSetAlg( unsigned char RstLvl )
{
 if( ! ( RstLvl & RstAlg ) ) {  // reset algorithmic parameters - - - - - - -
  ParIter = ParSS = 0;     // reset iterations count
  CSSCntr = CNSCntr = 0;   // ... comprised consecutive NS/SS count

  if( t != tInit ) {       // reset t
   t = tInit;
   tHasChgd = true;
   }

  CmptaBPX();  // reset the dynamic number of fetched items

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
       v_c05f[ Master->WComponent( i ) - 1 ]->delete_linearization( i ); // se e' di una linear??
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
    	v_c05f[ Master->WComponent( i ) - 1 ]->delete_linearization( i ); // se e' di una linear??
	    Delete( i );
        }
     }

 if( ! ( RstLvl & RstFiV ) )  // reset the current value of Fi( Lambda ) - - -
  UpFiLmb[NrFi] = Inf<double>();

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

 FreList.push_back( int(i) ); // ?? controllare ??
 push_heap(FreList.begin() , FreList.end());
 FreDim++;

 OOBase[ i ] = Inf<SIndex>();

 // compacting FreList[] if it's too big- - - - - - - - - - - - - - - - - - -
 // remove from FreList[] every name >= Master->MaxName(); note that every
 // ordered set *is* a Heap. apart from efficiency reasons, this is
 // needed because Master->MaxName() - FreDim is the only way in which the
 // Bundle can compute the number of "live" items

 cIndex MxNm = Master->MaxName();
 if( FreDim > MxNm ) {
  FreDim = 0;
  for( Index i = 0 ; i < MxNm ; i++ )
   if( OOBase[ i ] == Inf<SIndex>() )
    FreList[ FreDim++ ] = i;
  }

 assert( Master->MaxName() >= FreDim );

 }  // end( BundleSolver::Delete() ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdtaBP3( void )
{
 const Index tBP3 = ( BPar3 >= 0 ? Index( std::ceil( BPar3 ) ) :
		              Index( std::ceil( NrFi * ( - BPar3 ) ) ) );
 switch( BPar6 ) {
  case( 4 ):
   if( UpFiLmb[NrFi] > -Inf<double>() )
    aBP3 = ( BPar5 > 0 ? aBP4 : tBP3 ) +
           Index( BPar5 / std::log10( EpsU / RelAcc ) );
    break;
  case( 3 ):
   if( UpFiLmb[NrFi] > -Inf<double>() )
    aBP3 = ( BPar5 > 0 ? aBP4 : tBP3 ) +
                         Index( BPar5 / std::sqrt( EpsU / RelAcc ) );
   break;
  case( 2 ):
   if( UpFiLmb[NrFi] > -Inf<double>() )
    aBP3 = ( BPar5 > 0 ? aBP4 : tBP3 ) +
           Index( BPar5 * ( RelAcc / EpsU ) );
   break;
  case( 1 ):
   if( BPar5 && ( ! ( ParIter % Index( std::abs( BPar5 ) ) ) ) ) {
    if( BPar5 > 0 )
     aBP3++;
    else
     aBP3--;
    }
  }

 if( aBP3 > tBP3 )
  aBP3 = tBP3;
 else
  if( aBP3 < aBP4 )
   aBP3 = aBP4;

 }  // end( BundleSolver::UpdtaBP3 ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::CmptaBPX( void )
{
 aBP3 = ( BPar3 >= 0 ? Index( ceil( BPar3 ) ) :
                       Index( ceil( NrFi * ( - BPar3 ) ) ) );

 aBP4 = std::min( aBP3 , ( BPar4 >= 0 ? Index( ceil( BPar4 ) ) :
			                Index( ceil( NrFi * ( - BPar4 ) ) ) )
		  );

 if( BPar6 && ( BPar5 > 0 ) )
  aBP3 = aBP4;
 }  // end( BundleSolver::CmptaBPX() ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::IsOptimal( double eps ) const
{
 double FiL = UpFiLmb[0];

 if( FiL == Inf<double>() || vStar[ NrFi ] == -Inf<double>() )
  return( false );
 else {
  if( FiL < 0 ) FiL = - FiL;
  if( FiL < 1 ) FiL = 1;

  if( eps <= 0 )
   eps = RelAcc;

  if( tStar > 0 )
   return( DSTS + Sigma <= eps * FiL );
  else
   return( ( Master->ReadDStart( 1 ) <= AAccSol ) && ( Sigma <= eps * FiL ) );

  }
 } // end( BundleSolver::IsOptimal() ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::CheckAlfa( const bool All )
{
 return( Sigma >= - t * m3 * Master->ReadDStart( t ) );
 }  // end( CheckAlfa )  - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

BundleSolver::FakeFiOracle::FakeFiOracle( BundleSolver *solver ) : FiOracle()
{
 bslv = solver;

 GiNameVcblr.resize( bslv->BPar2 );
 auto it =  GiNameVcblr.begin();
 for( Index i = 0 ; i < bslv->v_c05f.size() ; ++i )
  for( Index j = 0 ; j < bslv->v_c05f[i]->get_int_par( C05Function::intGPMaxSz);
       j++ )
   *it = std::make_tuple( j , i , true );
 } // end ( FakeFiOracle::FakeFiOracle( ) )  - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetNDOSolver( NDOSolver *NwSlvr ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::SetNDOSolver() ) - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetFiLog( ostream *outs , const char lvl ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::SetFiLog() ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetFiTime( const bool TimeIt ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::SetFiTime() )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetMaxName( cIndex MxNme ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::SetMaxName() ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetNumVar( void ) const {
 return( bslv->NumVar );
 } // end ( FakeFiOracle::GetNumVar() )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetNrFi( void ) const {
 return( bslv->v_c05f.size( ) );
 } // end ( FakeFiOracle::GetNrFi() )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetMaxName( void ) const {
 return( bslv->BPar2 );
 } // end ( FakeFiOracle::GetMaxName() ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetMinusInfinity( void ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetMinusInfinity() ) - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetMaxNZ( cIndex wFi ) const {
 if( wFi != Inf<Index>() )
  throw( std::logic_error( "GetMaxNZ can be called with wFi = Inf only" ) );
 return( bslv->NumVar );
 } // end ( FakeFiOracle::GetMaxNZ() ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetMaxCNZ( cIndex wFi ) const {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetMaxCNZ() )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::FakeFiOracle::GetUC( cIndex i ) {

 double lb_value = bslv->LamVcblr[ i ]->get_lb();
 if( lb_value == -Inf<ColVariable::VarValue>() )
  return( true );

 if( lb_value != ColVariable::VarValue(0) )
  throw( std::logic_error( "any value different from zero is not allowed" ) );

 return( false );
 } // end ( FakeFiOracle::GetUC() )  - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

LMNum BundleSolver::FakeFiOracle::GetUB( cIndex i ) {
 return( bslv->LamVcblr[ i ]->get_ub() );
 } // end ( FakeFiOracle::GetUB() )  - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

LMNum BundleSolver::FakeFiOracle::GetBndEps(  ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetBndEps() )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetGlobalLipschitz( cIndex wFi ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetGlobalLipschitz() )   - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetBNC( cIndex wFi ) {

 return( bslv->MILP_s[ wFi - 1 ]->get_numcols() );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetBNR( cIndex wFi ) {

 return( bslv->MILP_s[ wFi -1 ]->get_numrows() );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetBNZ( cIndex wFi ) {

 return( bslv->MILP_s[ wFi -1 ]->get_nzelements() );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::GetBDesc( cIndex wFi , int *Bbeg , int *Bind , double *Bval ,
 			  double *lhs , double *rhs , double *cst ,
 			  double *lbd , double *ubd ) {

 auto MILPSlv = bslv->MILP_s[wFi-1];

 int num_col = MILPSlv->get_numcols();

 std::copy( MILPSlv->get_matbeg().begin() , MILPSlv->get_matbeg().end() , Bbeg );
 std::copy( MILPSlv->get_matind().begin() , MILPSlv->get_matind().end() , Bind );
 std::copy( MILPSlv->get_matval().begin() , MILPSlv->get_matval().end() , Bval );

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

NDOSolver * BundleSolver::FakeFiOracle::GetNDOSolver( void ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetNDOSolver() ) - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SETTING LAMBDA ------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::SetLambda( cLMRow Lmbd ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

void BundleSolver::FakeFiOracle::SetLamBase( cIndex_Set LmbdB  , cIndex LmbdBD ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

bool BundleSolver::FakeFiOracle::SetPrecision( HpNum Eps ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR COMPUTING Fi() ----------------------*/
/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::Fi( cIndex wFi ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::Fi( ) )  - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/

bool BundleSolver::FakeFiOracle::NewGi( cIndex wFi ) {
 if( wFi == 0 )
  throw( std::invalid_argument( "asking for the 0th component" ) );
 last_c05 =  wFi-1;
 return( true );
 } // end ( FakeFiOracle::NewGi( ) ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetGi( SgRow SubG , cIndex_Set &SGBse ,
			cIndex Name , cIndex strt , Index stp  ) {

 bslv->v_c05f[ std::get<1>(GiNameVcblr[Name]) ]->get_linearization_coefficients(
	 SubG , std::get<0>(GiNameVcblr[Name]) , {} , strt , stp );

 SGBse = nullptr;
 return( stp - strt );
 } // end ( FakeFiOracle::GetGi( ) )   - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetVal( cIndex Name )
{
 return( bslv->v_c05f[ std::get<1>(GiNameVcblr[Name]) ]->
 		 get_linearization_constant( std::get<0>(GiNameVcblr[Name]) ) );
 } // end ( FakeFiOracle::GetVal( ) )  - - - - - - - - - - - - - - - - - - - -

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

 bslv->v_c05f[ last_c05 ]->store_linearization( std::get<0>( *it ) );

 } // end ( FakeFiOracle::SetGiName( ) ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/

HpNum BundleSolver::FakeFiOracle::GetLowerBound( cIndex wFi ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }  // end ( FakeFiOracle::GetLowerBound( ) )  - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

FiOracle::FiStatus BundleSolver::FakeFiOracle::GetFiStatus( Index wFi ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }  // end ( FakeFiOracle::GetFiStatus( ) )  - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::Deleted( cIndex i ) {

 bslv->v_c05f[ std::get<1>(GiNameVcblr[i]) ]->
     delete_linearization( std::get<0>(GiNameVcblr[i]) );

 std::get<2>(GiNameVcblr[i]) = true;
 } // end ( FakeFiOracle::Deleted( ) ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::Aggregate( cHpRow Mlt , cIndex_Set NmSt ,
		cIndex Dm , cIndex NwNm ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }  // end ( FakeFiOracle::Aggregate( ) )  - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*----------------------- End File BundleSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
