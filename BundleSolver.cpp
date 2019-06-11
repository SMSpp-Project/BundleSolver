/*--------------------------------------------------------------------------*/
/*------------------------ File BundleSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BunldeSolver class.
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINES -----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BundleSolver.h"

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

#define NOISE_REDUCTION_FIRST 1

/* If NOISE_REDUCTION_FIRST > 0, then in case of trouble the Bundle will first
   try to raise t, then ask the FiOracle to increase the accuracy. Otherwise
   the order is reversed. */

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BundleSolver to the Solver factory
SMSpp_insert_in_factory_cpp_0( BundleSolver );

/*--------------------------------------------------------------------------*/
// define and initialize here the vector of int parameters names
const std::vector< std::string > BundleSolver::int_pars_str =
             { "intBPar1" , "intBPar2" , "intBPar6" ,
               "intEStps" , "intMnSSC" , "intMnNSC" ,
               "inttSPar1" , "intPPar1", "intPPar2" ,
			   "intSPar3" };

// define and initialize here the vector of double parameters names
const std::vector< std::string > BundleSolver::dbl_pars_str =
		     { "dbltStar"  , "dblEInit" , "dblEFnal"   ,
		       "dblEDcrs" , "dblBPar3" ,  "dblBPar4"  ,
		       "dblBPar5"  , "dblm1" , "dblm3" ,
			   "dblmxIncr" ,  "dblmnIncr" ,  "dblmxDecr" ,
			   "dblmnDecr" ,  "dbltMaior" ,  "dbltMinor" ,
			   "dbltInit" ,  "dbltSPar2" ,  "dblMPEFsb" ,
			   "dblMPEOpt"  };

// define and initialize here the map for int parameters names
const std::map< std::string , BundleSolver::idx_type > BundleSolver::int_pars_map =
                   { { "intBPar1"  , BundleSolver::intBPar1  } ,
		     // { "intBPar2" , BundleSolver::intBPar2 } ,
		     { "intBPar6" , BundleSolver::intBPar6 } ,
		     { "intEStps" , BundleSolver::intEStps } ,
		     { "intMnSSC" , BundleSolver::intMnSSC } ,
		     { "intMnNSC" , BundleSolver::intMnNSC } ,
		     { "inttSPar1", BundleSolver::inttSPar1 } ,
		   // { "intPPar1" , BundleSolver::intPPar1 } ,
		   // { "intPPar2" , BundleSolver::intPPar2 } ,
		   // { "intPPar3" , BundleSolver::intPPar3 }
			 };

// define and initialize here the map for double parameters names
const std::map< std::string , BundleSolver::idx_type > BundleSolver::dbl_pars_map =
                   { { "dbltStar" , BundleSolver::dbltStar } ,
		     { "dblEInit" , BundleSolver::dblEInit } ,
			 { "dblEFnal" , BundleSolver::dblEFnal } ,
			 { "dblEDcrs" , BundleSolver::dblEDcrs } ,
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
			 { "dblMPEFsb", BundleSolver::dblMPEFsb } ,
		     { "dblMPEOpt"  , BundleSolver::dblMPEOpt  } };

// define and initialize here the default int parameters
const std::vector<int> BundleSolver::dflt_int_par =
        {    10 ,  // intBPar1
		//	100 ,  // intBPar2
			  0 ,  // intBPar6
			  0 ,  // intEStps
			  0 ,  // intMnSSC
			  0 ,  // intMnNSC
			  0   // intSPar1
		//	 30 ,  // intPPar1
		//	 10 ,  // intPPar2
		//	  5    // intPPar3
			 };

// define and initialize here the default double parameters
const std::vector<double> BundleSolver::dflt_dbl_par =
           { 1e2 ,    // dbltStar
			 1e-2 ,   // dblEInit
			 1e6 ,    // dblEFnal
			 0.95 ,   // dblEDcrs
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
			  0.1 ,   // dbltSPar2
			 1e-6 ,   // dblMPEFsb
			 1e-6     // dblMPEOpt
               };

/*--------------------------------------------------------------------------*/

static const HpNum Nearly  = 1.01;
static const HpNum Nearly2 = 1.02;

static const HpNum DefMPEFsb = 1e-6;  // default value for MPEFsb
static const HpNum DefMPEOpt = 1e-6;  // default value for MPEOpt

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



int BundleSolver::compute( bool changedvars )
{
 // basic sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! Master )
  throw( std::logic_error( "Master not set yet" ) );

 if( v_c05f.empty() )
  throw( std::logic_error( "C05Function not set yet" ) );

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

  if( Result )  // problems in the Master Problem solver
   break;

  // a little bookkeeping - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // update out-of-base counters- - - - - - - - - - - - - - - - - - - - - - -

  UpdtCntrs();

  // some log - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log1();

  // contrast MPsolver::Alfa[] with Oracle::Alfa[]- - - - - - - - - - - - - -

  StrongCheckAlfa();

  // hook for derived classes - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int rs = EveryIteration();

  if( rs == kEIAbort ) {
   BLOG( 1 , " EveryIteration():STOP" << std::endl );
   Result = kStopIter;
   break;
   }

  if( rs == kEILoopNow ) {
   BLOG( 1 , " EveryIteration():loop" << std::endl );
   continue;
   }

  // check the status of the FiOracle and take the necessary action - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // FiOracle::FiStatus fs = Oracle->GetFiStatus(); // ?? cosa mettere ??

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( IsOptimal() )
   break;

  // Hard Long-Term t-strategy- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the hard long-term t-strategy requires t to increase if the step is too
  // small, and therefore has to be checked before the others

  if( ( ( tSPar1 & tSP1Msk ) == kHLTTS ) && ( FiLambda[0] < Inf<double>() ) ) {
   HpNum AFL = std::abs( FiLambda[0] );
   if( AFL < 1 )
    AFL = 1;

   if( vStar <= tSPar2 * EpsU * AFL ) {
    BLOG( 1 , "small v => increase t" << std::endl << "           " );

    // collect two numbers vc and vl such that v( tNew ) >= vc + tNew * vl
    // we require that v( tNew ) >= vc + tNew * vl = tSPar2 * EpsU * AFL
    // ==> tNew = ( tSPar2 * EpsU * AFL - vc ) / vl

    HpNum vl , vc;
    Master->SensitAnals( vl , vc );

    HpNum tt;
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
   }  // end if( Hard t-strategy )


  // a real iteration (iterations where Fi() is not evaluated do not count) -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ParIter++;

  // change the "precision" in computing Fi() - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is the "regular mechanism"; gradually decrease the precision along
  // the iterations

  if( EStps ) {
   cIndex nstps = ceil( double( EStps > 0 ? ParIter : ParSS ) /
			double( std::abs( EStps ) ) );
   EpsCurr = EDcrs >= 0 ? std::abs( EInit ) : EpsU;
   if( EFnal >= 0 )
    EpsCurr *= pow( std::abs( EDcrs ) , EFnal * nstps );
   else
    EpsCurr *= std::abs( EDcrs ) * ( nstps ? pow( nstps , EFnal ) : 1 );

   EpsCurr = std::max( RelAcc , std::min( EpsCurr , std::abs( EInit ) ) );
   }
  else
   EpsCurr = EDcrs >= 0 ? std::abs( EInit ) : ( std::abs( EDcrs ) * EpsU );

  if( ( EpsCurr < EpsFi ) ||
     ( ( EpsCurr > EpsFi ) && ( EInit < 0 ) ) ) {
   // only allow increasing EpsFi if EInit < 0, always allow decreasing it
   BLOG( 1 , " ~ changing precision to " << EpsCurr << std::endl
	 << "           " );
   // Oracle->SetPrecision( EpsFi = EpsCurr ); // bisgona cambiare precisione??
   }

  // calculate Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormLambda1( t );

  // update the number of items to be fetched from the oracle - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UpdtaBP3();

  // eliminate outdated info- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // This is done *after* the call to Oracle->GetFiStatus(), as well as after
  // the call to Master->SensitAnals() in the Hard Long-Term t-strategy and
  // to FormLambda1(), because elimination of items from the bundle may make
  // the current solution of the master problem invalid, and therefore all
  // solution information may be lost. In theory this should not happen, since
  // only items "out of base" are eliminated, and therefore the solution
  // remains optimal; however, not all MPSolvers may behave in this respect.

  SimpleBStrat();

  // calculate Fi( Lambda1 )- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool MPchgs;  // true if no cycling will occur
  for( ;; ) {   // ... possibly more than once due to precision issues

   // actually compute Fi and collect subgrads- - - - - - - - - - - - - - - -
   // meanwhile check if the new subgradients change the CP model enough
   MPchgs = FiAndGi();

   if( Result == kError )
    break;

   // if there are negative alphas, something has indeed changed
   if( FiLambda[0] < Inf<double>() )
    MPchgs |= CheckAlfa();

   MPchgs |= DoSS();  // doing a SS clearly changes the MP

   if( MPchgs )       // if something changes
    break;            // all done

   // check for running time - - - - - - - - - - - - - - - - - - - - - - - - -
   // if we get here there is something wrong with the FiOracle's precision;
   // the possible solution is to give it a few more resources to try to do
   // it, but this is not possible if we have ran out of time

   // ?? mettere tempo ??

   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // if neither a SS is done nor the Master Problem changes, there must be
   // something wrong with the FiOracle's precision

   #if( NOISE_REDUCTION_FIRST  )
    // we try to patch this by playing with t first, which corresponds to
    // saying that we think it faster to re-solve the Master Problem rather
    // than to get more precision from the FiOracle

    if( ( DSTS >= RelAcc * std::max( HpNum( 1 ) , std::abs( FiLambda[0] ) ) ) &&
  	( t < tMaior ) ) {
       t = std::min( t * mxIncr , tMaior );
       BLOG( 1 , " ~ noise reduction: t increased to " << t << std::endl );
       tHasChgd = true;
       break;
       }
   #endif

   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // Now check if any of the components has not been solved with the required
   // accuracy yet; in this case allow the component(s) to be re-computed
   // again (to increase accuracy). Note that you expect this to happen only
   // a limited number of times, both because eventually a "good enough"
   // solution will be found by the oracle, and because eventually time will
   // run out.

   FiOracle::FiStatus FiStt = FiOracle::kFiNorm;
   for( Index k = 0 ; k++ < NrFi ; )
    if( ( ( IsEasy.empty() ) || ( ! IsEasy[ k ] ) ) &&
      ( FiStatus[ k ] == FiOracle::kFiStop ) ) {
     FiStt = FiOracle::kFiStop;
     break;
     }

   if( FiStt == FiOracle::kFiStop )  // at least one non-easy component did
    continue;                        // not finish, go give them more time

   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // the oracle reports to have computed the function to the required
   // precision, but this is not enough: here goes the "emergency mechanism"
   // that tries to increase the accuracy, but of course this all depends on
   // if the oracle is actually available to do it

   EpsCurr = std::max( RelAcc , std::min( EpsCurr , EpsFi ) * std::abs( EDcrs ) );

   // ?? precisione oracolo da cambaire ??
   /* if( EpsCurr < EpsFi )  // unless we have already hit EpsLin
    if( Oracle->SetPrecision( EpsFi = EpsCurr ) ) {  // if it is
     BLOG( 1 , " ~ increasing precision to " << EpsCurr << std::endl );
     VectAssign( FiStatus + 1 , FiOracle::kFiStop , NrFi );
     // reset FiStatus[]: if a component has said "kFiOK" before, with a
     // coarser accuracy, this does not mean the computation is still OK
     // now that the accuracy has increased, so we must assume it is not
     continue;                                       // go compute Fi() again
     } */

   #if( ! NOISE_REDUCTION_FIRST  )
    // we try to patch this by playing with t last, which corresponds to
    // saying that we think it faster to get more precision from the
    // FiOracle rather than to re-solve the Master Problem

    if( ( DSTS >= RelAcc * std::max( HpNum( 1 ) , std::abs( FiLambda[0] ) ) ) &&
  	   ( t < tMaior ) ) {
     t = std::min( t * mxIncr , tMaior );
     BLOG( 1 , " ~ noise reduction: t increased to " << t << std::endl );
     tHasChgd = true;
     break;
     }
   #endif

   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // Now the Bundle is potentially in trouble: the new informations does not
   // change the model, changing t does not help, and the oracle is not
   // willing to provide any more precision right now. However, the oracle
   // may "know" this, this being an iteration where the Bundle would have
   // stopped already, but it was told not to by the oracle. This is clearly
   // at risk of cycling, so the oracle must have some internal mechanism to
   // avoid that. in this case, trust the oracle: just go to the next
   // iteration and hope that something eventually will change.

   // ?? vedere parte di sotto ??
   /*
   if( fs == FiOracle::kFiCont ) {  // if the oracle wants to rule
    MPchgs = true;                  // pretend to believe no cycling will
    break;                          // occur since the oracle says so
    }
   else {                           // really, nothing else to do but quit
    BLOG( 1 , " ~ too low precision in the FiOracle" << std::endl );
    Result = kLwPrcsn;
    break;
    } */

   }  // end( for( ever ) )

  // check whether the Lower Bounds have changed- - - - - - - - - - - - - - -

  UpdtLowerBound();

  // some log about the newly obtained information- - - - - - - - - - - - - -

  Log2();

  if( FiLambda1[0] == - Inf<double>() ) {
   Result = kUnbounded;
   break;
   }

  // check whether either any error has occurred or time has expired- - - - -

  if( ( Result == kError ) || ( Result == kStopTime ) )
   break;

  // if( Result == kLwPrcsn ) ??
  //  break;

  if( ( ~ MPchgs ) && tHasChgd )  // "noise reduction": t has changed,
   continue;                      // so go solve the master problem again
                                  // (no NS/SS decision can be made)

  // check the Lower Bound- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // note: TrueLB is true if *LowerBound > GetMinusInfinity(). Termination
  //       "by objective function value only" is only enabled if TrueLB is
  // true. In the Lagrangian case, if one sets as LowerBound the value of a
  // feasible solution it may stop here is that solution is EpsLin-optimal.
  // However, doing so might "disrupt the convexified solution", because the
  // Master Problem is not solved and therefore the optimal multipliers are
  // not computed. To avoid that, the Oracle can return the same value as
  // GetMinusInfinity(), thereby disabling this termination test and leaving
  // only the standard one using the Master Problem solution. However, if
  // unboundedness was to be declared when *FiBest <= *LowerBound, in this
  // case one could end up declaring the problem unbounded below. This is why
  // a value slightly smaller than *LowerBound is used instead.

  if( FiLambda[0] < Inf<double>() ) {  // .. but only if Fi( Lambda ) is defined
   if( TrueLB )
    if( FiBest[0] - RelAcc * std::abs( FiBest[0] ) <= LowerBound[0] )
     break;

   if( FiBest[0] <= LowerBound[0] *
                  ( 1 - ( LowerBound[0] > 0 ? RelAcc : - RelAcc ) ) ) {
    Result = kUnbounded;
    break;
    }
   }

  // avoid the t-changing phase if Lambda1 is unfeasible- - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // note: one possible alternative t-strategy would be to set t to the
  // largest value that would have produced a feasible point, i.e.
  // t := *Alfa1 / ( - *ScPr1 )

  if( FiLambda1[0] == Inf<double>() )
   continue;
  else
   if( FiLambda[0] == Inf<double>() ) {  // if reached feasibility  - - - - -
    GotoLambda1();             // go to the feasible point
    continue;                  // and start the actual minimization of Fi()
    }

  // the NS / SS decision - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SSDone = DoSS();

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
    if( Alfa1[0] <= m3 * Sigma ) {
     BLOG( 1 , " ~ small Alfa1" );
     tt = t;
     }
    else
     switch( tSPar1 & tSP1Msk ) {
      case( kSLTTS ):
      case( kHLTTS ):
       if( vStar <= tSPar2 * EpsU * std::max( std::abs( FiLambda[0] ) , HpNum( 1 ) ) ) {
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
   if( DSTS < RelAcc * std::max( std::abs( FiLambda[0] ) , HpNum( 1 ) ) / 10 ) {
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

  auto obj = boost::any_cast<FRealObjective *>( f_Block->get_objective() );
  if( obj == nullptr )
   throw( std::logic_error( "the objective is not a real function" ) );

  auto c05f = dynamic_cast<C05Function *>( (obj)->get_function() );
  if( c05f == nullptr )
   throw( std::logic_error( "the objective is not a C05Function" ) );

  v_c05f.push_back( c05f );
  linf = nullptr;

  }
 else {

  // the objective function of each block must be a LinearFunction - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( f_Block->get_objective().empty() )
   linf =  nullptr;
  else {
   auto obj = boost::any_cast<FRealObjective *>( f_Block->get_objective() );
   if( obj == nullptr )
    throw( std::logic_error( "the objective is not a real function" ) );

   linf = dynamic_cast<LinearFunction *>( (obj)->get_function() );
   if( linf == nullptr )
    throw( std::logic_error( "the objective is not a LinearFunction" ) );
   }


  auto sb = f_Block->get_nested_Blocks();
  v_c05f.resize( sb.size() );

  for( Index i = 0 ; i < sb.size() ; ++i ) { // for each sub-block

   // the objective function of each sub-block must be a C05Function - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto obj = boost::any_cast<FRealObjective *>( sb[ i ]->get_objective() );
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
  if( un_any_thing_0( ColVariable , el , [ & ]{ ++NumVar; } ) )
   break;
  if( un_any_thing_1( ColVariable , el , [ & ]{ NumVar += var.size(); } ) )
   break;
  if( un_any_thing_K( ColVariable , el , [ & ]{ NumVar += var.size(); } ) )
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

 PPar1 = PPar2 = PPar3 = 0;

 // read information about the function  - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector<Index> BNC( v_c05f.size() );

 if( v_c05f.size() > 1 ) {
  bool HasEasy = false;
  IsEasy.resize( v_c05f.size() );
  for( Index k = 0 ; k < v_c05f.size() ; ++k )
    if( BNC[ k ] )
     IsEasy[ k ] = HasEasy = true;
    else
     IsEasy[ k ] = false;

   if( ! HasEasy ) {
    IsEasy.clear();
    }
   }

 // allocate memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Lambda.resize( NumVar );    // the default starting point
 Lambda1.resize( NumVar );   // the tentative point

 if( KpBstL )  // best point found so far
  LmbdBst.resize( NumVar);
 // else
 // LmbdBst.clear(); ?? come si pulisce lo sparsevector

 OOBase.resize( BPar2 , Inf<SIndex>() );  // counter for eliminating outdated
                                          // items: Inf<SIndex>() means empty

 FreList.resize( BPar2 );       // list of free bundle slots
 whisZ.resize( v_c05f.size() ); // for each component, the name of its "Z" if it is
                                // in the bunlde

 NrFi = ( linf )?  v_c05f.size() + 1 : v_c05f.size();

 FiLambda.resize( NrFi );        // current, ...
 FiLambda1.resize( NrFi );       // tentative, ...
 FiBest.resize( NrFi );          // best, ...
 RfrncFi.resize( NrFi );         // and reference Fi() values

 FiStatus.resize( NrFi , kStopIter );
 LowerBound.resize( NrFi , -Inf<double>() ); // lower bounds
 TrueLB = false;

 RfrncFi.resize( NrFi , Inf<double>() );
 FiLambda1[0] = FiBest[0] = Inf<double>();       // Fi( Lambda ) is not known

 whisG1.resize( v_c05f.size() , Inf<Index>() );  // no representative yet

 ScPr1.resize( NrFi , 0 );
 Alfa1.resize( NrFi , 0 );
 DeltaAlfa.resize( v_c05f.size() );

 FreDim = 0;
 BHasChgd = true;  // ensure SetLamBase() is called at least once
 Result = kError;
 SSDone = false;

 // ReSetAlg( RstCrr | RstSbg | RstCnt );  ?? da fare ??   // Fi( Lambda ) is reset inside

 // Oracle->SetPrecision( EpsFi = ABS( EInit  ) );

 // warning: the following things can only be done *after* that
 // Oracle->SetMaxName() has been invoked, because they use methods of the
 // oracle which depends on knowledge of the MaxName to work properly
 // read b0- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // here one could initialize b0, if that was found to be of any use
 // b0 = Oracle->GetVal( BPar2 );

 // initialize the MP Solver, if any - - - - - - - - - - - - - - - - - - - -

 if( Master ) {
  InitMP();
  FakeFi = new FakeFiOracle( this );
  }

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
  // case( intBPar2 ):
  // BPar2 = value;
  // break;
  case( intBPar6 ):
   BPar6 = value;
   break;
  case( intEStps ):
   EStps = value;
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
 // case( intPPar1 ):
 //  PPar1 = value;
 //  break;
 // case( intPPar2 ):
 //  PPar2 = value;
 //  break;
 // case( intPPar3 ):
 //  PPar3 = value;
 //  break;
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
  case( dblUpCutOff ):
   UpCutOff = value;
   break;
  case( dblLwCutOff ):
   LwCutOff = value;
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
  case( dblEFnal ):
   EFnal = value;
   break;
  case( dblEDcrs ):
   EDcrs = value;
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
  case( dblMPEFsb ):
   MPEFsb = value;
   break;
  case( dblMPEOpt ):
   MPEOpt = value;
   break;
  default:
   CDASolver::set_par( par , value );
  }

 } // end (BundleSolver::set_par( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::SetMPSolver( MPSolver *MPS )
{
 if( Master )        // a MP solver is set ??? dove metterlo ???
  Master->SetDim();  // clear all its internal state

 if( FakeFi )
  delete[] FakeFi;

 // construct a FakeFiOracle to handle the MPSolver, which has to
 // interface with a FiOracle object, the FiOracle needs of all the
 // description of the Function including the 0th component

 Master = MPS;
 if( Master && f_Block ) {
  InitMP();
  FakeFi = new FakeFiOracle( this );
  }

 } // end( BundleSolver::SetMPSolver )  - - - - - - - - - - - - - - - - - - -

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
  return( dflt_dbl_par[ par - intBPar1 ] );
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
  // case( intBPar2 ):
  // return( BPar2 );
  // break;
  case( intBPar6 ):
   return( BPar6 );
   break;
  case( intEStps ):
   return( EStps );
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
 // case( intPPar1 ):
 //  return( PPar1 );
 //  break;
 // case( intPPar2 ):
 //  return( PPar2 );
 //  break;
 // case( intPPar3 ):
 //  return( PPar3 );
 //  break;
  default:
   return( get_dflt_int_par( par ) );
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
  case( dblUpCutOff ):
   return( UpCutOff );
   break;
  case( dblLwCutOff ):
   return( LwCutOff );
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
  case( dblEFnal ):
   return( EFnal );
   break;
  case( dblEDcrs ):
   return( EDcrs );
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
  case( dblMPEFsb ):
   return( MPEFsb );
   break;
  case( dblMPEOpt ):
   return( MPEOpt );
   break;
  default:
   return( get_dflt_dbl_par( par ) );
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

 if( LBHasChgd && ( FiLambda[0] < Inf<double>() ) ) {
  if( LowerBound[0] > - Inf<double>() )
   Master->SetLowerBound( LowerBound[0] - FiLambda[0] );
  else
   Master->SetLowerBound( - Inf<double>() );

  for( Index k = 0 ; k++ < NrFi ; ) {
   if( IsEasy.size() && IsEasy[ k ] )  // skip easy components
    continue;

   if( LowerBound[ k ] > - Inf<double>() )
    Master->SetLowerBound( LowerBound[ k ] - FiLambda[ k ] , k );
   else
    Master->SetLowerBound( - Inf<double>() , k );
   }

  LBHasChgd = false;
  }

 // set termination criterion - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLambda[0] < Inf<double>() )
  Master->SetPar( MPSolver::kZero ,
		  RelAcc * std::max( std::abs( FiLambda[0] ) , double( 1 ) )
		         / std::max( tStar / t , HpNum( 1 ) ) );

          //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for(;;)  // price-in loop- - - - - - - - - - - - - - - - - - - - - - - - - -
 {        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for(;;)  // error-handling loop - - - - - - - - - - - - - - - - - - - - - -
  {        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


   // ensure the MPSolver does not exceed the remaining time
   // Master->SetPar( MPSolver::kMaxTme , MaxTime - NDOt->Read() );


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
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Sigma = Master->ReadSigma();                  // read Sigma*

  vStar = - Master->ReadFiBLambda();            // read v*

  if( IsEasy.size() ) {                                // there are easy components
   for( Index k = 0 ; k++ < NrFi ; )            // read the *exact* Fi-value
    if( IsEasy[ k ] )                           // for all them
     FiLambda1[ k ] = Master->ReadFiBLambda( k );

   if( FiLambda[0] < Inf<double>() )
    for( Index k = 0 ; k++ < NrFi ; )
     if( IsEasy[ k ] )
      vStar += RfrncFi[ k ];
   }

  DSTS = Master->ReadDStart( tStar );           // D_{t*,\beta,x}
  Deltav = vStar;
  if( m1 < 0 )                                  // use - z( P_{t,\beta,x} )
   Deltav -= Master->ReadDt( t );

  // Sigma* + D*_{t*}( -z* ) is the "maximum expected increase" used in
  // the stopping criterion, EpsU is that relative to Fi( Lambda )
  if( FiLambda[0] < Inf<double>() )
   EpsU = ( DSTS + Sigma ) / std::max( std::abs( FiLambda[0] ) , double( 1 ) );
  else
   EpsU = 1;  // ensure EpsU is initialized somehow

  // the z[ i ] have changed, so in principle they are no longer in the
  // bundle: it may be the case that they actually are, but this is
  // taken care of in UpdtCntrs()
  whisZ.assign( NrFi , Inf<Index>() );

  // the scalar products have changed
  ScPr1.assign( NrFi , Inf<double>() );

  if( ! PPar2 )  // no L.V.G.
   return;       // nothing else to do

  // if( LamDim == NumVar )  // all variables are there
  // break;                 // no "price in" to do (but possibly "price out")

  // LVG: price in- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // do pricing only for the first PPar1 iterations and then once every PPar2
  // iterations: however, do it also if convergence is detected or if the
  // subproblem is primal unfeasible

  if( ( Result == kOK ) && ( ParIter > Index( PPar1  ) ) &&
      ( ( ParIter - Index( PPar1  ) ) % PPar2 ) && ( ! IsOptimal() ) )
   return;

  // ask for *all* the entries of d[] - - - - - - - - - - - - - - - - - - - -

  const double* tdir = Master->Readd( true );

  // now construct the new LamBase- - - - - - - - - - - - - - - - - - - - - -

  Index nBD = 0;
  Index oBD = 0;
  double epsDir = Master->EpsilonD() * t;

  /*
  Dynamic generation of variables not handled yet.

  Index_Set NewStuff = LamBase + LamDim;

  if( ! Master->NumNNVars() ) {  // there are no NN variables - - - - - - - -
   for( Index k = 0 ; k < NumVar ; k++ )
    if( LamBase[ oBD ] == k ) {
     oBD++;
     nBase[ nBD++ ] = k;
     }
    else
     if( ABS( tdir[ k ] ) > epsDir ) {
      nBase[ nBD++ ] = *(++NewStuff) = k;
      if( PPar3 )
       InctvCtr[ k ] = 0;
      BLOGb( LogVar , std::endl << " Created variable " << k );
      }
   }
  else
   if( Master->NumNNVars() == NumVar ) {  // there are only NN variables- - -
    for( Index k = 0 ; k < NumVar ; k++ )
     if( LamBase[ oBD ] == k ) {
      oBD++;
      nBase[ nBD++ ] = k;
      }
     else
      if( tdir[ k ] > epsDir ) {
       nBase[ nBD++ ] = *(++NewStuff) = k;
       if( PPar3 )
        InctvCtr[ k ] = 0;

       BLOGb( LogVar , std::endl << " Created variable " << k << " (>= 0)" );
       }
    }
   else {  // there are both NN and UC variables- - - - - - - - - - - - - - -
    for( Index k = 0 ; k < NumVar ; k++ )
     if( LamBase[ oBD ] == k ) {
      oBD++;
      nBase[ nBD++ ] = k;
      }
     else
      if( ( tdir[ k ] > epsDir ) ||
	  ( ( tdir[ k ] < - epsDir ) && ( ! Master->IsNN( k ) ) ) ) {
       nBase[ nBD++ ] = *(++NewStuff) = k;
       if( PPar3 )
        InctvCtr[ k ] = 0;

       BLOGb( LogVar , std::endl << " Created variable " << k );
       BLOG2b( LogVar , Master->IsNN( k ) , " (>= 0)" );
       }
    }  // end else( there are both NN and UC variables )- - - - - - - - - - -

  if( nBD == oBD )  // no changes in LamBase- - - - - - - - - - - - - - - - -
   break;
  else {  // LamBase has changed- - - - - - - - - - - - - - - - - - - - - - -
   BHasChgd = true;                // signal it
   *(++NewStuff) = InINF;          // and put termination marks to the
   nBase[ LamDim = nBD ] = InINF;  // vector of newly added stuff and
                                          // nBase
   // signal the changes to the MP solver
   Master->AddActvSt( LamBase + oBD + 1 , nBD - oBD , nBase );

   // add the entries to Lambda
   nBase--; LamBase--; Lambda--;

   for( ; nBD > oBD ; nBD-- )
    if( LamBase[ oBD ] == nBase[ nBD ] )
     Lambda[ nBD ] = Lambda[ oBD-- ];
    else
     Lambda[ nBD ] = 0;

   nBase++; LamBase++; Lambda++;

   // set the new LamBase
   Swap( LamBase , nBase );
   } */

  break; // ?? giusto ??

  // end price in - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  }  // end( for(;;) )

 /*
 // LVG: price out- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ParIter > Index( PPar1 ) ) {  // skip the first PPar1 iterations
  // (all of them if PPar3 == 0 ==> PPar1 == InINF )

  // ask for the entries of d[] corresponding to "active" variables - - - - -

  cLMRow tdir = Master->Readd( false );

  // now construct the new LamBase- - - - - - - - - - - - - - - - - - - - - -

  Index nBD = 0;
  Index oBD = 0;
  LMNum epsDir = Master->EpsilonD() * t;
  Index_Set DltdStuff = nBase + MaxNumVar;

  if( ! Master->NumNNVars() )  // there are no NN variables - - - - - - - - -
   for( Index k ; ( k = LamBase[ oBD ] ) < InINF ; oBD++ )
    if( ( ABS( Lambda[ oBD ] ) <= epsDir ) && ( ABS( tdir[ k ] ) <= epsDir )
	&& ( ++InctvCtr[ k ] >= Index( PPar3 ) ) ) {
     *(DltdStuff--) = k;
     BLOGb( LogVar , std::endl << " Eliminated variable " << k );
     }
    else {
     Lambda[ nBD ] = Lambda[ oBD ];
     InctvCtr[ nBase[ nBD++ ] = k ] = 0;
     }
  else
   if( Master->NumNNVars() == NumVar )  // there are only NN variables- - - -
    for( Index k ; ( k = LamBase[ oBD ] ) < InINF ; oBD++ )
     if( ( Lambda[ oBD ] <= epsDir ) && ( tdir[ k ] <= epsDir ) &&
         ( ++InctvCtr[ k ] >= Index( PPar3 ) ) ) {
      *(DltdStuff--) = k;
      BLOGb( LogVar , std::endl << " Eliminated variable " << k << " (>= 0)" );
      }
     else {
      Lambda[ nBD ] = Lambda[ oBD ];
      InctvCtr[ nBase[ nBD++ ] = k ] = 0;
      }
   else {  // there are both NN and UC variables- - - - - - - - - - - - - - -
    for( Index k ; ( k = LamBase[ oBD ] ) < InINF ; oBD++ )
     if( ( ABS( Lambda[ oBD ] ) <= epsDir )                   &&
	 ( ( tdir[ k ] <= epsDir ) &&
	   ( ( tdir[ k ] >= - epsDir ) || Master->IsNN( k ) ) ) &&
         ( ++InctvCtr[ k ] >= Index( PPar3 ) ) ) {
      *(DltdStuff--) = k;
      BLOGb( LogVar , std::endl << " Eliminated variable " << k << " (>= 0)" );
      BLOG2b( LogVar , Master->IsNN( k ) , " (>= 0)" );
      }
     else {
      Lambda[ nBD ] = Lambda[ oBD ];
      InctvCtr[ nBase[ nBD++ ] = k ] = 0;
      }
    }  // end else( there are both NN and UC variables )- - - - - - - - - - -

  if( nBD < oBD ) {  // some variables have been eliminated - - - - - - - - -
   BHasChgd = true;                // signal it
   nBase[ LamDim = nBD ] = InINF;  // and put the termination mark to nBase

   // invert the vector of deleted stuff, since it is in reverse order

   Index_Set tDSH = ++DltdStuff;
   for( Index_Set tDST = nBase + MaxNumVar ; tDSH < tDST ; )
    Swap( *(tDSH++) , *(tDST--) );

   // set the new LamBase

   Swap( LamBase , nBase );

   // signal the changes to the MP solver

   Master->RmvActvSt( DltdStuff , oBD - nBD , LamBase );
   }
  }  // end if( skip the first PPar1 iterations ) - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */

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
    if( ( *Mlt >= 1 - MPEOpt ) && Master->IsSubG( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ Master->WComponent( i ) ] = i;
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
    if( ( *Mlt >= 1 - MPEOpt ) && Master->IsSubG( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ Master->WComponent( i ) ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }

 }  // end( UpdtCntrs )

/*--------------------------------------------------------------------------*/

int BundleSolver::EveryIteration( void )
{
 return( kOK );

 }  // end( BundleSolver::EveryIteration ) - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FormLambda1( HpNum Tau )
{
 Master->MakeLambda1( Lambda.valuePtr() , Lambda1.valuePtr() , Tau );

  //  if( ( Oracle->GetBndEps() < MPEFsb ) && Master->NumBxdVars() ) ???

 if( Master->NumBxdVars() ) {
  // as the relative precision required to the MPSolver is not enough to
  // ensure that the bounds on the variables will be satisfied with the
  // precision required by the FiOracle, the (upper and lower) bounds are
  // strictly enforced here

  SparseVector tL1 = Lambda1;

  if( Master->NumNNVars() )             // there are NN vars and UB vars
   if( Master->NumNNVars() == NumVar )  // actually, all variables are NN
  /*  if( PPar2 ) {
     Index_Set tLB = LamBase;
     for( Index h ; ( h = *(tLB++) ) < InINF ; tL1++ ) {
      if( *tL1 < 0 )
       *tL1 = 0;

      cLMNum UBh = Oracle->GetUB( h );
      if( *tL1 > UBh )
       *tL1 = UBh;
      }
     }
    else */
     for( Index h = 0 ; h < NumVar ; h++ ) {
      if( tL1.coeffRef(h) < 0 )
       tL1.coeffRef(h) = 0;

      const double UBh = LamVcblr[h]->get_ub();
      if( tL1.coeffRef(h) > UBh )
       tL1.coeffRef(h) = UBh;
      }
   else                                 // not all variables are NN
   /* if( PPar2 ) {
     Index_Set tLB = LamBase;
     for( Index h ; ( h = *(tLB++) ) < InINF ; tL1++ ) {
      if( Master->IsNN( h ) && ( *tL1 < 0 ) )
       *tL1 = 0;

      cLMNum UBh = Oracle->GetUB( h );
      if( *tL1 > UBh )
       *tL1 = UBh;
      }
     }
    else */
     for( Index h = 0 ; h < NumVar ; h++ ) {
      if( Master->IsNN( h ) && ( tL1.coeffRef(h) < 0 ) )
       tL1.coeffRef(h) = 0;

      const double UBh = LamVcblr[h]->get_ub();
      if( tL1.coeffRef(h) > UBh )
       tL1.coeffRef(h) = UBh;
      }
  else  // there are only UB vars
  /* if( PPar2 ) {
    Index_Set tLB = LamBase;
    for( Index h ; ( h = *(tLB++) ) < InINF ; tL1++ ) {
     const double UBh = LamVcblr[h]->get_ub();
     if( tL1.coeffRef(h) > UBh )
      tL1.coeffRef(h) = UBh;
     }
    }
   else */
    for( Index h = 0 ; h < NumVar ; h++ ) {
     const double UBh = LamVcblr[h]->get_ub();
     if( tL1.coeffRef(h) > UBh )
      tL1.coeffRef(h) = UBh;
     }

  }  // end( if( the bounds have to be enforced ) )

 LHasChgd = true;  // signal that Lambda1 has changed

 FiStatus.assign( NrFi , kStopIter );

 // if( BHasChgd && LamBase )
 //  VectAssign( Lam1Bse , LamBase , LamDim + 1 );

 whisG1.assign( NrFi , Inf<Index>() );
 // reset the representatives

 }  // end( BundleSolver::FormLambda1 )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::FiAndGi( void )
{
 // set the new Lambda1[] and Lam1Bse[], if any - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 /* ??
 if( BHasChgd ) {  // if LamBase has changed, pass the new one to the Oracle
  if( LamDim < NumVar )
   Oracle->SetLamBase( Lam1Bse , LamDim );
  else
   Oracle->SetLamBase( NULL , NumVar );

  BHasChgd = false;
  } */

 if( LHasChgd ) {  // if Lambda has changed, pass the new one to the Oracle
  for( Index i = 0 ; i < NumVar ; i++ )
   LamVcblr[ i ]->set_value( Lambda1.coeff(i) );
  LHasChgd = false;
  for( auto c05 : v_c05f )
   c05->compute( true );
  if( linf )
   linf->compute( true );
  }

 // call Fi() - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiLambda1[ 0 ] = linf->get_value();  // the 0-th component

 for( Index k = 0 ; k++ < NrFi ; )
  // now all other components, one by one - - - - - - - - - - - - - - - - - -
  if( IsEasy.size() && IsEasy[ k ] )     // an easy component
   FiLambda1[ 0 ] += FiLambda1[ k ];  // nothing to do: Fi[ k ] is known already
  else {                          // an hard component
   /* FiStatus[ k ] == kFiStop it means that the function value of component
      k still needs to be computed "with enough accuracy" (this is set in
      FormLambda1() when the function has not been computed at all); here we
      make sure to compute again only the components that have not been
      "solved accurately enough yet" (at the first call after a FormLambda1(),
      all of them). */

   if( FiStatus[ k ] == kStopIter ) {
    FiLambda1[ k ] = v_c05f[ k ]->get_value();
    //FiStatus[ k ] = Oracle->GetFiStatus( k );   // ?? da aggiustare ??
    if( FiStatus[ k ] == kError ) {
     BLOG( 1 , " ~ Error in the FiOracle" << std::endl );
     return( false );
     }
    }

   if( FiLambda1[0] < Inf<double>() ) {
    if( FiLambda1[ k ] >= Inf<double>() )
     FiLambda1[0] = Inf<double>();
    else
     FiLambda1[0] += FiLambda1[ k ];
    }
   }

 if( FiLambda1[0] == Inf<double>() )    // Fi() is not defined in Lambda1
  DeltaFi = Inf<double>();
 else
  DeltaFi = RfrncFi[0] - FiLambda1[0];

 SSDone = false;

 if( FiLambda1[0] == - Inf<double>() )  // error in computing Fi()
  return( true );

 FiEvaltns++;

 // update FiBest, if necessary - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLambda1[0] < FiBest[0] ) {
  FiBest = FiLambda1; //?? va bene?? <--VectAssign( FiBest , FiLambda1 , NrFi + 1 );
  if( KpBstL ) {
  /*if( LamDim < NumVar )
    VectAssignB( LmbdBst , Lambda1 , Lam1Bse , NumVar );
   else */
   LmbdBst = Lambda1; // VectAssign( LmbdBst , Lambda1 , NumVar );
   }

  }

 // call (possibly many times) Gi() - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool dChanges = false;
 Index wFi = 0;  // which component of Fi the new item shall belong to; it
                 // *can* be 0, as the 0-th component (if HpINF) may
                 // generate global constraints

 for( Index Ftchd = 0 ; Ftchd < aBP3 ; ) {
  // check if the Oracle has any new items- - - - - - - - - - - - - - - - - -

  if( ! FindNextSG( wFi ) )
   break;

  // check if aggregation has to be performed - - - - - - - - - - - - - - - -
  // doing this now could occasionally result in useless aggregations, but it
  // avoids complications in the interface of MPSolver (inserting some
  // Z[ wFi ] while inserting the new item)

  Index wh = BStrategy( wFi );

  // get the space for the item from the MPSolver - - - - - - - - - - - - - -

  SparseVector G1;
  if( linf )
   linf->get_linearization_coefficients( G1 );

  // SgRow G1 = Master->GetItem( wFi );

  // fetch the item from the Oracle - - - - - - - - - - - - - - - - - - - - -

  v_c05f[ wFi ]->get_linearization_coefficients( G1 );

  cIndex_Set SGBse = (Index *)G1.innerNonZeroPtr();
  cIndex SGBDm = G1.nonZeros();

  HpNum Alfa1k = v_c05f[ wFi ]->get_linearization_constant();

  GiEvaltns++;

  // pass the base to the MP Solver - - - - - - - - - - - - - - - - - - - - -

  Master->SetItemBse( SGBse , SGBDm );

  // calculate ScPr1k and Alfa1k- - - - - - - - - - - - - - - - - - - - - - -

  Index cp;
  HpNum ScPr1k;

  if( FiLambda1[ wFi ] == Inf<double>() )  // it is a constraint
   cp = Master->CheckCnst( Alfa1k , ScPr1k , Lambda.valuePtr() );
  else                             // it is a subgradient
   cp = Master->CheckSubG( FiLambda1[ wFi ] - RfrncFi[ wFi ] ,
                           t , Alfa1k , ScPr1k );

  // check if the item changes the solution of the MP - - - - - - - - - - - -

  if( FiLambda[0] < Inf<double>() )  // ... but only if Fi( Lambda ) is defined
   dChanges |= Master->ChangesMPSol();
  else
   dChanges = true;               // any first item changes the solution

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
     dChanges = true;          // ensure that the outer Fi-cycle ends
     }
    else
     BLOG( 1 , std::endl << " WARNING: No space in the bundle" << std::endl );

    break;                     // the cycle ends
    }

   Master->SetItem( wh );      // insert the item in the MP Solver
   OOBase[ wh ] = -1;          // ensure it won't be touched again this round

   #if( LOG_BND )
    if( LogVerb ) {
     if( FiLambda1[ wFi ] == Inf<double>() )
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

   if( FiLambda1[ wFi ] < Inf<double>() ) {  // it is a subgradient
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
  }  // end of item-collecting loop - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // compute *Alfa1 and *ScPr1 - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Alfa1[0] = 0;
 ScPr1[0] = Master->ReadGid();

 for( wFi = 0 ; wFi++ < NrFi ; )
  if( whisG1[ wFi ] < InINF ) {
   if( Alfa1[ wFi ] == Inf<double>() )
    Alfa1[ wFi ] = (Master->ReadLinErr())[ whisG1[ wFi ] ];

   Alfa1[0] += Alfa1[ wFi ];

   if( ScPr1[ wFi ] == Inf<double>() )
    ScPr1[ wFi ] = Master->ReadGid( whisG1[ wFi ] );

   ScPr1[0] += ScPr1[ wFi ];
   }
  else
   Alfa1[ wFi ] = ScPr1[ wFi ] = 0;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( dChanges );

 }  // end( BundleSolver::FiAndGi() )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::GotoLambda1( void )
{
 // compute the DeltaFi[] vector- - - - - - - - - - - - - - - - - - - - - - -

 // VectDiff( FiLambda , FiLambda1 , RfrncFi , NrFi + 1 ); ??
 for( Index i = 0 ; i <= NrFi ; i++ )
  FiLambda[ i ] = FiLambda[ i ]	- RfrncFi[ i ];


 // do the move - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Lambda.swap(Lambda1);
 FiLambda.swap(FiLambda1);
 RfrncFi = FiLambda; //???

 // change the current point in the MP Solver - - - - - - - - - - - - - - - -

 Master->ChangeCurrPoint( t , FiLambda1.data() );

 // signal that Alfa1[] is not reliable - - - - - - - - - - - - - - - - - - -

 Alfa1.assign( NrFi + 1 , Inf<double>() );

 // check the new linearization errors- - - - - - - - - - - - - - - - - - - -

 CheckAlfa( true );

 // signal that the latest Fi() point is current- - - - - - - - - - - - - - -

 SSDone = true;

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
 /* ???
 cHpNum LwrBnd = Oracle->GetLowerBound();
 if( LwrBnd != LowerBound[ 0 ] ) {
  LowerBound[ 0 ] = LwrBnd;
  LBHasChgd = true;
  TrueLB = ( LwrBnd > Oracle->GetMinusInfinity() );
  }

 for( Index k = 0 ; k++ < NrFi ; ) {
  if( IsEasy && IsEasy[ k ] )  // skip easy components
   continue;

  cHpNum LwrBndk = Oracle->GetLowerBound( k );
  if( LwrBndk != LowerBound[ k ] ) {
   LowerBound[ k ] = LwrBndk;
   LBHasChgd = true;
   }
  }
  */
 } // end( BundleSolver::UpdtLowerBound )  - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::Log1( void )
{
 #if( LOG_BND )
  if( LogVerb > 1 ) {
   *f_log << std::endl << "{" << SCalls << "-" << ParIter << "-"
	   << Master->MaxName() - FreDim << "-" << MBDim << "} t = " << t
	   << " ~ D*_1( z* ) = " << Master->ReadDStart( 1 )
	   << " ~ Sigma = " << Sigma << std::endl << "           ";
  // if( PPar2 )
  //  *f_log << "LamDim = " << LamDim << " ~ ";

   *f_log <<  " Fi = ";

   if( FiLambda[0] == Inf<double>() )
    *f_log << " - INF";
   else
    *f_log <<  - FiLambda[0] << " ~ eU = " << EpsU;

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

   if( LowerBound[0] > - Inf<double>() )
    *f_log << "UB = " << - LowerBound[0] << " ~ ";

   *f_log << "Fi1 = ";

   if( FiLambda1[0] == - Inf<double>() )
    *f_log << "+ INF => STOP." << std::endl;
   else
    if( FiLambda1[0] == Inf<double>() )
     *f_log << " - INF" << std::endl;
    else
     *f_log << - FiLambda1[0] << " ~ Alfa1 = " << Alfa1[0]
	     << " ~ Gi1xd = " << - ScPr1[0] << std::endl;
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

 Master->SetDim( BPar2 , FakeFi , PPar2 ? true : false );

 Master->SetPar( MPSolver::kOptEps , MPEOpt );
 Master->SetPar( MPSolver::kFsbEps , MPEFsb );

 // insert the constant subgradient of the 0-th component - - - - - - - - - -

 SparseVector SubG;
 if( linf )
  linf->get_linearization_coefficients( SubG );

 // Master->GetItem( 0 ) ?? copy data di SubG in Master->GetItem( 0 ) ???

 Master->SetItemBse( (Index *)SubG.innerNonZeroPtr() , SubG.nonZeros() ); //???
 Master->SetItem( InINF );   //?? indice e vettore resta in memoria ??

 tHasChgd = LBHasChgd = true;

 }  // end( BundleSolver::InitMP( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::FindNextSG( Index &wFi )
{
 // checks if the Oracle is willing to give out any subgradient of some
 // component, beginning from the current value of wFi and going round-robin

 bool GiSet = false;
 for( Index i = NrFi + 1 ; i-- ; wFi = ( wFi < NrFi ? wFi + 1 : 0 ) ) {
  if( wFi ) {
   if( IsEasy.size() && IsEasy[ wFi ] )
    continue;
   }
  else
   if( FiLambda1[0] < Inf<double>() )    // Fi[ 0 ]( Lambda1 ) is defined
    continue;                         // Fi[ 0 ] cannot give anything new

  if( GiSet = v_c05f[ wFi ]->compute_new_linearization() ) // Oracle->NewGi( wFi ) )
   break;
  }

 return( GiSet );
 } // end( BundleSolver::FindNextSG() )  - - - - - - - - - - - - - - - - - - -

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

 if( FreDim ) {                               // there are deleted items
  // wh = HeapDel( FreList , --FreDim );         // pick one
  // ?? da aggiustare ??
  wh = FreList.front();
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

 // tell the Oracle what is going to happen - - - - - - - - - - - - - - - - -

 LinearCombination coefficients;

 v_c05f[ wFi ]->store_combination_of_linearizations( coefficients , whr );

 MBDm = coefficients.size();

 // ?? settare Mlt, MBse, MBDm
 //Oracle->Aggregate( Mlt , MBse , MBDm , whr );

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
  return( t * ( ( DeltaFi + Alfa1[0] ) / ( 2 * Alfa1[0] ) ) );
 } // end( BundleSolver::Heuristic1() ) - -  - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum BundleSolver::Heuristic2( void )
{
 if( std::abs( vStar - DeltaFi ) < Eps<double>() )
  return( tMaior );
 else
  return( t * ( vStar / ( 2 * ( vStar - DeltaFi ) ) ) );
 } // end( BundleSolver::Heuristic2() ) - -  - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::DoSS( void )
{
 if( vStar == - Inf<double>() )
  return( false );

 return( DeltaFi >= std::abs( m1 ) * Deltav );
 } // end( BundleSolver::DoSS() )  - - - - - - - - - - - - - - - - - - - - - -

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

 if( linf ) {
  int GPMaxSz = linf->get_int_par( C05Function::intGPMaxSz );
  for( Index i = 0 ; i < GPMaxSz ; i++ )
   linf->delete_linearization( i );
  }

 FreDim = 0;

 }  // end( BundleSolver::RemoveItems )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::guts_of_destructor( void )
{
 // LmbdBst.clear(); ?? clear per sparsevector ??
 // Lambda1.clear();
 // Lambda.clear();

 DeltaAlfa.clear();
 Alfa1.clear();
 ScPr1.clear();
 whisG1.clear();

 LowerBound.clear();
 FiStatus.clear();

 FiLambda1.clear();
 RfrncFi.clear();
 FiBest.clear();
 FiLambda.clear();

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

  // reset the precision of Fi computations
  EpsCurr = std::abs( EInit );
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
  FiLambda[0] = Inf<double>();

 }  // end( BundleSolver::ReSetAlg ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::Delete( cIndex i )
{
 if( Master ) {
  // check if this item was the "representative" for its component- - - - - -

  cIndex wFi = Master->WComponent( i );

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
   if( FiLambda[0] < Inf<double>() )
    aBP3 = ( BPar5 > 0 ? aBP4 : tBP3 ) +
           Index( BPar5 / std::log10( EpsU / RelAcc ) );
    break;
  case( 3 ):
   if( FiLambda[0] < Inf<double>() )
    aBP3 = ( BPar5 > 0 ? aBP4 : tBP3 ) +
                         Index( BPar5 / std::sqrt( EpsU / RelAcc ) );
   break;
  case( 2 ):
   if( FiLambda[0] < Inf<double>() )
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
 double FiL = FiLambda[0];

 if( FiL == Inf<double>() || vStar == -Inf<double>() )
  return( false );
 else {
  if( FiL < 0 ) FiL = - FiL;
  if( FiL < 1 ) FiL = 1;

  if( eps <= 0 )
   eps = RelAcc;

  // aggiungere alro test ??

  return( DSTS + Sigma <= eps * FiL );
  }
 } // end( BundleSolver::IsOptimal() ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool BundleSolver::CheckAlfa( const bool All )
{

 // ?? come bisgona comportarsi con alfa negativi ??

 }  // end( CheckAlfa )

/*--------------------------------------------------------------------------*/

void BundleSolver::StrongCheckAlfa( void )
{

 }  // end( StrongCheckAlfa )

/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*----------------------- End File BundleSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
