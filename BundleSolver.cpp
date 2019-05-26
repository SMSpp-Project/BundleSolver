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
#include "FakeFiOracle.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

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
		     { "intBPar2" , BundleSolver::intBPar2 } ,
		     { "intBPar6" , BundleSolver::intBPar6 } ,
		     { "intEStps" , BundleSolver::intEStps } ,
		     { "intMnSSC" , BundleSolver::intMnSSC } ,
		     { "intMnNSC" , BundleSolver::intMnNSC } ,
		     { "inttSPar1", BundleSolver::inttSPar1 } ,
		     { "intPPar1" , BundleSolver::intPPar1 } ,
		     { "intPPar2" , BundleSolver::intPPar2 } ,
			 { "intPPar3" , BundleSolver::intPPar3 } };

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
			100 ,  // intBPar2
			  0 ,  // intBPar6
			  0 ,  // intEStps
			  0 ,  // intMnSSC
			  0 ,  // intMnNSC
			  0 ,  // intSPar1
			 30 ,  // intPPar1
			 10 ,  // intPPar2
			  5    // intPPar3
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
static const double HpINF = SMSpp_di_unipi_it::Inf<double>();

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::set_Block( Block * block )
{
 Solver::set_Block( block );  // attach to the new Block

 /* Two types of block can be handled by the BundleSolver:

     1. Only one single non-smooth function
     2. a sum of non-smooth functions

    The algorithm here developed aims at solving non-constrained non-smooth
    optimization. The block can have box constraints at the most.
    It is expected the block to have in the first case a FRealObjective
    whose the function is a C05Function one and having no children, while in
    the second case a FRealObjective one whose the function is a LinearFunction
    and having as many sub-blocks as the number of components. In the latter case,
    each sub-block must not contain any Variable or Constraint. */

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

  }
 else {

  // the objective function of each block must be a LinearFunction - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto obj = boost::any_cast<FRealObjective *>( f_Block->get_objective() );
  if( obj == nullptr )
   throw( std::logic_error( "the objective is not a real function" ) );

  lf = dynamic_cast<LinearFunction *>( (obj)->get_function() );
  if( lf == nullptr )
   throw( std::logic_error( "the objective is not a LinearFunction" ) );

  for( auto & sb : f_Block->get_nested_Blocks() ) { // for each sub-block

   // the objective function of each sub-block must be a C05Function - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto obj = boost::any_cast<FRealObjective *>( sb->get_objective() );
   if( obj == nullptr )
    throw( std::logic_error( "the objective is not a real function" ) );

   auto c05f = dynamic_cast<C05Function *>( (obj)->get_function() );
   if( c05f == nullptr )
	throw( std::logic_error( "the objective is not a C05Function" ) );
   v_c05f.push_back( c05f );

   // nephew are not allowed - - - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( sb->get_nested_Blocks().size() )
	throw( std::logic_error( "nephew are not allowed" ) );

   }
  } // end decomposed case - - - - - - - - - - - - - - - - - - - - - - - - - -

 // construct a FakeFiOracle to handle the MPSolver, which has to
 // interface with a FiOracle object - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // FakeFiOracle * Fi = new FakeFiOracle( c05f );

 }  // end( BundleSolver::set_Block )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::SetMPSolver( MPSolver *MPS )
{

 } // end( BundleSolver::SetMPSolver )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

int BundleSolver::compute( bool changedvars )
{
 // basic sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if( ! Master )
 // throw( std::logic_error( "Master not set ye" ) );

 if( v_c05f.empty() )
  throw( std::logic_error( "C05Function not set yet" ) );

 // initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if( NDOt )  ?? tempo
 // NDOt->Start();

 Result = kOK;
 SCalls++;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 do {
  // construct the direction d- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // FormD();

 } while( ( ! MaxIter ) || ( ParIter < MaxIter ) );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if( MaxIter && ( ParIter >= MaxIter ) && ( ! Result ) )
 Result = kStopIter;

// if( NDOt )   tempo
// NDOt->Stop();

return( Result );

} // end( BundleSolver::compute ) - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::get_dual_solution( Configuration *solc ) {

 for( auto & zel : zA )
  v_c05f[ zel.first ]->set_important_linearization( std::move(zel.second) , zel.first );
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
  case( intBPar2 ):
   return( BPar2 );
   break;
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
  case( intPPar1 ):
   return( PPar1 );
   break;
  case( intPPar2 ):
   return( intPPar2 );
   break;
  case( intPPar3 ):
   return( intPPar3 );
   break;
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
 }

/*--------------------------------------------------------------------------*/

const std::string & BundleSolver::dbl_par_idx2str( const idx_type idx )
 const
{
 if( ( idx >= dbltStar ) && ( idx < dblLastBndSlvPar ) )
  return( dbl_pars_str[ idx - intBPar1 ] );
 else
  return( CDASolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File BundleSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
