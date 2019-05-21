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
               "intSPar1" , "intPPar1", "intPPar2" ,
			   "intSPar3" };

// define and initialize here the vector of double parameters names
const std::vector< std::string > BundleSolver::dbl_pars_str =
		     { "dbltStar"  , "dblEInit" , "dblEFnal"   ,
		       "dblEDcrs" , "dblBPar3" ,  "dblBPar4"  ,
		       "dblBPar5"  , "dblm1" , "dblm3" ,
			   "dblmxIncr" ,  "dblmnIncr" ,  "dblmxDecr" ,
			   "dblmmDecr" ,  "dbltMaior" ,  "dbltMinor" ,
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
		     { "intSPar1" , BundleSolver::intSPar1 } ,
		     { "intPPar1" , BundleSolver::intPPar1 } ,
			 { "intPPar2" , BundleSolver::intPPar3 } };

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
			 { "dblmmDecr" , BundleSolver::dblmmDecr } ,
			 { "dbltMaior" , BundleSolver::dbltMaior } ,
			 { "dbltMinor" , BundleSolver::dbltMinor } ,
			 { "dbltInit" , BundleSolver::dbltInit } ,
			 { "dbltSPar2" , BundleSolver::dbltSPar2 } ,
			 { "dblMPEFsb", BundleSolver::dblMPEFsb } ,
		     { "dblMPEOpt"  , BundleSolver::dblMPEOpt  } };

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

 // the objective function of the inner block must be linear - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto obj = boost::any_cast<FRealObjective *>( f_Block->get_objective() );
 if( obj == nullptr )
  throw( std::logic_error( "the objective is not a real function" ) );

 c05f = dynamic_cast<C05Function *>( (obj)->get_function() );
 if( c05f == nullptr )
  throw( std::logic_error( "the objective is not a C05Function" ) );

 // construct a FakeFiOracle to handle the MPSolver, which has to
 // interface with a FiOracle object - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FakeFiOracle * Fi = new FakeFiOracle( c05f );

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

 if( !c05f )
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

const std::string & BundleSolver::int_par_idx2str( const idx_type idx )
 const
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
