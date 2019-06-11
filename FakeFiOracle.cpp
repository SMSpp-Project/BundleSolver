/*--------------------------------------------------------------------------*/
/*------------------------ File FakeFiOracle.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the FakeFiOracle class.
 *
 * \version 0.01
 *
 * \date 19 - 05 - 2010
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

#include "FakeFiOracle.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

FakeFiOracle::FakeFiOracle( BundleSolver *solver ) : FiOracle()
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

void FakeFiOracle::SetNDOSolver( NDOSolver *NwSlvr ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::SetNDOSolver() ) - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void FakeFiOracle::SetFiLog( ostream *outs , const char lvl ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::SetFiLog() ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void FakeFiOracle::SetFiTime( const bool TimeIt ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::SetFiTime() )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void FakeFiOracle::SetMaxName( cIndex MxNme ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::SetMaxName() ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetNumVar( void ) const {
 return( bslv->NumVar );
 } // end ( FakeFiOracle::GetNumVar() )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetNrFi( void ) const {
 return( bslv->v_c05f.size( ) );
 } // end ( FakeFiOracle::GetNrFi() )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetMaxName( void ) const {
 return( bslv->BPar2 );
 } // end ( FakeFiOracle::GetMaxName() ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum FakeFiOracle::GetMinusInfinity( void ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetMinusInfinity() ) - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetMaxNZ( cIndex wFi ) const {
 if( wFi != Inf<Index>() )
  throw( std::logic_error( "GetMaxNZ can be called with wFi = Inf only" ) );
 return( bslv->NumVar );
 } // end ( FakeFiOracle::GetMaxNZ() ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetMaxCNZ( cIndex wFi ) const {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetMaxCNZ() )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

bool FakeFiOracle::GetUC( cIndex i ) {

 double lb_value = bslv->LamVcblr[ i ]->get_lb();
 if( lb_value == -Inf<ColVariable::VarValue>() )
  return( true );

 if( lb_value != ColVariable::VarValue(0) )
  throw( std::logic_error( "any value different from zero is not allowed" ) );

 return( false );
 } // end ( FakeFiOracle::GetUC() )  - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

LMNum FakeFiOracle::GetUB( cIndex i ) {
 return( bslv->LamVcblr[ i ]->get_ub() );
 } // end ( FakeFiOracle::GetUB() )  - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

LMNum FakeFiOracle::GetBndEps(  ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetBndEps() )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum FakeFiOracle::GetGlobalLipschitz( cIndex wFi ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetGlobalLipschitz() )   - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

NDOSolver * FakeFiOracle::GetNDOSolver( void ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::GetNDOSolver() ) - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SETTING LAMBDA ------------------------*/
/*--------------------------------------------------------------------------*/

void FakeFiOracle::SetLambda( cLMRow Lmbd ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

void FakeFiOracle::SetLamBase( cIndex_Set LmbdB  , cIndex LmbdBD ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

bool FakeFiOracle::SetPrecision( HpNum Eps ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR COMPUTING Fi() ----------------------*/
/*--------------------------------------------------------------------------*/

HpNum FakeFiOracle::Fi( cIndex wFi ) {
 throw( std::logic_error( "this method cannot be called" ) );
 } // end ( FakeFiOracle::Fi( ) )  - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/

bool FakeFiOracle::NewGi( cIndex wFi ) {
 if( wFi == 0 )
  throw( std::invalid_argument( "asking for the 0th component" ) );
 last_c05 =  wFi-1;
 return( true );
 } // end ( FakeFiOracle::NewGi( ) ) - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetGi( SgRow SubG , cIndex_Set &SGBse ,
			cIndex Name , cIndex strt , Index stp  ) {

 bslv->v_c05f[ std::get<1>(GiNameVcblr[Name]) ]->get_linearization_coefficients(
	 SubG , std::get<0>(GiNameVcblr[Name]) , {} , strt , stp );

 SGBse = nullptr;
 return( stp - strt );
 } // end ( FakeFiOracle::GetGi( ) )   - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

HpNum FakeFiOracle::GetVal( cIndex Name )
{
 return( bslv->v_c05f[ std::get<1>(GiNameVcblr[Name]) ]->
 		 get_linearization_constant( std::get<0>(GiNameVcblr[Name]) ) );
 } // end ( FakeFiOracle::GetVal( ) )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void FakeFiOracle::SetGiName( cIndex Name )
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

HpNum FakeFiOracle::GetLowerBound( cIndex wFi ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }  // end ( FakeFiOracle::GetLowerBound( ) )  - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

FiOracle::FiStatus FakeFiOracle::GetFiStatus( Index wFi ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }  // end ( FakeFiOracle::GetFiStatus( ) )  - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void FakeFiOracle::Deleted( cIndex i ) {

 bslv->v_c05f[ std::get<1>(GiNameVcblr[i]) ]->
     delete_linearization( std::get<0>(GiNameVcblr[i]) );

 std::get<2>(GiNameVcblr[i]) = true;
 } // end ( FakeFiOracle::Deleted( ) ) - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void FakeFiOracle::Aggregate( cHpRow Mlt , cIndex_Set NmSt ,
		cIndex Dm , cIndex NwNm ) {
 throw( std::logic_error( "this method cannot be called" ) );
 }  // end ( FakeFiOracle::Aggregate( ) )  - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*----------------------- End File FakeFiOracle.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
