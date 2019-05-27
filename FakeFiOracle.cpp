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
 SGBse1 = new Index[ NumVar + 1 ];
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void FakeFiOracle::SetMaxName( cIndex MxNme ) {
 int MaxName_ = MxNme / bslv->v_c05f.size();

 for( auto c05_it = bslv->v_c05f.begin() ; c05_it != bslv->v_c05f.end() ; ++c05_it )
  if( c05_it == bslv->v_c05f.begin() )
   (*c05_it)->set_par( C05Function::intGPMaxSz , MaxName_ + int(MxNme % bslv->v_c05f.size()) );
  else
   (*c05_it)->set_par( C05Function::intGPMaxSz , MaxName_ );
 MaxName = MxNme;
 } // end ( FakeFiOracle::SetMaxName() ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetNumVar( void ) const {
 return( bslv->NumVar );
 } // end ( FakeFiOracle::GetNumVar() )  - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetNrFi( void ) const {
 return( (bslv->lf == nullptr)? bslv->v_c05f.size() : bslv->v_c05f.size() + 1 );
 } // end ( FakeFiOracle::GetNrFi() )  - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetMaxName( void ) const {
 return( MaxName );
 } // end ( FakeFiOracle::GetMaxName() ) - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/

Index FakeFiOracle::GetGi( SgRow SubG , cIndex_Set &SGBse ,
			cIndex Name , cIndex strt , Index stp  ) {

 // c05f->get_linearization_coefficients( SubG , Name, {} , strt , stp );
 SGBse = nullptr;
 return( stp - strt );
 }

/*--------------------------------------------------------------------------*/

HpNum FakeFiOracle::GetVal( cIndex Name )
{
 if( Name < MaxName )
  throw( NDOException( "GetVal: past information is not recorded" ) );

 return( 0 );
 }


/*--------------------------------------------------------------------------*/
/*----------------------- End File FakeFiOracle.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
