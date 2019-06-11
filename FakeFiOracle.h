/*--------------------------------------------------------------------------*/
/*------------------------- File FakeFiOracle.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the FakeFiOracle class, which implements a "fake",
 * FiOracle whose role is not really that of creating an Object to be passed
 * to the MPsolver within the NDOSolver interface.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef _FakeFiOracle
 #define _FakeFiOracle  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/


#include "FiOracle.h"
#include "C05Function.h"
#include "LinearFunction.h"

#include "BundleSolver.h"

#include "OPTtypes.h"
#include "OPTUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/


#if( OPT_USE_NAMESPACES )
 using namespace NDO_di_unipi_it;
#endif

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

 class BundleSolver;     // forward declaration of class BundleSolver

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS FakeFiOracle  ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

class FakeFiOracle : public FiOracle
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
    @{ */


/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

/** Constructor of the class: takes no arguments, since everything that
    concerns the real evaluation of the function must be done in derived
    classes, which will have their parameters. */

   FakeFiOracle( BundleSolver *solver );

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

   virtual void SetNDOSolver( NDOSolver *NwSlvr = 0 ) override;

/*--------------------------------------------------------------------------*/

   virtual void SetFiLog( ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

   virtual void SetFiTime( const bool TimeIt = true ) override;

/*--------------------------------------------------------------------------*/

   virtual void SetMaxName( cIndex MxNme = 0 ) override;

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

/// get the number of Variable
/** Variable cannot be changed. This means that is used the default
 *  implementation of GetMaxNumVar(). The maximum number of variables is
 *  equal to the current number of variable*/

  virtual Index GetNumVar( void ) const override;

/*--------------------------------------------------------------------------*/

  virtual Index GetNrFi( void ) const override;

/*--------------------------------------------------------------------------*/

  virtual Index GetMaxName( void ) const override;

/*--------------------------------------------------------------------------*/

  virtual HpNum GetMinusInfinity( void ) override;

/*--------------------------------------------------------------------------*/

  virtual Index GetMaxNZ( cIndex wFi = Inf<Index>() ) const override;

/*--------------------------------------------------------------------------*/

  virtual Index GetMaxCNZ( cIndex wFi = Inf<Index>() ) const override;

/*--------------------------------------------------------------------------*/

  virtual bool GetUC( cIndex i ) override;

/*--------------------------------------------------------------------------*/

  virtual LMNum GetUB( cIndex i ) override;

/*--------------------------------------------------------------------------*/

  virtual LMNum GetBndEps( void ) override;

/*--------------------------------------------------------------------------*/

  virtual HpNum GetGlobalLipschitz( cIndex wFi = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

  virtual NDOSolver *GetNDOSolver( void ) override;

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SETTING LAMBDA ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Setting Lambda
   @{ */

  virtual void SetLambda( cLMRow Lmbd = 0 ) override;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  virtual void SetLamBase( cIndex_Set LmbdB = 0 , cIndex LmbdBD = 0 );

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  virtual bool SetPrecision( HpNum Eps );

/*@} -----------------------------------------------------------------------*/
/*------------------------ METHODS FOR COMPUTING Fi() ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Computing Fi()
   @{ */

   virtual HpNum Fi( cIndex wFi = Inf<Index>() );

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/


   virtual bool NewGi( cIndex wFi = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   virtual Index GetGi( SgRow SubG , cIndex_Set &SGBse ,
			cIndex Name = Inf<Index>() ,
			cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   virtual HpNum GetVal( cIndex Name = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   virtual void SetGiName( cIndex Name ) override;

/*@} -----------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading other results
   @{ */


   virtual HpNum GetLowerBound( cIndex wFi = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   virtual FiStatus GetFiStatus( Index wFi = Inf<Index>() ) override;

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
   @{ */

   virtual void Deleted( cIndex i = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   virtual void Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm ,
			   cIndex NwNm );

/*@} -----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~FakeFiOracle()
   {
	GiNameVcblr.clear();
    }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

  BundleSolver *bslv;
  std::vector< std::tuple< Index , Index , bool > >  GiNameVcblr; /* vocabulary
     for handling the items name; this is done to map the item name
     from the FiOracle to that of C05Function.  */

  Index last_c05;

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Standard fields

/*@} -----------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class FakeFiOracle )

/*--------------------------------------------------------------------------*/

 }
/*#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* FakeFiOracle.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File FakeFiOracle.h ------------------------*/
/*--------------------------------------------------------------------------*/
