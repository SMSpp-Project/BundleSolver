/*--------------------------------------------------------------------------*/
/*------------------- File test_zeroth_quadratic.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Self-contained internal test for the isotropic quadratic "0-th" component
 * of MasterProblemBlock [see MasterProblemBlock::set_zeroth_quadratic()].
 *
 * The component is not solved for here: what is checked is the arithmetic it
 * rests on, i.e. that the primal master problem carries
 * ( rho / 2 ) || x ||^2 by shrinking its own proximal parameter and shifting
 * its linear part,
 * \f[
 *   \frac{\rho}{2} \| \bar{x} + d \|^2 + \frac{1}{2t} \| d \|^2 =
 *   \frac{1}{2t'} \| d \|^2 + \rho \, \bar{x}^\top d + const \; , \qquad
 *   t' = \frac{t}{1 + \rho t} \; ,
 * \f]
 * and that with rho = 0 nothing whatsoever changes, which is the condition
 * the feature was accepted under. The coefficients of the Objective are read
 * back and compared with the closed form, at the first build and after a
 * change of t. The linear part is not here: it is installed into the
 * Objective once it exists, so the rho * x_bar that goes with it belongs to
 * the refresh done before every solve, and is exercised by a solve rather
 * than by this test.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>

#include "DQuadFunction.h"

#include "FRealObjective.h"

#include "MasterProblemBlock.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;

using namespace SMSpp_di_unipi_it;

#define RED( x )   "\x1B[31m" #x "\033[0m"
#define GREEN( x ) "\x1B[32m" #x "\033[0m"

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static int failed = 0;

static void check( bool ok , const std::string & what )
{
 if( ok )
  cout << GREEN( ok ) << "   " << what << endl;
 else {
  cout << RED( KO ) << "   " << what << endl;
  ++failed;
  }
 }

/*--------------------------------------------------------------------------*/

static void check_close( double got , double expected , double tol ,
                         const std::string & what )
{
 const bool ok = std::abs( got - expected ) <=
                 tol * std::max( 1.0 , std::abs( expected ) );
 if( ! ok )
  cout << "     got " << got << " , expected " << expected << endl;
 check( ok , what );
 }

/*--------------------------------------------------------------------------*/

/// the DQuadFunction of the Objective of @p mpb, which must be there

static DQuadFunction * objective_of( MasterProblemBlock * mpb )
{
 auto obj = dynamic_cast< FRealObjective * >( mpb->get_objective() );
 if( ! obj )
  return( nullptr );
 return( dynamic_cast< DQuadFunction * >( obj->get_function() ) );
 }

/*--------------------------------------------------------------------------*/

/// builds a primal, proximal MasterProblemBlock of @p n variables

static MasterProblemBlock * make_mpb( int n , double t ,
                                      const std::vector< double > & b ,
                                      const std::vector< double > & x_bar ,
                                      double rho )
{
 auto mpb = new MasterProblemBlock();
 mpb->configure( true ,    // primal MP
                 10 ,      // max bundle size
                 n ,       // number of variables
                 1 ,       // one hard component
                 {} ,      // no easy component
                 {} ,      // nothing ignored
                 MasterProblemBlock::kProximal );
 mpb->generate_abstract_variables();
 mpb->generate_abstract_constraints();
 mpb->set_t( t );
 mpb->set_linear_part( b );
 mpb->set_x_bar( x_bar );
 if( rho )
  mpb->set_zeroth_quadratic( rho );
 mpb->generate_objective();
 return( mpb );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 const int n = 3;
 const double t = 0.5;
 const double rho = 4.0;
 const std::vector< double > b = { 1.0 , -2.0 , 0.5 };
 const std::vector< double > x_bar = { 2.0 , -1.0 , 3.0 };

 // with no quadratic 0-th component nothing changes - - - - - - - - - - - -

 { auto mpb = make_mpb( n , t , b , x_bar , 0.0 );

   auto dqf = objective_of( mpb );
   check( dqf , "the primal objective is a DQuadFunction" );
   if( dqf ) {
    check_close( dqf->get_quadratic_coefficient( 0 ) , 1.0 / ( 2.0 * t ) ,
                 1e-12 , "rho = 0 leaves the proximal coefficient alone" );
    }
   }

 // with it, the master shrinks t and shifts the linear part - - - - - - - -

 { auto mpb = make_mpb( n , t , b , x_bar , rho );

   auto dqf = objective_of( mpb );
   const double t_prime = t / ( 1.0 + rho * t );

   check_close( dqf->get_quadratic_coefficient( 0 ) ,
                1.0 / ( 2.0 * t_prime ) , 1e-12 ,
                "the proximal coefficient is 1 / ( 2 t' )" );


   // a change of t is followed - - - - - - - - - - - - - - - - - - - - - - -

   const double t2 = 0.1;
   auto mpb2 = make_mpb( n , t2 , b , x_bar , rho );
   const double t2_prime = t2 / ( 1.0 + rho * t2 );
   check_close( objective_of( mpb2 )->get_quadratic_coefficient( 0 ) ,
                1.0 / ( 2.0 * t2_prime ) , 1e-12 ,
                "a smaller t gives a smaller t'" );

   // and so is a change of the stability centre, on which the linear part
   // depends only because of the quadratic component


   // the effective t is bounded by 1 / rho, i.e. a strongly convex objective
   // cannot be destabilized: t' -> 1 / rho as t -> infinity

   auto mpb4 = make_mpb( n , 1e+10 , b , x_bar , rho );
   check_close( objective_of( mpb4 )->get_quadratic_coefficient( 0 ) ,
                rho / 2.0 , 1e-6 ,
                "the effective t is bounded by 1 / rho" );

   }

 // what is not supported says so - - - - - - - - - - - - - - - - - - - - - -

 { auto mpb = new MasterProblemBlock();
   mpb->configure( true , 10 , n , 1 , {} , {} ,
                   MasterProblemBlock::kDoublyStabilized );
   bool thrown = false;
   try { mpb->set_zeroth_quadratic( rho ); }
   catch( const std::logic_error & ) { thrown = true; }
   check( thrown , "the level row cannot carry it, and it says so" );
   }

 cout << ( failed ? RED( Shit happened!! ) : GREEN( All tests passed!! ) )
      << endl;

 return( failed ? 1 : 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------- End File test_zeroth_quadratic.cpp -------------------*/
/*--------------------------------------------------------------------------*/
