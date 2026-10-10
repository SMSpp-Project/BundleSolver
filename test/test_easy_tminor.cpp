/*--------------------------------------------------------------------------*/
/*--------------------- File test_easy_tminor.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test of the verdict of the easy components at the first master problem,
 * which BundleSolver solves with an empty bundle and t at dbltMinor.
 *
 * The problem is
 * \f[
 *   \max_y \; L( y ) + h( y ) \; , \quad
 *   L( y ) = \min \{ ( 1 + y ) x \,:\, x \in X \} \; , \quad
 *   h( y ) = \min \{ y - y^* \, , \, y^* - y \} \; ,
 * \f]
 * where L() is a LagBFunction whose inner Block is a linear program, hence
 * an easy component that the master problem carries whole, and h() is a
 * concave PolyhedralFunction, a hard component whose bundle is empty at the
 * first master. The cases are
 *
 * - \f$ X = [ 0 , +\infty ) \f$ and the centre \f$ y = 0 \f$: the first
 *   master has a solution at tMinor, and the optimum is 0 at
 *   \f$ y = y^* \f$;
 *
 * - the same with the centre \f$ y = -2 \f$, outside the domain
 *   \f$ y \geq -1 \f$ of L(): the dual master at tMinor is
 *   \f$ \min \{ -x + ( t / 2 ) x^2 \,:\, x \geq 0 \} \f$, whose minimiser
 *   \f$ 1 / t \f$ is enormous, and with dbltMinor 1e-16 Gurobi
 *   declares it unbounded, i.e., the primal master empty, while the master
 *   solved with the previous t has a solution; BundleSolver has to find the
 *   optimum 0 rather than report the problem infeasible;
 *
 * - X empty (\f$ x \geq 0 \f$ and \f$ x \leq -1 \f$), whence
 *   \f$ L( y ) = +\infty \f$ and the problem is unbounded, which the
 *   master says at the first master and once more after t is restored:
 *   BundleSolver has to return kUnbounded, after solving the master again.
 *
 * The configuration of BundleSolver is read from BSPar-easy.txt, that of
 * its master problem from MPBCfg-easy.txt.
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
#include <sstream>
#include <string>
#include <vector>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "ColVariable.h"

#include "FRealObjective.h"

#include "LagBFunction.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "PolyhedralFunction.h"

#include "BundleSolver.h"

/*--------------------------------------------------------------------------*/

using namespace std;
using namespace SMSpp_di_unipi_it;

using Index = Block::Index;
using MultiVector = PolyhedralFunction::MultiVector;
using RealVector = PolyhedralFunction::RealVector;

#define RED( x )   "\x1B[31m" #x "\033[0m"
#define GREEN( x ) "\x1B[32m" #x "\033[0m"

const double BND = 1e+6;     ///< global upper bound of h()
const double YSTAR = -0.5;   ///< the maximiser of h()

static int failed = 0;

/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 cout << ( ok ? GREEN( ok ) : RED( KO ) ) << "   " << what << endl;
 if( ! ok )
  ++failed;
 }

/*--------------------------------------------------------------------------*/
/// solves the instance, returns the status and writes the value in val

static int solve( bool empty_X , double y0 , double & val ,
                  std::string & log )
{
 auto root = new AbstractBlock();
 auto y = new ColVariable();
 y->is_positive( false , eNoMod );
 y->is_unitary( false , eNoMod );
 root->add_static_variable( *y , "y" );

 // the easy component: L( y ) = min { ( 1 + y ) x : x in X }
 auto inner = new AbstractBlock();
 auto x = new ColVariable();
 x->is_positive( true , eNoMod );     // x >= 0
 inner->add_static_variable( *x , "x" );
 if( empty_X ) {
  auto ub = new BoxConstraint();
  ub->set_variable( x , eNoMod );
  ub->set_lhs( - Inf< double >() , eNoMod );
  ub->set_rhs( -1 , eNoMod );
  inner->add_static_constraint( *ub , "ub" );
  }
 auto iobj = new FRealObjective( inner , new LinearFunction(
                            LinearFunction::v_coeff_pair{ { x , 1 } } ) );
 iobj->set_sense( Objective::eMin , eNoMod );
 inner->set_objective( iobj , eNoMod );

 auto lbf = new LagBFunction( inner );
 LagBFunction::v_dual_pair dp;
 dp.emplace_back( y , new LinearFunction(
                            LinearFunction::v_coeff_pair{ { x , 1 } } ) );
 lbf->set_dual_pairs( std::move( dp ) );

 auto easy = new AbstractBlock( root );
 auto eobj = new FRealObjective( easy , lbf );
 eobj->set_sense( Objective::eMax , eNoMod );
 easy->set_objective( eobj , eNoMod );
 root->add_nested_Block( easy );

 // the hard component: h( y ) = min { y - y^* , y^* - y }
 auto hard = new AbstractBlock( root );
 auto PF = new PolyhedralFunction( PolyhedralFunction::VarVector{ y } ,
                                   MultiVector{ { 1 } , { -1 } } ,
                                   RealVector{ - YSTAR , YSTAR } , BND ,
                                   false );
 auto hobj = new FRealObjective( hard , PF );
 hobj->set_sense( Objective::eMax , eNoMod );
 hard->set_objective( hobj , eNoMod );
 root->add_nested_Block( hard );

 y->set_value( y0 );

 auto bsc = dynamic_cast< BlockSolverConfig * >(
                             Configuration::deserialize( "BSPar-easy.txt" ) );
 if( ! bsc ) {
  cerr << "Error: BSPar-easy.txt does not hold a BlockSolverConfig" << endl;
  exit( 1 );
  }
 bsc->apply( root );

 auto slvr = root->get_registered_solvers().front();
 std::ostringstream os;
 slvr->set_log( &os );
 const int status = slvr->compute( false );
 val = slvr->get_var_value();
 log = os.str();

 bsc->clear();
 bsc->apply( root );
 delete bsc;

 // the LagBFunction owns its inner Block; the Objective of easy only
 // clear()s its Function when the Block goes, hence it is deleted here
 eobj->set_function( nullptr , eNoMod , true );
 delete root;
 return( status );
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 const bool verbose = ( argc > 1 );
 const std::string again = "solving again with t";

 double val;
 std::string log;

 // the central case: the first master has a solution at tMinor
 int st = solve( false , 0 , val , log );
 if( verbose ) cout << log << endl;
 check( ( st == Solver::kOK ) && ( std::abs( val ) <= 1e-6 ) ,
        "centre in the domain: status " + std::to_string( st ) +
        ", value " + std::to_string( val ) );

 // the centre outside the domain of the easy component
 st = solve( false , -2 , val , log );
 if( verbose ) cout << log << endl;
 check( ( st == Solver::kOK ) && ( std::abs( val ) <= 1e-6 ) ,
        "centre outside the domain: status " + std::to_string( st ) +
        ", value " + std::to_string( val ) +
        ( log.find( again ) != std::string::npos ?
          " (master solved again)" : "" ) );

 // the easy region empty: the verdict survives t being restored
 st = solve( true , 0 , val , log );
 if( verbose ) cout << log << endl;
 check( st == Solver::kUnbounded ,
        "easy region empty: status " + std::to_string( st ) );
 check( log.find( again ) != std::string::npos ,
        "easy region empty: master solved again with the previous t" );

 cout << ( failed ? RED( Shit happened!! ) : GREEN( All tests passed!! ) )
      << endl;

 return( failed ? 1 : 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*--------------------- End File test_easy_tminor.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
