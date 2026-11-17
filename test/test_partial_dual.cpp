/*--------------------------------------------------------------------------*/
/*--------------------- File test_partial_dual.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test of the part of the dual solution that BundleSolver::
 * get_dual_solution() writes, i.e., which components get their "important
 * linearization", for each format of the Configuration it is given.
 *
 * The problem is
 * \f[
 *   \max_y \; h_0( y ) + h_1( y ) \; , \quad
 *   h_i( y ) = \min \{ y - a_i \, , \, a_i - y \} \; ,
 * \f]
 * with \f$ a_0 = -1 \f$ and \f$ a_1 = 1 \f$, two hard components, each a
 * concave PolyhedralFunction in the Objective of a sub-Block of the root,
 * an AbstractBlock: h_0 in a TwigBlock and h_1 in a LeafBlock, two
 * AbstractBlock with a classname() of their own defined here. Once solved,
 * the important linearization of both is emptied before each call, which
 * then has to give one back to the components it is asked for and only to
 * them:
 *
 * - { 0 } and { }, the positional format, which with no easy component at
 *   all writes the hard ones;
 *
 * - { LeafBlock : nullptr } and { TwigBlock : nullptr }, the "meta" format,
 *   a component chosen by the classname() of the Block its Objective is in;
 *
 * - { AbstractBlock : nullptr } and { * : nullptr }, the Block of the Solver
 *   chosen, hence all the components;
 *
 * - { NoSuchBlock : nullptr }, none of them;
 *
 * - nullptr, all of them.
 *
 * The configuration of BundleSolver is read from BSPar.txt, that of its
 * master problem from MPBCfg.txt.
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

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "ColVariable.h"

#include "FRealObjective.h"

#include "PolyhedralFunction.h"

#include "BundleSolver.h"

/*--------------------------------------------------------------------------*/

using namespace std;
using namespace SMSpp_di_unipi_it;

using MultiVector = PolyhedralFunction::MultiVector;
using RealVector = PolyhedralFunction::RealVector;
using Meta = SimpleConfiguration< std::map< std::string , Configuration * > >;

#define RED( x )   "\x1B[31m" #x "\033[0m"
#define GREEN( x ) "\x1B[32m" #x "\033[0m"

const double BND = 1e+6;     ///< global upper bound of each h_i()

static int failed = 0;

/*--------------------------------------------------------------------------*/
/// an AbstractBlock with a classname() of its own, for the "meta" format

class TwigBlock : public AbstractBlock {
 public:
 explicit TwigBlock( Block * father = nullptr ) : AbstractBlock( father ) {}
 SMSpp_insert_in_factory_h;
 };

SMSpp_insert_in_factory_cpp_1( TwigBlock );

/// another one, for the second component

class LeafBlock : public AbstractBlock {
 public:
 explicit LeafBlock( Block * father = nullptr ) : AbstractBlock( father ) {}
 SMSpp_insert_in_factory_h;
 };

SMSpp_insert_in_factory_cpp_1( LeafBlock );

/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 cout << ( ok ? GREEN( ok ) : RED( KO ) ) << "   " << what << endl;
 if( ! ok )
  ++failed;
 }

/*--------------------------------------------------------------------------*/
/// h( y ) = min { y - a , a - y } in the Objective of block

static PolyhedralFunction * add_peak( Block * block , ColVariable * y ,
				      double a )
{
 auto PF = new PolyhedralFunction( PolyhedralFunction::VarVector{ y } ,
                                   MultiVector{ { 1 } , { -1 } } ,
                                   RealVector{ - a , a } , BND , false );
 auto obj = new FRealObjective( block , PF );
 obj->set_sense( Objective::eMax , eNoMod );
 block->set_objective( obj , eNoMod );
 return( PF );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 auto root = new AbstractBlock();
 auto y = new ColVariable();
 y->is_positive( false , eNoMod );
 y->is_unitary( false , eNoMod );
 root->add_static_variable( *y , "y" );

 auto twig = new TwigBlock( root );
 auto leaf = new LeafBlock( root );
 root->add_nested_Block( twig );
 root->add_nested_Block( leaf );
 std::vector< PolyhedralFunction * > h = { add_peak( twig , y , -1 ) ,
                                           add_peak( leaf , y , 1 ) };

 auto bsc = dynamic_cast< BlockSolverConfig * >(
                                  Configuration::deserialize( "BSPar.txt" ) );
 if( ! bsc ) {
  cerr << "Error: BSPar.txt does not hold a BlockSolverConfig" << endl;
  return( 1 );
  }
 bsc->apply( root );

 auto slvr = dynamic_cast< BundleSolver * >(
                                    root->get_registered_solvers().front() );
 const int status = slvr->compute( false );
 check( status == Solver::kOK , "solved, status " + std::to_string( status ) );

 // the components whose important linearization the call gives back, as
 // a string of their indices
 auto written = [ & ]( Configuration * solc ) {
  for( auto f : h )
   f->set_important_linearization( C05Function::LinearCombination() );
  slvr->get_dual_solution( solc );
  std::string w;
  for( std::size_t i = 0 ; i < h.size() ; ++i )
   if( ! h[ i ]->get_important_linearization_coefficients().empty() )
    w += std::to_string( i );
  return( w );
  };

 if( status == Solver::kOK ) {
  SimpleConfiguration< std::vector< int > > first( std::vector< int >{ 0 } );
  SimpleConfiguration< std::vector< int > > none( std::vector< int >{} );
  check( written( & first ) == "0" , "positional { 0 }" );
  check( written( & none ) == "" , "positional { }" );

  const std::vector< std::pair< std::string , std::string > > metas = {
   { "LeafBlock" , "1" } , { "TwigBlock" , "0" } ,
   { "AbstractBlock" , "01" } , { "*" , "01" } , { "NoSuchBlock" , "" } };
  for( const auto & [ key , w ] : metas ) {
   Meta solc( { { key , nullptr } } );
   check( written( & solc ) == w , "meta { " + key + " : nullptr }" );
   }

  check( written( nullptr ) == "01" , "nullptr" );
  }

 bsc->clear();
 bsc->apply( root );
 delete bsc;
 delete root;

 if( failed )
  cout << failed << " checks failed" << endl;
 else
  cout << "All tests passed!!" << endl;
 return( failed ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------- End File test_partial_dual.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
