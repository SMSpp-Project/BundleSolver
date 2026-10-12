/*--------------------------------------------------------------------------*/
/*----------------------- File test_mp_recovery.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test of the recovery of BundleSolver from a failure of the master problem
 * right after the new items have entered the bundle: these are not
 * removable at this iteration, and when they are all the bundle has (as at
 * the first master problem of a call) the recovery has to remove them as
 * well, rather than declaring the failure unrecoverable.
 *
 * The problem is
 * \f[
 *   \min_x \{ \, \sum_{k = 1}^K \max_{i} \{ a_{ki}^\top x + c_{ki} \} \, :
 *            \, -1 \leq x \leq 1 \, \} \; ,
 * \f]
 * with K components, solved by the two BundleSolver of BSPar-recovery.txt:
 * the first has its master problem given to TransientMPSolver, which fails
 * at its first FAIL calls and then has all of them solved by HiGHS [see
 * MPBCfg-transient.txt], the second has HiGHS for the master problem from
 * the start. At the first master problem the bundle only has the
 * linearizations of the first evaluation, so the first has to remove some
 * of them and go on: it has to end in kOK with the optimal value of the
 * second.
 *
 * The test is skipped (exit status 77) if HiGHSMILPSolver is not in the
 * build.
 *
 *       ./mp_recovery_test [seed]
 *       seed: pseudo-random generator seed of the instance [0]
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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "CDASolver.h"

#include "ColVariable.h"

#include "FRealObjective.h"

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

const Index N = 5;           ///< number of variables
const Index M = 5;           ///< number of pieces of each component
const Index K = 10;          ///< number of components
const double BND = 1e+6;     ///< bound of each component
const long FAIL = 3;         ///< master problems that fail at the start

static int failed = 0;

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS TransientMPSolver -------------------------*/
/*--------------------------------------------------------------------------*/
/// the Solver of the master problem: the first FAIL calls fail, the others
/// are answered by the Solver registered after it [see
/// MPBCfg-transient.txt]

class TransientMPSolver : public CDASolver {
 public:
 TransientMPSolver( void ) : CDASolver() {}
 ~TransientMPSolver() override = default;

 int compute( bool changedvars = true ) override {
  if( ++f_calls <= FAIL ) {
   ++Failures;
   return( kError );
   }
  return( real()->compute( changedvars ) );
  }

 bool has_var_solution( void ) override {
  return( real()->has_var_solution() );
  }
 bool is_var_feasible( void ) override {
  return( real()->is_var_feasible() );
  }
 void get_var_solution( Configuration * solc = nullptr ) override {
  real()->get_var_solution( solc );
  }
 bool has_dual_solution( void ) override {
  return( real()->has_dual_solution() );
  }
 bool is_dual_feasible( void ) override {
  return( real()->is_dual_feasible() );
  }
 void get_dual_solution( Configuration * solc = nullptr ) override {
  real()->get_dual_solution( solc );
  }
 OFValue get_var_value( void ) override { return( real()->get_var_value() ); }
 OFValue get_lb( void ) override { return( real()->get_lb() ); }
 OFValue get_ub( void ) override { return( real()->get_ub() ); }

 void add_Modification( sp_Mod & mod ) override {}

 /// how many calls have failed, in all the instances
 static inline long Failures = 0;

 SMSpp_insert_in_factory_h;

 private:

 CDASolver * real( void ) {
  return( static_cast< CDASolver * >(
	      *std::next( f_Block->get_registered_solvers().begin() ) ) );
  }

 long f_calls = 0;
 };

SMSpp_insert_in_factory_cpp_0( TransientMPSolver );

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 cout << ( ok ? GREEN( ok ) : RED( KO ) ) << "   " << what << endl;
 if( ! ok )
  ++failed;
 }

/*--------------------------------------------------------------------------*/

/// the instance: K components of M pieces over N variables in [ -1 , 1 ]

static AbstractBlock * make_block( unsigned seed )
{
 std::mt19937 rg( seed );
 std::uniform_real_distribution<> dis( -1.0 , 1.0 );

 auto root = new AbstractBlock();

 auto x = new std::vector< ColVariable >( N );
 root->add_static_variable( *x , "x" );

 auto box = new std::vector< BoxConstraint >( N );
 for( Index i = 0 ; i < N ; ++i ) {
  ( *box )[ i ].set_variable( &(*x)[ i ] );
  ( *box )[ i ].set_lhs( -1 );
  ( *box )[ i ].set_rhs( 1 );
  }
 root->add_static_constraint( *box , "box" );

 for( Index k = 0 ; k < K ; ++k ) {
  MultiVector A( M , RealVector( N ) );
  for( auto & Ai : A )
   for( auto & aij : Ai )
    aij = 10 * dis( rg );
  RealVector c( M );
  for( auto & ci : c )
   ci = 10 * dis( rg );

  auto sub = new AbstractBlock( root );
  PolyhedralFunction::VarVector vars( N );
  for( Index i = 0 ; i < N ; ++i )
   vars[ i ] = &(*x)[ i ];

  auto PF = new PolyhedralFunction( std::move( vars ) , std::move( A ) ,
                                    std::move( c ) , - BND , true );
  auto obj = new FRealObjective( sub , PF );
  obj->set_sense( Objective::eMin , eNoMod );
  sub->set_objective( obj , eNoMod );
  root->add_nested_Block( sub );
  }

 return( root );
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 unsigned seed = 0;
 if( argc > 1 )
  seed = std::stoul( argv[ 1 ] );

 // the master problems are HiGHS's [see MPBCfg-transient.txt]
 if( ! Solver::has_Solver( "HiGHSMILPSolver" ) ) {
  cout << "HiGHSMILPSolver is not in this build: nothing to do" << endl;
  return( 77 );
  }

 auto root = make_block( seed );
 auto bsc = dynamic_cast< BlockSolverConfig * >(
                      Configuration::deserialize( "BSPar-recovery.txt" ) );
 if( ! bsc ) {
  cerr << "Error: BSPar-recovery.txt does not hold a BlockSolverConfig"
       << endl;
  return( 1 );
  }
 bsc->apply( root );

 const auto & slvs = root->get_registered_solvers();
 auto failing = dynamic_cast< BundleSolver * >( *slvs.begin() );
 auto plain = dynamic_cast< BundleSolver * >( *std::next( slvs.begin() ) );

 const auto st0 = failing->compute( false );
 const double v0 = failing->get_ub();
 cout << "failing master: status " << st0 << ", "
      << TransientMPSolver::Failures << " failed master problems, value "
      << v0 << endl;

 const auto st1 = plain->compute( false );
 const double v1 = plain->get_ub();
 cout << "plain master: status " << st1 << ", value " << v1 << endl;

 check( TransientMPSolver::Failures == FAIL ,
	"the first master problems fail" );
 check( st0 == Solver::kOK , "kOK after the failures of the master" );
 check( st1 == Solver::kOK , "kOK with the plain master" );
 check( std::abs( v0 - v1 ) <=
	1e-5 * std::max( { 1.0 , std::abs( v0 ) , std::abs( v1 ) } ) ,
	"the same optimal value" );

 bsc->clear();
 bsc->apply( root );
 delete bsc;
 delete root;

 if( failed )
  cout << RED( KO ) << "   " << failed << " checks failed" << endl;
 else
  cout << GREEN( ok ) << "   all checks passed" << endl;

 return( failed ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File test_mp_recovery.cpp -------------------*/
/*--------------------------------------------------------------------------*/
