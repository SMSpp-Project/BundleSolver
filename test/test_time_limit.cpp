/*--------------------------------------------------------------------------*/
/*----------------------- File test_time_limit.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test of the time limit of BundleSolver [see dblMaxTime] while the master
 * problem fails: the recovery from a failure of the master solves it again
 * after each item it removes from the bundle, which with a large bundle may
 * take much longer than the time that is left, and the time limit has to
 * hold there as well.
 *
 * The problem is
 * \f[
 *   \min_x \{ \, \sum_{k = 1}^K \max_{i} \{ a_{ki}^\top x + c_{ki} \} \, :
 *            \, -1 \leq x \leq 1 \, \} \; ,
 * \f]
 * with K components. The master problem is given to FailingMPSolver,
 * which has the first few of them solved by HiGHS and then takes a while and
 * returns kError at each call [see MPBCfg-fail.txt], when the bundle is
 * large enough for the recovery to take much longer than the time limit.
 * The instance is solved by the two BundleSolver of BSPar-time.txt:
 *
 * - the first, with no time limit, has to end in kLowPrecision at the
 *   stability centre: the recovery removes the items one at a time, those
 *   of the current iteration as well, the master problem with an empty
 *   bundle is the centre itself, which then the noise reduction cannot get
 *   away from, until t is at its maximum;
 *
 * - the second, with dblMaxTime well below the time the recovery takes,
 *   has to end in kStopTime soon after the time limit, with the value of
 *   the stability centre (get_ub()), which is finite.
 *
 * The test is skipped (exit status 77) if HiGHSMILPSolver is not in the
 * build.
 *
 *       ./time_limit_test [seed]
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

#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <thread>
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
const Index K = 60;          ///< number of components
const double BND = 1e+6;     ///< bound of each component
const double WAIT = 0.05;    ///< seconds each master problem takes
const double MAXT = 1;       ///< dblMaxTime of the 2nd BundleSolver
const long GOOD = 5;         ///< master problems solved before the failures

static int failed = 0;

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS FailingMPSolver --------------------------*/
/*--------------------------------------------------------------------------*/
/// the Solver of the master problem: the first GOOD calls are answered by
/// the Solver registered after it [see MPBCfg-fail.txt], the following ones
/// take WAIT seconds and fail

class FailingMPSolver : public CDASolver {
 public:
 FailingMPSolver( void ) : CDASolver() {}
 ~FailingMPSolver() override = default;

 int compute( bool changedvars = true ) override {
  if( ++f_calls > GOOD ) {
   ++Failures;
   std::this_thread::sleep_for( std::chrono::duration< double >( WAIT ) );
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

SMSpp_insert_in_factory_cpp_0( FailingMPSolver );

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

 // the master problems that do not fail are HiGHS's [see MPBCfg-fail.txt]
 if( ! Solver::has_Solver( "HiGHSMILPSolver" ) ) {
  cout << "HiGHSMILPSolver is not in this build: nothing to do" << endl;
  return( 77 );
  }

 auto root = make_block( seed );
 auto bsc = dynamic_cast< BlockSolverConfig * >(
                          Configuration::deserialize( "BSPar-time.txt" ) );
 if( ! bsc ) {
  cerr << "Error: BSPar-time.txt does not hold a BlockSolverConfig" << endl;
  return( 1 );
  }
 bsc->apply( root );

 const auto & slvs = root->get_registered_solvers();
 auto nolimit = dynamic_cast< BundleSolver * >( *slvs.begin() );
 auto limit = dynamic_cast< BundleSolver * >( *std::next( slvs.begin() ) );

 // with no time limit the recovery removes the items one at a time, at WAIT
 // seconds each, until the bundle is empty and the master problem gives the
 // centre, which the noise reduction does not change: the call ends with
 // the centre, as not optimal, after much longer than the time limit
 const auto st0 = nolimit->compute( false );
 const double t0 = nolimit->get_elapsed_time();
 const double v0 = nolimit->get_ub();  // the value at the stability centre
 const long f0 = FailingMPSolver::Failures;
 cout << "no time limit: status " << st0 << " in " << t0 << " s, " << f0
      << " failed master problems" << endl;
 check( st0 == Solver::kLowPrecision , "kLowPrecision with no time limit" );
 check( std::isfinite( v0 ) , "the value of the stability centre with no "
	                      "time limit" );
 check( t0 > MAXT + 10 * WAIT , "the recovery outlasts the time limit" );

 // with the time limit it stops in the recovery, soon after the limit, with
 // the stability centre
 const auto st1 = limit->compute( false );
 const double t1 = limit->get_elapsed_time();
 const double v1 = limit->get_ub();  // the value at the stability centre
 const long f1 = FailingMPSolver::Failures - f0;
 cout << "dblMaxTime " << MAXT << ": status " << st1 << " in " << t1
      << " s, " << f1 << " failed master problems, value " << v1 << endl;
 check( f1 > 0 , "the master fails before the time limit" );
 check( st1 == Solver::kStopTime , "kStopTime with the time limit" );
 check( t1 < MAXT + 10 * WAIT , "stopped soon after the time limit" );
 check( std::isfinite( v1 ) , "the value of the stability centre" );

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
/*------------------------ End File test_time_limit.cpp --------------------*/
/*--------------------------------------------------------------------------*/
