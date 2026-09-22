/*--------------------------------------------------------------------------*/
/*--------------------------- File test_kbm.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The large-scale nonsmooth test problems of
 *
 *   N. Karmitsa, A. Bagirov, M.M. Makela "Comparing different nonsmooth
 *   minimization methods and software", Optimization Methods and Software
 *   27(1):131-153, 2012
 *
 * as coded in tnsunc.f, the file of test problems distributed with the
 * limited memory bundle method LMBM of Karmitsa: this tester is linked with
 * that very file, and it takes the starting point from its STARTX and the
 * value and a subgradient from its FUNC, so that BundleSolver solves the
 * same functions from the same points as the Fortran codes. The function is
 * a single component (an OracleFunction, i.e., a C05Function whose oracle is
 * a callback, with a global pool) of an unconstrained Block; the value and
 * the number of evaluations are printed, and if the optimal value is given
 * the test fails if the relative error is larger than tol.
 *
 * Usage: kbm_test problem n BSPar [ fstar [ tol ] ]
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/

#include "AbstractBlock.h"
#include "BlockSolverConfig.h"
#include "C05Function.h"
#include "FRealObjective.h"
#include "LinearFunction.h"
#include "Solver.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <unordered_map>

using namespace SMSpp_di_unipi_it;

// the two routines of tnsunc.f
extern "C" {
 void startx_( int * n , double * x , int * next );
 void func_( int * n , double * x , double * f , double * g , int * next );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS OracleFunction --------------------------*/
/*--------------------------------------------------------------------------*/
/// a C05Function given by a callback, with a global pool
/** The callback returns f( x ) and writes a subgradient g of f at x; the
 * linearization is ( g , alpha ), alpha = f( x ) - g x. The global pool
 * keeps the ( g , alpha ) that are stored and their combinations. */

class OracleFunction : public C05Function {

/*--------------------------------------------------------------------------*/

 using v_col_var = std::vector< ColVariable * >;

 class v_iterator : public ThinVarDepInterface::v_iterator {
  public:
   explicit v_iterator( v_col_var::iterator itr ) : itr_( itr ) {}
   v_iterator * clone( void ) override { return( new v_iterator( itr_ ) ); }
   void operator++( void ) override { ++itr_; }
   reference operator*( void ) const override { return( **itr_ ); }
   pointer operator->( void ) const override { return( *itr_ ); }
   bool operator==( const ThinVarDepInterface::v_iterator & rhs )
    const override {
    auto t = dynamic_cast< const v_iterator * >( & rhs );
    return( t ? itr_ == t->itr_ : false );
    }
   bool operator!=( const ThinVarDepInterface::v_iterator & rhs )
    const override { return( ! ( *this == rhs ) ); }
  private:
   v_col_var::iterator itr_;
  };

 class v_const_iterator : public ThinVarDepInterface::v_const_iterator {
  public:
   explicit v_const_iterator( v_col_var::const_iterator itr ) : itr_( itr ) {}
   v_const_iterator * clone( void ) override {
    return( new v_const_iterator( itr_ ) );
    }
   void operator++( void ) override { ++itr_; }
   reference operator*( void ) const override { return( **itr_ ); }
   pointer operator->( void ) const override { return( *itr_ ); }
   bool operator==( const ThinVarDepInterface::v_const_iterator & rhs )
    const override {
    auto t = dynamic_cast< const v_const_iterator * >( & rhs );
    return( t ? itr_ == t->itr_ : false );
    }
   bool operator!=( const ThinVarDepInterface::v_const_iterator & rhs )
    const override { return( ! ( *this == rhs ) ); }
  private:
   v_col_var::const_iterator itr_;
  };

/*--------------------------------------------------------------------------*/

 public:

 /// the oracle gives f( x ) and writes the gradient in g
 using Oracle = std::function< double( const double * x , double * g ) >;

 OracleFunction( v_col_var && vars , Oracle && oracle )
  : C05Function() , v_vars( std::move( vars ) ) ,
    f_oracle( std::move( oracle ) ) , f_val( Inf< FunctionValue >() ) ,
    f_g( v_vars.size() , 0 ) {
  for( Index i = 0 ; i < v_vars.size() ; ++i )
   f_pos.emplace( v_vars[ i ] , i );
  }

 ~OracleFunction() override = default;

 void clear( void ) override {}

/*--------------------------------------------------------------------------*/

 void set_par( idx_type par , int value ) override {
  if( par == intGPMaxSz ) {
   if( value < 0 )
    throw( std::invalid_argument( "OracleFunction::set_par: negative "
                                  "intGPMaxSz" ) );
   if( Index( value ) < v_pool.size() )
    v_pool.resize( value );
   f_gp_max = value;
   return;
   }
  if( par == intLPMaxSz )
   return;  // there is only one linearization anyway
  C05Function::set_par( par , value );
  }

 using C05Function::set_par;

 int get_int_par( idx_type par ) const override {
  if( par == intGPMaxSz )
   return( f_gp_max );
  return( C05Function::get_int_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 void remove_variable( Index i , ModParam issueMod = eModBlck ) override {
  throw( std::logic_error( "OracleFunction::remove_variable: not "
                           "allowed" ) );
  }

/*--------------------------------------------------------------------------*/

 int compute( bool changedvars = true ) override {
  const Index p = v_vars.size();
  f_x.resize( p );
  for( Index j = 0 ; j < p ; ++j )
   f_x[ j ] = v_vars[ j ]->get_value();
  f_val = f_oracle( f_x.data() , f_g.data() );
  ++f_nevals;
  f_alpha = f_val;
  for( Index j = 0 ; j < p ; ++j )
   f_alpha -= f_g[ j ] * f_x[ j ];
  return( kOK );
  }

 /// the number of calls of compute() so far
 long get_nevals( void ) const { return( f_nevals ); }

 FunctionValue get_value( void ) override { return( f_val ); }

 bool is_convex( void ) override { return( true ); }


/*--------------------------------------------------------------------------*/

 void store_linearization( Index name , ModParam issueMod = eModBlck )
  override {
  if( name >= f_gp_max )
   throw( std::invalid_argument( "OracleFunction::store_linearization: "
                                 "invalid name" ) );
  if( name >= v_pool.size() )
   v_pool.resize( name + 1 );
  v_pool[ name ].g = f_g;
  v_pool[ name ].alpha = f_alpha;
  v_pool[ name ].there = true;
  }

 bool is_linearization_there( Index name ) const override {
  return( ( name < v_pool.size() ) && v_pool[ name ].there );
  }

 void store_combination_of_linearizations(
                    c_LinearCombination & coefficients , Index name ,
                    ModParam issueMod = eModBlck ) override {
  if( name >= f_gp_max )
   throw( std::invalid_argument( "OracleFunction::store_combination_of_"
                                 "linearizations: invalid name" ) );
  std::vector< double > g( v_vars.size() , 0 );
  double alpha = 0;
  for( auto & [ i , c ] : coefficients ) {
   const auto & e = entry( i );
   for( Index j = 0 ; j < g.size() ; ++j )
    g[ j ] += c * e.g[ j ];
   alpha += c * e.alpha;
   }
  if( name >= v_pool.size() )
   v_pool.resize( name + 1 );
  v_pool[ name ].g = std::move( g );
  v_pool[ name ].alpha = alpha;
  v_pool[ name ].there = true;
  }

 void delete_linearization( Index name , ModParam issueMod = eModBlck )
  override {
  if( name < v_pool.size() ) {
   v_pool[ name ].there = false;
   v_pool[ name ].g.clear();
   }
  }

/*--------------------------------------------------------------------------*/

 void get_linearization_coefficients( FunctionValue * g ,
                                      Range range = INFRange ,
                                      Index name = Inf< Index >() )
  override {
  const auto & src = ( name == Inf< Index >() ) ? f_g : entry( name ).g;
  range.second = std::min( range.second , Index( src.size() ) );
  for( Index j = range.first ; j < range.second ; ++j )
   *(g++) = src[ j ];
  }

 void get_linearization_coefficients( FunctionValue * g ,
                                      c_Subset & subset ,
                                      bool ordered = false ,
                                      Index name = Inf< Index >() )
  override {
  const auto & src = ( name == Inf< Index >() ) ? f_g : entry( name ).g;
  for( auto j : subset )
   *(g++) = src[ j ];
  }

 FunctionValue get_linearization_constant( Index name = Inf< Index >() )
  override {
  return( ( name == Inf< Index >() ) ? f_alpha : entry( name ).alpha );
  }

/*--------------------------------------------------------------------------*/

 Index get_num_active_var( void ) const override { return( v_vars.size() ); }

 Index is_active( const Variable * var ) const override {
  auto it = f_pos.find( var );
  return( it == f_pos.end() ? Inf< Index >() : it->second );
  }

 void map_active( c_Vec_p_Var & vars , Subset & map , bool ordered = false )
  const override {
  map.resize( vars.size() );
  for( Index i = 0 ; i < vars.size() ; ++i ) {
   map[ i ] = is_active( vars[ i ] );
   if( map[ i ] == Inf< Index >() )
    throw( std::invalid_argument( "OracleFunction::map_active: not an "
                                  "active Variable" ) );
   }
  }

 Variable * get_active_var( Index i ) const override { return( v_vars[ i ] ); }

 ThinVarDepInterface::v_iterator * v_begin( void ) override {
  return( new v_iterator( v_vars.begin() ) );
  }
 ThinVarDepInterface::v_const_iterator * v_begin( void ) const override {
  return( new v_const_iterator( v_vars.begin() ) );
  }
 ThinVarDepInterface::v_iterator * v_end( void ) override {
  return( new v_iterator( v_vars.end() ) );
  }
 ThinVarDepInterface::v_const_iterator * v_end( void ) const override {
  return( new v_const_iterator( v_vars.end() ) );
  }

/*--------------------------------------------------------------------------*/

 protected:

 void print( std::ostream & output ) override {
  output << "OracleFunction [" << this << "] with " << v_vars.size()
         << " Variable";
  }

/*--------------------------------------------------------------------------*/

 private:

 struct Entry {
  std::vector< double > g;
  double alpha = 0;
  bool there = false;
  };

 const Entry & entry( Index name ) const {
  if( ! is_linearization_there( name ) )
   throw( std::invalid_argument( "OracleFunction: no linearization " +
                                 std::to_string( name ) + " in the pool" ) );
  return( v_pool[ name ] );
  }

 v_col_var v_vars;
 std::unordered_map< const Variable * , Index > f_pos;
 Oracle f_oracle;
 std::vector< double > f_x;
 long f_nevals = 0;
 FunctionValue f_val , f_alpha = 0;
 std::vector< double > f_g;
 Index f_gp_max = 0;
 std::vector< Entry > v_pool;

 };  // end( class OracleFunction )

/*--------------------------------------------------------------------------*/

static Solver * attach( Block * block , const std::string & fn )
{
 auto c = Configuration::deserialize( fn );
 auto bsc = dynamic_cast< BlockSolverConfig * >( c );
 if( ! bsc )
  throw( std::invalid_argument( "attach: not a BlockSolverConfig: " + fn ) );
 bsc->apply( block );
 delete( bsc );
 if( block->get_registered_solvers().empty() )
  throw( std::logic_error( "attach: no Solver from " + fn ) );
 return( block->get_registered_solvers().front() );
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 if( argc < 4 ) {
  std::cerr << "Usage: " << argv[ 0 ] << " problem n BSPar [ fstar [ tol ] ]"
            << std::endl;
  return( 1 );
  }

 int next = std::atoi( argv[ 1 ] );
 int n = std::atoi( argv[ 2 ] );
 std::vector< double > x0( n );
 startx_( & n , x0.data() , & next );
 if( next < 0 )
  throw( std::invalid_argument( "no such problem" ) );

 auto root = new AbstractBlock();
 auto x = new std::vector< ColVariable >( n );
 root->add_static_variable( *x );
 for( int j = 0 ; j < n ; ++j )
  ( *x )[ j ].set_value( x0[ j ] );
 LinearFunction::v_coeff_pair zero;
 for( auto & xj : *x )
  zero.emplace_back( & xj , 0.0 );
 auto robj = new FRealObjective( root , new LinearFunction( std::move( zero ) ) );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 std::vector< ColVariable * > vars( n );
 for( int j = 0 ; j < n ; ++j )
  vars[ j ] = & ( *x )[ j ];
 std::vector< double > xx( n );
 auto of = new OracleFunction( std::move( vars ) ,
                               [ n , next , &xx ]( const double * p ,
                                                   double * g ) mutable {
  std::copy( p , p + n , xx.begin() );
  double f;
  func_( & n , xx.data() , & f , g , & next );
  return( f );
  } );
 auto sub = new AbstractBlock( root );
 auto sobj = new FRealObjective( sub , of );
 sobj->set_sense( Objective::eMin );
 sub->set_objective( sobj );
 root->access_nested_Blocks().push_back( sub );

 auto s = attach( root , argv[ 3 ] );
 s->set_log( & std::cout );
 using clk = std::chrono::steady_clock;
 auto t0 = clk::now();
 auto st = s->compute();
 double tt = std::chrono::duration< double >( clk::now() - t0 ).count();
 double v = s->get_var_value();
 std::cout << std::setprecision( 12 ) << "problem " << next << ", n = " << n
           << ": status " << st << ", value " << v << ", iterations "
           << s->get_elapsed_iterations() << ", evaluations "
           << of->get_nevals() << ", time " << tt << " s" << std::endl;
 bool ok = ( st == Solver::kOK );
 if( argc > 4 ) {
  const double fs = std::atof( argv[ 4 ] );
  const double tol = ( argc > 5 ) ? std::atof( argv[ 5 ] ) : 1e-4;
  const double err = ( v - fs ) / ( 1 + std::abs( fs ) );
  std::cout << "error     : " << err << std::endl;
  ok = ok && ( err <= tol );
  }
 std::cout << ( ok ? "OK" : "KO" ) << std::endl;
 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- End File test_kbm.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
