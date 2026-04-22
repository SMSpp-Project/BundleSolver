/*--------------------------------------------------------------------------*/
/*-------------------- File MasterProblemBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/

#include "MasterProblemBlock.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <ostream>
#include <string>
#include <vector>

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- utility constants -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace
{
 static const std::string empty_string;
 static const std::vector< int > empty_vint;
 static const std::vector< std::string > empty_vstr;
}

/*--------------------------------------------------------------------------*/
/*------------------------- clear / reinitialize ---------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::clear()
{
 /* Reset the structural information owned by the master problem block. */

 MaxBSize = 0;

 NoEasyCmps = 0;
 NoHardCmps = 0;

 EasyCmps.clear();
 HardCmps.clear();

 InnerSolver = CDASolver{};
}

void MasterProblemBlock::SetDim( int MxBSz , int NVars , 
    int NrFi , int NrFiEasy ){
 
 //TBD
}

/*--------------------------------------------------------------------------*/
/*--------------------------- set parameters -------------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par( idx_type par , int value )
{
 /* No MPBlock-specific integer parameters are currently defined in the
  * header enum, so forward to the base Block interface.
  */
 Block::set_par( par , value );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par( idx_type par , double value )
{
 /* No MPBlock-specific double parameters are currently defined. */
 Block::set_par( par , value );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par( idx_type par , std::string && value )
{
 /* No MPBlock-specific string parameters are currently defined. */
 Block::set_par( par , std::move( value ) );
}

/*--------------------------------------------------------------------------*/
/* NOTE:
 * The header contains a duplicated signature:
 *
 *   void set_par( idx_type par , std::string && value ) override;
 *
 * twice in a row, while the comment says "vector-of-int".
 * The implementation below assumes the intended declaration is:
 *
 *   void set_par( idx_type par , std::vector< int > && value ) override;
 *
 * Therefore there is no second std::string&& overload here.
 */
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par( idx_type par , std::vector< int > && value )
{
 /* No MPBlock-specific vector<int> parameters are currently defined. */
 Block::set_par( par , std::move( value ) );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par(
 idx_type par , std::vector< std::string > && value )
{
 /* No MPBlock-specific vector<string> parameters are currently defined. */
 Block::set_par( par , std::move( value ) );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_log( std::ostream * log_stream )
{
 Block::set_log( log_stream );
}

/*--------------------------------------------------------------------------*/
/*---------------------------- get parameters ------------------------------*/
/*--------------------------------------------------------------------------*/

int MasterProblemBlock::get_int_par( idx_type par ) const
{
 /* No MPBlock-specific integer parameters currently implemented. */
 return Block::get_int_par( par );
}

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_dbl_par( idx_type par ) const
{
 /* No MPBlock-specific double parameters currently implemented. */
 return Block::get_dbl_par( par );
}

/*--------------------------------------------------------------------------*/

const std::string & MasterProblemBlock::get_str_par( idx_type par ) const
{
 /* No MPBlock-specific string parameters currently implemented. */
 return Block::get_str_par( par );
}

/*--------------------------------------------------------------------------*/

const std::vector< int > & MasterProblemBlock::get_vint_par(
 idx_type par ) const
{
 /* No MPBlock-specific vector<int> parameters currently implemented. */
 return Block::get_vint_par( par );
}

/*--------------------------------------------------------------------------*/

const std::vector< std::string > & MasterProblemBlock::get_vstr_par(
 idx_type par ) const
{
 /* No MPBlock-specific vector<string> parameters currently implemented. */
 return Block::get_vstr_par( par );
}

}  // end( namespace SMSpp_di_unipi_it )