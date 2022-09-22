#ifndef BIN_HASH_TABLE_H
#define BIN_HASH_TABLE_H

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

#define EXPANSION (6)

//#include <map>

class BinHashTable
{
    
    protected :
        
    long* indexArr;
    int* dataArr;
    
    // -- capacity of the has table --
    int capacity;
    
    // -- count of the elements in the hash table
    int count;
    
    inline int getHash( const long& i ) const
    {
        return int( i % capacity );
    }    

    
    int hasElement( const long i ) const
    {
        int hashCode = getHash(i);
       
        while( indexArr[hashCode] != i && indexArr[hashCode] >= 0 )
        {
           ++hashCode; 
	   hashCode = hashCode%capacity;
        }	
	
	return indexArr[hashCode] >= 0 ? hashCode : -1;
    }
    
    
    //std::map<long,int> binMap; 
    
    public :
        
    BinHashTable() 
    {
       
       capacity = 0;
       count = 0;
       
       indexArr = (long*)malloc( sizeof( long ) * capacity );
       dataArr = (int*)malloc( sizeof( int ) * capacity );
       
       for( int i = 0; i < capacity; ++i )
       {
          dataArr[i] = -1;
          indexArr[i] = -1;
       }
       
    }
    
    ~BinHashTable()
    {
        free( indexArr );
	free( dataArr );
    }
    
    inline int size() const
    {
        return count;
    }
    
    // -- constant access --
    inline int operator[]( const long& i ) const
    {
              
       const int hashCode = hasElement( i );
       
       return hashCode >= 0 ? dataArr[hashCode] : -1;
       
       /*if( binMap.count(i) == 0 )
	   return -1;
       else
	   return binMap.at(i);*/

    }
    
    
    // -- non-constant access -- 
    void set( const long& i, const int& value )
    {
	
	const int hashCode = hasElement(i);
	
	if( hashCode >= 0 )
	{ 
	    dataArr[hashCode] = value;
	}else
	{
	   
	    // -- create a new element, but do not allow memory increase --
	    ++count;
	    
	    if( EXPANSION*count > capacity )
	    {
		// incorrectly reallocated -> throw error
		std::cout<<"Count = "<<count<<" Capacity = "<<capacity<<std::endl;
		throw "Insufficient capacity";
	    }
	    
	    int newHashCode = getHash(i);
	    
	    while( indexArr[newHashCode] >= 0 )
	    {
	        ++newHashCode;
		newHashCode = newHashCode%capacity;
	    }

	    indexArr[newHashCode] = i;
	    dataArr[newHashCode] = value;
	    
	}
	
    }

    // -- overwrite all the elements in table --
    void clear()
    {
    
       count = 0;
       
       for( int i = 0; i < capacity; ++i )
       {
          dataArr[i] = -1;
          indexArr[i] = -1;       
       }
       
       /*
       const int& capacity_ = capacity;
        
       int* dataArr_ = (int*)__builtin_assume_aligned(dataArr, sizeof( int ) * EXPANSION );
       long* indexArr_ = (long*)__builtin_assume_aligned(indexArr, sizeof( long ) * EXPANSION );
       
       for( int i = 0; i < capacity_; ++i )
       {
          dataArr_[i] = -1;
          indexArr_[i] = -1;
       }
       */
    }

    // -- allocated more memory (if necessary) and overwrite all the current data
    void reserve( int n )
    {
         	 
	 if( capacity < EXPANSION * n )
	 {	 
	    capacity = EXPANSION * n;
	    
            indexArr = (long*)realloc( indexArr, capacity * sizeof( long ) );
            dataArr  = (int*) realloc( dataArr,  capacity * sizeof( int ) );	 
	 }
	 
	 clear();
	 
    }

};

#endif
