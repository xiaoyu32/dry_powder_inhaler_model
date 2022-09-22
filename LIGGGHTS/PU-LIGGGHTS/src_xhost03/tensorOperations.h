   
  class Tensor
  {
  
      private :
      
      bool trans;
      
      mutable bool isCpy;
      
      // - data -
      double* data;
      
      // - number of rows - 
      int n;
      
      // - number of columns - 
      int m;
      
      public :
      
      Tensor() :
      n( 0 ),
      m( 0 ),
      isCpy( false ),
      trans( false ),
      data( NULL )
      {}
      
      
      Tensor( double* v, int n_, bool isCpy_ = true ) :
      n( n_ ),
      m( 1 ),
      isCpy( isCpy_ ),
      trans( false )
      {
      
          if( !isCpy )
	  {
             data = (double*)malloc( sizeof( double ) * n );

	     for( int i = 0; i < n_; ++i )
	     {
		 data[i] = v[i];
	     }
	  }else
	  {
	     data = v;
	  }
	  
      }
      
      //constructor for 
      Tensor( double* v, bool isCpy_ = true ) :
      n( 3 ),
      m( 1 ),
      isCpy( isCpy_ ),
      trans( false )
      {
          if( !isCpy )
	  {
             data = (double*)malloc( sizeof( double ) * 3 );

	     for( int i = 0; i < 3; ++i )
	     {
		 data[i] = v[i];
	     }
	  }else
	  {
	     data = v;
	  }      
      }
      
      Tensor( int n_, int m_ ) :
      n( n_ ),
      m( m_ ),
      isCpy( false ),
      trans( false )
      {
	  data = (double*)malloc( sizeof( double ) * n * m );  
	  
	  for( int i = 0; i < n_*m_; ++i )
	     data[i] = 0.0;
      }
      
      Tensor( double* v, int n_, int m_, bool isCpy_ = true ) :
      n( n_ ),
      m( m_ ),
      isCpy( isCpy_ ),
      trans( false )
      {
          
	  if( !isCpy )
	  {
             
	     data = (double*)malloc( sizeof( double ) * n * m );

	     for( int i = 0; i < n_*m_; ++i )
	     {
		 data[i] = v[i];
	     }
	     
	  }else
	  {
	     data = v;
	  }      
	  
      }
      
      Tensor( const Tensor& t ) :
      n( t.n ),
      m( t.m ),
      isCpy( false ),
      data( NULL ),
      trans( t.trans )
      {
	  data = (double*)malloc( sizeof( double ) * n * m );

	  for( int i = 0; i < n*m; ++i )
	  {
	      data[i] = t.data[i];
	  }          
      } 
      
      ~Tensor()
      {
          if( !isCpy && data ) free( data );
      }
      
      // -- creator for identity tensor --
      static Tensor eye()
      {
         
	 Tensor t( 3, 3 );
	 
	 t.get(0,0) = 1;
	 t.get(1,1) = 1;
	 t.get(2,2) = 1;
	 
	 return t;
	 
      }
      
      void print() const
      {
         
	 for( int i = 0; i < n; ++i )
	 {
	    std::cout<<"[ ";
	    for( int j = 0; j < m; ++j )
	        std::cout<<get(i,j)<<" ";
            std::cout<<" ]"<<std::endl;
	 }
      }
      
      const double& operator[](int i ) const
      {
          return data[i];
      }
            
      inline const double& get( int i, int j ) const
      {
          if( !trans )
	     return data[i+n*j];
	  else
	     return data[i+m*j];
	  
      }

      inline double& get( int i, int j )
      {
          if( !trans )
	     return data[i+n*j];
	  else
	     return data[i+m*j];
	  
      }
      
      inline int getN() const
      {
          
	  if( !trans )
	      return n;
	  else 
	      return m;
	  
      }
      
      inline int getM() const
      {
          
	  if( !trans )
	      return m;
	  else 
	      return n;
	  
      }
      
      Tensor transpose() const
      {
	   Tensor t( *this );
	   t.trans = !trans;
	   return t;
      }

      Tensor operator-(const Tensor& t ) const
      {
	  Tensor tt( this->n, this->m );
	  
	  if( this->getM() != t.getM() || this->getN() != t.getN() ) return tt;
	  
	  for( int i = 0; i < this->getN(); ++i )
	  {
	     for( int j = 0; j < this->getM(); ++j )
	     {
		double& element = tt.get(i,j);
		element = this->get(i,j) - t.get(i,j);
	     } 
	  }
	  
	  return tt;      
      }
      
      Tensor operator+(const Tensor& t ) const
      {
	  Tensor tt( this->n, this->m );
	  
	  if( this->getM() != t.getM() || this->getN() != t.getN() ) return tt;
	  
	  for( int i = 0; i < this->getN(); ++i )
	  {
	     for( int j = 0; j < this->getM(); ++j )
	     {
		double& element = tt.get(i,j);
		element = this->get(i,j) + t.get(i,j);
	     } 
	  }
	  
	  return tt;      
      }

      Tensor& operator -=( const Tensor& t )
      {
          
	  if( this->getM() != t.getM() || this->getN() != t.getN() ) return *this;

	  for( int i = 0; i < this->getN(); ++i )
	  {
	     for( int j = 0; j < this->getM(); ++j )
	     {
		double& element = this->get(i,j);
		element -= t.get(i,j);
	     } 
	  }	  
	  
	  return *this;
      }
      
      Tensor& operator +=( const Tensor& t )
      {
          
	  if( this->getM() != t.getM() || this->getN() != t.getN() ) return *this;

	  for( int i = 0; i < this->getN(); ++i )
	  {
	     for( int j = 0; j < this->getM(); ++j )
	     {
		double& element = this->get(i,j);
		element += t.get(i,j);
	     } 
	  }	  
	  
	  return *this;
      }
      
      Tensor operator*( const double a ) const
      {
          
	  Tensor t( this->data, this->n, this->m, false );
	  
	  for( int i = 0; i < n*m; ++i )
	  {
	      t.data[i] *= a;
	  }
	  
	  return t;
	  
      }

      Tensor& operator*=( const double a )
      {
          	  
	  for( int i = 0; i < n*m; ++i )
	  {
	      this->data[i] *= a;
	  }
	  
	  return *this;
	  
      }
      
      Tensor operator*( const Tensor& t ) const
      {
          
	  Tensor tt( this->getN(), t.getM() );
	  
	  if( this->getM() != t.getN() ) return tt;
	  
	  for( int i = 0; i < this->getN(); ++i )
	  {
	     for( int j = 0; j < t.getM(); ++j )
	     {

		double& element = tt.get(i,j);

		for( int k = 0; k < t.getN(); ++k )
		{
		    element += this->get(i,k) * t.get(k,j);
		}

	     } 
	  }
	  
	  return tt;
	  
      }
      
      Tensor crossproduct( const Tensor& other ) const
      {
          
          Tensor t( 3, 1 );
	  
	  t.data[0] = this->data[1]*other.data[2] - this->data[2]*other.data[1];
	  t.data[1] = this->data[2]*other.data[0] - this->data[0]*other.data[2];
	  t.data[2] = this->data[0]*other.data[1] - this->data[1]*other.data[0];
	  
	  return t;	  
	  
      }
      
      double mag() const
      {
         double rr = 0;
	 
	 for( int i = 0; i < n*m; ++i )
	     rr += data[i] * data[i];
	 
	 return sqrt( rr );
	 
      }
      
      inline double value() const
      {
          return data[0];
      }
      
      void operator=( Tensor& other )
      {
          
	  if( this != &other )
	  {
	      
	      this->n = other.n;
	      this->m = other.m;
	      
	      this->trans = other.trans;
	      
	      // - if tensor is reserving memory, free the memory -
	      if( this->data && !isCpy ) free( this->data );	 
		 
	      if( other.isCpy ) // -- other is a copy of another tensor -> pass the cpy pointer --
	      {

		  this->isCpy = true;
		  this->data = other.data;

	      }else // -- other is not a copy of another tensor -> change ownership of the pointer --
	      { 

		  other.isCpy = true;
		  this->isCpy = false;
		  this->data = other.data;

	      }
 	      
	  }
	  	  
      }
  };
