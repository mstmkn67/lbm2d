///// E-comment /////
//===== Array1d, Array2d, Array3d ===========================================================
//
//		A family of template Array classes for 1d 2d, 3d
//
//				ver.0.5		Aug.22, 1998	M.Doi
//===========================================================================

#ifndef _ARRAY_H_
#define _ARRAY_H_

//============================================================================
//
//  class Array1d: handle the 1d array a[i] with arbitrary bounds min <= i <=max;
//                 The bounds can be set by the constructor or the setBounds function
//                 The bounds are given by a.begin(), a.end();
//                 The operation a=b is supported (the size of a and b may differ.)
//  class Array2d  handle the 2d array a[i][j] with arbitrary bounds
//                    min1 <= i <= max1,  min2<= j <= max2
//                 the bounds are given by a.begin, a[i].begin() ...
//  class Array3d  handle the 3d array a[i][j][k] with arbitrary bounds
//                    min1 <= i <= max1,  min2<= j <= max2, min3<= k <= max3
//                 the bounds are given by a.begin, a[i].begin(), a[i][j].begin() etc..
//
//     Japanese manual is attatched at the end of this file
//============================================================================


template<class T> class Array1d;
template<class T> class Array2d;
template<class T> class Array3d;


class BaseArray {
   public:	
  	int begin() { return min;}
	int end()	{ return max+1;}
	int size()	{ return sizeData;}
   protected:
    BaseArray(){}
    void setMinMax(int min_, int max_){ 
	   min=min_; max =max_; 
	   max=larger(min, max);
	   sizeData =max-min+1;
    }
	int larger ( int n, int m){return ( n > m ? n : m);}
	int smaller( int n, int m){return ( n < m ? n : m);}
    int min, max, sizeData;
};


template <class T> 
class Array1d: public BaseArray {
   public:
     Array1d(int min_=0, int max_=0)
	   {setMinMax(min_, max_);
	    head = new T[sizeData];
		head_minus_min =head - min;
	   }
     Array1d (const Array1d<T>& array) {
		setMinMax(array.min, array.max);
		head = new T [array.sizeData];
		head_minus_min =head - min;
        *this = array;
	 }
	 ~Array1d() {delete[]  head;}
	 T& operator[] (int i)
	 	 {return *(head_minus_min + i); }

     Array1d<T>& operator= ( const Array1d<T>& rhs) {
		  int copyStart = larger(min, rhs.min);
		  int copyEnd   = smaller (max, rhs.max);
		  for (int i = copyStart; i <=copyEnd; i++)
		     head_minus_min[i]=  rhs.head_minus_min[i];
		  return *this;
	 }

	 Array1d<T>& setBounds(int newMin, int newMax) 
	    { Array1d<T>* pArray1d = new Array1d<T> (newMin, newMax);
		  *pArray1d = *this;	
          delete[] head;
		  head = pArray1d->head;
		  head_minus_min = head - min;
		  setMinMax(newMin, newMax);
		  return *this;
	   }

//    friend class Array1d<T>;
    friend class Array2d<T>;

   private:
     T* head; 
	 T* head_minus_min;
};

template <class T> 
class Array2d: public BaseArray  {
   public:
     Array2d(int min1_=0, int max1_=0, int min2_=0, int max2_=0)
	 { setMinMax(min1_, max1_);
	   array1d = new Array1d<T>*[sizeData];
	   array1d_minus_min = array1d - min;
	   for (int i=0; i< sizeData; i++){
			array1d[i] = new Array1d<T>(min2_, max2_);
       }
	 }
	 Array2d (const Array2d<T>& array){
		setMinMax(array.min, array.max);
	    array1d = new Array1d<T>*[sizeData];
	    array1d_minus_min = array1d - min;
	    for (int i=0; i< sizeData; i++){
			array1d[i] = new Array1d<T>(*(array.array1d[i]));
       }
	 }
	 ~Array2d() { for (int i=0; i< sizeData; i++) delete array1d[i]; 
	              delete array1d;}

	 Array1d<T>& operator[] (int i)
	       { return **(array1d_minus_min + i); }
     Array2d<T>& operator= ( const Array2d<T>& rhs) {
		  int copyStart = larger(min, rhs.min);
		  int copyEnd   = smaller (max, rhs.max);
		  for (int i = copyStart; i <=copyEnd; i++){
			 **(array1d_minus_min +i) = **(rhs.array1d_minus_min + i);
		  } 
		  return *this;
	 }

	  Array2d<T>& setBounds(int newMin1, int newMax1, int newMin2, int newMax2) 
	  {  Array2d<T>* pArray2d 
	              = new Array2d<T>(newMin1, newMax1, newMin2, newMax2);
	     *pArray2d = *this;
		 for (int i=0; i< sizeData; i++){ delete array1d[i];}
		 delete array1d;
		 array1d = pArray2d->array1d;
         array1d_minus_min = array1d - newMin1;
		 setMinMax(newMin1, newMax1);
		 return *this;
      }		  

   friend class Array3d<T>;
   private:
     Array1d<T>** array1d;
     Array1d<T>** array1d_minus_min;
};


template <class T> 
class Array3d: public BaseArray  {
   public:
     Array3d(int min1_=0, int max1_=0, int min2_=0, int max2_=0,
			 int min3_=0, int max3_=0)
	 { setMinMax(min1_, max1_);
	   array2d = new Array2d<T>*[sizeData];
	   array2d_minus_min = array2d - min;
	   for (int i=0; i< sizeData; i++){
			array2d[i] = new Array2d<T>(min2_, max2_, min3_, max3_);
       }
	 }

	 Array3d (const Array3d<T>& array)  {
		setMinMax(array.min, array.max);
  	    array2d = new Array2d<T>*[sizeData];
	    array2d_minus_min = array2d - min;
	    for (int i=0; i< sizeData; i++){
			array2d[i] = new Array2d<T>(*(array.array2d[i]));
        }
    }

	 ~Array3d() { for (int i=0; i< sizeData; i++) delete array2d[i]; 
	              delete array2d;}

	 Array2d<T>& operator[] (int i)
	       { return **(array2d_minus_min + i); }
     Array3d<T> & operator= (const Array3d<T>& rhs)  {
		  int copyStart = larger(min, rhs.min);
		  int copyEnd   = smaller (max, rhs.max);
		  for (int i = copyStart; i <=copyEnd; i++){
			 **(array2d_minus_min +i) = **(rhs.array2d_minus_min + i);
		  } 
		  return *this;
	 }

	 Array3d<T> & setBounds(int newMin1, int newMax1, int newMin2, int newMax2,
	                int newMin3, int newMax3) 
	  {  Array3d<T>* pArray3d 
	              = new Array3d<T>(newMin1, newMax1, newMin2, newMax2,
				                   newMin3, newMax3);
	     *pArray3d = *this;
		 for (int i=0; i< sizeData; i++){ delete array2d[i];}
		 delete array2d;
		 array2d = pArray3d->array2d;
         array2d_minus_min = array2d - newMin1;
		 setMinMax(newMin1, newMax1);
		 return *this;
       }		  
   private:
     Array2d<T>** array2d;
     Array2d<T>** array2d_minus_min;
};

#endif // _ARRAY_H_
