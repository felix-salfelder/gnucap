
#include <stdio.h>
#include <complex>
#include <iostream>
// #include "io_error.h"
#include "m_matrix_extra.h"
#include "m_matrix_ev.h"
#include "io_matrix.h"
#include "m_matrix.h"

static unsigned nm[10]={0,1,2,3,4,5,6,7,8,9};

void error(int badness, const char* fmt, ...)
{
  if (badness >= 0) {
    char buffer[BIGBUFLEN] = "";
    va_list arg_ptr;
    va_start(arg_ptr,fmt);
    vsprintf(buffer,fmt,arg_ptr);
    va_end(arg_ptr);
	 std::cerr << buffer;
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void error(int badness, const std::string& message)
{
  if (badness >= 0) {
	 std::cerr << message;
  }else{
  }
}

using namespace std;

int main(){
   double gmin=1e-9;

	BSMATRIX<double> a;
	BSMATRIX<double> lu;

	trace0("have bsmatrices.");

	a.reinit(3);
	lu.reinit(3);
	trace0("done init.");

	a.iwant(1,1);
	a.iwant(1,3);
	a.iwant(3,1);
	a.iwant(1,2);
	a.iwant(2,1);
	a.iwant(2,2);
	a.iwant(3,3);

	lu.iwant(1,1);
	lu.iwant(3,1);
	lu.iwant(1,3);
	lu.iwant(1,2);
	lu.iwant(2,1);
	lu.iwant(2,2);
	lu.iwant(3,3);
	trace0("done iwant");

	a.iwant(nm);
	a.allocate();
	lu.iwant(nm);
	lu.allocate();
	trace0("done alloc");

	a.load_point(1,2, 0.1);
	a.load_point(2,1, 0.1);
	a.load_point(1,1, 2.0);
	a.load_point(2,2, 1.0);
	a.load_point(3,3, .000000);
	a.load_point(1,3, 5);
	a.load_point(3,1, 5);
	trace0("done load");

	cout << "a matrix:\n" << a << "\n";

	a.set_min_pivot(0.0);
	//a.lu_decomp();

	double v[4] = {0,1,1,1};
	double w[4] = {1,1,1,0};

	cout << "v: "<< v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << "\n";
	lu.lu_decomp(a, false);
	cout << "lu decomp:\n" << lu << "\n";

	lu.fbsub(w,v,w);
	//lu.fbsub(v-1);

	cout << "v: " << v[0] << ", " << v[1] << ", " << v[2]<< ", " << v[3] << "\n";
	cout << "w: " << w[0] << ", " << w[1] << ", " << w[2]<< ", " << w[3] << "\n";

	BSMATRIX<double>* C = a.copy();

	cout << "orig\n" << a << "\n";
	cout << "copy\n" << *C << "\n";


	BSMATRIX<complex<double> > c;

	c.reinit(4);
	c.iwant(1,1);
	c.iwant(1,2);
	c.iwant(2,1);
	c.iwant(2,2);
	c.iwant(3,3);
	c.iwant(3,4);
	c.iwant(1,4);
	c.iwant(4,2);
	c.iwant(4,3);
	c.iwant(4,4);
	c.iwant(nm);
	c.allocate();
	c.load_point(3,3, 0.1);
	c.load_point(1,2, 1.0);
	c.load_point(1,1, 1.0);
	c.load_point(4,4, 1.0);
	complex<double> mycomplex (2.0,2.0);
	c.load_point(1,4, mycomplex);
	mycomplex = complex<double> (0.1,0.5);
	c.load_point(2,2, mycomplex);
	c.load_point(4,2, mycomplex);

	BSMATRIX<double> R = c.real();
	// BSMATRIX<double> R( BSMATRIX<double>::_REAL, c);

	cout << "complex matrix\n" << c << "\n";
	cout << "real part\n" << c.real() << "\n";

	double* b = new double[5];
	double x[5] = {0,2,1,5,5};

	cout << "x=    " << x[1] << ", " << x[2] << ", " << x[3]<< ", " << x[4] << "\n";
	b = R.rmul(b,x);
	cout << "Rx=b= " << b[1] << ", " << b[2] << ", " << b[3]<< ", " << b[4] << "\n";

	R.dezero(gmin);
	R.lu_decomp();
	R.fbsub(b);
	cout << "x: " << b[1] << ", " << b[2] << ", " << b[3]<< ", " << b[4] << "\n";

	R.zero();
	R.load_point(3,3, 5);
	R.load_point(2,2, 1.0);
	R.load_point(1,1, 1.0);
	R.load_point(1,2, -2.0);
	R.load_point(2,1, -2.0);
	R.load_point(4,4, 1.0);
	R.load_point(3,4, -1.0);
	R.load_point(4,3, -1.0);

	double _aug[2*R.size()-1];
	double _aug2[2*R.size()-1];
//	double* _aug2 = new double[2*R.size()-1];

	BSMATRIX<double>S(R);

	R.augment(_aug, _aug+R.size()-1);
	S.augment(_aug2, _aug2+S.size()-1);
	S.zero();

	for( unsigned i=0; i<2*R.size()-1; i++) _aug[i]=i;
	for( unsigned i=0; i<2*R.size()-1; i++) _aug2[i]=0;
#if 0
	_aug[2]=1;
	_aug[3]=.1;
	_aug[5]=.1;
	_aug[6]=1;
	_aug[4]=0;
#endif

	cout << "augmented " << R.size() << "\n" << R << "\n";
//	cout << "not lu'd \n" << S << "\n";
   // R.dezero(gmin);
	S.lu_decomp(R, false);
	cout << "lu'd \n" << S << "\n";

	R.zero();
	R.load_point(5,5, 1.);
	R.load_point(4,4, 1.);
	R.load_point(3,3, 1.);
	R.load_point(2,2, 1.);
	R.load_point(1,1, 0.);
	R.load_point(1,2, -2.0);
	R.load_point(2,1, -2.0);
	cout << "non-LU example\n";
	cout << R << "\n";
	R.lu_decomp();
	cout << "decomposed\n" << R << "\n";
}
