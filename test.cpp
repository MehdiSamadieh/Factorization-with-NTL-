#include <bpas.h>
#include <NTL/ZZ.h>
#include <NTL/ZZXFactoring.h>
#include <fstream>      // std::fstream, to write to file

using namespace std;
using namespace NTL;

extern void testRationalNumber();
extern void testComplexRationalNumber();
extern void testUQPoperations(int);
extern void testSubresultant();
extern void testGCD();
extern void testFactoring();
extern SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> 
          exchanging (SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> first,
                      string x, string y);
extern std::vector <SparseUnivariatePolynomial<RN>> rationalFactorize(SparseUnivariatePolynomial<RN>& a);

int main(int argc, char *argv[])
 {

	int s = 1;
	if (argc > 1)
	s = atoi(argv[1]);
	testFactoring();
	return 0;
}

void testFactoring() {


/////////////////////// counter of shofting to square free Norm /////
    long s=0;
////////////////////   temporary SUP for a //////////////////////////

	SparseUnivariatePolynomial<RN> C;
	C.setVariableName("a");
    SparseUnivariatePolynomial<RN> CC;
	CC.setVariableName("z");
/////////////////////// minimal polynomila for a ///////////////////
    

	SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> m;

	m.setVariableName("x");
	C.setVariableName("z"); 
	C.setCoefficient(0,RN(-3,1));
	m.setCoefficient(0,C);
    C.zero();
	C.setCoefficient(0,RN(1,1));
	m.setCoefficient(4,C);
	cout<<"\n m="<<m<<endl;

	
////////////////////// main polynomial //////////////////////////////
 	 SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> PP;

	 PP.setVariableName("z");
	 C.setVariableName("a");
	 C.zero();
     C.setCoefficient(0,RN(-2,1));
     PP.setCoefficient(0,C);
     C.zero();
     C.setCoefficient(0,RN(1,1));
     C.setCoefficient(2,RN(1,1));
     C.setCoefficient(3,RN(-2,1));
     PP.setCoefficient(1,C);
      C.zero();
     C.setCoefficient(0,RN(2,1));
     C.setCoefficient(1,RN(1,1));
     C.setCoefficient(2,RN(-1,1));
     PP.setCoefficient(2,C);
      C.zero();
	C.setCoefficient(0,RN(1,1));
    PP.setCoefficient(3,C);
    C.zero();
	C.setCoefficient(0,RN(1,1));
    PP.setCoefficient(4,C);

    cout<<"\n PP= --->     "<<PP<<"exchanging------>"<<exchanging(PP,"a","z")<<"\n";

/////////////////////////// set PP to PP_s /////////////////////////

SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> PP_s;
PP_s.setVariableName("x");
C.zero();
 C.setVariableName("z");

C.setCoefficient(0,RN(-2,1));
C.setCoefficient(1,RN(1,1));
C.setCoefficient(2,RN(2,1));
C.setCoefficient(3,RN(1,1));
C.setCoefficient(4,RN(1,1));
PP_s.setCoefficient(0,C);
C.zero();
 

C.setCoefficient(2,RN(1,1));

PP_s.setCoefficient(1,C);
C.zero();
 

C.setCoefficient(1,RN(1,1));
C.setCoefficient(2,RN(-1,1));

PP_s.setCoefficient(2,C);
C.zero();
 

C.setCoefficient(1,RN(-2,1));


PP_s.setCoefficient(3,C);


cout<<"\nPP_s=----->    "<<PP_s<<"\n";


/////////////////////////  Resultant & Norm ////////////////////////////////////////
SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> norm;
norm.setVariableName("x");

SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> dnorm(norm);
dnorm.setVariableName("x");

SparseUnivariatePolynomial<RN> supTemp;
supTemp.setVariableName("z");



norm=m.resultant(PP_s);
 cout<< "\n norm ---->  =" << norm<<"\n";

supTemp = norm.coefficient(0);
cout<< " \n\nnorm.coefficient(0);--->"<< norm.coefficient(0);

supTemp.differentiate(1);
cout<<"\n\nsuptemp---->"<<supTemp;

dnorm.setCoefficient(0,supTemp);
//dnorm=norm;

 //dnorm.differentiate(1);
 cout<< "\n dnorm ---->  =" << dnorm<<"\n";
cout <<" \n GCD(Norm, dnorm)-------->"<<dnorm.gcd(norm)<<"\n";
//cout<< "\n degree(GCD(Norm, dnorm))------>"<<(dnorm.gcd(norm)).degree();

/////////////////////////////// shifting ///////////////////////////////////////////
//SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> temp_m;
//SparseUnivariatePolynomial<RN> C1;
///////////////////////////////////////////////////////////////////////////
SparseUnivariatePolynomial<RN> temp_m;
       temp_m.setVariableName("a");
temp_m.setCoefficient(0,RN(-3,1));
temp_m.setCoefficient(4,RN(1,1));


	/*temp_m.setVariableName("z");
	C.setVariableName("a");
	C.setCoefficient(0,RN(-3,1));
	temp_m.setCoefficient(0,C);
    C.zero();
	C.setCoefficient(0,RN(1,1));
	temp_m.setCoefficient(4,C);
	*/

//////////////////////////////////////////////////////////////

/*SparseUnivariatePolynomial<RN> temp_m;
temp_m.setVariableName("a");
temp_m.setCoefficient(0,RN(-3,1));
temp_m.setCoefficient(4,RN(1,1));*/



    SparseUnivariatePolynomial<RN>temp_shift ;
    SparseUnivariatePolynomial<RN>temp1_shift ;
   temp_shift.setVariableName("a");
   temp1_shift.setVariableName("a");


    SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>>temp_PP ;
    SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>>temp1_PP ;
    SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>>temp2_PP ;
     temp1_PP.setVariableName("z");
     temp_PP.setVariableName("z");
     temp2_PP.setVariableName("z");

    temp_shift.setCoefficient(1,RN(1,1));//a
    
    temp1_shift.setCoefficient(1,RN(1,1));
   // cout<<"\ntemp1_PP----------->"<<temp_PP.monicDivide(m);
   //  temp1_PP=(temp_PP+temp_shift)^5;


        temp1_PP.setCoefficient(1,RN(1,1));//z
        temp2_PP.setCoefficient(0,temp1_shift);//a
		temp_PP.setCoefficient(1,RN(1,1));//z

       SparseUnivariatePolynomial<RN> coeftemp;
       coeftemp.setVariableName("a");

         int S=0;

        while((dnorm.gcd(norm)).degree()) 

         {
	 
	 	temp2_PP=temp_PP- S*temp1_shift;//z-Sa
        temp1_PP=temp2_PP;//z-Sa
   		temp_shift=PP.coefficient(PP.degree());
   		temp2_PP=temp2_PP * temp_shift;
   		temp_shift=PP.coefficient(PP.degree()-1);
   		temp2_PP=temp2_PP + temp_shift;

/////////////////////////////////////new  mod for  the first time ///////////////////////////////////
   		coeftemp=temp2_PP.coefficient(0).monicDivide(temp_m);
   		temp2_PP.setCoefficient(0, coeftemp);

  int coefcount=1; //// coef counter

///////////////////////////// horner /////////////////////////////////////////////////////
	    for(int i=PP.degree()-1;i>=1;i-- )
	   {
		temp_shift=PP.coefficient(i-1);
		temp2_PP =(temp2_PP*temp1_PP)+temp_shift;
///////////////////////////new mod for the rest ////////////////////////////////////////////
      coeftemp=temp2_PP.coefficient(coefcount).monicDivide(temp_m);
       temp2_PP.setCoefficient(coefcount, coeftemp);

       coefcount=coefcount+1;
       
       
      coeftemp=temp2_PP.coefficient(coefcount).monicDivide(temp_m);
       temp2_PP.setCoefficient(coefcount, coeftemp);


////////////////////////////////////////old mod /////////////////////////////////////////////////////

        //temp2_PP=exchanging(temp2_PP,"z","a");
        //cout<<"\n temp2_PP------->\n"<<temp2_PP;
   		//temp2_PP.monicDivide(temp_m);
        //temp2_PP=exchanging(temp2_PP,"z","a");
        //cout<<"temp2_PP.monicDivide(temp_m);" <<
       
        }
   ////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"\n\ntemp2_PP------------->"<<temp2_PP; S=S+1;

    temp2_PP=exchanging(temp2_PP,"x","z");
    norm.zero();
    norm=m.resultant(temp2_PP);
    cout<< "\n norm ---->  =" << norm<<"\n";
    supTemp.zero();
    supTemp = norm.coefficient(0);
    cout<< " \n\nnorm.coefficient(0);--->"<< norm.coefficient(0);
           
    supTemp.differentiate(1);
    cout<<"\n\nsuptemp---->"<<supTemp;
    dnorm.zero();
    dnorm.setCoefficient(0,supTemp);


    }
/////////////////////////// factoring on z ////////////////////////////////////////
SparseUnivariatePolynomial<RN> b;
b.setVariableName("z");
b=norm.coefficient(0);

cout<<"\n \nnorm______________>>>"<<b;

std::vector <SparseUnivariatePolynomial<RN>> factors;
factors = rationalFactorize(b);


    ///////////////////////////////// sub_resultant ///////////////////////////////////////////
   	if (factors.size()==2)
   	{

   // return b;

   	}

SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> temp_gcd;
    temp_gcd.setVariableName("z");
    SparseUnivariatePolynomial<RN> temp1_c;
    temp1_c.setVariableName("a");

    
    if (S==0)
    {
    	         temp2_PP=PP;
            
    }

    /////////////////////
	
	std::cout << "temp2_PP = " << temp2_PP << std::endl;

    for (int j=1;j<factors.size();++j)
    {
    	for (int i=0; i<=factors.at(j).degree();++i)
		{
			temp1_c.zero();
			temp1_c.setCoefficient(0,factors[j].coefficient(i));
			temp_gcd.setCoefficient(i,temp1_c);

		}
		std::cout << "temp_gcd = " << temp_gcd << std::endl;
	    vector<  SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> > subresult;
	    subresult = temp2_PP.subresultantChain(temp_gcd);
	    
	    for (int k=0;k<subresult.size();k++) {
	    	std::cout << "subresult[" << k << "] = " << subresult[k] << std::endl;
	    	std::cout << "content(subresult[" << k << "]) = " << subresult[k].content() << std::endl;
	    }

	    // if modularRemainder(subresult[0],m) != 0
	    	// GCD(temp2_PP,temp_gcd) = 1
	    // else
	    	// i=1
	    	// while (true)
	    		// compute sContent = content(subresult[i])
	    		// if modularRemainder(content,m) != 0
	    			// GCD = subresult[i];
	    			// break the while loop
	    		// compute sPrimpart = subresult[i]/sContent;
	    		// for (int l=0;l<subresult[i].degree();l++){
	    			// if modularRemainder(sPrimpart.coefficient(l),m) != 0
	    				// GCD = subresult[i];
	    				// break the while loop
				//	}
	    	// end while
    }
  


    //////////////////////////////  add variable  ////////////////////////////////////////////
    

   	 /////////////////////// minimal polynomila for a ///////////////////
    

	
/*
cout<<"\n\ntemp_m=="<<temp_m;
cout<<"\n \n subresult[1]=="<<subresult[1].coefficient(1);

cout<<"\n monicDivide=="<<subresult[1].coefficient(1).monicDivide(temp_m)
<<"  \n \n subresult[1].coefficient(1)"<<
subresult[1].coefficient(1);
*/



//cout<<"\nR===========>"<<subresult.size()<<subresult[1];

/*
temp_gcd.setCoefficient(0,Qfactored[1]);

vector<  SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> > R = 
temp2_PP.subresultantChain(temp_gcd);
*/




    
   //  temp2_PP= exchanging(temp2_PP,"z","a");
         
         //temp2_PP=exchanging(PP,"a","z");
     //cout<<"\n \n\ntemp2_PP+++++++++++++>"<<subresult[1]<<endl; 

    



//////////////////////////////////////////////////////////////////////////////////
/*

       SparseUnivariatePolynomial <RN> P;
	   P.setVariableName("x");

	   SparseUnivariatePolynomial<RN> Q(P);
     	


	P.setCoefficient(0,RN(1,1));
	P.setCoefficient(2,RN(2,-1));
	Q.setCoefficient(0,RN(-2,3));
	Q.setCoefficient(1,RN(5,2));
	

	cout<<"coefficient  "<<Q.coefficient(1).get_den()<<"     ";
    

	mpz_class big=Q.coefficient(1).get_den();
	cout<<Q<<"first"<<endl;
    

    Q*=P;
    cout<<Q;
    Q.resultant(P);
    Q.differentiate(1);


    double big1=big.get_d();
    cout<<"big1        "<<big1;
    //Q/big;

	// Q.setCoefficient(2,RN(2,1));
	// Q.setCoefficient(3,RN(3,1));
	// Q.setCoefficient(4,RN(4,1));
	Q.setCoefficient(5,RN(5,1));
	std::cout << "P = " << P << ", Q = " << Q << std::endl;
	cout<<Q.degree()<<endl;
	cout<<P.degree()<<endl;
	// cout<<Q.numberOfTerms()<<endl;
	// cout<<Q.coefficient(0)<<endl;
	// cout<<Q.coefficient(1)<<endl;
	// Taking the LCM and making the Q-poly to Z-poly
	// int Denom_Prod = 1;
	// for(int i=0;i<Q.degree();i++)
	// {
	// 	Denom_Prod = Denom_Prod*Q.coefficient(i).get_den();
	// }
	// Q = Q*Denom_Prod;

	DenseUnivariateRationalPolynomial Q_dense;
	Q_dense=Q.convertToDUQP();
	// cout<<Q_dense<<endl;
	// cout<<Q_dense.degree()<<endl;
	// FILE * dense_Q_file=fopen("dense_Q_file.txt","w");
	// char tmp[1024];
	ZZ tmp_zz,c;
	ZZX Q_zzx;
	Vec< Pair< ZZX, long > > factors;
	// cout<<tmp_zz<<endl;	

	
	for (int i=0;i<=Q_dense.degree();i++)
	{
		// cout<<"Q["<<i<<"]="<<Q_dense.coefficient(i)<<endl;
		
		// ZZX.setCoefficient
		// sprintf(tmp,"%s",());
		// cout<<tmp<<endl;

		conv(tmp_zz,Q_dense.coefficient(i).get_str().c_str());
		SetCoeff(Q_zzx,i, tmp_zz);
		// cout<<Q_dense.coefficient(i).get_str()<<endl;
		// sprintf(tmp, "%s", Q_dense.coefficient(i).get_str().c_str());
		// cout<<tmp_zz<<endl;
		// fprintf(dense_Q_file,"%s \n",Q_dense.coefficient(i).get_str().c_str());

	}
	cout<<Q_zzx<<endl;
	// fclose(dense_Q_file);



*/
///////////////////////////////////////////////
/*
factor(c, factors, Q_zzx);

   cout <<"c=" << c << "\n";
   cout <<"factors = "<< factors << "\n";


   // The factored polynomial 
   SparseUnivariatePolynomial<RN> Qfactored;

   // Iterate over all the factors in the factored polynomial
   for(int i=0;i<factors.length();i++)
   {
   		// Iterate over all the coefficients in each of the factors
   		for(int j=0;j<deg(factors[i].a)+1;j++)
   		{
   			// Convert polynomial to long from ZZ 
   			long temp_coeff = 0;
   			conv(temp_coeff,factors[i].a[j]);

   			// Assign the coefficient to the BPAS polymial
   			Qfactored.setCoefficient(j,RN(temp_coeff,1));
   		}

   		// Multiply the resulting polynomial p times to get the
   		// p^th power 
   		for(int j=1;j<factors[i].b;j++)
   		{
   			Qfactored = Qfactored*Qfactored;
   		}
   }
   	// Multiply by c 
   	long temp_coeff = 0;
   	conv(temp_coeff,c);
	Qfactored = Qfactored*temp_coeff;
*/

}
//////////////////////////////////////////////////////////////////////////////////
SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> 
exchanging (SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> first,
 string x, string y)
{
SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> temp_main;
				temp_main.setVariableName(x);
		SparseUnivariatePolynomial<SparseUnivariatePolynomial<RN>> temp_2;
			temp_2.setVariableName(x);
			SparseUnivariatePolynomial<RN> temp_1;
            temp_1.setVariableName(y);
            SparseUnivariatePolynomial<RN>  Q_first;
	        Q_first.setVariableName(x);

	for (int i = 0; i <= first.degree(); ++i)
	{
		Q_first.zero();
		Q_first=first.coefficient(i);

		for (int j = 0; j <= Q_first.degree(); ++j)
		{
			temp_1.setCoefficient(i,Q_first.coefficient(j));///////cheng double to gmp.
			temp_2.setCoefficient(j,temp_1);
			
			temp_1.zero();

            temp_main=temp_main+temp_2;
            temp_2.zero();
		}
	}
	return temp_main;
};

////////////////////////////////////////////////////////////////////////////////////

std::vector <SparseUnivariatePolynomial<RN>> rationalFactorize(SparseUnivariatePolynomial<RN>& a) {

	RN rnTemp(0,1);
	int iTemp;
	long lTemp(0);
	std::string variable(a.variable());
	SparseUnivariatePolynomial<RN> supTemp;
	supTemp.setVariableName(variable);
	std::vector <SparseUnivariatePolynomial<RN>> output;

	//SparseUnivariatePolynomial<RN> b;
	//b.setVariableName("z");
	//b=norm.coefficient(0);

	//cout<<"\n \nnorm______________>>>"<<b;

	//DenseUnivariateRationalPolynomial b_dense;
	//b_dense=b.convertToDUQP();

	ZZ intNTL,intContentNTL;
	ZZX polyNTL;
	Vec < Pair< ZZX, long > > aFactors;
		//cout<<"\n \nb_factors==========>"<<b_dense<<"\n\n";

	for (int i=0;i<=a.degree();i++)
	{
		if (!a.coefficient(i).isZero())
			conv(intNTL,a.coefficient(i).get_str().c_str());
		else 
			conv(intNTL,rnTemp.get_str().c_str());
		SetCoeff(polyNTL,i,intNTL);
	}
	cout<<intNTL<<endl;
		

	factor(intContentNTL, aFactors, polyNTL);

	std::cout << "content=" << intContentNTL << std::endl;
	std::cout << "factors = " << aFactors << std::endl;

   	// Input the content as the zeroth element of aFactored
   	conv(lTemp,intContentNTL);
   	iTemp = (int)lTemp;
   	rnTemp = RN(iTemp,1);
   	supTemp.setCoefficient(0,rnTemp);
   	output.push_back(supTemp);
    std::cout << "factors[0]===>" << supTemp << std::endl;

	// Iterate over all the factors in the factored polynomial
	for(int i=0;i<aFactors.length();i++)
	{   
   		supTemp.zero();
   		// Iterate over all the coefficients in each of the factors
   		for(int j=0;j<=deg(aFactors[i].a);j++)
   		{
   			// Convert polynomial to long from ZZ 
   			conv(lTemp,aFactors[i].a[j]);
   			iTemp = (int)lTemp;
   			rnTemp = RN(iTemp,1);

   			// Assign the coefficient to the BPAS polymial
   			supTemp.setCoefficient(j,rnTemp);
   		}
   		output.push_back(supTemp);
        std::cout << "factors[" << i+1 << "]===>" << supTemp << std::endl;
    }
    return output;
}

SparseUnivariatePolynomial<RN> modularRemainder(SparseUnivariatePolynomial<RN>& a, SparseUnivariatePolynomial<RN>& m) {

	SparseUnivariatePolynomial<RN> output;

	// compute the remainder of a modulo m and return the result


	return output;
}