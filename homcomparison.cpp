#include <cstddef>
#include <sys/time.h>
#include "FHE.h"
#include "EncryptedArray.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <gmp.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include "Ctxt.h"
#include "polyEval.h"
#include <algorithm>
#include <math.h>



//define how many numbers you want to compare at VECTOR_COUNT
#define VECTOR_COUNT 2
#define VECTOR_SIZE 1

// Simple class to measure time for each method 
class Timer
{
public:
    void start() { m_start = my_clock(); }
    void stop() { m_stop = my_clock(); }
    double elapsed_time() const {
        return m_stop - m_start;
    }

private:
    double m_start, m_stop;
    double my_clock() const {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
};

//Equality recursive function at page 4 
Ctxt z(int i,int j, std::vector<Ctxt> x,std::vector<Ctxt> y, Ctxt enc1){
    if (j==1){
        Ctxt g=x[i];
        g+=y[i];
        g+=enc1;
        return g;    
    }
    else {
        int l=ceil(j/2);
        Ctxt f= z(i+l,j-l,x,y, enc1);
        f*=z(i,l,x,y,enc1);
        return f;                     
    }
       
}

//Inequality recursive function at page 4 
Ctxt t(int i,int j, std::vector<Ctxt> x,std::vector<Ctxt> y,  Ctxt enc1){
    if (j==1){
    Ctxt  H=x[i];
    H*=y[i];
    H+=x[i];
    return H;
    }
    else {
        int l1=ceil(j/2);
        Ctxt we=z(i+l1,j-l1,x,y,enc1); 
        we*=t(i,l1,x,y,enc1); 
        we+=t(i+l1,j-l1,x,y,enc1);
        return we;                                    
    }
 }

//selection  function with input two numbers (extracted to bits) with output the minimum of two numbers
std::vector<Ctxt> sel(std::vector<Ctxt> digits, std::vector<Ctxt> digits1,std::vector<Ctxt> digitsenc1,std::vector<Ctxt> digitsqw){
    for (int i=0;i<digits1.size();i++){
        Ctxt reset=digitsenc1[0];
        digits1[i]*=digitsqw[0];
        digitsenc1[0]-=digitsqw[0];
        digits[i]*=digitsenc1[0];
        digits1[i]+=digits[i];
        digitsenc1[0]=reset;
    }
    return digits1;

}

//function with input a list of encrypted numbers and output the minimum of them
std::vector<Ctxt> minimum(std::vector<Ctxt> list, Ctxt enc1){
    Ctxt mini=list[0];

    std::vector<Ctxt> digits1;
    extractDigits(digits1, mini);

   
    
    for (int i=1;i<list.size();i++){

         std::vector<Ctxt> digitsenc1;
        extractDigits(digitsenc1,enc1);
        
        std::vector<Ctxt> digitsi;
        extractDigits(digitsi, list[i]); 
         
        Ctxt compute_t=t(0,digits1.size(),digits1,digitsi, enc1);   

        std::vector<Ctxt> digitsqw;
        extractDigits(digitsqw, compute_t);
        std::vector<Ctxt> mi=sel( digits1, digitsi, digitsenc1, digitsqw);
        digits1=mi;
    }
    return digits1;

}





int main(int argc, char **argv)
{
    /*** BEGIN INITIALIZATION ***/
    long m = 0;                   // Specific modulus
    long p = 2;                 // Plaintext base [default=2], should be a prime number
    long r = 14;                   // Lifting [default=1]
    long L = 35;                  // Number of levels in the modulus chain [default=heuristic]
    long c = 2;                   // Number of columns in key-switching matrix [default=2]
    long w = 5;                  // Hamming weight of secret key
    long d = 1;                   // Degree of the field extension [default=1]
    long k = 80;                  // Security parameter [default=80] 
    long s = 0;                   // Minimum number of slots [default=0]
    
    Timer tInit;
    tInit.start();
	
    std::cout << "Finding m... " << std::flush;
    m = FindM(k, L, c, p, d, s, 0);           // Find a value for m given the specified values
    
    std::cout << "m = " << m << std::endl;
	
    std::cout << "Initializing context... " << std::flush;
    FHEcontext context(m, p, r); 	          // Initialize context
    buildModChain(context, L, c);             // Modify the context, adding primes to the modulus chain
    std::cout << "OK!" << std::endl;

    std::cout << "Generating keys... " << std::flush;
    
    fstream pubKeyFile("pk.txt", fstream::out|fstream::trunc);  
    assert(pubKeyFile.is_open());
    writeContextBase(pubKeyFile,context);
    pubKeyFile << context << std::endl;



    FHESecKey sk(context);                    // Construct a secret key structure
    const FHEPubKey& pk = sk;                 // An "upcast": FHESecKey is a subclass of FHEPubKey
    sk.GenSecKey(w);                          // Actually generate a secret key with Hamming weight
    //addSome1DMatrices(sk);                    // Extra information for relinearization
    std::cout << "OK!" << std::endl;


    pubKeyFile << pk << std::endl; 
    pubKeyFile.close();

    /****INITIALIZATION END****/


    std::ifstream infile("message1.txt");
   
   //open the message.txt file each line of this txt is a vector


    std::vector< std::vector<int> > e;
    e.resize(VECTOR_COUNT);
    for (int i=0; i<VECTOR_COUNT; i++){
        e[i].resize(VECTOR_SIZE);
    }
    
    for (int i=0; i<VECTOR_COUNT; i++){
        for (int j=0; j<VECTOR_SIZE ;j++){
            infile >> e[i][j]; 
        }
    }
    
    std::cout << "starting"<< std::endl;

   
    //put the first line vector to u and the second line vector to v

    /*******************************/
    /************CLIENT*************/
    /*******************************/
    long int u,v,y,z; 
    u=e[0][0];
    v=e[1][0];
    //y=e[2][0];
    //z=e[3][0];
    std::cout << "u:" << u << std::endl; 
    std::cout << "v:" << v << std::endl;
    //std::cout << "y:" << y << std::endl;
    //std::cout << "z:" << z << std::endl;
    std::cout << "encryption of two number from the file message1.txt" << std::endl;

    Ctxt encU(pk),encV(pk),encG(pk),enc1(pk),enc11(pk),enc0(pk),encY(pk),encZ(pk);


    
    

    pk.Encrypt(encU,to_ZZX(u));
    pk.Encrypt(encV,to_ZZX(v));
    pk.Encrypt(enc0,to_ZZX(0));
     pk.Encrypt(enc1,to_ZZX(1));
     //pk.Encrypt(encY,to_ZZX(y));
     pk.Encrypt(enc11,to_ZZX(1));
     //pk.Encrypt(encZ,to_ZZX(z));
     std::vector<Ctxt> digits0;
     extractDigits(digits0, enc0);

    
    //extractdigits of encU and store them at vector:digits
    std::vector<Ctxt> digitsU;
    extractDigits(digitsU, encU);
    
     

  	 /*************************************************************/
    /***decrypt each digit of encU to see if extractDigits work*****/
    /***************************************************************/

   

    long res[digitsU.size()];
    for (int i=0;i<digitsU.size();i++){
    ZZX result;
    sk.Decrypt(result,digitsU[i]);
    if (result[0]>(pow(p,r))/2){ result[0]=result[0]-pow(p,r);}
    res[i]=conv<long>(result[0]);
    }

    std::cout<< "U:";
    size_t res_size = sizeof(res)/sizeof(res[0]);
    std::reverse(res, res + res_size);
    for (int i=0;i<digitsU.size();i++){
        std::cout << res[i] << "," ;
    }
    std::cout << std::endl;

    /****CORRECT******/


    //extractdigits of encV and store them at vector:digitsV
    std::vector<Ctxt> digitsV;
    extractDigits(digitsV, encV);

     
    std::vector<Ctxt> digitsenc1;
    extractDigits(digitsenc1, enc1);

    /*************************************************************/
    /***decrypt each digit of encV to see if extractDigits work*****/
    /***************************************************************/


    
    long res1[digitsV.size()];
    for (int i=0;i<digitsV.size();i++){
    ZZX result;
    sk.Decrypt(result,digitsV[i]);
    if (result[0]>(pow(p,r))/2){ result[0]=result[0]-pow(p,r);}
    res1[i]=conv<long>(result[0]);
    }
    size_t res1_size = sizeof(res1)/sizeof(res1[0]);
    std::reverse(res1, res1 + res1_size);
    std::cout<< "V:";
    for (int i=0;i<digitsV.size();i++){
        std::cout << res1[i] << ",";
    }
    std::cout<< std::endl;
    /*******CORRECT******/

     //extractdigits of encV and store them at vector:digitsY
    //td::vector<Ctxt> digitsY;
    //extractDigits(digitsY, encY);


    /*************************************************************/
    /***decrypt each digit of encY to see if extractDigits work*****/
    /***************************************************************/

    /*long dig[digitsY.size()];
    for (int i=0;i<digitsY.size();i++){
     ZZX resultY;
    sk.Decrypt(resultY,digitsY[i]);
    if (resultY[0]>(pow(p,r))/2){ resultY[0]=resultY[0]-pow(p,r);}
    dig[i]=conv<long>(resultY[0]);
    }
    size_t dig_size = sizeof(dig)/sizeof(dig[0]);
    std::reverse(dig, dig + dig_size);
    std::cout<< "Y:";
    for (int i=0;i<digitsY.size();i++){
        std::cout << dig[i] << ",";
    }
    std::cout<< std::endl;
    */

    //input all encrypted numbers on list 
    std::vector<Ctxt> list;
    
        list.push_back(encU);
        list.push_back(encV);
        //list.push_back(encY);
        //list.push_back(encZ);
    

    Timer timecompare;
    timecompare.start();
    
   	//calculate the encrypted minimum of the list
    std::vector<Ctxt> kappa=minimum( list, enc1);
     
     timecompare.stop();
    

   


     //decryption of the minimum 
    long res2[digitsV.size()];
    for (int i=0;i<digitsV.size();i++){
    
        ZZX result2;
        sk.Decrypt(result2,kappa[i]);
        
        res2[i]=conv<long>(result2[0]);
    }
    //the number
    std::cout<< "the minimum number is:";
    int min=0,flow;
    for (int i=0;i<digitsV.size();i++){
        flow=pow(2,i)*res2[i];
        min+=flow;
    }
    std::cout << min << std::endl;
    std::cout << "and his binary form is:" ;
    //the binary form of the number
    size_t res2_size = sizeof(res2)/sizeof(res2[0]);
    std::reverse(res2, res2 + res2_size);
    for (int i=0;i<digitsU.size();i++){
        std::cout << res2[i] << "," ;
    }
    std::cout<< std::endl;
    std::cout << "the time to compute minimum is:" << timecompare.elapsed_time() << std::endl;

    

    std::cout << std::endl;

     return 0;




}