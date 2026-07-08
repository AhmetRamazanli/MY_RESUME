#include <iostream>
#include <string>
#include <fstream>
#include <tensors.h>
#include <cmath>

using namespace std;

// faktoryel fonksiyonu//

double factorial(int x)
{
        double faktor=1;
        int i;
        for(i=x; i>=1; i--)
        {
                faktor=faktor*i;
        }

        return faktor;
}

// kac farkli olasilik oldugunu hesaplama//

double combination(int a,int b)
{
     double comb = factorial(a) / ( factorial(b)*factorial(a-b) );
          return comb;
}


int main ()
{
        int nactmo, nelec;
        double nconfig;
        SharedTensor2i Wvertex;
        SharedTensor2i matrixvector;

        cout << "NMO degerini giriniz = ";
        cin >> nactmo;
        cout << endl <<"AEC degerini giriniz = ";
        cin >> nelec;

        const char *outfile;
        FILE *output;
        outfile = "output.dat";
        output = fopen(outfile, "w");

        nconfig = combination(nactmo,nelec);

        //cikti << "Number of configurations = " << nconfig << endl;
        fprintf(output,"\n Number of configurations = %3.0f \n", nconfig);
        fflush(output);


        // compute vertex weights
        Wvertex = SharedTensor2i(new Tensor2i("Wvertex", nactmo, nelec));
        for(int k=0; k < nactmo; k++){
                for(int i=0; i < nelec; i++){
                        //if (k-1 >= i) Wvertex[k][i] = (int)combination(k,i+1);
                        if (k-1 >= i) Wvertex->set(k,i,(int)combination(k,i+1));
                        else Wvertex->set(k,i,0);
                }
        }
        Wvertex->print(output);

        //malloc for matrixvector
        //double **matrixvector=init2d(nelec,nconfig) ;
        //matrixvector = SharedTensor2d(new Tensor2d("matrixvector", nelec,nconfig));
        matrixvector = SharedTensor2i(new Tensor2i("matrixvector", nelec,nconfig));

        for(int j=0; j<nconfig; j++){
                int I = j+1;
                for(int i=nelec-1 ;i>=0; i--){

                        for(int k=nactmo-1;k>=0; k--){
                                if(I>Wvertex->get(k,i) ){

                                        matrixvector->set(i,j,k+1);
                                        I = I - Wvertex->get(k,i);
                                        break;
                                }
                        }
                }
        }

        /*
        for(int j=0; j<nconfig; j++){
                for(int i=0; i<nelec; i++){
                        //cikti <<matrixvector->get(i,j) << " ";
                        fprintf(output, "%2.f   \n", <matrixvector->get(i,j));
                }
                //cikti <<endl;
        }
        */
        matrixvector->print(output);

        fclose(output);

        Wvertex.reset();
        matrixvector.reset();
        

	return 0;
}
