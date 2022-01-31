#include "NTL/ZZ.h"
#include "math.h"
#include <fstream>

NTL_CLIENT
using namespace NTL;
using std::cout;
using std::endl;

int aesSbox[] = {99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254, 215, 171, 118, 202, 130, 201, 125, 250,
                       89, 71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 183, 253, 147, 38, 54, 63, 247, 204, 52, 165,
                       229, 241, 113, 216, 49, 21, 4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39, 178, 117,
                       9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41, 227, 47, 132, 83, 209, 0, 237, 32, 252,
                       177, 91, 106, 203, 190, 57, 74, 76, 88, 207, 208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2, 127,
                       80, 60, 159, 168, 81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16, 255, 243, 210, 205,
                       12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61, 100, 93, 25, 115, 96, 129, 79, 220, 34, 42, 144,
                       136, 70, 238, 184, 20, 222, 94, 11, 219, 224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98, 145,
                       149, 228, 121, 231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234, 101, 122, 174, 8, 186, 120,
                       37, 46, 28, 166, 180, 198, 232, 221, 116, 31, 75, 189, 139, 138, 112, 62, 181, 102, 72, 3, 246, 14,
                       97, 53, 87, 185, 134, 193, 29, 158, 225, 248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233, 206,
                       85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104, 65, 153, 45, 15, 176, 84, 187, 22};
                       
const int N = 10000;
const int size = 256;
int population[N][size];
int pBest[N][size]; 
int gBest[size];
int pBestNL [N];
int gBestNL [1];   
    ZZ p [N];
    ZZ g [1];
    ZZ v [N];
    ZZ x [N]; 
long long pBestCost [N];
long long gBestCost [1];  

ifstream file ("NL-To-Search.txt"); 
   int NLToSearch;           

int *numberToFactorial(ZZ x, int n);

int *numberToSubstitution(int *number, int size);

int *substitutionToFactorial(int *sub, int size);

NTL::ZZ factorialCounting(int x);

NTL::ZZ factorialNumberToNumber(int *number, int size);

int raiseToPower(int num, int pow);

void FisherYates(int *arr, int n);

int *SBoxGeneratingDec(int n, int m, int counter);

int *valueToBinary(int i, int rank);

int *binaryElementsApprox(int *arr, int size, int count);

int *SBoxApprox(int *sbox, int size, int count);

int *elemsForN(int size);

int LATMax(int *sbox, int size, int count);

int myModulusDec(int number, int mod);

int *binaryElements(int *arr, int size, int count);

int *SBoxToBooleanFunc(int *sbox, int size, int count);

int *linearCombinations(const int *arr, int size, int count);

int myModulus(int number, int mod);

void bubbleSort(int *data, int size);

int *WHTSpectrumForLinearComb(const int *arr, int size, int count);

int *toPolarTable(const int *function, int size);

int *autoCorrelation(int *func, int size, int count);

int *ACForLinearComb(const int *arr, int size, int count);

int linearRedundancy(int *sbox, int size, int count);

void buildOneRow(int *arr, int *monomials);

int rankCalculation(int arr[256][697]);

int numOfCombinations(int n, int d);

int algebraicImmunity(const int *sbox, int size, int count);

int deltaUniformity(const int *arr, int size, int count);

int particleSwarmOptimization(int size, int count, int N);

long long costFunctionWHS(int *sbox_d, int size, int count);

long raiseToPowerLong(int num, int pow);

long long costFunctionCuba(int *sbox_d, int size, int count);

int *HadamardCoefficients(const int *func, int size, int count);

int HadamardMax(int *arr, int size);

int NLOfLinearCombinations(const int *arr, int size, int count);

int NLOfSBoxDec(int *sbox, int size, int count);

int HadamardNLinearity(int max, int count);

long long costFunctionCubaLinComb(int *sbox_d, int size, int count);

int main()
{

   file >> NLToSearch;  
   int ar = particleSwarmOptimization(256, 8, N);
   	
     	
   return 0;
}

NTL::ZZ factorialCounting(int x){
	ZZ result = conv<ZZ>(1);
	for (int i = 1; i <= x; ++i){
	result = result*i;
	}
    return result;
}

int *numberToFactorial(ZZ x, int n){
    int q = 2;
    int counter = n-2;
    int *positions = (int*) calloc (n,sizeof(int));
    positions[n-1] = 0;
    while (counter > 0){
        ZZ temp = x/q;
        //cout << "q =" << q;
        //cout << "temp =" << temp;
        ZZ val = (x - (temp*q));
        //cout << "val =" << val;
        int t;
        conv(t,val);
        positions[counter] = t;
        x = temp;
        ++q;
        counter--;
    }
    int u;
    conv(u,x);
    int lastp = u;
    positions[counter] = lastp;
    cout << "\nPOSITIONS VECTOR\n";
    for (int i = 0; i < n ; ++i){
    	cout << positions[i] << ",";
    }
    int *positionsRev = (int*) calloc (n,sizeof(ZZ));
    for (int r = 0,t = n-1; r < n, t>=0; ++r, t--){
        positionsRev[t] = positions[r];
    }
    /*cout << "\n\nPOSITIONS VECTOR\n";
    for (int i = 0; i < n ; ++i){
        cout << positionsRev[i];
    }*/
    //printf("\n");
    free(positions);
    return positionsRev;
}

int *numberToSubstitution(int *number, int size){
    int *S = (int*) calloc (size,sizeof(int));
    int *result = (int*) calloc (size,sizeof(int));
    int newSize = size;
    int counter = 0;
    int coeffNum = 0;
    int innerCounter = 0;
    while(counter < size){
        int *emptyPos = (int*) calloc (newSize,sizeof(int));
        for (int i = 0; i < size; ++i){
            if (S[i]==0){
                emptyPos[innerCounter] = i;
                innerCounter++;
            }
        }
        for (int k =0; k <newSize; ++k){
            //printf("%d, ",emptyPos[k]);
        }
        for (int q = 0; q < newSize; ++q){
            if (q==number[newSize-1]){
                //printf("q = %d ",q);
                //printf("MT pos = %d ",emptyPos[q]);
                S[emptyPos[q]] = 1;
                result[coeffNum] = emptyPos[q];
                ++coeffNum;
            }
        }
        innerCounter = 0;
        //printf("\n");
        for (int k =0; k <size; ++k){
            //printf("%d, ",S[k]);
        }
        //printf("\n");
        ++counter;
        newSize--;
        free(emptyPos);
    }
    int numberInArr = size;
    int *sub = (int*) calloc (size,sizeof(int));
    for (int u = 0; u < size; ++u){
        sub[result[u]] = numberInArr;
        numberInArr--;
    }
    int *subRev = (int*) calloc (size,sizeof(int));
    for (int r = 0,t = size-1; r < size, t>=0; ++r, t--){
        subRev[t] = sub[r];
    }
    free(sub);
    free(S);
    free(result);
    return subRev;
}

int *substitutionToFactorial(int *sub, int size) {
    int *result = (int*) calloc(size,sizeof(int));
    int value = size;
    int flag = 0;
    int innerCounter = 0;
    int counter = 0;
    while (counter<size) {
        for (int i = 0; i < size; ++i) {
            if (flag == 1 && sub[i] < value) {
                ++innerCounter;
            }
            if (sub[i] == value) {
                flag = 1;
            }
        }
        //printf("IC = %d ",innerCounter);
        result[counter] = innerCounter;
        innerCounter = 0;
        flag = 0;
        ++counter;
        value--;
    }
    return result;
}

NTL::ZZ factorialNumberToNumber(int *number, int size){
    int *numberRev = (int *) calloc (size,sizeof(int));
    ZZ result = conv<ZZ>(0);
    for (int r = 0,t = size-1; r < size, t>=0; ++r, t--){
        numberRev[t] = number[r];
    }
    for (int i = 0; i < size; ++i){
        result += numberRev[i]*factorialCounting(i);
        //cout << "\n" << factorialCounting(numberRev[i]);
        //cout << "\n" << result;
        //cout << "\n";
    }
    free(numberRev);
    return result;
}

//Функція визначення нелінійності функції через коефіцієнти Уолдша-Адамара

int HadamardNLinearity(int max, int count) {
    int nl = (raiseToPower(2, count) - max) / 2;
    return nl;
}

//Функція знаходження мінімальної нелінійності серед лінійних комбінацій

int NLOfLinearCombinations(const int *arr, int size, int count) {
    int *minimalNL = (int *)  calloc(size - 1, sizeof(int));
    int result;
    int *temp = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size - 1; ++i) {
        for (int j = 0; j < size; ++j) {
            temp[j] = arr[i * size + j];
        }
        int *fxarr = HadamardCoefficients(temp, size, count);
        int max1 = HadamardMax(fxarr, size);
        int nl2 = HadamardNLinearity(max1, count);
        minimalNL[i] = nl2;
        free(fxarr);
    }
    int min = 0;
    min = minimalNL[0];
    //printf("\nNON-LINEARITIES ARRAY");
    //printf("\n");
    for (int r = 0; r < size - 1; ++r) {
        //printf("%d ", minimalNL[r]);
        if (minimalNL[r] < min) {
            min = minimalNL[r];
        }
    }
    result = min;
    free(minimalNL);
    free(temp);
    return result;
}

//Функція знаходження нелінійності S-Box'у в десятковому вигляді

int NLOfSBoxDec(int *sbox, int size, int count) {
    int result;
    int *ar1 = SBoxToBooleanFunc(sbox,size,count);
    int *ar2 = linearCombinations(ar1, size, count);
    result = NLOfLinearCombinations(ar2, size, count);
    free(ar2);
    return result;
}

//Функція зведення до ступеня

int raiseToPower(int num, int pow) {
    int res = 1;
    for (int i = 0; i < pow; ++i) {
        res *= num;
    }
    return res;
}

//Функція зведення до ступеня long

long raiseToPowerLong(int num, int pow) {
    long res = 1;
    for (int i = 0; i < pow; ++i) {
        res *= num;
    }
    return res;
}

//Перемішування Фішера-Йейтса

void FisherYates(int *arr, int n) {
    int i, j, tmp;

    for (i = n - 1; i > 0; i--) {
        j = rand() % (i + 1);
        tmp = arr[j];
        arr[j] = arr[i];
        arr[i] = tmp;
    }
}

//Функція визначення найбільшого коефіцієнта перетворення Уолдша-Адамара

int HadamardMax(int *arr, int size) {
    int maxCoefficient = abs(arr[0]);
    for (int i = 0; i < size; ++i) {
        if (abs(arr[i]) > abs(maxCoefficient)) {
            maxCoefficient = abs(arr[i]);
        }
    }
    //printf("max coef %d", maxCoefficient);
    return maxCoefficient;
}

//Функція генерації десяткового S-Box'у

int *SBoxGeneratingDec(int n, int m, int counter) {
    int size = raiseToPower(2, n);
    int *dec = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size;) {
        dec[i] = rand() % size;
        int contains = 0;
        for (int j = 0; j < i; ++j) {
            if (dec[i] == dec[j]) {
                contains = 1;
                break;
            }
        }
        if (!contains) {
            i++;
        }
    }
    /*printf("Generated s-box: ");
    for (int i = 0; i < size; ++i) {
        printf("%d, ", dec[i]);
    }
    printf("\n");*/
    FisherYates(dec,size);
    return dec;
}

//Функція перетворення числа з десяткової СЧ у двійкову СЧ

int *valueToBinary(int i, int rank) {
    int *res = (int *) calloc(rank, sizeof(int));
    for (int j = 0; j < rank; ++j) {
        res[rank - 1 - j] = i >> j & 1;
    }
    return res;
}

//Функція генерації двійкових елементів у зворотньому порядку

int *binaryElementsApprox(int *arr, int size, int count) {
    int *result = (int *) calloc(size * count, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin = valueToBinary(arr[i], count);
        for (int j = count - 1; j >= 0; j--) {
            result[j * size + i] = bin[j];
        }
        free(bin);
    }
    return result;
}

//Функція генерації S-Box'у у зворотньому порядку

int *SBoxApprox(int *sbox, int size, int count) {
    int *result = binaryElementsApprox(sbox, size, count);
    return result;
}

//Функція генерації чисел для вхідних векторів ступеня N

int *elemsForN(int size) {
    int *result = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        result[i] = i;
    }
    return result;
}

//Функція знаходження LAT та її максимуму

int LATMax(int *sbox, int size, int count) {
    int *ar = SBoxApprox(sbox, size, count);
    int *elems = elemsForN(size);
    int *binelems = binaryElementsApprox(elems, size, count);
    int *temp = (int *) calloc(size, sizeof(int));
    int *temp2 = (int *) calloc(size, sizeof(int));
    int *coefficients = (int *) calloc(size * size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin1 = valueToBinary(i, count);
        for (int k = count - 1; k >= 0; k--) {
            if (bin1[k]) {
                //printf("K===%d ",k);
                //printf("X == \n ");
                for (int l = 0; l < size; ++l) {
                    temp[l] = temp[l] ^ binelems[k * size + l];
                    //printf("%d ",temp[l]);
                }
                //printf("\n ");
            }
        }
        //printf("\n ");
        for (int j = 0; j < size; ++j) {
            int *bin2 = valueToBinary(j, count);
            for (int q = count - 1; q >= 0; q--) {
                if (bin2[q]) {
                    //printf("K===%d ",k);
                    //printf("\nY [%d]== \n ", j);
                    for (int w = 0; w < size; ++w) {
                        temp2[w] = temp2[w] ^ ar[q * size + w];
                        //printf("%d ",temp2[l]);
                    }
                }
            }
            //printf("\n ");
            int calc = 0;
            for (int r = 0; r < size; ++r) {
                temp2[r] = temp2[r] ^ temp[r];
                //printf("%d ", temp2[l]);
                if (temp2[r] == 0) {
                    ++calc;
                }
                temp2[r] = 0;
            }
            int result = 0;
            result = calc - (size / 2);
            //printf("COEFFS = %d ", result);
            coefficients[i * size + j] = result;
            free(bin2);
        }
        for (int t = 0; t < size; ++t) {
            temp[t] = 0;
        }
        //printf("\n ");
        free(bin1);
    }
    for (int n = 0; n < size; ++n) {
        for (int m = 0; m < size; ++m) {
            //printf("%d ", coefficients[n*size+m]);
        }
        //printf("\n");
    }
    int result = 0;
    for (int p = 1; p < size * size; p++) {
        if (abs(coefficients[p]) > result)
            result = abs(coefficients[p]);
    }
    free(ar);
    free(elems);
    free(binelems);
    free(temp);
    free(temp2);
    free(coefficients);
    return result;
}

//Функція приведення числа за модулем дуякого числа

int myModulusDec(int number, int mod) {
    if (number < 0) {
        while (number < 0) {
            number = number + mod;
        }
    }
    return number % mod;
}

//Функція перетворення елементів з десяткової СЧ у двійкову СЧ, для певного ступеня N

int *binaryElements(int *arr, int size, int count) {
    int *result = (int *) calloc(size * count, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin = valueToBinary(arr[i], count);
        for (int j = 0, k = count - 1; j < count; ++j, k--) {
            result[j * size + i] = bin[k];
        }
        free(bin);
    }
    return result;
}

//Функція перетворення S-Box'у на набір булевих функцій

int *SBoxToBooleanFunc(int *sbox, int size, int count) {
    //printf("\nS-BOX\n");
    /*for (int i = 0; i < size; ++i) {
        printf("%d ", sbox[i]);
    }*/
    //printf("\n");
    //printf("\nS-BOX IN BOOLEAN FUNCTIONS REPRESENTATION\n");
    int *result = binaryElements(sbox, size, count);
    /*for (int i = 0; i < count; ++i) {
        printf("Function %d = ", i + 1);
        for (int j = 0; j < size; ++j) {
            printf("%d ", result[i * size + j]);
        }
        printf("\n");
    }*/

    //printf("\n");

    /*for (int i = 0; i < count; ++i) {
         int *temp = calloc(size, sizeof(int));
         //printf("Function %d", i);
         for (int j = 0; j < size; ++j) {
             temp[j] = result[i * size + j];
         }
         int weight = HammingWeight(temp, size);
         int flag = funcIsBalanced(weight, count);
         //printf("\n");
         free(temp);
     }*/
    return result;
}

//Функція знаходження лінійних комбінацій для булевих функцій S-Box'у

int *linearCombinations(const int *arr, int size, int count) {
    int *result = (int *) calloc(size*(size-1), sizeof(int));
    int *calc = (int *) calloc(size, sizeof(int));
    for (int i = 1; i < size; ++i) {
        int *bin = valueToBinary(i, count);
        for (int j = 0, k = count - 1; j < count, k >= 0; ++j, k--) {
            if (bin[k] == 1) {
                for (int w = 0; w < size; ++w) {
                    calc[w] = calc[w] ^ arr[j * size + w];
                    //printf(" %d", arr[j*size]);
                    //printf(" %d", j * size + k);
                    //printf("calc =  %d", calc[k]);
                    //result[(i-1)*size+k] = calc[k];
                    //printf("result  =  %d", (i-1)*size+k);
                }
                //printf("\n");
            }
            for (int r = 0; r < size; ++r) {
                result[(i - 1) * size + r] = calc[r];
            }
            //printf(" %d", bin[j]);
        }
        for (int l = 0; l < size; ++l) {
            //printf("calc =  %d", calc[l]);
            //result[(i-1) * size + l] = calc[l];
            calc[l] = 0;
        }
        //printf("\n");
        free(bin);
    }
    free(calc);
    return result;
}

//Функція приведення числа за модулем дуякого числа

int myModulus(int number, int mod) {
    if (number < 0) {
        while (number < 0) {
            number = number + mod;
        }
    }
    return number % 2;
}

//Функція визначення коефіцієнтів перетворення Уолдша-Адамара

int *HadamardCoefficients(const int *func, int size, int count) {
    int *result = (int *) calloc(size, sizeof(int));
    int *test = (int *) calloc(size * count, sizeof(int));
    int *functions2 = elemsForN(size);
    /*for (int i = 0; i < size; ++i) {
        printf(" %d",functions2 [i]);
    }*/
    //printf("\n");
    for (int i = 0; i < size; ++i) {
        int *bin = valueToBinary(functions2[i], count);
        for (int j = 0; j < count; ++j) {
            //printf(" bin j = %d", bin[j]);
            //*(functions + i * cols + j) = (i >> cols - j - 1) & 1u;
            test[i * count + j] = bin[j];
            //printf(" %d",test [i * count + j]);
        }
        //printf("\n");
        free(bin);
    }
    int *w = (int *) calloc(count, sizeof(int));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < count; ++j) {
            w[j] = test[i * count + j];
            //printf("w = %d", w[j]);
        }
        int res = 0;
        for (int j = 0; j < size; ++j) {
            int r = 0;
            for (int k = 0; k < count; ++k) {
                r += myModulus(w[k] * test[j * count + k], 2);
            }
            res += raiseToPower(-1, myModulus(func[j] + r, 2));
        }
        result[i] = res;
    }
    free(test);
    free(functions2);
    free(w);
    return result;
}

//Функція бульбашкового сортування

void bubbleSort(int *data, int size) {
    int i, j;
    for (i = 0; i < size; ++i) {
        for (j = size - 1; j > i; --j) {
            if (data[j] < data[j-1]) {
                int t = data[j - 1];
                data[j - 1] = data[j];
                data[j] = t;
            }
        }
    }
}

//Функція знаходження спектру Уолдша-Адамара для кожної з лінійних комбінацій

int *WHTSpectrumForLinearComb(const int *arr, int size, int count) {
    int *result = (int *) calloc(size*(size-1), sizeof(int));
    int *temp = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size - 1; ++i) {
        //printf("\nCombination %d", i + 1);
        for (int j = 0; j < size; ++j) {
            temp[j] = arr[i * size + j];
            //printf("%d ", temp[j]);
        }
        int *fxarr = HadamardCoefficients(temp, size, count);
        for (int g = 0; g< size; ++g){
            fxarr[g] = abs(fxarr[g]);
        }
        bubbleSort(fxarr,size);
        for (int g = 0; g< size; ++g){
            result[i*size+g] = fxarr[g];
        }
        //printf("\nHADAMARD COEFFICIENTS");
        //printf("\n");
        /*for (int q = 0; q < size; ++q) {
            //printf("%d ", fxarr[q]);
            result[i*size+q] = fxarr[q];
        }*/
        free(fxarr);
    }
    free(temp);
    return result;
}

//Функція представлення таблиці істиності в полярному вигляді

int *toPolarTable(const int *function, int size) {
    int *res = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        res[i] = raiseToPower(-1, function[i]);
    }
    return res;
}

//Функція обчислення автокореляційної функції

int *autoCorrelation(int *func, int size, int count) {
    int temp = 0;
    int *acFunc = (int *) calloc(size, sizeof(int));
    int *polFunc = toPolarTable(func, size);
    int *polFunc2 = (int *) calloc(size, sizeof(int));
    for (int i =0, k = size-1; i<size, k>=0; ++i, k--) {
        polFunc2[i] = polFunc[k];
        //printf("\npf2= %d", polFunc2[i]);
    }
    acFunc[0] = raiseToPower(2, count);
    for (int f = 1; f < size; ++f) {
        for (int j = 0; j < size; ++j) {
            temp = polFunc2[j] * polFunc2[j ^ f];
            //printf("\nj = %d", j);
            //printf("\nj^i = %d", j^i);
            //printf("\ntemp= %d", temp);
            acFunc[f] = acFunc[f] + temp;
        }
        //printf("ac i = %d", acFunc[i]);
    }
    free(polFunc);
    free(polFunc2);
    return acFunc;
}

//Функція знаходження автокореляційної функції для кожної з лінійних комбінацій

int *ACForLinearComb(const int *arr, int size, int count) {
    int *result = (int *) calloc(size*(size-1), sizeof(int));
    int *temp = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size - 1; ++i) {
        //printf("\nCombination %d", i + 1);
        for (int j = 0; j < size; ++j) {
            temp[j] = arr[i * size + j];
            //printf("%d ", temp[j]);
        }
        int *ar = autoCorrelation(temp, size, count);

        for (int g = 0; g< size; ++g){
            ar[g] = abs(ar[g]);
        }

        bubbleSort(ar,size);
        for (int q = 0; q < size; ++q) {
            //printf("%d ", fxarr[q]);
            result[i*size+q] = ar[q];
        }
        free(ar);
    }
    free(temp);
    return result;
}

//Функція знаходження лінійної збитковості S-Box'у

int linearRedundancy(int *sbox, int size, int count) {
    int result;
    int sp [size-1][size];
    int ac [size-1][size];
    int *ar1 = SBoxToBooleanFunc(sbox, size, count);
    //printf("\n1");
    int *ar2 = linearCombinations(ar1, size, count);
    /*for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", sp[i][j]);
        }
        printf("%d ",counter);
        printf("\n");
        counter++;
    }*/
    free(ar1);
    //printf("\n2");
    //int *ar5 = calloc((size-1),sizeof(int));
    //int ar33[size-1][size];
    int *ar3 = WHTSpectrumForLinearComb(ar2, size, count);
    //printf("\n3");
    /*for (int i = 0; i < size*(size-1); ++i){
        printf("%d ", ar3[i]);
    }*/
    int *ar4 = ACForLinearComb(ar2, size, count);
    free(ar2);
    //printf("\n4");
    //ar5 = DegForLinearComb(ar2,size,count);
    for (int b = 0; b < size - 1; ++b) {
        for (int q = 0; q < size; ++q) {
            sp[b][q] = ar3[b * size + q];
        }
    }
    for (int z = 0; z < size - 1; ++z) {
        for (int n = 0; n < size; ++n) {
            ac[z][n] = ar4[z * size + n];
        }
    }
    free(ar3);
    free(ar4);
    /*int dg [size-1];
    for (int i = 0; i < size - 1; ++i) {
        dg[i] = ar5[i];
    }*/
    /*FILE *file;
    fopen_s(&file, "Linear comb and WHT Spectrum.txt", "w");
    if (file == NULL) {
        printf("ERROR: Can't save sbox to file!\n");
        for (;;);
    }
    fprintf(file, "\n");
    for (int i = 0; i < size-1; ++i){
        fprintf(file,"\nLINEAR COMBINATION\n");
        for(int j = 0; j<size; ++j){
            ar33[i][j] = ar2[i*size+j];
            fprintf(file,"%d, ", ar33[i][j]);;
        }
        fprintf(file,"\n");
        fprintf(file,"\nHADAMARD SPECTRUM\n");
        for(int k = 0; k<size; ++k){
            fprintf(file,"%d ", sp[i][k]);
        }
        fprintf(file,"\n");
    }
    fprintf(file, "\n");
    fclose(file);*/

    /*printf("\nHADAMARD SPECTRUM\n");
    for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", sp[i][j]);
        }
        printf("\n");
    }
    printf("\nAUTO CORRELATION FUNCTIONS\n");
    for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", ac[i][j]);
        }
        printf("\n");
    }
    printf("\nDEGREES\n");
    for (int i = 0; i < size-1; ++i){
        printf("%d ", dg[i]);
    }
    printf("\n");*/
    int innerCounter = 0;
    int OuterCounter = 0;
    int finalCounter = 0;
    for (int i = 0; i < size - 1; ++i) {
        for (int h = i + 1; h < size - 1; ++h) {
            for (int j = 0; j < size; ++j) {
                if (i != h) {
                    if (sp[i][j] == sp[h][j] && sp[h][j] != -999 &&
                        (ac[i][j] == ac[h][j] && ac[h][j] != -999) /*(dg[i] == dg[h] && dg[h] != -999)*/) {
                        innerCounter++;
                    }
                }
            }
            //printf("\nInner counter = %d %d %d", i, h,innerCounter);
            /*printf("\n");
            printf("\n");
            for (int j = 0; j < size; ++j){
                printf("%d ", sp[i][j]);
            }
            printf("\n");
            printf("\n");
            for (int j = 0; j < size; ++j) {
                printf("%d ", sp[h][j]);
            }*/
            if (innerCounter == size) {
                //if (i!=h) {
                //if (ar4[i]==ar4[h] && ar5[i] == ar5[h]) {
                //printf("\nSTRINGS ARE EQUAL");
                //dg[h] = -999;
                OuterCounter++;
                for (int d = 0; d < size; ++d) {
                    sp[h][d] = -999;
                    //ac[h][j] = -999;
                    //dg[h] = -999;
                }
                //}
                //}
            }
            innerCounter = 0;
        }
        /*for (int j = 0; j < size; ++j){
            sp[i][j] = -999;
            ac[i][j] = -999;
            dg[i] = -999;
        }*/
        //printf("\nOC ==%d ", OuterCounter);
    }
    /*int counter = 0;
    printf("\nHADAMARD SPECTRUM AFTER\n");
    for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", sp[i][j]);
        }
        printf("%d ",counter);
        printf("\n");
        counter++;
    }*/
    /*printf("\nAUTO CORRELATION FUNCTIONS AFTER\n");
    for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", ac[i][j]);
        }
        printf("\n");
    }
    printf("\nDEGREES AFTER\n");
    for (int i = 0; i < size-1; ++i){
        printf("%d ", dg[i]);
    }*/
    finalCounter = finalCounter + OuterCounter;
    //printf("\nFINAL ==%d ", finalCounter);
    result = (size - 1) - finalCounter;
    return result;
}

//Функція побудування одного рядка матриці, що описує S-Box

void buildOneRow(int *arr, int *monomials) {
    monomials[0] = 1;
    //monomials x1,x8,y1,...,y8
    for (int i = 1; i <= 16; i++)
        monomials[i] = arr[i - 1];
    int pos = 17;
    //monomials x1x2
    for (int i = 1; i < 16; i++) {
        for (int j = i + 1; j <= 16; j++) {
            monomials[pos] = monomials[i] & monomials[j];
            pos++;
        }
    }
    //monomials x1x2x3
    for (int i = 1; i < 15; i++) {
        for (int j = i + 1; j <= 16; j++) {
            for (int k = j + 1; k <= 16; k++) {
                monomials[pos] = monomials[i] & monomials[j] & monomials[k];
                pos++;
            }
        }
    }
}

//Функція обчислення рангу матриці

int rankCalculation(int arr[256][697]) {
    int m = 697;
    int n = 256;

    int rank = 697;
    int  line_used[697] = { 0, };
    for (int i = 0;i < 697;i++)
        line_used[i] = 0;
    for (int i = 0; i < m; ++i) {
        int j;
        for (j = 0; j < n; ++j)
            if (!line_used[j] && arr[j][i])
                break;
        if (j == n)
            --rank;
        else {
            line_used[j] = 1;
            for (int k = 0; k < n; ++k)
                if (k != j && arr[k][i])
                    for (int p = i + 1; p < m; ++p)
                        arr[k][p] ^= arr[j][p] & arr[k][i];
        }
    }
    return rank;
}

//Функція знаходження кількості сполучень

int numOfCombinations(int n, int d) {
    if (n == d)
        return 1;
    if (d == 1)
        return n;
    if (d == 0)
        return 1;
    return numOfCombinations(n - 1, d - 1) + numOfCombinations(n - 1, d);
}

//Функція обчислення алгебраїчного імунітету

int algebraicImmunity(const int *sbox, int size, int count) {
    int mat[256][697];
    int tmp[697];
    int values[16];
    int *input_values = (int *) calloc(size * count, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin = valueToBinary(i, count);
        for (int j = 0; j < count; ++j) {
            input_values[i * count + j] = bin[j];
            //printf("%d ", input_values[i*count+j]);
        }
        //printf("\n");
        free(bin);
    }
    for (int i = 0; i < 256; i++) {
        int y = sbox[i];
        values[0] = input_values[i * count + 0];
        values[1] = input_values[i * count + 1];
        values[2] = input_values[i * count + 2];
        values[3] = input_values[i * count + 3];
        values[4] = input_values[i * count + 4];
        values[5] = input_values[i * count + 5];
        values[6] = input_values[i * count + 6];
        values[7] = input_values[i * count + 7];
        values[8] = input_values[y * count + 0];
        values[9] = input_values[y * count + 1];
        values[10] = input_values[y * count + 2];
        values[11] = input_values[y * count + 3];
        values[12] = input_values[y * count + 4];
        values[13] = input_values[y * count + 5];
        values[14] = input_values[y * count + 6];
        values[15] = input_values[y * count + 7];
        buildOneRow((int *) &values, (int *) &mat[i]);
    }
    int rank = rankCalculation(mat);
    free(input_values);
    //printf("%d", rank);
    return rank == 256 ? 3 : 2;
}

//Функція знаходження дельта-рівномірності

int deltaUniformity(const int *arr, int size, int count) {
    int result;
    int max = 0;
    for (int a = 1; a < size; ++a) {
        for (int b = 0; b < size; ++b) {
            result = 0;
            for (int x = 0; x < size; ++x) {
                if ((arr[x] ^ arr[x ^ a]) == b) {
                    ++result;
                }
            }
            if (result > max) {
                max = result;
            }
        }

    }
    return max;
}

//Функція "вартості" WHS

long long costFunctionWHS(int *sbox_d, int size, int count) {
    long long res = 0;
    int X1 = 21;
    int R1 = 7;
    long long calc = 1;
    int *sbox = SBoxToBooleanFunc(sbox_d,size,count);
    for (int i = 0; i < count; ++i) {
        int *fxarr = HadamardCoefficients(sbox + i * size, size, count);
        /*printf("\nHADAMARD COEFFICIENTS");
        printf("\n");
        for (int q = 0; q < size; ++q) {
            printf("%d ", fxarr[q]);
        }*/
        for (int j = 0; j < size; ++j) {
            //printf (" %d ", fxarr[j]);
            int w = abs(abs(fxarr[j]) - X1);
            //printf(" %d ", w);
            res += raiseToPowerLong(w, R1);
            //printf("res = %lld ", res);
        }
        //printf("\n");
        free(fxarr);
    }
    /*int cost;
    cost = costArray[0];
    //printf("\n");
    //printf("\nCOST ARRAY");
    //printf("\n");
    for (int t = 0; t < size - 1; ++t) {
        //printf("%d ", costArray[t]);
        if (costArray[t] > cost) {
            cost = costArray[t];
        }
    }
    free(costArray);*/
    free(sbox);
    return res;
}

//Функція "вартості" S-Box'у з Cuba

long long costFunctionCuba(int *sbox_d, int size, int count) {
    long long res = 0;
    int X1 = 21;
    int R1 = 7;
    long long calc = 1;
    int *sbox_b = SBoxToBooleanFunc(sbox_d, size, count);
    int C [] = {0,4,8,12,16,20,24,28,32};
    for (int i = 0; i < count; ++i) {
        int *fxarr = HadamardCoefficients(sbox_b + i * size, size, count);
        /*printf("\nHADAMARD COEFFICIENTS");
        printf("\n");
        for (int q = 0; q < size; ++q) {
            printf("%d ", fxarr[q]);
        }*/
        for (int j = 0; j < size; ++j) {
            //printf (" %d ", fxarr[j]);
            for (int k = 0; k < 9; ++k){
                long long w = abs(abs(fxarr[j]) - C[k]);
                //printf("w=  %d ", w);
                calc = calc*w;
                //printf("calc=  %d ", calc);
            }
            //printf(" %d ", w);
            res += calc;
            calc = 1;
            //printf("res = %lld ", res);
        }
        //printf("\n");
        free(fxarr);
    }
    /*int cost;
    cost = costArray[0];
    //printf("\n");
    //printf("\nCOST ARRAY");
    //printf("\n");
    for (int t = 0; t < size - 1; ++t) {
        //printf("%d ", costArray[t]);
        if (costArray[t] > cost) {
            cost = costArray[t];
        }
    }
    free(costArray);*/
    free(sbox_b);
    return res;
}

//Функція "вартості" S-Box'у з Cuba (розрахунок за лінійними комбінаціями)

long long costFunctionCubaLinComb(int *sbox_d, int size, int count) {
    long long res = 0;
    long long calc = 1;
    int C [] = {0,4,8,12,16,20,24,28,32};
    int *sbox_b = SBoxToBooleanFunc(sbox_d,size,count);
    int *ar1 = linearCombinations(sbox_b, size, count);
    for (int i = 0; i < size - 1; ++i) {
        int *temp = (int *) calloc(size, sizeof(int));
        for (int j = 0; j < size; ++j) {
            temp[j] = ar1[i * size + j];
        }
        int *fxarr = HadamardCoefficients(temp, size, count);
        /*printf("\nHADAMARD COEFFICIENTS");
        printf("\n");
        for (int q = 0; q < size; ++q) {
            printf("%d ", fxarr[q]);
        }*/
        for (int j = 0; j < size; ++j) {
            for (int k = 0; k < 9; ++k){
                long long w = abs(abs(fxarr[j]) - C[k]);
                calc = calc*w;
            }
            res += calc;
            calc = 1;
        }
        free(fxarr);
        free(temp);
    }
    free(ar1);
    return res;
}

//Функція генерації S-Box'у за допомогою методу Рою Часток

int particleSwarmOptimization(int size, int count, int N){
    cout << "Required NL = "<< NLToSearch;
    cout << "\nPopulation Initialization...\n";
    int flag102 = 0;
    clock_t start = clock();
    srand(time(NULL));
    int iter = 0;
    int flag = rand()%size;
    /*for (int i = 0; i < size; ++i){
        population[0][i] = aesSbox[i];
    }*/
    for (int q = 0; q < N; ++q){
        int *ar1 = SBoxGeneratingDec(count,count,q+flag);
        for(int w = 0; w < size; ++w) {
            population[q][w] = ar1[w];
        }
        free(ar1);
    }
    long long minCost = costFunctionCubaLinComb(population[0],size,count);
    for (int q = 0; q < N; ++q){
      	for(int w = 0; w < size; ++w){
            printf("%d ",population[q][w]);
        }
        int LAT = LATMax(population[q],size,count);
        int NL = raiseToPower(2, count - 1) - LAT;
        long long cost = costFunctionCubaLinComb(population[q],size,count);
        printf( "\nNon-linearity from LAT = %d \n", NL);
        printf("\n");
        pBestNL[q] = NL;
        pBestCost [q] = cost;
        if (cost <= minCost){
        minCost = cost;
            for (int m = 0; m < size; ++m){
        	gBest[m] = population[q][m];
    		}
    		int LAT4gB = LATMax(gBest,size,count);
    		int maxNL = raiseToPower(2, count - 1) - LAT4gB;
    		printf( "\nNL of GBest = %d \n", maxNL);
    		printf("\n");
  		gBestNL[0] = maxNL;
  		gBestCost[0] = cost;
        }
    }
    for (int y = 0; y < N; ++y){
     	for(int w = 0; w < size; ++w){
     	    pBest[y][w] = population[y][w];  
        }
    }
    printf("\ngBest\n");
    for (int m = 0; m < size; ++m){
        printf("%d, ",gBest[m]);
    }
    cout << "\nGBest Cost" << gBestCost[0] << "\n";
    printf("\n\npBest\n");
    for (int y = 0; y < N; ++y){
     	for(int w = 0; w < size; ++w){
            printf("%d ",population[y][w]);
        }
        cout << "\nPBest Cost" << pBestCost[y] << "\n";
        printf("\n\n");
    }
    
    
    //gBest
    
    int *sub = (int *) calloc (256, sizeof(int));
        for (int r = 0; r < 256; ++r){
            sub[r] = gBest[r]+1;
        }
        int *factorialNumber = substitutionToFactorial(sub,256);
   	//cout << "\n\nFACTORIAL NUMBER\n";
      	for (int i = 0; i < 256; ++i){
        //cout << factorialNumber[i] << ",";
   	}
   	ZZ gRes = factorialNumberToNumber(factorialNumber,256);
   	g[0] = gRes;
   	cout << "g = " <<g[0];
   	free(sub);
   	free(factorialNumber);
    
    //

    cout << "\n";
    cout << "\n";
    
    //pBest
    
    for (int i = 0; i < N; ++i){
    int *sub = (int *) calloc (256, sizeof(int));
        for (int r = 0; r < 256; ++r){
            sub[r] = pBest[i][r]+1;
        }
        int *factorialNumber = substitutionToFactorial(sub,256);
   	//cout << "\n\nFACTORIAL NUMBER\n";
      	for (int i = 0; i < 256; ++i){
        //cout << factorialNumber[i] << ",";
   	}
   	ZZ pRes = factorialNumberToNumber(factorialNumber,256);
   	p[i] = pRes;
   	cout << "\n";
   	cout << "p  = " <<p[i];
    	cout << "\n";
   	free(sub);
   	free(factorialNumber);
    }
    
    //x
    
    for (int i = 0; i < N; ++i){
    int *sub = (int *) calloc (256, sizeof(int));
        for (int r = 0; r < 256; ++r){
            sub[r] = population[i][r]+1;
        }
        int *factorialNumber = substitutionToFactorial(sub,256);
   	//cout << "\n\nFACTORIAL NUMBER\n";
      	for (int i = 0; i < 256; ++i){
        //cout << factorialNumber[i] << ",";
   	}
   	ZZ pRes = factorialNumberToNumber(factorialNumber,256);
   	x[i] = pRes;
   	cout << "\n";
   	cout << "x  = " <<x[i];
    	cout << "\n";
   	free(sub);
   	free(factorialNumber);
   	cout << "\n";
    }
    
    //
    
    //v
    
    for (int i = 0; i < N; ++i){
    v[i] = conv<ZZ>(0);
    cout << "v  = " <<v[i];
    	cout << "\n";
    	}
    
    //
    
    int chack = 1;
    ofstream PBResults;
    ofstream GBest;
   
    PBResults.open("PBResults.txt");
    PBResults << "PB of S-Boxes in Population\n";
    PBResults.close();

    while (chack){
    PBResults.open("PBResults.txt",std::ios::app);
    PBResults << "Iteration № " << chack << "\n";
    PBResults.close();
    ++chack;
    int Q = 1000000;
        int alpha = rand() % (Q);
        cout << "\n" << alpha;
        int beta = Q - alpha;
        cout << "\n" << beta;
        cout << "\n";
        cout << "NL of G-Best on iterration = " << gBestNL[0] << "\n";
        cout << "COST of G-Best on iterration = " << gBestCost[0] << "\n";
        cout << g[0]<< "\n";
        int speed = Q;
        /*if (gBestNL[0] >= 98){
        	--speed;
        }*/
    for (int i = 0; i < N; ++i){
        v[i] = (v[i] + ((p[i]-x[i])*alpha)/Q + ((g[0]-x[i])*beta))/speed;
        cout << "\n";
        cout << "vi = " << v[i];
        cout << "\n";
        cout << "\n";
        x[i] = (x[i] + v[i]);
        cout << "xi = " << x[i];
        cout << "\n";
        cout << "\n";
       int *factNum = numberToFactorial(x[i],256);
	//cout << "\n\nFACTORIAL NUMBER\n";
	int *subFin = numberToSubstitution(factNum,256);
	cout << "\n\nSUBSTITUTION\n";
	for (int r = 0; r < 256; ++r){
	  subFin[r] = subFin[r]-1;
	  cout << subFin[r] << ", ";
	}
	int LAT3 = LATMax(subFin, size, count);
        int NL3 = raiseToPower(2, count - 1) - LAT3;
        printf("\nNon-linearity = %d \n", NL3);
        long long cost = costFunctionCubaLinComb(subFin,size,count);
        cout << "\nCost = "<< cost << "\n";
        if (cost <= pBestCost[i]){
           if (cost <= gBestCost[0]){
           	gBestNL[0] = NL3;
           	g[0] = x[i];
           	gBestCost [0] = cost;
           }
        pBestNL[i] = NL3;
        p[i] = x[i];
        pBestCost[i] = cost; 
        }
        int ai = algebraicImmunity(subFin, 256, 8);
        cout << "NL3 = " << NL3;
        cout << "ai = " << ai;
        free(factNum);
    	free(subFin);
    	if (NL3 >= NLToSearch && ai == 3){
    	GBest.open("GBest.txt");
    	GBest << "Best S-box\n";
    	GBest << "Number: ";
   	GBest << p[i] << "\n";
   	GBest << "Cost: ";
   	GBest << pBestCost[i] << "\n";
	GBest << "\nNon-Linearity: ";
   	GBest << pBestNL[i] << "\n";
   	GBest.close();
    	GBest.close();
    	chack = 0;
    	}
    	PBResults.open("PBResults.txt",std::ios::app);
    	PBResults << "Number in population: " << i;
    	PBResults << "Number: ";
   	PBResults << p[i] << "\n";
   	PBResults  << "Cost: ";
   	PBResults  << pBestCost[i] << "\n";
	PBResults << "\nNon-Linearity: ";
   	PBResults << pBestNL[i] << "\n";
   	PBResults.close();
    }
    PBResults.open("PBResults.txt",std::ios::app);
   	PBResults << "\n\n";
   	PBResults.close();  
    }
    
    
    cout << "\n";
    cout << "\nPBEST NLs";
    for (int i = 0; i < N; ++i){
    cout << pBestNL[i];
    }
    
     cout << "\n";
    cout << "\n";
    
    cout << "\n";
    cout << "\nGBEST NL";
    cout << gBestNL[0];
    
        cout << "\n";
    cout << "\n";
    
    for (int i = 0; i < N; ++i){
        int *arr = numberToFactorial(p[i],256);
	int *factNum = (int *) calloc (256,sizeof(int));
	//cout << "\n\nFACTORIAL NUMBER\n";
	for (int i = 0; i < 256 ; ++i){
	  conv(factNum[i],arr[i]);
	  //cout << arr[i] << ",";
	}
	int *subFin = numberToSubstitution(factNum,256);
	cout << "\n\nSUBSTITUTION\n";
	for (int i = 0; i < 256; ++i){
	  subFin[i] = subFin[i]-1;
	  cout << subFin[i] << ", ";
	}
	int LAT3 = LATMax(subFin, size, count);
        int NL3 = raiseToPower(2, count - 1) - LAT3;
        printf("\nNon-linearity = %d \n", NL3);
        free(factNum);
    	free(subFin);
    	free(arr);
    }
    return 0;
}

