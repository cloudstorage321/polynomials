#include<stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include<string.h>

void print_sigma_single( int i, int max_deg, int alpha_deg,uint32_t ***sigma, const char *varname) {
    printf("%s^{(%d)}(x, a) = ", varname, i);
    int first = 1;

    for (int j = 0; j < max_deg; j++) {
        for (int k = 0; k < alpha_deg; k++) {
            if (sigma[i][j][k]) {
                if (!first) printf(" + ");
                first = 0;

                // Print a part
                if (k == 0)
                    ; // skip a^0
                else if (k == 1)
                    printf("a");
                else
                    printf("a^%d", k);

                // Print x part
                if (j == 1 && k != 0)
                    printf("x");
                else if (j == 1 && k == 0)
                    printf("x");
                else if (j > 1)
                    printf("x^%d", j);
                else if (j > 0 && k != 0)
                    printf("x^%d", j);
                else if (j > 0 && k == 0)
                    printf("x^%d", j);

                // If it's the constant 1
                if (j == 0 && k == 0)
                    printf("1");
            }
        }
    }

    if (first)
        printf("0");

    printf("\n");
}

// Function to find the degree of x for each sigma[i] polynomial
int find_degree_of_x_for_i(int x_deg, int alpha_deg, uint32_t **sigma) {
    int max_degree_x = -1;  // Degree of x (initialize to -1)

    // Iterate over all x degrees (index j)
    for (int j = 0; j <= x_deg; j++) {
        for (int k = 0; k <= alpha_deg; k++) {
            if (sigma[j][k] != 0) {
                // If sigma[j][k] is non-zero, update the degree of x
                if (j > max_degree_x) {
                    max_degree_x = j;
                }
            }
        }
    }

    return max_degree_x;  // Return the highest degree of x found
}

int find_best_rho1(int mu, int *d, int *l) {
    int best_rho = -1;
    int max_value = -1;

    for (int rho = 0; rho < mu; rho++) {
        // Check if d_rho is non-zero
        if (d[rho] != 0) {
            int value = 2 * rho - l[rho];
            
            // Maximize 2*rho - l_rho
            if (value > max_value) {
                max_value = value;
                best_rho = rho;
            }
        }
    }

    // Print the row from which rho was selected
    if (best_rho != -1) {
        printf("Best rho: %d selected from row: %d\n", best_rho, best_rho);
    } else {
        printf("No valid rho found.\n");
    }

    return best_rho;
}

// Multiply two field elements modulo the irreducible polynomial
uint32_t gf_multiply(uint32_t a, uint32_t b, uint32_t mod_poly, int m) {
    uint32_t result = 0;
    for (int i = 0; i < m * 2; i++) {
        if (b & 1)
            result ^= a;

        b >>= 1;
        int overflow = a & (1 << (m - 1));
        a <<= 1;
        if (overflow)
            a ^= mod_poly;
    }
    return result ;//& ((1 << m) - 1); // Keep result in field
}

// Exponentiation by squaring in GF(2^m)
uint32_t gf_pow(uint32_t base, uint32_t exp, uint32_t mod_poly, int m) {
    uint32_t result = 1;
    while (exp > 0) {
        if (exp & 1)
            result = gf_multiply(result, base, mod_poly, m);
        base = gf_multiply(base, base, mod_poly, m);
        exp >>= 1;
    }
    return result;
}

// Compute a^(-1) in GF(2^m)
uint32_t gf_inverse(uint32_t a, uint32_t mod_poly, int m) {
    uint32_t exp = (1 << m) - 2; // 2^m - 2
    return gf_pow(a, exp, mod_poly, m);
}

uint32_t poly_multiply(uint32_t a, uint32_t b) {
    uint32_t result = 0;
    while (b) {
        if (b & 1) result ^= a;
        a <<= 1;
        b >>= 1;
    }
    return result;
}

void bivariate_mult(uint32_t term1, uint32_t x_pow, int x, int alpha_deg, uint32_t **result) {
    for (int i = 0; i <= x; i++) {
        if ((x_pow >> i) & 1) {  // If x^i is present in x_pow
            for (int j = 0; j <= alpha_deg; j++) {
                if ((term1 >> j) & 1) {  // If y^j is present in term1
                    result[i][j] ^= 1;   // Flip bit in result for x^i y^j
                }
            }
        }
    }
}


void print_bivariate2d(int x_deg, int alpha_deg, uint32_t **result)
 {
 	int need_plus=0;
    for(int i=x_deg;i>=0;i--)
    {
    	for(int j=alpha_deg;j>=0;j--)
    	{
    		if(result[i][j])
    		{
    			if(need_plus) printf(" + ");
    			if(i==0 && j==0)  printf(" 1 ");
    			else if(i==1 && j==1) printf(" x * a");
    			else if(i==0 && j!=0) printf("a^%d",j);
    			else if(i!=0 && j==0) printf("x^%d",i);
    			else printf("x^%d a^%d",i,j);
    			need_plus=1;
			}
		}
	}
	if(!need_plus) printf(" 0 ");
}

void print_bivariate3d(int index,int x_deg, int alpha_deg, uint32_t ***result)
 {
 	int need_plus=0;
    for(int i=x_deg;i>=0;i--)
    {
    	for(int j=alpha_deg;j>=0;j--)
    	{
    		if(result[index][i][j])
    		{
    			if(need_plus) printf(" + ");
    			if(i==0 && j==0)  printf(" 1 ");
    			else if(i==1 && j==1) printf(" x * a");
    			else if(i==0 && j!=0) printf("a^%d",j);
    			else if(i!=0 && j==0) printf("x^%d",i);
    			else printf("x^%d a^%d",i,j);
    			need_plus=1;
			}
		}
	}
	if(!need_plus) printf(" 0 ");
}

void decimal_to_polynomial(uint32_t poly)
{
	int need_plus=0;
	for(int i=32;i>=0;i--)
	{
		if((poly>>i)&1)
		{
			if(need_plus) printf(" + ");
			if(i==0) printf(" 1 ");
			else if(i==1) printf(" a ");
			else printf("a^%d",i);
			need_plus=1;
		}
	}
	if(!need_plus) printf(" 0 ");
}

// Function to calculate the degree of the polynomial (highest set bit)
int degree(uint32_t poly) {
    int deg = -1;
    while (poly) {
        deg++;
        poly >>= 1;
    }
    return deg;
}
/*
void multiply_2d_with_3d(int x1,int alpha1,uint32_t **result,int index,int x_deg,int alpha_deg,uint32_t ***sigma,uint32_t ***temp,uint32_t irred_poly,int m)
{
	printf("\n=================================================\n");
	int index1=0;
	for(int i=0;i<=x1;i++)
	{
		for(int j=0;j<=alpha1;j++)
		{
			
			if(result[i][j]==0) continue;
			for(int k=0;k<=x_deg;k++)
			{
				for(int l=0;l<=alpha_deg;l++)
				{
					printf("\n result[%d][%d]=%u \n sigma[%d][%d][%d]=%u \n",i,j,result[i][j],index,k,l,sigma[index][k][l]);
					uint32_t ele=sigma[index][k][l];
					if(ele==0) continue;
					uint32_t prod=result[i][j]*ele;
					printf("prod=result*sigma-->%u\n",prod);
					temp[index1][i+k][j+l]^=prod;
					printf("temp[%d][%d][%d]=%u\n",index1,i+k,j+l,temp[index1][i+k][j+l]);
					//x ^i a^ j)·(x ^k a ^l)=x ^i+k (a ^j+l)
				}
			}
		}
	}
}*/
// Addition in GF(2^m) is XOR (since char=2)
uint32_t gf_add(uint32_t a, uint32_t b) {
    return a ^ b;
}
// Multiply a 2D poly with a 3D poly and store result in temp
void multiply_2d_with_3d(
    int x1, int alpha1,
    uint32_t **result,           // result[x][a]
    int index,
    int x_deg, int alpha_deg,
    uint32_t ***sigma,           // sigma[index][x][a]
    uint32_t ***temp,            // temp[index1][x][a]
    uint32_t irred_poly,
    int m
) {
    //printf("\n================ MULTIPLICATION START ================\n");
    int index1 = 0;  // Storing in temp[index1] for now

    for (int i = 0; i <= x1; i++) {
        for (int j = 0; j <= alpha1; j++) {
            if (result[i][j] == 0) continue;

            for (int k = 0; k <= x_deg; k++) {
                for (int l = 0; l <= alpha_deg; l++) {
                    uint32_t ele = sigma[index][k][l];
                    if (ele == 0) continue;

                    //printf("\nresult[%d][%d]=%u \nsigma[%d][%d][%d]=%u\n", i, j, result[i][j], index, k, l, ele);

                    uint32_t prod = gf_multiply(result[i][j], ele, irred_poly, m);
                    //printf("prod = result * sigma --> %u\n", prod);

                    temp[index1][i + k][j + l] = gf_add(temp[index1][i + k][j + l], prod);
                    //printf("temp[%d][%d][%d] = %u\n", index1, i + k, j + l, temp[index1][i + k][j + l]);
                }
            }
        }
    }

    //printf("================ MULTIPLICATION END ==================\n");
}
/*
void multiply_polynomials(
    int arr1[MAX_X_DEG][MAX_ALPHA_DEG],
    int arr2[MAX_X_DEG][MAX_ALPHA_DEG],
    int result[2 * MAX_X_DEG][2 * MAX_ALPHA_DEG]
) {
    // Initialize result to zero
    for (int i = 0; i < 2 * MAX_X_DEG; i++) {
        for (int j = 0; j < 2 * MAX_ALPHA_DEG; j++) {
            result[i][j] = 0;
        }
    }

    // Perform multiplication
    for (int i1 = 0; i1 < MAX_X_DEG; i1++) {
        for (int j1 = 0; j1 < MAX_ALPHA_DEG; j1++) {
            for (int i2 = 0; i2 < MAX_X_DEG; i2++) {
                for (int j2 = 0; j2 < MAX_ALPHA_DEG; j2++) {
                    // Multiply coefficients and add to the appropriate position
                    result[i1 + i2][j1 + j2] += arr1[i1][j1] * arr2[i2][j2];
                    // If working in F2:
                    result[i1 + i2][j1 + j2] %= 2;
                }
            }
        }
    }
}
*/

int main()
{
	int t;
	printf("Enter Error Correcting Capability t=");
	scanf("%d",&t);
	
	int m;
	printf("Enter m degree of extension field =");
	scanf("%d",&m);
	
	char input[32];
	printf("\n Enter irreducible polynomial of degree %d in Binary=",m);
	scanf("%s",input);
	uint32_t irred_poly=strtoul(input,NULL,2);
	printf("\n");
	decimal_to_polynomial(irred_poly);
	
	int max_iter=t+2;
	int x_deg=t;
	int syn_count=2*t;
	
	uint32_t syndromes[syn_count];
	//maximum degree of syndrome is m-1
	
	printf("\n Enter the %d syndromes :\n",syn_count);
	for(int i=0;i<syn_count;i++)
	{
		printf("\nSyndromes[%d]=",i);
		scanf("%u",&syndromes[i]);
		
	}
	for(int i=0;i<syn_count;i++)
	{
		printf("\n Syndromes[%d]=%u  ---->",i,syndromes[i]);
		decimal_to_polynomial(syndromes[i]);
	}
	printf("\n--------------------------\n");
	// constructing table
	double mue[max_iter];
	mue[0]=-0.5;
	for(int i=1;i<max_iter;i++)
		mue[i]=i-1;
	for(int i=0;i<max_iter;i++)
	{
		printf("\n mue[%d]=%.1f\n",i,mue[i]);
	}
	printf("\n--------------------------\n");
	int alpha_deg=m-1; //alpha
	//sigma
	// Allocate 3D array
	uint32_t ***sigma = (uint32_t ***)malloc(max_iter * sizeof(uint32_t **));
	for (int i = 0; i < max_iter; i++)
	 {
    	sigma[i] = (uint32_t **)malloc(x_deg * sizeof(uint32_t *));
    	for (int j = 0; j <=x_deg; j++)
		 {
        	sigma[i][j] = (uint32_t *)malloc(alpha_deg * sizeof(uint32_t));
       		 for (int k = 0; k <= alpha_deg; k++) 
				{
            	sigma[i][j][k] = 0;
        		}	
    	}
	}
	
	/*for (int i = 0; i < max_iter; i++)
	 {
    	for (int j = 0; j <=x_deg; j++) 
		{
        for (int k = 0; k <= alpha_deg; k++) {
            printf("sigma[%d][%d][%d] = %u\n", i, j, k, sigma[i][j][k]);
        }
    }
}*/
printf("\n--------------------------\n");
	sigma[0][0][0]=1;
	sigma[1][0][0]=1;
	
	/*for (int i = 0; i < max_iter; i++)
	 {
    	for (int j = 0; j <=x_deg; j++) 
		{
        for (int k = 0; k <= alpha_deg; k++) {
            printf("sigma[%d][%d][%d] = %u\n", i, j, k, sigma[i][j][k]);
        }
    }
}*/
	printf("\n sigma polynomial 0---->");
	 //print_sigma_single(0,x_deg, alpha_deg,sigma,"Sigma");
	 print_bivariate3d(0,x_deg,alpha_deg,sigma);
	 printf("\n Sigma polynomial 1 --->");
	 //print_sigma_single(1,x_deg,alpha_deg,sigma,"Sigma");
	 print_bivariate3d(1,x_deg,alpha_deg,sigma);
	 //-----d_mue
	 uint32_t d[max_iter];
	 for(int i=0;i<max_iter;i++)
	 {
	 	d[i]=0;
	 }
	 /*printf("\n--------------------\n Intializing d to 0 \n\n");
	 for(int i=0;i<max_iter;i++)
	 {
	 	printf("d[%d]=%u \n",i,d[i]);
	   }  
	 	  */
    d[0]=1;
    d[1]=syndromes[0];
    printf("\n--------------------------\n");
    /*for(int i=0;i<max_iter;i++)
    {
    	printf("d[%d]=%u\n",i,d[i]);
	}
    printf("\n--------------------------\n");*/
    //l_mue
    int l[max_iter];
    memset(l,0,sizeof(l));
    //for(int i=0;i<max_iter;i++)
    //	printf("l[%d]=%d \n",i,l[i]);
    //printf("\n--------------------------\n");
    for(int i=0;i<max_iter;i++)
    {
    	l[i]=find_degree_of_x_for_i(x_deg, alpha_deg, sigma[i]);
    	printf("\n l[%d]=%d",i,l[i]);
	}
    printf("\n---------------------------------\n");
	int update_cri[max_iter];
	memset(update_cri,0,sizeof(update_cri));
	
	
	for(int i=0;i<max_iter;i++)
	{
		update_cri[i]=2*mue[i]-l[i];
		printf("\n 2mu-l_mue[%d]=%d \n",i,update_cri[i]);
	}
	for(int mu=1;mu<=max_iter;mu++)
	{
		printf("\n -----------------------  mu=%d  -----------------------------\n",mu);
		int rho=find_best_rho1(mu,d,l);
		uint32_t d_inv=gf_inverse(d[rho],irred_poly,m);
		printf("\n Inverse of d[rho]--->%u",d_inv);
		uint32_t term1=gf_multiply(d[mu],d_inv,irred_poly,m);
		//uint32_t term1=d[mu]*d_inv;
		printf("\n d[%d]=%u \n",mu,d[mu]);
		printf("\n term1=d[mu]*d_inv---> %u\n",term1);
		int exponent=2*(mue[mu]-mue[rho]);  // eg: exponent = 4
   		 // Step 2: Create x^exponent as monomial
   		uint32_t x_pow = 1 << exponent;  // eg: x^4 = 0b00010000
   		int x_deg1=exponent;
    	int alpha_deg1=degree(term1);
    	printf("\n x_pow degree --->%d",x_deg1);
    	printf("\n alpha_deg1 degree ---->%d\n",alpha_deg1);
    	if(alpha_deg1==-1) alpha_deg1=1;
		printf("\n alpha_deg1 degree ---->%d\n",alpha_deg1);
		uint32_t **result = malloc(x_deg1 * sizeof(uint32_t *));
		for (int i = 0; i <=x_deg1; i++)
    		result[i] = malloc(alpha_deg1 * sizeof(uint32_t));
    	
    	//initiliazing result[][] to zeroes
    	for(int i=0;i<=x_deg1;i++)
    	{
    		for(int j=0;j<=alpha_deg1;j++)
    		{
    			result[i][j]=0;
			}
		}
		printf("\n term1 --->%u \n",term1);
		printf("\n x_pow ---->%u \n",x_pow);
		bivariate_mult(term1,x_pow,x_deg1,alpha_deg1,result);
		printf("\n ----------------------- \n Result \n");
		for(int i=0;i<=x_deg1;i++)
		{
			for(int j=0;j<=alpha_deg1;j++)
			{
				printf("Result[%d][%d]=%u \n",i,j,result[i][j]);
			}
		}
		print_bivariate2d(x_deg1, alpha_deg1,result);
		int depth=2;
		uint32_t ***temp = (uint32_t ***)malloc(depth * sizeof(uint32_t **));
		for (int i = 0; i < depth; i++)
		 {
    		temp[i] = (uint32_t **)malloc((2*x_deg1 + 1) * sizeof(uint32_t *));
    		for (int j = 0; j <= 2*x_deg1; j++)
			 {
       			 temp[i][j] = (uint32_t *)malloc((2*alpha_deg1 + 1) * sizeof(uint32_t));
        		for (int k = 0; k <= 2*alpha_deg1; k++)
				 {
            			temp[i][j][k] = 0;
        		}
    		}
		}
		int index=0;
		multiply_2d_with_3d(x_deg1,alpha_deg1,result,index,x_deg,alpha_deg,sigma,temp,irred_poly,m);
		printf("\n-------------------------------------------\n Temp of sigm*result =\n");
		for(int j=0;j<=2*x_deg1;j++)
		{
			for(int k=0;k<=2*alpha_deg1;k++)
			{
				if(temp[0][j][k])
					printf("\n temp[%d][%d][%d]=%u",j,k,temp[0][j][k]);
			}
		}
		printf("\n Polynomial form of above--->");
		print_bivariate3d(0,x_deg1,alpha_deg1,temp);
		printf("\n-------------------------\n Adding \n");
		int max_x_deg = (x_deg > x_deg1) ? x_deg : x_deg1;
		int max_alpha_deg = (alpha_deg > alpha_deg1) ? alpha_deg : alpha_deg1;
	
		uint32_t ***temp1 = (uint32_t ***)malloc(depth * sizeof(uint32_t **));
		for (int i = 0; i < depth; i++)
		 {
    		temp1[i] = (uint32_t **)malloc((max_x_deg + 1) * sizeof(uint32_t *));
    		for (int j = 0; j <= max_x_deg; j++) 
			{
        		temp1[i][j] = (uint32_t *)malloc((max_alpha_deg + 1) * sizeof(uint32_t));
        		for (int k = 0; k <= max_alpha_deg; k++)
				 {
				 	temp1[i][j][k] = 0;
        		}
    		}
		}
		

		for (int j = 0; j <= max_x_deg; j++) {
    		for (int k = 0; k <= max_alpha_deg; k++) {
        	uint32_t a = 0, b = 0;
        	if (j <= x_deg && k <= alpha_deg)
        	{
            	a = sigma[1][j][k];
            	//printf("\nSigma[1][%d][%d]=%u",j,k,sigma[1][j][k]);
        }
        	if (j <= x_deg1 && k <= alpha_deg1){
            	b = temp[0][j][k];
            	//printf("\nTemp[0][%d][%d]=%u",j,k,temp[0][j][k]);
            }
        temp1[0][j][k] = a ^ b;
        //printf("\n temp1[0][%d][%d]=%u",j,k,temp1[0][j][k]);
    	}
	}	
	printf("\n final polynomial -->");
	print_bivariate3d(0,max_x_deg,max_alpha_deg,temp1);

	// Final print to check temp1 after loops
    printf("\nFinal temp1 values:\n");
    for (int j = 0; j <= max_x_deg; j++) {
        for (int k = 0; k <= max_alpha_deg; k++) {
        	if(temp1[0][j][k])
           		 printf("temp1[0][%d][%d] = %u\n", j, k, temp1[0][j][k]);
        }
    }
    printf("\nResult polynomial ---------->");
    print_bivariate3d(0,max_x_deg,max_alpha_deg,temp1);
    printf("\n Copying to sigma[2] \n");
	for (int j = 0; j <= max_x_deg; j++) {
    for (int k = 0; k <= max_alpha_deg; k++) {
        sigma[mu + 1][j][k] = temp1[0][j][k];
    	}
	}
	print_bivariate3d(mu+1,max_x_deg,max_alpha_deg,sigma);
	l[mu+1]=find_degree_of_x_for_i(max_x_deg, max_alpha_deg, sigma[mu+1]);
	printf("\n l[mu+1]=l[%d]-->%d",mu+1,l[mu+1]);
	
	/*
	// Loop through the coefficients of sigma^(mu+1)
	for (int i = 0; i <= l[mu+1]; i++) 
	{
    	// Get the coefficient of x^i in sigma^(mu+1), i.e., sigma_i^(mu+1)
    	uint32_t sigma_i_mu_plus_1 = sigma[mu + 1][i];  // Adjust if needed based on array indexing
    	printf("\n Sigma_i_mu_plus_1=>%u",sigma_i_mu_plus_1)

    	// Calculate the corresponding syndrome index (2*mu + 3 - i)
    	int syndrome_index = 2*mue[mu] + 3 - i;

    	// Get the corresponding syndrome value, S[syndrome_index]
    	uint32_t syndrome_value = syndromes[syndrome_index];  // Make sure S is indexed correctly

    	// Add to d_mu_plus_1 (d_{mu+1})
    	d_mu_plus_1 ^= gf_multiply(sigma_i_mu_plus_1, syndrome_value, irred_poly, m);
	}

// Now d_mu_plus_1 holds the calculated value
printf("d_(%d) = %u\n", mu+1, d_mu_plus_1);	*/
	update_cri[mu+1]=2*mue[mu+1]-l[mu+1];



	
	for (int i = 0; i < depth; i++) {
    	for (int j = 0; j <= 2*x_deg1; j++) {
        	free(temp[i][j]);
    	}
    	free(temp[i]);
	}
	free(temp);

	for (int i = 0; i < depth; i++) {
    	for (int j = 0; j <= 2*x_deg1; j++) {
        	free(temp1[i][j]);
    	}
    	free(temp1[i]);
	}
	free(temp1);
// Assuming result is a dynamically allocated 2D array:
// result = malloc(x_deg1 * sizeof(uint32_t *));
// and each row is allocated as: result[i] = malloc(alpha_deg1 * sizeof(uint32_t));

// Free the memory allocated for result
for (int i = 0; i <= x_deg1; i++) {
    free(result[i]);  // Free each row (result[i] is a uint32_t*)
}

free(result);  // Free the top-level array (result is a uint32_t**)







	}
//freeing the memory
	 for (int i = 0; i < max_iter; i++) {
    for (int j = 0; j <=x_deg; j++) {
        free(sigma[i][j]);
    }
    free(sigma[i]);
	}
	free(sigma);



	
}
